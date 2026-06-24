"""Tab widgets that visualize the reorganized Sage/Percolator search output.

These views render the tables that reorganize-results.py writes
(files.tsv, peptides.tsv, proteins.tsv and the per-file quant tables) plus
the LFQ quantities Sage produced. Nothing here computes spectra, ions, or
isotope distributions; it only surfaces values already on disk.
"""

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .plots import plot_bars
    from .session import safe_float
except ImportError:
    from plots import plot_bars
    from session import safe_float


FRACTION_COLORS = {
    "GuHCl": "#4c72b0",
    "NaCl": "#dd8452",
    "": "#888888",
}


def make_item(value):
    """Table item that sorts numerically when the value looks like a number."""
    text = "" if value is None else str(value)
    item = QTableWidgetItem(text)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    number = safe_float(text)

    if number is not None and text.strip() != "":
        item.setData(Qt.EditRole, number)
        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

    return item


def fill_table(table, columns, rows):
    table.setSortingEnabled(False)
    table.setColumnCount(len(columns))
    table.setHorizontalHeaderLabels([label for _, label in columns])
    table.setRowCount(len(rows))

    for row_i, row in enumerate(rows):
        for col_i, (field, _) in enumerate(columns):
            item = make_item(row.get(field, ""))
            item.setData(Qt.UserRole, row_i)
            table.setItem(row_i, col_i, item)

    table.resizeColumnsToContents()
    table.setSortingEnabled(True)


def configure_table(table):
    table.setSelectionBehavior(QAbstractItemView.SelectRows)
    table.setSelectionMode(QAbstractItemView.SingleSelection)
    table.setEditTriggers(QAbstractItemView.NoEditTriggers)
    table.setAlternatingRowColors(True)


def filter_rows(rows, text, fields):
    text = text.strip().lower()

    if not text:
        return rows

    out = []

    for row in rows:
        for field in fields:
            if text in str(row.get(field, "")).lower():
                out.append(row)
                break

    return out


class OverviewView(QWidget):
    """Per-file run summary from files.tsv, with a bar chart of counts."""

    COLUMNS = [
        ("filename", "file"),
        ("fraction", "fraction"),
        ("replicate", "replicate"),
        ("n_psms", "PSMs"),
        ("n_peptides", "peptides"),
        ("n_proteins", "proteins"),
        ("n_quant_rows", "quant rows"),
        ("n_scans", "scans"),
    ]

    def __init__(self, session):
        super().__init__()
        self.session = session

        self.summary_label = QLabel()
        self.summary_label.setWordWrap(True)

        self.table = QTableWidget()
        configure_table(self.table)

        self.plot = pg.PlotWidget()

        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.table)
        splitter.addWidget(self.plot)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.summary_label)
        layout.addWidget(splitter, stretch=1)

        self.load()

    def load(self):
        rows = self.session.files()
        fill_table(self.table, self.COLUMNS, rows)

        if not rows:
            self.summary_label.setText(
                "No folder open. Double-click anywhere or use File ▸ Open reorganized folder… (Ctrl+O)."
            )
            return

        summary = self.session.summary() or {}
        counts = summary.get("rows", {})

        parts = [f"reorganized: {self.session.reorganized}"]

        if counts:
            parts.append(
                "files {files} · psms {psms} · peptides {global_peptides} · "
                "proteins {global_proteins} · quant {quant}".format(
                    files=counts.get("files", "?"),
                    psms=counts.get("psms", "?"),
                    global_peptides=counts.get("global_peptides", "?"),
                    global_proteins=counts.get("global_proteins", "?"),
                    quant=counts.get("quant", "?"),
                )
            )

        if summary.get("q_max") is not None:
            parts.append(f"q ≤ {summary.get('q_max')}")

        self.summary_label.setText("\n".join(parts))

        labels = []
        peptide_counts = []
        colors = []

        for row in rows:
            labels.append(row.get("filename", ""))
            peptide_counts.append(safe_float(row.get("n_peptides"), 0.0) or 0.0)
            colors.append(FRACTION_COLORS.get(row.get("fraction", ""), "#4c72b0"))

        self.plot.clear()
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setTitle("identified peptides per file (colored by fraction)")
        self.plot.setLabel("left", "peptides")

        if labels:
            x = list(range(len(labels)))
            bars = pg.BarGraphItem(
                x=x,
                height=peptide_counts,
                width=0.6,
                brushes=colors,
                pen=pg.mkPen("#22222288"),
            )
            self.plot.addItem(bars)
            self.plot.getAxis("bottom").setTicks(
                [[(i, label) for i, label in enumerate(labels)]]
            )


class PeptidesView(QWidget):
    """Global peptide table with per-file LFQ quantity across runs."""

    COLUMNS = [
        ("peptide", "peptide"),
        ("proteins", "proteins"),
        ("n_psms", "PSMs"),
        ("n_files", "files"),
        ("percolator_q", "perc q"),
        ("percolator_score", "perc score"),
        ("has_lfq", "lfq"),
    ]

    def __init__(self, session):
        super().__init__()
        self.session = session
        self.rows = []

        self.search = QLineEdit()
        self.search.setPlaceholderText("filter peptides or proteins…")
        self.search.textChanged.connect(self.apply_filter)

        self.table = QTableWidget()
        configure_table(self.table)
        self.table.itemSelectionChanged.connect(self.on_selected)

        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)

        self.plot = pg.PlotWidget()

        right = QSplitter(Qt.Vertical)
        right.addWidget(self.detail)
        right.addWidget(self.plot)
        right.setStretchFactor(0, 1)
        right.setStretchFactor(1, 2)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self.search)
        left_layout.addWidget(self.table, stretch=1)

        main = QSplitter(Qt.Horizontal)
        main.addWidget(left)
        main.addWidget(right)
        main.setStretchFactor(0, 3)
        main.setStretchFactor(1, 2)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(main)

        self.load()

    def load(self):
        self.rows = self.session.global_peptides()
        self.apply_filter()

    def apply_filter(self):
        filtered = filter_rows(self.rows, self.search.text(), ["peptide", "proteins"])
        self.filtered = filtered
        fill_table(self.table, self.COLUMNS, filtered)

    def selected_row(self):
        items = self.table.selectedItems()

        if not items:
            return None

        row_i = items[0].data(Qt.UserRole)

        if row_i is None or row_i >= len(self.filtered):
            return None

        return self.filtered[row_i]

    def on_selected(self):
        row = self.selected_row()

        if row is None:
            return

        peptide_plain = row.get("peptide_plain", "") or row.get("peptide", "")

        lines = [
            f"peptide: {row.get('peptide', '')}",
            f"plain: {peptide_plain}",
            f"proteins: {row.get('proteins', '')}",
            f"best PSM: {row.get('best_psm_id', '')} in {row.get('best_filename', '')}",
            f"percolator q: {row.get('percolator_q', '')}  score: {row.get('percolator_score', '')}",
            f"PSMs: {row.get('n_psms', '')}  files: {row.get('n_files', '')}",
            f"files seen: {row.get('files', '')}",
        ]

        totals = self.session.quant_for_peptide(peptide_plain)
        file_order = [f.get("filename", "") for f in self.session.files()]

        labels = []
        values = []

        for filename in file_order:
            labels.append(filename)
            values.append(totals.get(filename, 0.0))

        if any(values):
            lines.append("")
            lines.append("LFQ quantity by file (summed over charge):")

            for filename in file_order:
                if totals.get(filename):
                    lines.append(f"  {filename}: {totals[filename]:.4g}")
        else:
            lines.append("")
            lines.append("no LFQ quantity recorded for this peptide")

        self.detail.setPlainText("\n".join(lines))

        plot_bars(
            self.plot,
            labels,
            values,
            title=f"LFQ quantity: {peptide_plain}",
            y_label="quantity",
        )


class ProteinsView(QWidget):
    """Global protein table with its peptide members and supporting files."""

    COLUMNS = [
        ("protein_id", "protein"),
        ("n_percolator_peptides", "peptides"),
        ("n_linked_psms", "PSMs"),
        ("n_files", "files"),
        ("percolator_q", "perc q"),
    ]

    PEPTIDE_COLUMNS = [
        ("peptide", "peptide"),
    ]

    def __init__(self, session):
        super().__init__()
        self.session = session
        self.rows = []

        self.search = QLineEdit()
        self.search.setPlaceholderText("filter proteins…")
        self.search.textChanged.connect(self.apply_filter)

        self.table = QTableWidget()
        configure_table(self.table)
        self.table.itemSelectionChanged.connect(self.on_selected)

        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)

        self.peptide_table = QTableWidget()
        configure_table(self.peptide_table)

        right = QSplitter(Qt.Vertical)
        right.addWidget(self.detail)
        right.addWidget(self.peptide_table)
        right.setStretchFactor(0, 1)
        right.setStretchFactor(1, 2)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self.search)
        left_layout.addWidget(self.table, stretch=1)

        main = QSplitter(Qt.Horizontal)
        main.addWidget(left)
        main.addWidget(right)
        main.setStretchFactor(0, 3)
        main.setStretchFactor(1, 2)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(main)

        self.load()

    def load(self):
        self.rows = self.session.global_proteins()
        self.apply_filter()

    def apply_filter(self):
        filtered = filter_rows(self.rows, self.search.text(), ["protein_id", "peptides"])
        self.filtered = filtered
        fill_table(self.table, self.COLUMNS, filtered)

    def selected_row(self):
        items = self.table.selectedItems()

        if not items:
            return None

        row_i = items[0].data(Qt.UserRole)

        if row_i is None or row_i >= len(self.filtered):
            return None

        return self.filtered[row_i]

    def on_selected(self):
        row = self.selected_row()

        if row is None:
            return

        lines = [
            f"protein: {row.get('protein_id', '')}",
            f"group: {row.get('protein_group_id', '')}",
            f"percolator q: {row.get('percolator_q', '')}  pep: {row.get('percolator_pep', '')}",
            f"peptides: {row.get('n_percolator_peptides', '')}  linked PSMs: {row.get('n_linked_psms', '')}",
            f"files: {row.get('n_files', '')}",
            f"files seen: {row.get('files', '')}",
        ]

        self.detail.setPlainText("\n".join(lines))

        peptides = [p for p in str(row.get("peptides", "")).split(";") if p]
        peptide_rows = [{"peptide": p} for p in peptides]
        fill_table(self.peptide_table, self.PEPTIDE_COLUMNS, peptide_rows)
