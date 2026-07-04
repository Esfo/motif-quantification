"""Tab 2 — Protein viewing.

Shows identified protein sequences residue-by-residue (N→C), with the tryptic
peptides the search *attempted* drawn as outlined rectangles whose background is
coloured by q-value (best/green → worst/red on the shared q-value gradient).

Layout mirrors the MS Data tab: a file selector at the top-left, an FDR spin box
just below it (the accepted-FDR that governs the protein list, defaulting to the
search's ``q_max``), and a single protein list on the left. The right side is a
vertical split — panel 1 (top) over panel 2 (bottom), a plain horizontal divider
between them, exactly like MS Data's panels.

* **Panel 1** renders the selected protein horizontally, wrapping to new lines.
  Each base tryptic peptide (the protease's cut segments) is an outlined
  rectangle; its fill is the best q-value of any peptide identified over those
  residues in the current file. Segments too short/long to be searched carry no
  outline or fill (plain letters). An **All** button combines the identifications
  across every file into one consensus colouring.
* **Panel 2** renders the same protein *vertically* (N at top → C at bottom), one
  column per file involved in the search, so a protein's identification can be
  compared across files at a glance. Clicking a column loads that file into
  panel 1. File names are drawn on a stylish 45° slant.

Sequences come from the project FASTA (``ViewerSession.protein_sequence``); the
reorganized tables carry none. The digestion parameters live in ``session.py``
and mirror the Sage enzyme config in execution.xsh.
"""

import math

from PySide6.QtCore import Qt, QRectF, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

try:
    from .session import digest_peptides, cleavage_sites, DIGEST_MIN_LEN, DIGEST_MAX_LEN
except ImportError:
    from session import digest_peptides, cleavage_sites, DIGEST_MIN_LEN, DIGEST_MAX_LEN


# q-value → colour on the shared gradient: best (low q) is green, worst (high q)
# is red, through yellow. Mapped on log10(q) between 1e-4 and 0.05 so the useful
# FDR band spans the full green→red sweep. q is None (identified, no q reported)
# counts as the best (green).
_Q_LOG_LO = -4.0
_Q_LOG_HI = math.log10(0.05)


def q_color(q, alpha=190):
    if q is None:
        t = 0.0
    else:
        q = max(float(q), 1e-6)
        t = (math.log10(q) - _Q_LOG_LO) / (_Q_LOG_HI - _Q_LOG_LO)
        t = min(1.0, max(0.0, t))
    # green (t=0) → yellow (t=0.5) → red (t=1)
    r = int(255 * min(1.0, t * 2.0))
    g = int(255 * min(1.0, (1.0 - t) * 2.0))
    return QColor(r, g, 40, alpha)


def base_segments(seq):
    """The fully-cleaved (0 missed cleavage) tryptic segments tiling ``seq``.

    Returns ``[(start, end, searchable), ...]`` covering the whole protein;
    ``searchable`` is False for segments outside the search's length bounds (the
    peptides the search never attempted → drawn plain, no rectangle)."""
    if not seq:
        return []
    sites = cleavage_sites(seq)
    out = []
    for i in range(len(sites) - 1):
        a, b = sites[i], sites[i + 1]
        searchable = DIGEST_MIN_LEN <= (b - a) <= DIGEST_MAX_LEN
        out.append((a, b, searchable))
    return out


def residue_q(seq, qmap):
    """Per-residue best (min) q-value from the identified peptides in ``qmap``.

    ``qmap`` maps plain peptide sequence → q (or None). Every in-silico peptide
    identified in a file is located on the protein and paints its residues with
    its q; overlaps keep the best. Residues with no identification stay None.
    A separate boolean array marks residues touched by *any* identification so a
    None-q (identified but unscored) peptide is still distinguishable from a miss.
    """
    n = len(seq)
    qs = [None] * n
    hit = [False] * n
    if not qmap:
        return qs, hit
    for _start, _end, pep in digest_peptides(seq):
        if pep not in qmap:
            continue
        q = qmap[pep]
        # locate every occurrence (a peptide can repeat within a protein)
        start = 0
        while True:
            idx = seq.find(pep, start)
            if idx < 0:
                break
            for i in range(idx, idx + len(pep)):
                hit[i] = True
                if q is not None and (qs[i] is None or q < qs[i]):
                    qs[i] = q
            start = idx + 1
    return qs, hit


def _row_runs(a, b, cols):
    """Split residue range [a, b) into per-wrap-row contiguous runs.

    Yields ``(row, col_start, count)`` so a peptide that wraps across the flow's
    line breaks is drawn as one rectangle per visual line."""
    i = a
    while i < b:
        row = i // cols
        col = i % cols
        count = min(b, (row + 1) * cols) - i
        yield row, col, count
        i += count


class HorizontalSequenceView(QWidget):
    """Panel 1: the protein spelled out horizontally, wrapping, with peptide
    rectangles coloured by q-value."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._seq = ""
        self._segments = []
        self._qs = []
        self._hit = []
        self._fg = QColor("#e6e6e6")
        self._panel_bg = QColor("#101216")
        self.setMinimumHeight(80)
        self._metrics()

    def _metrics(self):
        font = QFont("monospace")
        font.setStyleHint(QFont.Monospace)
        font.setPointSize(11)
        self._font = font
        fm = QFontMetrics(font)
        self.cw = max(12, fm.horizontalAdvance("W") + 6)
        self.ch = fm.height() + 8
        self.margin = 10

    def set_theme(self, fg, bg):
        self._fg = QColor(fg)
        self._panel_bg = QColor(bg)
        self.update()

    def set_protein(self, seq, qmap):
        self._seq = seq or ""
        self._segments = base_segments(self._seq)
        self._qs, self._hit = residue_q(self._seq, qmap or {})
        self._recompute_height()
        self.update()

    def _cols(self):
        usable = max(1, self.width() - 2 * self.margin)
        return max(1, usable // self.cw)

    def _recompute_height(self):
        n = len(self._seq)
        if n == 0:
            self.setMinimumHeight(80)
            return
        rows = math.ceil(n / self._cols())
        self.setMinimumHeight(rows * self.ch + 2 * self.margin)

    def resizeEvent(self, event):
        self._recompute_height()
        super().resizeEvent(event)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, False)
        painter.fillRect(self.rect(), self._panel_bg)
        if not self._seq:
            painter.setPen(self._fg)
            painter.drawText(self.rect(), Qt.AlignCenter,
                             "select a protein to view its sequence")
            painter.end()
            return

        cols = self._cols()
        painter.setFont(self._font)

        def cell(i):
            row = i // cols
            col = i % cols
            return self.margin + col * self.cw, self.margin + row * self.ch

        # 1) peptide rectangles (background fill + outline), behind the letters.
        for a, b, searchable in self._segments:
            if not searchable:
                continue
            # best q over the segment; None if identified-but-unscored; skip fill
            # entirely when no residue of the segment was identified.
            seg_q = None
            covered = False
            for i in range(a, b):
                if self._hit[i]:
                    covered = True
                    if self._qs[i] is not None and (seg_q is None or self._qs[i] < seg_q):
                        seg_q = self._qs[i]
            fill = q_color(seg_q) if covered else None
            for row, col, count in _row_runs(a, b, cols):
                x = self.margin + col * self.cw
                y = self.margin + row * self.ch
                rect = QRectF(x, y + 2, count * self.cw, self.ch - 4)
                if fill is not None:
                    painter.fillRect(rect, fill)
                painter.setPen(QPen(self._fg, 1))
                painter.drawRect(rect)

        # 2) residue letters on top (clip to the exposed region for long proteins).
        exposed = event.rect()
        painter.setPen(self._fg)
        for i, chr_ in enumerate(self._seq):
            x, y = cell(i)
            if y > exposed.bottom() or y + self.ch < exposed.top():
                continue
            painter.drawText(QRectF(x, y, self.cw, self.ch), Qt.AlignCenter, chr_)
        painter.end()


class VerticalMultiFileView(QWidget):
    """Panel 2: the protein spelled top→bottom (N→C), one column per file, each
    column coloured by that file's identification. Clicking a column emits its
    filename so panel 1 can switch to it."""

    file_clicked = Signal(str)

    HEADER_H = 116

    def __init__(self, parent=None):
        super().__init__(parent)
        self._seq = ""
        self._segments = []
        self._files = []          # [(filename, qs, hit)]
        self._fg = QColor("#e6e6e6")
        self._panel_bg = QColor("#101216")
        self._current = None
        self._metrics()

    def _metrics(self):
        font = QFont("monospace")
        font.setStyleHint(QFont.Monospace)
        font.setPointSize(10)
        self._font = font
        fm = QFontMetrics(font)
        self.ch = fm.height() + 2
        self.col_w = max(20, fm.horizontalAdvance("W") + 14)
        self.gap = 10
        self.margin = 10

    def set_theme(self, fg, bg):
        self._fg = QColor(fg)
        self._panel_bg = QColor(bg)
        self.update()

    def set_protein(self, seq, files_qmaps, current=None):
        """``files_qmaps`` = ``[(filename, qmap), ...]`` over every file."""
        self._seq = seq or ""
        self._segments = base_segments(self._seq)
        self._current = current
        self._files = []
        for filename, qmap in files_qmaps:
            qs, hit = residue_q(self._seq, qmap or {})
            self._files.append((filename, qs, hit))
        self._recompute_size()
        self.update()

    def set_current(self, current):
        self._current = current
        self.update()

    def _col_pitch(self):
        return self.col_w + self.gap

    def _recompute_size(self):
        n = len(self._seq)
        width = 2 * self.margin + max(1, len(self._files)) * self._col_pitch()
        height = self.HEADER_H + n * self.ch + self.margin
        self.setMinimumSize(int(width), int(height))

    def mousePressEvent(self, event):
        if not self._files:
            return
        x = event.position().x() if hasattr(event, "position") else event.x()
        idx = int((x - self.margin) // self._col_pitch())
        if 0 <= idx < len(self._files):
            self.file_clicked.emit(self._files[idx][0])

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), self._panel_bg)
        if not self._seq or not self._files:
            painter.setPen(self._fg)
            painter.drawText(self.rect(), Qt.AlignCenter,
                             "select a protein to compare across files")
            painter.end()
            return

        exposed = event.rect()
        pitch = self._col_pitch()
        for fidx, (filename, qs, hit) in enumerate(self._files):
            x0 = self.margin + fidx * pitch

            # highlight the column currently shown in panel 1
            if filename == self._current:
                painter.fillRect(QRectF(x0 - 3, 0, self.col_w + 6, self.height()),
                                 QColor(120, 150, 210, 40))

            # 45°-slanted file label
            painter.save()
            painter.translate(x0 + self.col_w / 2.0, self.HEADER_H - 8)
            painter.rotate(-45)
            painter.setPen(self._fg)
            label_font = QFont(self._font)
            label_font.setPointSize(9)
            painter.setFont(label_font)
            painter.drawText(0, 0, _short_name(filename))
            painter.restore()

            # per-residue background from base segments + letters
            painter.setFont(self._font)
            for a, b, searchable in self._segments:
                if not searchable:
                    continue
                seg_q = None
                covered = False
                for i in range(a, b):
                    if hit[i]:
                        covered = True
                        if qs[i] is not None and (seg_q is None or qs[i] < seg_q):
                            seg_q = qs[i]
                y_top = self.HEADER_H + a * self.ch
                y_bot = self.HEADER_H + b * self.ch
                rect = QRectF(x0, y_top + 1, self.col_w, (b - a) * self.ch - 2)
                if y_bot < exposed.top() or y_top > exposed.bottom():
                    continue
                if covered:
                    painter.fillRect(rect, q_color(seg_q))
                painter.setPen(QPen(self._fg, 1))
                painter.drawRect(rect)

            for i, chr_ in enumerate(self._seq):
                y = self.HEADER_H + i * self.ch
                if y > exposed.bottom() or y + self.ch < exposed.top():
                    continue
                painter.setPen(self._fg)
                painter.drawText(QRectF(x0, y, self.col_w, self.ch),
                                 Qt.AlignCenter, chr_)
        painter.end()


def _short_name(filename):
    name = filename
    for suffix in (".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name if len(name) <= 28 else name[:25] + "…"


class ProteinsTab(QWidget):
    """Tab 2 assembled: file selector + FDR + protein list on the left, panel 1
    over panel 2 on the right."""

    def __init__(self, session, theme="dark"):
        super().__init__()
        self.session = session
        self.theme = theme
        self.current_file = None
        self.current_protein = None
        self._combined = False   # panel 1 "All" mode

        self._build()
        self._populate_files()
        self._apply_theme(theme)
        self._refresh_protein_list()

    # ---- construction ----------------------------------------------------

    def _build(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)

        # left column: file selector, FDR, protein list
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(2, 2, 2, 2)

        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self._on_file_changed)
        left_layout.addWidget(self.file_combo)

        # FDR acceptance criteria (percent). Governs the protein list; defaults
        # to the search's q_max. Placed near the left, like MS Data's spin box.
        fdr_row = QHBoxLayout()
        self.fdr_edit = QDoubleSpinBox()
        self.fdr_edit.setDecimals(2)
        self.fdr_edit.setRange(0.0, 100.0)
        self.fdr_edit.setSingleStep(0.1)
        self.fdr_edit.setValue(self._default_fdr_percent())
        self.fdr_edit.setFixedWidth(70)
        self.fdr_edit.setToolTip("Proteins are listed when identified at or below "
                                 "this FDR (percent)")
        self.fdr_edit.valueChanged.connect(self._on_fdr_changed)
        fdr_row.addWidget(self.fdr_edit)
        fdr_row.addWidget(QLabel("% FDR"))
        fdr_row.addStretch(1)
        left_layout.addLayout(fdr_row)

        left_layout.addWidget(QLabel("proteins"))
        self.protein_list = QListWidget()
        self.protein_list.setStyleSheet(
            "QListWidget::item:selected,"
            "QListWidget::item:selected:!active"
            " { background-color: #2f6fb3; color: white; }"
        )
        self.protein_list.itemSelectionChanged.connect(self._on_protein_selected)
        left_layout.addWidget(self.protein_list, stretch=1)
        left.setMaximumWidth(320)

        # right side: panel 1 over panel 2
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)

        # panel 1 header bar with the All button
        bar = QHBoxLayout()
        bar.setContentsMargins(4, 2, 4, 2)
        self.all_btn = QPushButton("All")
        self.all_btn.setToolTip("Combine identifications across every file into a "
                                "single consensus colouring")
        self.all_btn.setCheckable(True)
        self.all_btn.setFixedWidth(48)
        self.all_btn.clicked.connect(self._on_all_clicked)
        self.p1_title = QLabel("")
        bar.addWidget(self.all_btn)
        bar.addWidget(self.p1_title)
        bar.addStretch(1)
        right_layout.addLayout(bar)

        self.panel1 = HorizontalSequenceView()
        p1_scroll = QScrollArea()
        p1_scroll.setWidgetResizable(True)
        p1_scroll.setWidget(self.panel1)
        p1_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.panel2 = VerticalMultiFileView()
        self.panel2.file_clicked.connect(self._on_panel2_file)
        p2_scroll = QScrollArea()
        p2_scroll.setWidgetResizable(True)
        p2_scroll.setWidget(self.panel2)

        splitter = QSplitter(Qt.Vertical)
        splitter.setHandleWidth(2)   # a plain horizontal divider, like MS Data
        splitter.addWidget(p1_scroll)
        splitter.addWidget(p2_scroll)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        right_layout.addWidget(splitter, stretch=1)

        root.addWidget(left)
        divider = QFrame()
        divider.setFrameShape(QFrame.VLine)
        root.addWidget(divider)
        root.addWidget(right, stretch=1)

    def _default_fdr_percent(self):
        q_max = self.session.summary().get("q_max") if self.session else None
        try:
            return max(0.0, float(q_max) * 100.0)
        except (TypeError, ValueError):
            return 1.0

    # ---- data ------------------------------------------------------------

    def _all_files(self):
        return [row.get("filename", "") for row in self.session.files()
                if row.get("filename")]

    def _populate_files(self):
        self.file_combo.blockSignals(True)
        self.file_combo.clear()
        for name in self._all_files():
            self.file_combo.addItem(name, name)
        self.file_combo.blockSignals(False)
        self.current_file = self.file_combo.currentData()

    def _fdr_threshold(self):
        return max(0.0, self.fdr_edit.value()) / 100.0

    def _identified_proteins(self, filename):
        """Proteins with a peptide identified at/below the FDR in ``filename``,
        as ``[(protein_id, row), ...]`` sorted by best protein q."""
        thr = self._fdr_threshold()
        out = []
        for row in self.session.file_proteins(filename or ""):
            pid = row.get("protein_id", "")
            if not pid:
                continue
            q = _safe_float(row.get("protein_q"))
            if q is not None and q > thr:
                continue
            out.append((q if q is not None else 1.0, pid, row))
        out.sort(key=lambda t: t[0])
        return [(pid, row) for _q, pid, row in out]

    # ---- events ----------------------------------------------------------

    def _on_file_changed(self):
        self.current_file = self.file_combo.currentData()
        self._refresh_protein_list()

    def _on_fdr_changed(self):
        self._refresh_protein_list()

    def _refresh_protein_list(self):
        prev = self.current_protein
        self.protein_list.blockSignals(True)
        self.protein_list.clear()
        for pid, _row in self._identified_proteins(self.current_file):
            item = QListWidgetItem(pid)
            item.setData(Qt.UserRole, pid)
            self.protein_list.addItem(item)
        self.protein_list.blockSignals(False)
        # keep the selection if the protein is still listed
        if prev is not None:
            for i in range(self.protein_list.count()):
                if self.protein_list.item(i).data(Qt.UserRole) == prev:
                    self.protein_list.setCurrentRow(i)
                    break
        self._update_panels()

    def _on_protein_selected(self):
        items = self.protein_list.selectedItems()
        self.current_protein = items[0].data(Qt.UserRole) if items else None
        self._update_panels()

    def _on_all_clicked(self):
        self._combined = self.all_btn.isChecked()
        self._update_panel1()

    def _on_panel2_file(self, filename):
        # clicking a panel-2 column brings that file up in panel 1
        self._combined = False
        self.all_btn.setChecked(False)
        idx = self.file_combo.findData(filename)
        if idx >= 0:
            self.file_combo.setCurrentIndex(idx)   # triggers list refresh
        else:
            self.current_file = filename
            self._update_panels()

    # ---- rendering -------------------------------------------------------

    def _combined_qmap(self, seq):
        """Best q per identified peptide across every file (the All view)."""
        combined = {}
        for filename in self._all_files():
            for pep, q in self.session.peptide_q_for_file(filename).items():
                if pep not in combined:
                    combined[pep] = q
                elif q is not None and (combined[pep] is None or q < combined[pep]):
                    combined[pep] = q
        return combined

    def _update_panels(self):
        self._update_panel1()
        self._update_panel2()

    def _update_panel1(self):
        seq = self.session.protein_sequence(self.current_protein) if self.current_protein else None
        if not seq:
            self.panel1.set_protein("", {})
            if self.current_protein and not self.session.has_sequences:
                self.p1_title.setText("no FASTA found — cannot show sequences")
            else:
                self.p1_title.setText("")
            return
        if self._combined:
            qmap = self._combined_qmap(seq)
            self.p1_title.setText(f"{self.current_protein}   ·   All files "
                                  f"({len(seq)} aa)")
        else:
            qmap = self.session.peptide_q_for_file(self.current_file or "")
            self.p1_title.setText(f"{self.current_protein}   ·   "
                                  f"{_short_name(self.current_file or '')} ({len(seq)} aa)")
        self.panel1.set_protein(seq, qmap)

    def _update_panel2(self):
        seq = self.session.protein_sequence(self.current_protein) if self.current_protein else None
        if not seq:
            self.panel2.set_protein("", [])
            return
        files_qmaps = [(f, self.session.peptide_q_for_file(f)) for f in self._all_files()]
        self.panel2.set_protein(seq, files_qmaps, current=self.current_file)

    # ---- theme -----------------------------------------------------------

    def apply_theme(self, theme):
        self._apply_theme(theme)

    def _apply_theme(self, theme):
        self.theme = theme
        if theme == "light":
            fg, bg = "#202020", "#fafafa"
        else:
            fg, bg = "#e6e6e6", "#101216"
        self.panel1.set_theme(fg, bg)
        self.panel2.set_theme(fg, bg)


def _safe_float(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None
