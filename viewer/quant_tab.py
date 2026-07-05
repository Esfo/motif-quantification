"""Tab 3 — Quantitative Comparisons.

Compare peptide / protein quantities across files, grouped **entirely** by the
project ``experimental-setup`` file. Nothing about the design is hard-coded: the
tab reads whatever columns the setup file has and treats every one of them as a
generic category. The only special designation is optional and made *by the
user*: marking one column as the **replicate** column (so its runs may be
averaged; every other column is compared, never averaged).

Layout
------
Top half (horizontal split):
  * **left** — a *faceted* view of the selected feature's quantities. You choose,
    by depth, which categories nest the plot (Level 1 splits the panel into
    side-by-side sub-panels, Level 2 splits each of those again, …) and which
    category forms the x-axis of the leaf plots. This is the "organize the axes
    by depth" view: pick conditions first → the panel splits per condition, then
    a time-series/replicate axis inside each, and so on.
  * **right** — a scatter of **every** feature: x = log2 fold change between two
    chosen category values (a plain log difference, no statistical test),
    y = mean log2 abundance. Click a point to select that feature.

Bottom half:
  * a **Peptides ⇄ Proteins** switch and a **unique-only** filter, over a table
    whose columns are the design's files (with their category metadata) showing
    each feature's per-file quantity, plus a unique/non-unique column.

Everything is reactive — changing a category, the replicate column, the contrast,
or the level re-computes immediately; there is no run button.
"""

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .quant_model import QuantModel
    from .theming import palette, style_plot
except ImportError:
    from quant_model import QuantModel
    from theming import palette, style_plot


NONE_LABEL = "(none)"
FILE_LABEL = "(file)"
MAX_FACET_LEVELS = 3


# Local widget-chrome palette (mirrors proteins_tab) so list/label backgrounds
# and outlines adapt to the theme too, not just the pyqtgraph plots.
THEMES = {
    "dark": {"fg": "#e6e6e6", "bg": "#101216", "muted": "#8a8f98",
             "panel": "#16181d", "line": "#2a2d33"},
    "light": {"fg": "#202020", "bg": "#fafafa", "muted": "#606060",
              "panel": "#ffffff", "line": "#c8ccd2"},
}


def _sorted_values(values):
    """Numeric-aware sort of category values."""
    values = list(values)

    def as_num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None
    nums = [as_num(v) for v in values]
    if values and all(n is not None for n in nums):
        return [v for _, v in sorted(zip(nums, values))]
    return sorted(values)


class NumericItem(QTableWidgetItem):
    """Table item sorting by an underlying float (blank/NaN sort last)."""

    def __init__(self, text, value):
        super().__init__(text)
        self.setData(Qt.UserRole, value)

    def __lt__(self, other):
        a = self.data(Qt.UserRole)
        b = other.data(Qt.UserRole) if isinstance(other, NumericItem) else None
        a_bad = a is None or (isinstance(a, float) and math.isnan(a))
        b_bad = b is None or (isinstance(b, float) and math.isnan(b))
        if a_bad:
            return False
        if b_bad:
            return True
        return a < b


class QuantTab(QWidget):
    feature_selected = Signal(str)

    def __init__(self, session, experimental, theme="dark", parent=None):
        super().__init__(parent)
        self.session = session
        self.experimental = experimental
        self.theme = theme

        self.model = QuantModel(session)
        self.level = "peptide"
        self.unique_only = False
        self.selected_feature = None
        self._active = not experimental.is_empty()

        self._categories = [c for c in (experimental.columns() if self._active else [])
                            if c != "filename"]

        self._build_ui()
        if self._active:
            self._refresh_table()
            self._refresh_fold_change()
            self._auto_select_first()
        self.apply_theme(theme)

    # ---- design helpers --------------------------------------------------

    def _replicate_column(self):
        col = self.replicate_combo.currentText()
        return col if col and col != NONE_LABEL else None

    def _row_for(self, filename):
        return self.experimental.by_filename.get(filename, {})

    def _files_with(self, feature):
        """Files where this feature has a positive quantity."""
        per_file = self.model.matrix(self.level).get(feature, {})
        return [f for f, q in per_file.items() if q and q > 0]

    def _visible_features(self):
        """Feature keys currently in scope (respects the unique-only filter for
        peptides; proteins are already unique-quantified)."""
        matrix = self.model.matrix(self.level)
        feats = list(matrix.keys())
        if self.level == "peptide" and self.unique_only:
            feats = [f for f in feats if self.model.peptide_is_unique(f)]
        return sorted(feats)

    # ---- UI construction -------------------------------------------------

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        if not self._active:
            warn = QLabel(
                "No experimental-setup file found for this project.\n\n"
                "Quantitative Comparisons groups runs entirely from the project's "
                "experimental-setup csv (its first column is the mzML filename; "
                "every other column is a category you can compare or facet by). "
                "Add one beside distributions/ and searches/, then reload.")
            warn.setWordWrap(True)
            warn.setAlignment(Qt.AlignCenter)
            outer.addWidget(warn)
            return

        self.v_split = QSplitter(Qt.Vertical)
        outer.addWidget(self.v_split, 1)

        top = QSplitter(Qt.Horizontal)
        top.addWidget(self._build_facet_panel())
        top.addWidget(self._build_fold_panel())
        top.setStretchFactor(0, 3)
        top.setStretchFactor(1, 2)
        top.setSizes([760, 460])
        self.v_split.addWidget(top)
        self.v_split.addWidget(self._build_table_panel())
        self.v_split.setStretchFactor(0, 3)
        self.v_split.setStretchFactor(1, 2)

    def _build_facet_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)

        self.facet_title = QLabel("Select a feature below")
        self.facet_title.setStyleSheet("font-weight: bold; font-size: 13px;")
        lay.addWidget(self.facet_title)

        # organize-by-depth controls
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Split by depth:"))
        self.facet_combos = []
        cat_options = [NONE_LABEL] + self._categories
        for i in range(MAX_FACET_LEVELS):
            combo = QComboBox()
            combo.addItems(cat_options)
            # sensible default: first level = first category, rest none
            if i == 0 and self._categories:
                combo.setCurrentText(self._categories[0])
            combo.currentTextChanged.connect(self._rebuild_facets)
            controls.addWidget(QLabel(f"L{i + 1}"))
            controls.addWidget(combo)
            self.facet_combos.append(combo)

        controls.addSpacing(12)
        controls.addWidget(QLabel("X axis:"))
        self.xaxis_combo = QComboBox()
        self.xaxis_combo.addItems([FILE_LABEL] + self._categories)
        if len(self._categories) > 1:
            self.xaxis_combo.setCurrentText(self._categories[-1])
        self.xaxis_combo.currentTextChanged.connect(self._rebuild_facets)
        controls.addWidget(self.xaxis_combo)

        self.logy_check = QCheckBox("log2 Y")
        self.logy_check.setChecked(True)
        self.logy_check.stateChanged.connect(self._rebuild_facets)
        controls.addWidget(self.logy_check)
        controls.addStretch(1)
        lay.addLayout(controls)

        # replicate designation (the only special role, user-chosen)
        rep_row = QHBoxLayout()
        rep_row.addWidget(QLabel("Replicate column (averaged; others compared):"))
        self.replicate_combo = QComboBox()
        self.replicate_combo.addItems([NONE_LABEL] + self._categories)
        self.replicate_combo.currentTextChanged.connect(self._on_replicate_changed)
        rep_row.addWidget(self.replicate_combo)
        rep_row.addStretch(1)
        lay.addLayout(rep_row)

        self.facet_area = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.facet_area)
        lay.addWidget(holder, 1)
        return panel

    def _build_fold_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)
        title = QLabel("All features — log2 fold change vs abundance")
        title.setStyleSheet("font-weight: bold; font-size: 13px;")
        lay.addWidget(title)

        form = QHBoxLayout()
        form.addWidget(QLabel("Compare"))
        self.compare_combo = QComboBox()
        self.compare_combo.addItems(self._categories)
        self.compare_combo.currentTextChanged.connect(self._on_compare_changed)
        form.addWidget(self.compare_combo)
        form.addWidget(QLabel("A"))
        self.a_combo = QComboBox()
        self.a_combo.currentTextChanged.connect(self._refresh_fold_change)
        form.addWidget(self.a_combo)
        form.addWidget(QLabel("B"))
        self.b_combo = QComboBox()
        self.b_combo.currentTextChanged.connect(self._refresh_fold_change)
        form.addWidget(self.b_combo)
        form.addStretch(1)
        lay.addLayout(form)

        self.fold_plot = pg.PlotWidget()
        self.fold_plot.setLabel("bottom", "log2 fold change (A − B)")
        self.fold_plot.setLabel("left", "mean log2 abundance")
        for name in ("left", "bottom"):
            self.fold_plot.getAxis(name).enableAutoSIPrefix(False)
        lay.addWidget(self.fold_plot, 1)
        self.fold_status = QLabel("")
        self.fold_status.setStyleSheet("color: #8a8f98;")
        lay.addWidget(self.fold_status)

        # populate A/B for the initial compare column
        if self._categories:
            self._on_compare_changed(self.compare_combo.currentText())
        return panel

    def _build_table_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)

        bar = QHBoxLayout()
        self.pep_button = QPushButton("Peptides")
        self.prot_button = QPushButton("Proteins")
        for b in (self.pep_button, self.prot_button):
            b.setCheckable(True)
            b.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.pep_button.setChecked(True)
        self.pep_button.clicked.connect(lambda: self._set_level("peptide"))
        self.prot_button.clicked.connect(lambda: self._set_level("protein"))
        bar.addWidget(self.pep_button)
        bar.addWidget(self.prot_button)
        bar.addSpacing(16)
        self.unique_check = QCheckBox("Unique peptides only")
        self.unique_check.stateChanged.connect(self._on_unique_toggled)
        bar.addWidget(self.unique_check)
        bar.addStretch(1)
        self.count_label = QLabel("")
        self.count_label.setStyleSheet("color: #8a8f98;")
        bar.addWidget(self.count_label)
        lay.addLayout(bar)

        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.verticalHeader().setVisible(False)
        self.table.itemSelectionChanged.connect(self._on_table_selection)
        lay.addWidget(self.table, 1)
        return panel

    # ---- reactive handlers ----------------------------------------------

    def _on_replicate_changed(self, _):
        self._rebuild_facets()
        self._refresh_fold_change()

    def _on_compare_changed(self, column):
        vals = _sorted_values(self.experimental.values(column)) if column else []
        for combo, idx in ((self.a_combo, 0), (self.b_combo, 1)):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(vals)
            if len(vals) > idx:
                combo.setCurrentIndex(idx)
            combo.blockSignals(False)
        self._refresh_fold_change()

    def _on_unique_toggled(self, _):
        self.unique_only = self.unique_check.isChecked()
        self._refresh_table()
        self._refresh_fold_change()

    def _set_level(self, level):
        self.pep_button.setChecked(level == "peptide")
        self.prot_button.setChecked(level == "protein")
        self.unique_check.setEnabled(level == "peptide")
        if level == self.level:
            return
        self.level = level
        self.selected_feature = None
        self._refresh_table()
        self._refresh_fold_change()
        self._auto_select_first()

    # ---- feature table ---------------------------------------------------

    def _refresh_table(self):
        matrix = self.model.matrix(self.level)
        files = [f for f in self.model.filenames()]
        feats = self._visible_features()

        # columns: Feature | Unique | one per file
        headers = ["Feature", "Unique"] + files
        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        # annotate file headers with their category metadata as a tooltip
        for ci, f in enumerate(files):
            item = self.table.horizontalHeaderItem(ci + 2)
            if item is not None:
                row = self._row_for(f)
                meta = ", ".join(f"{c}={row.get(c, '')}" for c in self._categories)
                item.setToolTip(meta)

        self.table.setRowCount(len(feats))
        for r, feat in enumerate(feats):
            self.table.setItem(r, 0, QTableWidgetItem(feat))
            if self.level == "peptide":
                uniq = "yes" if self.model.peptide_is_unique(feat) else "no"
            else:
                uniq = "yes"
            self.table.setItem(r, 1, QTableWidgetItem(uniq))
            per_file = matrix.get(feat, {})
            for ci, f in enumerate(files):
                q = per_file.get(f)
                if q and q > 0:
                    self.table.setItem(r, ci + 2, NumericItem(f"{q:.4g}", float(q)))
                else:
                    self.table.setItem(r, ci + 2, NumericItem("", float("nan")))
        self.table.setSortingEnabled(True)

        label = "peptides" if self.level == "peptide" else "proteins"
        self.count_label.setText(f"{len(feats)} {label} × {len(files)} files")

    def _auto_select_first(self):
        if self.table.rowCount() > 0:
            self.table.selectRow(0)

    def _on_table_selection(self):
        items = self.table.selectedItems()
        if not items:
            return
        feat_item = self.table.item(items[0].row(), 0)
        if feat_item is None:
            return
        self.selected_feature = feat_item.text()
        self.feature_selected.emit(self.selected_feature)
        self._rebuild_facets()

    def select_feature(self, feature):
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and item.text() == feature:
                self.table.selectRow(r)
                self.table.scrollToItem(item)
                return

    # ---- fold-change scatter (plain log difference) ---------------------

    def _group_mean_log2(self, per_file, files):
        """Mean log2 quantity over ``files`` for one feature.

        When a replicate column is set, replicates (files differing only in the
        replicate column) are first averaged in linear space, then logged — so
        replicate runs count once, not N times. Otherwise every file counts."""
        rep = self._replicate_column()
        vals = []
        if rep:
            groups = {}
            for f in files:
                q = per_file.get(f)
                if q and q > 0:
                    row = self._row_for(f)
                    key = tuple((c, row.get(c, "")) for c in self._categories if c != rep)
                    groups.setdefault(key, []).append(q)
            for qs in groups.values():
                vals.append(math.log2(sum(qs) / len(qs)))
        else:
            for f in files:
                q = per_file.get(f)
                if q and q > 0:
                    vals.append(math.log2(q))
        if not vals:
            return None
        return sum(vals) / len(vals)

    def _refresh_fold_change(self):
        if not self._active:
            return
        self.fold_plot.clear()
        pal = palette(self.theme)
        col = self.compare_combo.currentText()
        a_val = self.a_combo.currentText()
        b_val = self.b_combo.currentText()
        if not col or not a_val or not b_val or a_val == b_val:
            self.fold_status.setText("Pick a compare column with two distinct values.")
            return

        files_a = set(self.experimental.filenames_for(**{col: a_val}))
        files_b = set(self.experimental.filenames_for(**{col: b_val}))
        matrix = self.model.matrix(self.level)

        xs, ys, feats = [], [], []
        for feat in self._visible_features():
            per_file = matrix.get(feat, {})
            ma = self._group_mean_log2(per_file, files_a)
            mb = self._group_mean_log2(per_file, files_b)
            if ma is None or mb is None:
                continue
            xs.append(ma - mb)
            ys.append((ma + mb) / 2.0)
            feats.append(feat)

        if not xs:
            self.fold_status.setText(
                f"No features quantified in both {a_val} and {b_val}.")
            return

        point = pg.mkBrush(pal["points"][0], pal["points"][1], pal["points"][2], 180)
        outline = pg.mkPen(pal["fg"], width=0.6)
        scatter = pg.ScatterPlotItem(x=xs, y=ys, brush=point, pen=outline, size=6)
        scatter.feats = feats
        scatter.sigClicked.connect(self._on_fold_click)
        self.fold_plot.addItem(scatter)
        self._fold_scatter = scatter

        guide = pg.mkPen(pal["fg"], width=1, style=Qt.DashLine)
        self.fold_plot.addLine(x=0.0, pen=guide)
        self.fold_plot.setLabel("bottom", f"log2 fold change ({a_val} − {b_val})")
        self.fold_status.setText(
            f"{len(xs)} features · {a_val} (n={len(files_a)}) vs "
            f"{b_val} (n={len(files_b)})")

    def _on_fold_click(self, scatter, points):
        if not points:
            return
        idx = points[0].index()
        feats = getattr(scatter, "feats", [])
        if 0 <= idx < len(feats):
            self.select_feature(feats[idx])

    # ---- faceted feature view -------------------------------------------

    def _clear_facets(self):
        while self.facet_area.count():
            item = self.facet_area.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

    def _rebuild_facets(self, *args):
        if not self._active:
            return
        self._clear_facets()
        feat = self.selected_feature
        if not feat:
            self.facet_title.setText("Select a feature below")
            return
        self.facet_title.setText(feat)

        levels = []
        for combo in self.facet_combos:
            c = combo.currentText()
            if c != NONE_LABEL and c not in levels:
                levels.append(c)

        matrix = self.model.matrix(self.level)
        per_file = matrix.get(feat, {})
        files = [f for f in self.model.filenames() if per_file.get(f)]
        widget = self._facet_node(feat, per_file, files, levels)
        self.facet_area.addWidget(widget)

    def _facet_node(self, feat, per_file, files, levels):
        if not levels or not files:
            return self._leaf_plot(feat, per_file, files)
        col = levels[0]
        values = _sorted_values({self._row_for(f).get(col, "") for f in files})
        split = QSplitter(Qt.Horizontal)
        pal = THEMES["dark" if self.theme != "light" else "light"]
        for v in values:
            sub = [f for f in files if self._row_for(f).get(col, "") == v]
            wrap = QWidget()
            wlay = QVBoxLayout(wrap)
            wlay.setContentsMargins(2, 2, 2, 2)
            head = QLabel(f"{col} = {v}")
            head.setAlignment(Qt.AlignCenter)
            head.setStyleSheet(
                f"color: {pal['fg']}; font-weight: bold; "
                f"border-bottom: 1px solid {pal['line']};")
            wlay.addWidget(head)
            wlay.addWidget(self._facet_node(feat, per_file, sub, levels[1:]), 1)
            split.addWidget(wrap)
        return split

    def _leaf_plot(self, feat, per_file, files):
        plot = pg.PlotWidget()
        style_plot(plot, palette(self.theme))
        for name in ("left", "bottom"):
            plot.getAxis(name).enableAutoSIPrefix(False)
        logy = self.logy_check.isChecked()
        plot.setLabel("left", "log2 quantity" if logy else "quantity")
        xcol = self.xaxis_combo.currentText()
        pal = palette(self.theme)
        rep = self._replicate_column()

        def yval(q):
            return math.log2(q) if logy else q

        if xcol == FILE_LABEL:
            xcol = None

        # bucket files by x value
        buckets = {}
        for f in files:
            q = per_file.get(f)
            if not q or q <= 0:
                continue
            xv = self._row_for(f).get(xcol, "") if xcol else f
            buckets.setdefault(xv, []).append(q)

        xvals = _sorted_values(list(buckets.keys()))
        xindex = {v: i for i, v in enumerate(xvals)}
        point = (pal["points"][0], pal["points"][1], pal["points"][2])

        sx, sy, mx, my = [], [], [], []
        for xv in xvals:
            qs = buckets[xv]
            i = xindex[xv]
            for q in qs:
                sx.append(i)
                sy.append(yval(q))
            # mean marker (linear-mean then transform, so replicates average)
            mean_q = sum(qs) / len(qs)
            mx.append(i)
            my.append(yval(mean_q))

        plot.plot(sx, sy, pen=None, symbol="o", symbolSize=7,
                  symbolBrush=pg.mkBrush(*point, 220),
                  symbolPen=pg.mkPen(pal["fg"], width=0.5))
        if rep and len(mx) > 1:
            plot.plot(mx, my, pen=pg.mkPen(pal["accent"], width=2))
        plot.getAxis("bottom").setTicks([list(enumerate(xvals))])
        plot.setLabel("bottom", xcol if xcol else "file")
        return plot

    # ---- theming ---------------------------------------------------------

    def apply_theme(self, theme):
        self.theme = theme
        if not self._active:
            return
        pal = palette(theme)
        if getattr(self, "fold_plot", None) is not None:
            style_plot(self.fold_plot, pal)
        self._refresh_fold_change()
        self._rebuild_facets()
