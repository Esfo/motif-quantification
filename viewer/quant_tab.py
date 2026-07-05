"""Tab 3 — Quantitative Comparisons.

Compare peptide / protein quantities across files, grouped **entirely** by the
project ``experimental-setup`` file. Nothing about the design is hard-coded: the
tab reads whatever columns the setup file has and treats every one of them as a
generic category. The only special designation is optional and made *by the
user*: marking one column as the **replicate** column (its runs may be averaged;
every other column is compared, never averaged).

Layout
------
Top half (horizontal split):
  * **left** — a *faceted* view of the selected feature's quantities. Above the
    plot is an **organizer** (a small growable pseudo-table): each row picks a
    category and how it organizes the view — *split into columns*, *split into
    rows*, or *use as the x-axis*. Add as many layers as you like, in any order,
    to nest the panels by depth (e.g. split by condition → then a time-series
    x-axis inside each). There is no fixed number of levels.
  * **right** — a scatter of **every** feature: x = log2 fold change between two
    chosen category values (a plain log difference, no statistical test),
    y = mean log2 abundance. Click a point to select that feature.

Bottom half:
  * a **Peptides ⇄ Proteins** switch and a **unique-only** filter, over a
    long-form table whose columns are **the experimental-setup columns**
    (filename + every design category) plus the feature id, a unique flag, and
    the quantity — one row per (feature, file), so you can see and sort every
    quantity across conditions and files.

Filenames from the search tables (which may carry a ``.centroid.mzML`` suffix)
are matched to the design's ``filename`` column by stripped stem, so the join is
robust to extension differences.

Everything is reactive — no run button.
"""

import json
import math

from PySide6.QtCore import Qt, QSettings, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
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
    from .session import strip_ms_suffix
except ImportError:
    from quant_model import QuantModel
    from theming import palette, style_plot
    from session import strip_ms_suffix


FILE_LABEL = "(file)"
MODE_COLS = "Split → columns"
MODE_ROWS = "Split ↓ rows"
MODE_XAXIS = "X-axis"
MODES = [MODE_COLS, MODE_ROWS, MODE_XAXIS]

THEMES = {
    "dark": {"fg": "#e6e6e6", "bg": "#101216", "muted": "#8a8f98",
             "panel": "#16181d", "line": "#2a2d33"},
    "light": {"fg": "#202020", "bg": "#fafafa", "muted": "#606060",
              "panel": "#ffffff", "line": "#c8ccd2"},
}


def _norm(name):
    """Filename → comparable stem (drop mzML/centroid/raw suffixes, lowercase)."""
    return strip_ms_suffix(name).lower()


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
        self.settings = QSettings("motif-quantification", "viewer")
        self._saved = self._load_state()
        self._restoring = True

        self.level = self._saved.get("level", "peptide")
        self.unique_only = bool(self._saved.get("unique", False))
        self.selected_feature = None
        self._active = not experimental.is_empty()

        self._categories = [c for c in (experimental.columns() if self._active else [])
                            if c != "filename"]
        # robust filename join: model filename -> experimental row
        self._exp_by_norm = {}
        if self._active:
            for row in experimental.rows:
                self._exp_by_norm[_norm(row.get("filename", ""))] = row

        self._layers = []  # list of {"cat": str, "mode": str} organizer rows

        self._build_ui()
        if self._active:
            self._refresh_table()
            self._refresh_fold_change()
            self._auto_select_first()
        self._restoring = False
        self.apply_theme(theme)

    # ---- persisted state (organizer + contrast, like the panel layouts) --

    def _load_state(self):
        raw = self.settings.value("quant_state")
        if not raw:
            return {}
        try:
            return json.loads(raw)
        except Exception:
            return {}

    def _save_state(self):
        if getattr(self, "_restoring", False) or not self._active:
            return
        state = {
            "level": self.level,
            "unique": self.unique_only,
            "logy": self.logy_check.isChecked(),
            "replicate": self.replicate_combo.currentText(),
            "layers": self._layers,
            "compare": self.compare_combo.currentText(),
            "a": self.a_combo.currentText(),
            "b": self.b_combo.currentText(),
        }
        self.settings.setValue("quant_state", json.dumps(state))

    # ---- design helpers --------------------------------------------------

    def _row_for(self, filename):
        return self._exp_by_norm.get(_norm(filename), {})

    def _cat_value(self, filename, column):
        return self._row_for(filename).get(column, "")

    def _replicate_column(self):
        col = self.replicate_combo.currentText()
        return col if col in self._categories else None

    def _files_in(self, column, value):
        """Model filenames whose design row has column == value."""
        return [f for f in self.model.filenames() if self._cat_value(f, column) == value]

    def _visible_features(self):
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

        # organizer pseudo-table: growable list of layer rows
        org_header = QHBoxLayout()
        org_header.addWidget(QLabel("Organize by (top → bottom = outer → inner):"))
        org_header.addStretch(1)
        self.logy_check = QCheckBox("log2 Y")
        self.logy_check.setChecked(bool(self._saved.get("logy", True)))
        self.logy_check.stateChanged.connect(self._on_logy_changed)
        org_header.addWidget(self.logy_check)
        rep_lbl = QLabel("Replicate:")
        org_header.addWidget(rep_lbl)
        self.replicate_combo = QComboBox()
        self.replicate_combo.addItems(["(none)"] + self._categories)
        saved_rep = self._saved.get("replicate")
        if saved_rep and saved_rep in self._categories:
            self.replicate_combo.setCurrentText(saved_rep)
        self.replicate_combo.currentTextChanged.connect(self._on_replicate_changed)
        org_header.addWidget(self.replicate_combo)
        lay.addLayout(org_header)

        self.layer_area = QVBoxLayout()
        self.layer_area.setSpacing(2)
        layer_holder = QWidget()
        layer_holder.setLayout(self.layer_area)
        lay.addWidget(layer_holder)

        add_row = QHBoxLayout()
        self.add_layer_btn = QPushButton("+ Add layer")
        self.add_layer_btn.clicked.connect(lambda: self._add_layer())
        add_row.addWidget(self.add_layer_btn)
        add_row.addStretch(1)
        lay.addLayout(add_row)

        self.facet_area = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.facet_area)
        lay.addWidget(holder, 1)

        # No auto-fill — the organizer starts however the user last left it
        # (persisted), otherwise empty. It is the user's choice what to add.
        for layer in self._saved.get("layers", []):
            cat = layer.get("cat")
            mode = layer.get("mode")
            if cat in self._categories and mode in MODES:
                self._layers.append({"cat": cat, "mode": mode})
        self._rebuild_layer_rows()
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
        saved_compare = self._saved.get("compare")
        if saved_compare and saved_compare in self._categories:
            self.compare_combo.setCurrentText(saved_compare)
        self.compare_combo.currentTextChanged.connect(self._on_compare_changed)
        form.addWidget(self.compare_combo)
        form.addWidget(QLabel("A"))
        self.a_combo = QComboBox()
        self.a_combo.currentTextChanged.connect(lambda _: self._on_ab_changed())
        form.addWidget(self.a_combo)
        form.addWidget(QLabel("B"))
        self.b_combo = QComboBox()
        self.b_combo.currentTextChanged.connect(lambda _: self._on_ab_changed())
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

        if self._categories:
            self._on_compare_changed(self.compare_combo.currentText())
            # restore the saved A/B values now that the combos are populated
            sa, sb = self._saved.get("a"), self._saved.get("b")
            if sa and self.a_combo.findText(sa) >= 0:
                self.a_combo.setCurrentText(sa)
            if sb and self.b_combo.findText(sb) >= 0:
                self.b_combo.setCurrentText(sb)
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
        self.pep_button.setChecked(self.level == "peptide")
        self.prot_button.setChecked(self.level == "protein")
        self.pep_button.clicked.connect(lambda: self._set_level("peptide"))
        self.prot_button.clicked.connect(lambda: self._set_level("protein"))
        bar.addWidget(self.pep_button)
        bar.addWidget(self.prot_button)
        bar.addSpacing(16)
        self.unique_check = QCheckBox("Unique peptides only")
        self.unique_check.setChecked(self.unique_only)
        self.unique_check.setEnabled(self.level == "peptide")
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

    # ---- organizer (dynamic layers) -------------------------------------

    def _add_layer(self, cat=None, mode=None, rebuild=True):
        if not self._categories:
            return
        self._layers.append({
            "cat": cat or self._categories[0],
            "mode": mode or MODE_COLS,
        })
        if rebuild:
            self._rebuild_layer_rows()
            self._rebuild_facets()
            self._save_state()

    def _remove_layer(self, index):
        if 0 <= index < len(self._layers):
            self._layers.pop(index)
            self._rebuild_layer_rows()
            self._rebuild_facets()
            self._save_state()

    def _rebuild_layer_rows(self):
        # clear existing widgets
        while self.layer_area.count():
            item = self.layer_area.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

        for i, layer in enumerate(self._layers):
            row = QWidget()
            rl = QHBoxLayout(row)
            rl.setContentsMargins(0, 0, 0, 0)
            rl.addWidget(QLabel(f"{i + 1}."))
            cat = QComboBox()
            cat.addItems(self._categories)
            cat.setCurrentText(layer["cat"])
            cat.currentTextChanged.connect(lambda v, idx=i: self._set_layer(idx, "cat", v))
            rl.addWidget(cat, 1)
            mode = QComboBox()
            mode.addItems(MODES)
            mode.setCurrentText(layer["mode"])
            mode.currentTextChanged.connect(lambda v, idx=i: self._set_layer(idx, "mode", v))
            rl.addWidget(mode)
            rm = QPushButton("✕")
            rm.setFixedWidth(26)
            rm.clicked.connect(lambda _=False, idx=i: self._remove_layer(idx))
            rl.addWidget(rm)
            self.layer_area.addWidget(row)

    def _set_layer(self, index, key, value):
        if 0 <= index < len(self._layers):
            self._layers[index][key] = value
            self._rebuild_facets()
            self._save_state()

    # ---- reactive handlers ----------------------------------------------

    def _on_logy_changed(self, _):
        self._rebuild_facets()
        self._save_state()

    def _on_replicate_changed(self, _):
        self._rebuild_facets()
        self._refresh_fold_change()
        self._save_state()

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
        self._save_state()

    def _on_ab_changed(self):
        self._refresh_fold_change()
        self._save_state()

    def _on_unique_toggled(self, _):
        self.unique_only = self.unique_check.isChecked()
        self._refresh_table()
        self._refresh_fold_change()
        self._save_state()

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
        self._save_state()

    # ---- feature table (long form; columns = experimental-setup cols) ----

    def _refresh_table(self):
        matrix = self.model.matrix(self.level)
        feats = self._visible_features()

        feat_label = "peptide" if self.level == "peptide" else "protein"
        # Columns are the experimental-setup categories only — NOT the filename
        # (each file is already fully described by its category values).
        cat_cols = [c for c in self.experimental.columns() if c != "filename"]
        headers = [feat_label] + cat_cols + ["unique", "quantity"]

        self.table.setSortingEnabled(False)
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        self.table.setRowCount(0)

        rows = []
        for feat in feats:
            uniq = ("yes" if (self.level == "protein"
                              or self.model.peptide_is_unique(feat)) else "no")
            for fname, q in matrix.get(feat, {}).items():
                if not q or q <= 0:
                    continue
                rows.append((feat, fname, uniq, float(q)))

        self.table.setRowCount(len(rows))
        for r, (feat, fname, uniq, q) in enumerate(rows):
            design = self._row_for(fname)
            self.table.setItem(r, 0, QTableWidgetItem(feat))
            for ci, col in enumerate(cat_cols):
                self.table.setItem(r, ci + 1, QTableWidgetItem(design.get(col, "")))
            self.table.setItem(r, 1 + len(cat_cols), QTableWidgetItem(uniq))
            self.table.setItem(r, 2 + len(cat_cols), NumericItem(f"{q:.4g}", q))

        self.table.setSortingEnabled(True)
        label = "peptides" if self.level == "peptide" else "proteins"
        self.count_label.setText(f"{len(feats)} {label} · {len(rows)} rows")

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
        """Mean log2 quantity over ``files`` for one feature. With a replicate
        column set, replicates are averaged in linear space first."""
        rep = self._replicate_column()
        vals = []
        if rep:
            groups = {}
            for f in files:
                q = per_file.get(f)
                if q and q > 0:
                    key = tuple((c, self._cat_value(f, c))
                                for c in self._categories if c != rep)
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

        files_a = self._files_in(col, a_val)
        files_b = self._files_in(col, b_val)
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
                f"No features quantified in both {a_val} (n={len(files_a)}) and "
                f"{b_val} (n={len(files_b)}).")
            return

        pts = pal["points"]
        scatter = pg.ScatterPlotItem(
            x=xs, y=ys, size=6,
            brush=pg.mkBrush(pts[0], pts[1], pts[2], 180),
            pen=pg.mkPen(pal["fg"], width=0.6))
        scatter.feats = feats
        scatter.sigClicked.connect(self._on_fold_click)
        self.fold_plot.addItem(scatter)
        self.fold_plot.addLine(x=0.0, pen=pg.mkPen(pal["fg"], width=1, style=Qt.DashLine))
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

        splits = [l for l in self._layers if l["mode"] in (MODE_COLS, MODE_ROWS)]
        xaxis = next((l["cat"] for l in self._layers if l["mode"] == MODE_XAXIS), None)

        per_file = self.model.matrix(self.level).get(feat, {})
        files = [f for f in self.model.filenames() if per_file.get(f)]
        widget = self._facet_node(feat, per_file, files, splits, 0, xaxis)
        self.facet_area.addWidget(widget)

    def _facet_node(self, feat, per_file, files, splits, depth, xaxis):
        if depth >= len(splits) or not files:
            return self._leaf_plot(feat, per_file, files, xaxis)
        layer = splits[depth]
        col = layer["cat"]
        orient = Qt.Horizontal if layer["mode"] == MODE_COLS else Qt.Vertical
        values = _sorted_values({self._cat_value(f, col) for f in files})
        pal = THEMES["dark" if self.theme != "light" else "light"]
        split = QSplitter(orient)
        for v in values:
            sub = [f for f in files if self._cat_value(f, col) == v]
            wrap = QWidget()
            wlay = QVBoxLayout(wrap)
            wlay.setContentsMargins(2, 2, 2, 2)
            head = QLabel(f"{col} = {v}" if v != "" else f"{col} = (blank)")
            head.setAlignment(Qt.AlignCenter)
            head.setStyleSheet(
                f"color: {pal['fg']}; font-weight: bold; "
                f"border-bottom: 1px solid {pal['line']};")
            wlay.addWidget(head)
            wlay.addWidget(self._facet_node(feat, per_file, sub, splits, depth + 1, xaxis), 1)
            split.addWidget(wrap)
        return split

    def _leaf_plot(self, feat, per_file, files, xaxis):
        plot = pg.PlotWidget()
        style_plot(plot, palette(self.theme))
        for name in ("left", "bottom"):
            plot.getAxis(name).enableAutoSIPrefix(False)
        logy = self.logy_check.isChecked()
        plot.setLabel("left", "log2 quantity" if logy else "quantity")
        pal = palette(self.theme)
        rep = self._replicate_column()
        xcol = xaxis if xaxis in self._categories else None

        def yval(q):
            return math.log2(q) if logy else q

        buckets = {}
        for f in files:
            q = per_file.get(f)
            if not q or q <= 0:
                continue
            xv = self._cat_value(f, xcol) if xcol else f
            buckets.setdefault(xv, []).append(q)

        xvals = _sorted_values(list(buckets.keys()))
        xindex = {v: i for i, v in enumerate(xvals)}
        pts = pal["points"]

        sx, sy, mx, my = [], [], [], []
        for xv in xvals:
            qs = buckets[xv]
            i = xindex[xv]
            for q in qs:
                sx.append(i)
                sy.append(yval(q))
            mean_q = sum(qs) / len(qs)
            mx.append(i)
            my.append(yval(mean_q))

        plot.plot(sx, sy, pen=None, symbol="o", symbolSize=7,
                  symbolBrush=pg.mkBrush(pts[0], pts[1], pts[2], 220),
                  symbolPen=pg.mkPen(pal["fg"], width=0.5))
        if rep and len(mx) > 1:
            # mean-across-replicates line, drawn in the theme fg (white on dark)
            plot.plot(mx, my, pen=pg.mkPen(pal["fg"], width=2))
        # short tick labels (avoid long filenames overflowing)
        ticks = [(i, (v if len(str(v)) <= 14 else str(v)[:12] + "…"))
                 for i, v in enumerate(xvals)]
        plot.getAxis("bottom").setTicks([ticks])
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
