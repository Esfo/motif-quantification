"""Tab 3 — Quantitative Comparisons.

Differential-expression workbench over the reorganized per-file quant tables,
driven by the project ``experimental-setup`` design file.

Layout
------
Top half (horizontal split):
  * **Design** (left) — every experimental-setup column gets a *role*:
        Group        — a categorical factor whose values define the samples
                        being compared (the DE contrast is between two of its
                        values);
        Series-axis  — an ordered factor (a titration / calibration / time
                        course); the visualization plots quantity *along* it;
        Replicate    — repeats aggregated over (the spread within a group);
        Pair         — matched-sample id for the paired test;
        Ignore       — not used.
    Below the roles: the contrast (which Group column, value A vs value B),
    optional "restrict to" filters on any other Group columns, the test
    (Welch / paired), a min-replicates spin box, and **Run DE**.
  * **Visualization** (center) — for the currently-selected feature: a titration
    line plot when a Series-axis is assigned (one line per group, x = the ordered
    series values, points = replicates), otherwise a per-group strip/mean plot.
  * **Differential expression** (right) — a volcano plot (log2 fold change vs
    -log10 FDR) of the last DE run; click a point to select that feature.

Bottom half:
  * a **Peptides ⇄ Proteins** switch (+ protein roll-up selector) and a sortable
    table of every feature with its DE stats. Selecting a row drives the
    visualization above; switching Peptides/Proteins re-analyzes everything.

The flexible column-role model is the point: nothing here is hard-coded to
condition/fraction/replicate. Assign *any* column as the Group being contrasted
and *any* column as the ordered Series-axis — e.g. "fraction is the group,
replicate is the series-axis" groups by fraction then orders quantities along
replicate, exactly the titration/calibration framing requested.
"""

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .quant_model import QuantModel
    from .de_stats import differential_expression, log2_values
    from .theming import palette, style_plot
except ImportError:
    from quant_model import QuantModel
    from de_stats import differential_expression, log2_values
    from theming import palette, style_plot


ROLES = ["Ignore", "Group", "Series-axis", "Replicate", "Pair"]

# A stable, colour-blind-friendly categorical palette for group lines/points.
GROUP_COLORS = [
    (76, 114, 176), (221, 132, 82), (85, 168, 104), (196, 78, 82),
    (129, 114, 179), (147, 120, 96), (218, 139, 195), (140, 140, 140),
    (204, 185, 116), (100, 181, 205),
]


def _sorted_series_values(values):
    """Order series-axis values numerically when they all parse as numbers,
    else lexically. Returns a new sorted list."""
    def as_num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None
    nums = [as_num(v) for v in values]
    if all(n is not None for n in nums):
        return [v for _, v in sorted(zip(nums, values))]
    return sorted(values)


class NumericItem(QTableWidgetItem):
    """Table item that sorts by an underlying float (NaN/None sort last)."""

    def __init__(self, text, value):
        super().__init__(text)
        self.setData(Qt.UserRole, value)

    def __lt__(self, other):
        a = self.data(Qt.UserRole)
        b = other.data(Qt.UserRole) if isinstance(other, NumericItem) else None
        a_bad = a is None or (isinstance(a, float) and math.isnan(a))
        b_bad = b is None or (isinstance(b, float) and math.isnan(b))
        if a_bad and b_bad:
            return False
        if a_bad:
            return False
        if b_bad:
            return True
        return a < b


class QuantTab(QWidget):
    """Quantitative Comparisons / differential-expression tab."""

    feature_selected = Signal(str)

    def __init__(self, session, experimental, theme="dark", parent=None):
        super().__init__(parent)
        self.session = session
        self.experimental = experimental
        self.theme = theme

        self.model = QuantModel(session)
        self.level = "peptide"      # or "protein"
        self.rollup = "sum"
        self.roles = {}             # column -> role
        self.selected_feature = None
        self.de_records = []
        self._de_by_feature = {}
        self._restrict_combos = {}  # column -> QComboBox
        self._legend = None
        self._active = not experimental.is_empty()

        self._seed_default_roles()
        self._build_ui()
        if self._active:
            self._refresh_features()
        self.apply_theme(theme)

    # ---- role defaults ---------------------------------------------------

    def _seed_default_roles(self):
        """Sensible defaults from the canonical column names, if present."""
        cols = self.experimental.columns() if not self.experimental.is_empty() else []
        for col in cols:
            if col == "filename":
                continue
            low = col.lower()
            if low == "condition":
                self.roles[col] = "Group"
            elif low == "fraction":
                self.roles[col] = "Group"
            elif low == "replicate":
                self.roles[col] = "Replicate"
            elif low == "pair_id":
                self.roles[col] = "Pair"
            else:
                self.roles[col] = "Ignore"

    def _group_columns(self):
        return [c for c, r in self.roles.items() if r == "Group"]

    def _series_column(self):
        for c, r in self.roles.items():
            if r == "Series-axis":
                return c
        return None

    def _pair_column(self):
        for c, r in self.roles.items():
            if r == "Pair":
                return c
        return None

    # ---- UI construction -------------------------------------------------

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        if self.experimental.is_empty():
            warn = QLabel(
                "No experimental-setup file found for this project.\n\n"
                "Quantitative Comparisons reads the project's experimental-setup "
                "csv (filename,condition,fraction,replicate,pair_id,…) to group "
                "runs. Add one beside distributions/ and searches/ and reload.")
            warn.setWordWrap(True)
            warn.setAlignment(Qt.AlignCenter)
            outer.addWidget(warn)
            return

        self.v_split = QSplitter(Qt.Vertical)
        outer.addWidget(self.v_split, 1)

        # --- top half ---
        top = QSplitter(Qt.Horizontal)
        top.addWidget(self._build_design_panel())
        top.addWidget(self._build_viz_panel())
        top.addWidget(self._build_de_panel())
        top.setStretchFactor(0, 0)
        top.setStretchFactor(1, 1)
        top.setStretchFactor(2, 1)
        top.setSizes([320, 520, 520])
        self.h_split = top
        self.v_split.addWidget(top)

        # --- bottom half ---
        self.v_split.addWidget(self._build_table_panel())
        self.v_split.setStretchFactor(0, 3)
        self.v_split.setStretchFactor(1, 2)

    def _build_design_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)

        title = QLabel("Design")
        title.setStyleSheet("font-weight: bold; font-size: 14px;")
        lay.addWidget(title)

        hint = QLabel("Assign a role to each design column:")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #8a8f98;")
        lay.addWidget(hint)

        # role rows
        self.role_combos = {}
        roles_form = QFormLayout()
        roles_form.setLabelAlignment(Qt.AlignRight)
        for col in self.experimental.columns():
            if col == "filename":
                continue
            combo = QComboBox()
            combo.addItems(ROLES)
            combo.setCurrentText(self.roles.get(col, "Ignore"))
            combo.currentTextChanged.connect(
                lambda role, c=col: self._on_role_changed(c, role))
            self.role_combos[col] = combo
            n_vals = len(self.experimental.values(col))
            roles_form.addRow(f"{col} ({n_vals})", combo)
        lay.addLayout(roles_form)

        divider = QFrame()
        divider.setFrameShape(QFrame.HLine)
        lay.addWidget(divider)

        # contrast
        contrast_title = QLabel("Contrast (A vs B)")
        contrast_title.setStyleSheet("font-weight: bold;")
        lay.addWidget(contrast_title)

        cform = QFormLayout()
        self.contrast_combo = QComboBox()
        self.contrast_combo.currentTextChanged.connect(self._on_contrast_changed)
        cform.addRow("Group column", self.contrast_combo)
        self.group_a_combo = QComboBox()
        self.group_b_combo = QComboBox()
        cform.addRow("Group A", self.group_a_combo)
        cform.addRow("Group B", self.group_b_combo)
        lay.addLayout(cform)

        # restrict-to filters (rebuilt dynamically)
        self.restrict_label = QLabel("Restrict to")
        self.restrict_label.setStyleSheet("font-weight: bold;")
        lay.addWidget(self.restrict_label)
        self.restrict_container = QVBoxLayout()
        lay.addLayout(self.restrict_container)

        # test + params
        pform = QFormLayout()
        self.test_combo = QComboBox()
        self.test_combo.addItems(["Welch's t-test", "Paired t-test"])
        pform.addRow("Test", self.test_combo)
        self.min_rep_spin = QSpinBox()
        self.min_rep_spin.setRange(2, 100)
        self.min_rep_spin.setValue(2)
        pform.addRow("Min replicates", self.min_rep_spin)
        lay.addLayout(pform)

        self.run_button = QPushButton("Run DE")
        self.run_button.clicked.connect(self.run_de)
        lay.addWidget(self.run_button)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color: #8a8f98;")
        lay.addWidget(self.status_label)

        lay.addStretch(1)

        self._rebuild_contrast_options()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        scroll.setMinimumWidth(280)
        return scroll

    def _build_viz_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)
        self.viz_title = QLabel("Select a feature below")
        self.viz_title.setStyleSheet("font-weight: bold; font-size: 13px;")
        lay.addWidget(self.viz_title)
        self.viz_plot = pg.PlotWidget()
        self.viz_plot.setLabel("left", "log2 quantity")
        lay.addWidget(self.viz_plot, 1)
        return panel

    def _build_de_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(6, 6, 6, 6)
        title = QLabel("Differential expression — volcano")
        title.setStyleSheet("font-weight: bold; font-size: 13px;")
        lay.addWidget(title)
        self.volcano_plot = pg.PlotWidget()
        self.volcano_plot.setLabel("bottom", "log2 fold change (A − B)")
        self.volcano_plot.setLabel("left", "-log10 FDR")
        lay.addWidget(self.volcano_plot, 1)
        self.volcano_scatter = None
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
        bar.addWidget(QLabel("Protein roll-up:"))
        self.rollup_combo = QComboBox()
        self.rollup_combo.addItems(["sum", "median", "unique"])
        self.rollup_combo.currentTextChanged.connect(self._on_rollup_changed)
        self.rollup_combo.setEnabled(False)
        bar.addWidget(self.rollup_combo)

        bar.addStretch(1)
        self.count_label = QLabel("")
        self.count_label.setStyleSheet("color: #8a8f98;")
        bar.addWidget(self.count_label)
        lay.addLayout(bar)

        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels(
            ["Feature", "n files", "mean A", "mean B", "log2FC", "p", "FDR"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.verticalHeader().setVisible(False)
        self.table.itemSelectionChanged.connect(self._on_table_selection)
        lay.addWidget(self.table, 1)
        return panel

    # ---- role / contrast wiring -----------------------------------------

    def _on_role_changed(self, column, role):
        # a column may hold only one Series-axis / Pair role at a time
        if role in ("Series-axis", "Pair"):
            for other, combo in self.role_combos.items():
                if other != column and combo.currentText() == role:
                    combo.blockSignals(True)
                    combo.setCurrentText("Ignore")
                    combo.blockSignals(False)
                    self.roles[other] = "Ignore"
        self.roles[column] = role
        self._rebuild_contrast_options()
        self._update_viz()

    def _rebuild_contrast_options(self):
        groups = self._group_columns()
        cur = self.contrast_combo.currentText()
        self.contrast_combo.blockSignals(True)
        self.contrast_combo.clear()
        self.contrast_combo.addItems(groups)
        if cur in groups:
            self.contrast_combo.setCurrentText(cur)
        self.contrast_combo.blockSignals(False)
        self._on_contrast_changed(self.contrast_combo.currentText())

    def _on_contrast_changed(self, column):
        vals = self.experimental.values(column) if column else []
        for combo, default_idx in ((self.group_a_combo, 0), (self.group_b_combo, 1)):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(vals)
            if len(vals) > default_idx:
                combo.setCurrentIndex(default_idx)
            combo.blockSignals(False)
        self._rebuild_restrict_filters(exclude=column)

    def _rebuild_restrict_filters(self, exclude):
        # clear existing
        while self.restrict_container.count():
            item = self.restrict_container.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._restrict_combos = {}

        others = [c for c in self._group_columns() if c != exclude]
        self.restrict_label.setVisible(bool(others))
        for col in others:
            row = QWidget()
            rlay = QHBoxLayout(row)
            rlay.setContentsMargins(0, 0, 0, 0)
            rlay.addWidget(QLabel(col))
            combo = QComboBox()
            combo.addItem("(any)")
            combo.addItems(self.experimental.values(col))
            rlay.addWidget(combo, 1)
            self._restrict_combos[col] = combo
            self.restrict_container.addWidget(row)

    # ---- level / rollup switching ---------------------------------------

    def _set_level(self, level):
        if level == self.level:
            # keep the pressed button checked
            self.pep_button.setChecked(self.level == "peptide")
            self.prot_button.setChecked(self.level == "protein")
            return
        self.level = level
        self.pep_button.setChecked(level == "peptide")
        self.prot_button.setChecked(level == "protein")
        self.rollup_combo.setEnabled(level == "protein")
        self.selected_feature = None
        self.de_records = []
        self._de_by_feature = {}
        self._refresh_features()
        self._draw_volcano()
        self._update_viz()

    def _on_rollup_changed(self, rollup):
        self.rollup = rollup
        if self.level == "protein":
            self.de_records = []
            self._de_by_feature = {}
            self._refresh_features()
            self._draw_volcano()
            self._update_viz()

    # ---- feature table ---------------------------------------------------

    def _current_matrix(self):
        return self.model.matrix(self.level, self.rollup)

    def _refresh_features(self):
        """Populate the table from the current matrix; DE columns blank until a
        run has been made (or filled from cached DE records)."""
        matrix = self._current_matrix()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)

        features = sorted(matrix.keys())
        self.table.setRowCount(len(features))
        for r, feat in enumerate(features):
            per_file = matrix[feat]
            n_files = sum(1 for q in per_file.values() if q and q > 0)
            rec = self._de_by_feature.get(feat)

            self.table.setItem(r, 0, QTableWidgetItem(feat))
            self.table.setItem(r, 1, NumericItem(str(n_files), float(n_files)))
            if rec:
                self.table.setItem(r, 2, NumericItem(f"{rec['mean_a']:.2f}", rec["mean_a"]))
                self.table.setItem(r, 3, NumericItem(f"{rec['mean_b']:.2f}", rec["mean_b"]))
                self.table.setItem(r, 4, NumericItem(f"{rec['log2fc']:+.2f}", rec["log2fc"]))
                p = rec["p"]
                fdr = rec["fdr"]
                self.table.setItem(r, 5, NumericItem(_fmt_p(p), p))
                self.table.setItem(r, 6, NumericItem(_fmt_p(fdr), fdr))
            else:
                for c in range(2, 7):
                    self.table.setItem(r, c, NumericItem("", float("nan")))

        self.table.setSortingEnabled(True)
        self.rollup_combo.setEnabled(self.level == "protein")
        label = "peptides" if self.level == "peptide" else "proteins"
        self.count_label.setText(f"{len(features)} {label}")

    def _on_table_selection(self):
        items = self.table.selectedItems()
        if not items:
            return
        row = items[0].row()
        feat_item = self.table.item(row, 0)
        if feat_item is None:
            return
        self.selected_feature = feat_item.text()
        self.feature_selected.emit(self.selected_feature)
        self._update_viz()

    def select_feature(self, feature):
        """Programmatically select a feature (e.g. from a volcano click)."""
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and item.text() == feature:
                self.table.selectRow(r)
                self.table.scrollToItem(item)
                return

    # ---- DE run ----------------------------------------------------------

    def _restrict_filters(self):
        filters = {}
        for col, combo in self._restrict_combos.items():
            val = combo.currentText()
            if val and val != "(any)":
                filters[col] = val
        return filters

    def _samples_for(self, contrast_col, value, extra_filters):
        filters = dict(extra_filters)
        filters[contrast_col] = value
        return self.experimental.filenames_for(**filters)

    def _paired_pairs(self, samples_a, samples_b):
        """Match A↔B filenames by their Pair column value."""
        pair_col = self._pair_column()
        if not pair_col:
            return None
        by_pair_a = {}
        for f in samples_a:
            row = self.experimental.by_filename.get(f, {})
            by_pair_a.setdefault(row.get(pair_col, ""), f)
        pairs = []
        for f in samples_b:
            row = self.experimental.by_filename.get(f, {})
            key = row.get(pair_col, "")
            if key in by_pair_a:
                pairs.append((by_pair_a[key], f))
        return pairs

    def run_de(self):
        contrast_col = self.contrast_combo.currentText()
        group_a = self.group_a_combo.currentText()
        group_b = self.group_b_combo.currentText()
        if not contrast_col or not group_a or not group_b:
            self.status_label.setText("Pick a Group column and two values first.")
            return
        if group_a == group_b:
            self.status_label.setText("Group A and B must differ.")
            return

        extra = self._restrict_filters()
        samples_a = self._samples_for(contrast_col, group_a, extra)
        samples_b = self._samples_for(contrast_col, group_b, extra)
        if not samples_a or not samples_b:
            self.status_label.setText(
                f"No files for {group_a} ({len(samples_a)}) / {group_b} "
                f"({len(samples_b)}) under the current filters.")
            return

        paired = None
        if self.test_combo.currentText().startswith("Paired"):
            paired = self._paired_pairs(samples_a, samples_b)
            if not paired:
                self.status_label.setText(
                    "Paired test needs a Pair column with matching ids across "
                    "the two groups — none found.")
                return

        matrix = self._current_matrix()
        records = differential_expression(
            matrix.keys(), samples_a, samples_b, matrix,
            paired_pairs=paired, min_replicates=self.min_rep_spin.value())

        self.de_records = records
        self._de_by_feature = {r["feature"]: r for r in records}
        self._contrast_labels = (group_a, group_b)

        n_sig = sum(1 for r in records
                    if r["fdr"] is not None and not math.isnan(r["fdr"])
                    and r["fdr"] < 0.05)
        test_name = "paired" if paired else "Welch"
        pair_note = f", {len(paired)} pairs" if paired else ""
        self.status_label.setText(
            f"{group_a} (n={len(samples_a)}) vs {group_b} (n={len(samples_b)}) — "
            f"{test_name} test{pair_note}. {len(records)} features tested, "
            f"{n_sig} at FDR<0.05.")

        self._refresh_features()
        self._sort_table_by_p()
        self._draw_volcano()

    def _sort_table_by_p(self):
        # column 5 = p
        self.table.sortItems(5, Qt.AscendingOrder)

    # ---- visualization ---------------------------------------------------

    def _update_viz(self):
        if not self._active:
            return
        self.viz_plot.clear()
        # clear() leaves a stale legend that would accumulate duplicate labels
        if self._legend is not None:
            try:
                self.viz_plot.getPlotItem().removeItem(self._legend)
            except Exception:
                pass
            self._legend = None
        pal = palette(self.theme)
        feat = self.selected_feature
        if not feat:
            self.viz_title.setText("Select a feature below")
            return

        matrix = self._current_matrix()
        per_file = matrix.get(feat, {})
        self.viz_title.setText(feat)

        series_col = self._series_column()
        contrast_col = self.contrast_combo.currentText()
        if series_col:
            self._draw_series(feat, per_file, series_col, contrast_col, pal)
        else:
            self._draw_groups(feat, per_file, contrast_col, pal)

    def _group_key_for_file(self, filename, columns):
        row = self.experimental.by_filename.get(filename, {})
        return tuple(row.get(c, "") for c in columns)

    def _draw_series(self, feat, per_file, series_col, contrast_col, pal):
        """Titration/time-course: x = ordered series values, one line per group
        (all Group columns combined), y = mean log2 quantity with replicate dots."""
        self.viz_plot.setLabel("bottom", series_col)
        group_cols = [c for c in self._group_columns()]
        series_vals = _sorted_series_values(self.experimental.values(series_col))
        x_index = {v: i for i, v in enumerate(series_vals)}

        # group_key -> {series_val: [log2 quantities]}
        grouped = {}
        for fname, q in per_file.items():
            if not q or q <= 0:
                continue
            row = self.experimental.by_filename.get(fname, {})
            sval = row.get(series_col, "")
            if sval not in x_index:
                continue
            gkey = tuple(row.get(c, "") for c in group_cols) if group_cols else ("all",)
            grouped.setdefault(gkey, {}).setdefault(sval, []).append(math.log2(q))

        self._legend = self.viz_plot.addLegend(offset=(-10, 10))
        for gi, (gkey, series_map) in enumerate(sorted(grouped.items())):
            color = GROUP_COLORS[gi % len(GROUP_COLORS)]
            xs, ys, sx, sy = [], [], [], []
            for sval in series_vals:
                vals = series_map.get(sval)
                if not vals:
                    continue
                xi = x_index[sval]
                xs.append(xi)
                ys.append(sum(vals) / len(vals))
                for v in vals:
                    sx.append(xi)
                    sy.append(v)
            name = "/".join(gkey) if group_cols else "all"
            if xs:
                self.viz_plot.plot(xs, ys, pen=pg.mkPen(color, width=2), name=name)
            if sx:
                self.viz_plot.plot(sx, sy, pen=None, symbol="o", symbolSize=6,
                                   symbolBrush=color, symbolPen=None)
        ax = self.viz_plot.getAxis("bottom")
        ax.setTicks([list(enumerate(series_vals))])

    def _draw_groups(self, feat, per_file, contrast_col, pal):
        """Per-group strip plot with a mean bar for each value of the contrast
        column (falls back to all Group columns combined if no contrast set)."""
        self.viz_plot.setLabel("bottom", contrast_col or "group")
        group_cols = [contrast_col] if contrast_col else self._group_columns()
        if not group_cols:
            self.viz_plot.setLabel("bottom", "file")
            group_cols = None

        buckets = {}  # label -> [log2 quantities]
        for fname, q in per_file.items():
            if not q or q <= 0:
                continue
            if group_cols:
                row = self.experimental.by_filename.get(fname, {})
                label = "/".join(row.get(c, "") for c in group_cols)
            else:
                label = fname
            buckets.setdefault(label, []).append(math.log2(q))

        labels = sorted(buckets.keys())
        for gi, label in enumerate(labels):
            color = GROUP_COLORS[gi % len(GROUP_COLORS)]
            vals = buckets[label]
            jitter = [gi + (0.12 * ((k % 5) - 2) / 2.0) for k in range(len(vals))]
            self.viz_plot.plot(jitter, vals, pen=None, symbol="o", symbolSize=7,
                               symbolBrush=color, symbolPen=None)
            mean = sum(vals) / len(vals)
            self.viz_plot.plot([gi - 0.25, gi + 0.25], [mean, mean],
                               pen=pg.mkPen(color, width=3))
        ax = self.viz_plot.getAxis("bottom")
        ax.setTicks([list(enumerate(labels))])

    # ---- volcano ---------------------------------------------------------

    def _draw_volcano(self):
        if not self._active:
            return
        self.volcano_plot.clear()
        if not self.de_records:
            return
        pal = palette(self.theme)
        xs, ys, brushes, feats = [], [], [], []
        for rec in self.de_records:
            p = rec.get("fdr")
            fc = rec.get("log2fc")
            if p is None or math.isnan(p) or fc is None or math.isnan(fc):
                continue
            y = -math.log10(max(p, 1e-300))
            xs.append(fc)
            ys.append(y)
            feats.append(rec["feature"])
            sig = p < 0.05 and abs(fc) >= 1.0
            if sig and fc > 0:
                brushes.append(pg.mkBrush(196, 78, 82, 200))    # up in A: red
            elif sig and fc < 0:
                brushes.append(pg.mkBrush(76, 114, 176, 200))   # up in B: blue
            else:
                brushes.append(pg.mkBrush(140, 140, 140, 120))

        scatter = pg.ScatterPlotItem(x=xs, y=ys, brush=brushes, pen=None, size=7)
        scatter.feats = feats
        scatter.sigClicked.connect(self._on_volcano_click)
        self.volcano_plot.addItem(scatter)
        self.volcano_scatter = scatter

        # significance guides
        fg = pal["fg"]
        line_pen = pg.mkPen(fg, width=1, style=Qt.DashLine)
        self.volcano_plot.addLine(y=-math.log10(0.05), pen=line_pen)
        self.volcano_plot.addLine(x=1.0, pen=line_pen)
        self.volcano_plot.addLine(x=-1.0, pen=line_pen)

        if getattr(self, "_contrast_labels", None):
            a, b = self._contrast_labels
            self.volcano_plot.setLabel("bottom", f"log2 fold change ({a} − {b})")

    def _on_volcano_click(self, scatter, points):
        if not points:
            return
        idx = points[0].index()
        feats = getattr(scatter, "feats", [])
        if 0 <= idx < len(feats):
            self.select_feature(feats[idx])

    # ---- theming ---------------------------------------------------------

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)
        for plot in (getattr(self, "viz_plot", None),
                     getattr(self, "volcano_plot", None)):
            if plot is not None:
                style_plot(plot, pal)
        self._update_viz()
        self._draw_volcano()


def _fmt_p(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return ""
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}"
