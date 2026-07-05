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
  * a **Peptides ⇄ Proteins** switch and a **unique-only** filter, over a table
    that is **flat by default and pivots on demand**. With **no Nested Layers**
    (the default) it shows the full flat view: one column per distinct value of
    every category, side by side, all visible immediately. Adding **Nested Layers**
    (ordered dropdowns — layer 1 outermost, "+ Add layer" appends the next) *pivots*
    that same data into nested combinations with a spanning multi-level header (e.g.
    condition over replicate). Fixed ``feature`` + ``unique`` columns lead. Each cell
    is the feature's quantity **averaged across every file in that column** (zeros/
    missing never averaged in); a cell that averaged ≥2 files is shaded grey.
    Independent of the panel-1 organizer and of the replicate designation.

An optional **Normalize** mode (median-center) shifts each file so its median
log2 quantity matches the grand median, correcting systematic per-run intensity
differences before anything is compared.

Filenames from the search tables (which may carry a ``.centroid.mzML`` suffix)
are matched to the design's ``filename`` column by stripped stem, so the join is
robust to extension differences.

Everything is reactive — no run button.
"""

import json
import math
import statistics

from PySide6.QtCore import Qt, QRect, QSettings, QSize, Signal
from PySide6.QtGui import QColor, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
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
CHOOSE_LABEL = "(choose column)"
MODE_COLS = "Split → columns"
MODE_ROWS = "Split ↓ rows"
MODE_XAXIS = "X-axis"
MODES = [MODE_COLS, MODE_ROWS, MODE_XAXIS]


def _norm(name):
    """Filename → comparable stem (drop mzML/centroid/raw suffixes, lowercase)."""
    return strip_ms_suffix(name).lower()


class RotatedAxisItem(pg.AxisItem):
    """Bottom axis that draws its tick labels slanted by ``angle`` degrees, so
    long or numerous category labels stay legible instead of overlapping (or
    being silently dropped, as pyqtgraph does when horizontal labels collide).
    ``angle=0`` falls back to the normal horizontal rendering."""

    def __init__(self, *args, angle=45, **kwargs):
        super().__init__(*args, **kwargs)
        self._angle = angle

    def set_angle(self, angle):
        self._angle = angle
        self.picture = None
        self.update()

    def drawPicture(self, p, axisSpec, tickSpecs, textSpecs):
        if not self._angle:
            return super().drawPicture(p, axisSpec, tickSpecs, textSpecs)
        p.setRenderHint(p.RenderHint.TextAntialiasing, True)
        pen, p1, p2 = axisSpec
        p.setPen(pen)
        p.drawLine(p1, p2)
        for tpen, tp1, tp2 in tickSpecs:
            p.setPen(tpen)
            p.drawLine(tp1, tp2)
        if self.style["tickFont"] is not None:
            p.setFont(self.style["tickFont"])
        p.setPen(self.textPen())
        for rect, flags, text in textSpecs:
            # Anchor at the top-centre of each label slot (just under its tick)
            # and rotate; point-based drawText doesn't clip, so long labels show
            # in full running down-right.
            p.save()
            p.translate(rect.center().x(), rect.top())
            p.rotate(self._angle)
            p.drawText(0, 0, text)
            p.restore()


class GroupedHeaderView(QHeaderView):
    """Horizontal header painting a multi-level, spanning column header: each data
    column carries a path of category values (outer→inner) and the header draws one
    row per nesting level, merging adjacent columns that share a prefix. Leading
    fixed columns (feature, unique) are drawn full-height with a corner label."""

    def __init__(self, parent=None):
        super().__init__(Qt.Horizontal, parent)
        self._paths = []      # per logical column: tuple(values) or None (fixed)
        self._nlevels = 1
        self._corners = []    # labels for the leading fixed columns
        self._fg = QColor("#e6e6e6")
        self._bg = QColor("#16181d")
        self._line = QColor("#2a2d33")
        self.setSectionsClickable(True)

    def _level_h(self):
        return self.fontMetrics().height() + 8

    def set_structure(self, paths, nlevels, corners, fg, bg, line):
        self._paths = paths
        self._nlevels = max(1, nlevels)
        self._corners = list(corners)
        self._fg, self._bg, self._line = QColor(fg), QColor(bg), QColor(line)
        self.setFixedHeight(self._level_h() * self._nlevels)
        self.updateGeometry()
        self.viewport().update()

    def sizeHint(self):
        return QSize(super().sizeHint().width(), self._level_h() * self._nlevels)

    def paintEvent(self, event):
        p = QPainter(self.viewport())
        lh = self._level_h()
        total_h = lh * self._nlevels
        n = self.count()
        nfixed = len(self._corners)
        p.fillRect(self.viewport().rect(), self._bg)

        for c in range(min(nfixed, n)):
            x0 = self.sectionViewportPosition(c)
            w0 = self.sectionSize(c)
            rect = QRect(int(x0), 0, int(w0), int(total_h))
            p.setPen(QPen(self._line))
            p.drawRect(rect)
            p.setPen(QPen(self._fg))
            p.drawText(rect, Qt.AlignCenter, self._corners[c])

        fm = self.fontMetrics()

        # Flat view (single level): draw each data column with its own value, no
        # merging (values from different categories must not span together).
        if self._nlevels == 1:
            for c in range(nfixed, n):
                path = self._paths[c] if c < len(self._paths) else None
                if path is None:
                    continue
                x0 = self.sectionViewportPosition(c)
                rect = QRect(int(x0), 0, int(self.sectionSize(c)), int(lh))
                p.setPen(QPen(self._line))
                p.drawRect(rect)
                val = path[0] if path else ""
                label = str(val) if val != "" else "(blank)"
                p.setPen(QPen(self._fg))
                p.drawText(rect, Qt.AlignCenter,
                           fm.elidedText(label, Qt.ElideRight, rect.width() - 4))
            p.end()
            return

        for level in range(self._nlevels):
            c = nfixed
            while c < n:
                path = self._paths[c] if c < len(self._paths) else None
                if path is None:
                    c += 1
                    continue
                key = path[:level + 1]
                start = c
                while (c < n and c < len(self._paths)
                       and self._paths[c] is not None
                       and self._paths[c][:level + 1] == key):
                    c += 1
                last = c - 1
                x0 = self.sectionViewportPosition(start)
                xend = self.sectionViewportPosition(last) + self.sectionSize(last)
                rect = QRect(int(x0), int(level * lh), int(xend - x0), int(lh))
                p.setPen(QPen(self._line))
                p.drawRect(rect)
                val = path[level] if level < len(path) else ""
                label = str(val) if val != "" else "(blank)"
                p.setPen(QPen(self._fg))
                p.drawText(rect, Qt.AlignCenter,
                           fm.elidedText(label, Qt.ElideRight, rect.width() - 4))
        p.end()


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
        self.on_theme_toggle = None   # set by MainWindow; wired to the theme button
        self.settings = QSettings("motif-quantification", "viewer")
        self._saved = self._load_state()
        self._restoring = True

        self.level = self._saved.get("level", "peptide")
        self.unique_only = bool(self._saved.get("unique", False))
        self.normalize = self._saved.get("normalize", "none")
        self._matrix_cache = {}
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

        # Table nesting layers — an ordered list of categories (layer 1 = outermost).
        # INDEPENDENT of the panel-1 organizer and of the replicate designation.
        # DEFAULT IS EMPTY: with no layers the table shows the full flat view (every
        # category value as a column); adding layers PIVOTS that same data.
        self._nest_layers = [c for c in self._saved.get("nest_layers", [])
                             if c in self._categories]

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
            "normalize": self.normalize,
            "logy": self.logy_check.isChecked(),
            "replicate": self.replicate_combo.currentText(),
            "layers": self._layers,
            "compare": self.compare_combo.currentText(),
            "a": self.a_combo.currentText(),
            "b": self.b_combo.currentText(),
            "nest_layers": self._nest_layers,
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

    def _matrix(self):
        """Current-level quantity matrix, per-file normalized if requested.

        ``median-center`` shifts every file so its median log2 quantity matches
        the grand median across files — the standard label-free correction for
        differing per-run loading/intensity (what makes a systematic GuHCl−NaCl
        offset go away). Works in log2 space, then maps back to linear so the
        rest of the tab is unchanged. Cached per (level, mode)."""
        key = (self.level, self.normalize)
        cached = self._matrix_cache.get(key)
        if cached is not None:
            return cached

        raw = self.model.matrix(self.level)
        if self.normalize != "median-center":
            self._matrix_cache[key] = raw
            return raw

        file_logs = {}
        for per_file in raw.values():
            for f, q in per_file.items():
                if q and q > 0:
                    file_logs.setdefault(f, []).append(math.log2(q))
        med = {f: statistics.median(v) for f, v in file_logs.items() if v}
        grand = statistics.median(list(med.values())) if med else 0.0

        out = {}
        for feat, per_file in raw.items():
            d = {}
            for f, q in per_file.items():
                if q and q > 0:
                    d[f] = 2.0 ** (math.log2(q) - med.get(f, 0.0) + grand)
            out[feat] = d
        self._matrix_cache[key] = out
        return out

    def _visible_features(self):
        matrix = self._matrix()
        feats = list(matrix.keys())
        if self.level == "peptide" and self.unique_only:
            feats = [f for f in feats if self.model.peptide_is_unique(f)]
        return sorted(feats)

    # ---- UI construction -------------------------------------------------

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        if not self._active:
            warn = QLabel("Double-click to open a folder")
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

        # Fold-change contrast summary — full-contrast, above the organizer
        # dropdowns (recoloured in apply_theme).
        self.fold_status = QLabel("")
        lay.addWidget(self.fold_status)

        # organizer pseudo-table: growable list of layer rows
        org_header = QHBoxLayout()
        org_header.addWidget(QLabel("Organize by (top → bottom = outer → inner):"))
        org_header.addStretch(1)
        self.theme_btn = QPushButton("Light Mode" if self.theme != "light" else "Dark Mode")
        self.theme_btn.setFixedWidth(90)
        self.theme_btn.clicked.connect(lambda: self.on_theme_toggle and self.on_theme_toggle())
        org_header.addWidget(self.theme_btn)
        org_header.addWidget(QLabel("Normalize:"))
        self.normalize_combo = QComboBox()
        self.normalize_combo.addItems(["none", "median-center"])
        if self.normalize in ("none", "median-center"):
            self.normalize_combo.setCurrentText(self.normalize)
        self.normalize_combo.setToolTip(
            "median-center: shift each file so its median log2 quantity matches "
            "the grand median — corrects systematic per-run loading/intensity "
            "differences before comparing across files.")
        self.normalize_combo.currentTextChanged.connect(self._on_normalize_changed)
        org_header.addWidget(self.normalize_combo)

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
        lay.addLayout(bar)

        # Nested Layers: ordered dropdowns (layer 1 = outermost). Choosing a column
        # does NOT change the table — the table re-pivots only when "+ Add layer" is
        # pressed (which applies the current layers and opens the next one).
        nest_row = QHBoxLayout()
        nest_row.addWidget(QLabel("Nested Layers:"))
        self.nest_area = QHBoxLayout()
        nest_row.addLayout(self.nest_area)
        self.add_nest_btn = QPushButton("+ Add layer")
        self.add_nest_btn.clicked.connect(self._add_nest)
        nest_row.addWidget(self.add_nest_btn)
        nest_row.addStretch(1)
        lay.addLayout(nest_row)
        self._rebuild_nest_rows()

        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.verticalHeader().setVisible(False)
        self._header = GroupedHeaderView(self.table)
        self.table.setHorizontalHeader(self._header)
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

    # ---- table nested layers --------------------------------------------

    def _rebuild_nest_rows(self):
        # Always keep at least one layer dropdown visible; an unset "(choose
        # column)" layer means the table stays flat until the user picks.
        if not self._nest_layers:
            self._nest_layers = [CHOOSE_LABEL]
        while self.nest_area.count():
            item = self.nest_area.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        for i, cat in enumerate(self._nest_layers):
            chip = QWidget()
            cl = QHBoxLayout(chip)
            cl.setContentsMargins(2, 0, 2, 0)
            cl.setSpacing(2)
            cl.addWidget(QLabel(f"{i + 1}."))     # hard-coded layer number, += 1
            combo = QComboBox()
            combo.addItems([CHOOSE_LABEL] + self._categories)
            combo.setCurrentText(cat if cat in self._categories else CHOOSE_LABEL)
            combo.setMinimumWidth(130)
            combo.currentTextChanged.connect(lambda v, idx=i: self._set_nest(idx, v))
            cl.addWidget(combo)
            rm = QPushButton("✕")
            rm.setFixedWidth(24)
            rm.clicked.connect(lambda _=False, idx=i: self._remove_nest(idx))
            cl.addWidget(rm)
            self.nest_area.addWidget(chip)

    # Choosing a column only updates the pending config; the table re-pivots when
    # "+ Add layer" (or removing a layer) is pressed — never on a mere selection.

    def _add_nest(self):
        if not self._categories:
            return
        self._refresh_table()                    # apply the current layers now
        self._nest_layers.append(CHOOSE_LABEL)   # then open the next (unset) layer
        self._rebuild_nest_rows()
        self._save_state()

    def _remove_nest(self, index):
        if 0 <= index < len(self._nest_layers):
            self._nest_layers.pop(index)
            self._rebuild_nest_rows()
            self._refresh_table()
            self._save_state()

    def _set_nest(self, index, cat):
        if 0 <= index < len(self._nest_layers):
            self._nest_layers[index] = cat
            self._save_state()

    # ---- reactive handlers ----------------------------------------------

    def _on_logy_changed(self, _):
        self._rebuild_facets()
        self._save_state()

    def _on_normalize_changed(self, mode):
        self.normalize = mode
        self._refresh_table()
        self._refresh_fold_change()
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

    # ---- feature table (nested pivot driven by the Nested Layers) --------

    def _columns_and_paths(self):
        """Return ``(columns, paths, nlevels, tooltips)`` for the current table.

        No Nested Layers → the FLAT full view: one column per distinct value of
        every category, side by side (paths are single-element, header stays flat).
        With Nested Layers → PIVOT: leaf columns = unique combinations of the
        layers' values, with a spanning multi-level header. Either way each column
        carries the list of files that feed it (averaged, non-zero, in the cells)."""
        files = self.model.filenames()
        levels = []
        for c in self._nest_layers:
            if c in self._categories and c not in levels:
                levels.append(c)

        def key(vals):
            out = []
            for v in vals:
                try:
                    out.append((0, float(v)))
                except (TypeError, ValueError):
                    out.append((1, str(v)))
            return tuple(out)

        if not levels:
            columns, paths, tips = [], [], []
            for cat in self._categories:
                for v in _sorted_values({self._cat_value(f, cat) for f in files
                                         if self._cat_value(f, cat) != ""}):
                    columns.append([f for f in files if self._cat_value(f, cat) == v])
                    paths.append((v,))
                    tips.append(f"{cat} = {v}")
            return columns, paths, 1, tips

        leaf_files = {}
        for f in files:
            leaf_files.setdefault(tuple(self._cat_value(f, c) for c in levels), []).append(f)
        leaves = sorted(leaf_files.keys(), key=key)
        columns = [leaf_files[p] for p in leaves]
        tips = [" / ".join(f"{lvl} = {val}" for lvl, val in zip(levels, p)) for p in leaves]
        return columns, list(leaves), len(levels), tips

    def _refresh_table(self):
        matrix = self._matrix()
        feats = self._visible_features()
        feat_label = "peptide" if self.level == "peptide" else "protein"
        columns, paths, nlevels, tips = self._columns_and_paths()

        grey = QColor(150, 150, 150, 235)   # opaque grey = an average of ≥2 files
        dark = QColor("#101216")
        center = Qt.AlignCenter
        fixed = [feat_label, "unique"]
        nfixed = len(fixed)

        self.table.setSortingEnabled(False)
        self.table.setColumnCount(nfixed + len(columns))
        self.table.setHorizontalHeaderLabels(fixed + [str(p[-1]) for p in paths])
        for ci, tip in enumerate(tips):
            hi = self.table.horizontalHeaderItem(nfixed + ci)
            if hi is not None:
                hi.setToolTip(tip)
        self.table.setRowCount(len(feats))

        for r, feat in enumerate(feats):
            fitem = QTableWidgetItem(feat)
            fitem.setTextAlignment(center)
            self.table.setItem(r, 0, fitem)
            uniq = ("yes" if (self.level == "protein"
                              or self.model.peptide_is_unique(feat)) else "no")
            uitem = QTableWidgetItem(uniq)
            uitem.setTextAlignment(center)
            self.table.setItem(r, 1, uitem)
            per = matrix.get(feat, {})
            for ci, colfiles in enumerate(columns):
                qs = [per.get(f) for f in colfiles]
                qs = [q for q in qs if q and q > 0]   # never average in zeros
                if not qs:
                    it = NumericItem("", float("nan"))
                else:
                    it = NumericItem(f"{sum(qs) / len(qs):.4g}", sum(qs) / len(qs))
                    if len(qs) > 1:                    # averaged across ≥2 files
                        it.setBackground(grey)
                        it.setForeground(dark)
                it.setTextAlignment(center)
                self.table.setItem(r, nfixed + ci, it)

        self.table.setSortingEnabled(True)

        fg, bg, line = self._chrome()
        self._header.set_structure([None] * nfixed + list(paths), nlevels,
                                   fixed, fg, bg, line)
        self._size_columns(feats, paths)

    def _chrome(self):
        if self.theme == "light":
            return ("#202020", "#ffffff", "#c8ccd2")
        return ("#e6e6e6", "#16181d", "#2a2d33")

    def _size_columns(self, feats, leaves):
        """Size each column to its own content — feature to the widest sampled
        sequence, the rest to their leaf-label / quantity text. No stretch-to-fill,
        so columns stay tight and re-tighten whenever layers change."""
        fm = QFontMetrics(self.table.font())
        wfeat = fm.horizontalAdvance("peptide") + 24
        for feat in feats[:400]:
            wfeat = max(wfeat, fm.horizontalAdvance(feat) + 20)
        self.table.setColumnWidth(0, wfeat)
        self.table.setColumnWidth(1, fm.horizontalAdvance("unique") + 20)
        qwidth = fm.horizontalAdvance("0.000e+00") + 16
        for ci, path in enumerate(leaves):
            self.table.setColumnWidth(
                2 + ci, max(qwidth, fm.horizontalAdvance(str(path[-1])) + 16))

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
        matrix = self._matrix()

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

        # Track leaf plots + global data extents so axes can be shared across the
        # whole facet grid (comparable heights/positions everywhere).
        self._leaf_plots = []
        self._y_lo, self._y_hi, self._x_max = math.inf, -math.inf, 0

        per_file = self._matrix().get(feat, {})
        files = [f for f in self.model.filenames() if per_file.get(f)]

        # One global x ordering shared by every leaf, so a given x position means
        # the same category everywhere (required for a shared x-axis to be valid).
        xcol = xaxis if xaxis in self._categories else None
        if xcol:
            self._facet_xvals = _sorted_values({self._cat_value(f, xcol) for f in files})
        else:
            self._facet_xvals = _sorted_values(files)
        self._facet_xcol = xcol

        widget = self._facet_node(feat, per_file, files, splits, 0, xaxis, [])
        self.facet_area.addWidget(widget)

        self._link_facet_axes()

    def _link_facet_axes(self):
        """Share X and Y across every leaf plot in the grid: link their views and
        set one common range spanning all panels' data, so a split on the x-axis
        shares x and a split on the y-axis shares y (and quantities stay directly
        comparable across the whole grid)."""
        plots = self._leaf_plots
        if len(plots) < 1:
            return
        base = plots[0]
        for p in plots[1:]:
            p.setXLink(base)
            p.setYLink(base)
        if math.isfinite(self._y_lo) and math.isfinite(self._y_hi):
            span = self._y_hi - self._y_lo or 1.0
            base.setYRange(self._y_lo - 0.05 * span, self._y_hi + 0.05 * span,
                           padding=0)
        base.setXRange(-0.5, self._x_max + 0.5, padding=0)

    def _facet_node(self, feat, per_file, files, splits, depth, xaxis, path):
        if depth >= len(splits) or not files:
            return self._leaf_plot(feat, per_file, files, xaxis, path)
        layer = splits[depth]
        col = layer["cat"]
        orient = Qt.Horizontal if layer["mode"] == MODE_COLS else Qt.Vertical
        values = _sorted_values({self._cat_value(f, col) for f in files})
        split = QSplitter(orient)
        for v in values:
            sub = [f for f in files if self._cat_value(f, col) == v]
            child = self._facet_node(feat, per_file, sub, splits, depth + 1,
                                     xaxis, path + [(col, v)])
            split.addWidget(child)
        return split

    def _leaf_plot(self, feat, per_file, files, xaxis, path):
        logy = self.logy_check.isChecked()
        pal = palette(self.theme)
        # Global (shared) x ordering — same category at the same position in every
        # leaf, so the linked x-axis is meaningful across the grid.
        xcol = self._facet_xcol
        xvals = self._facet_xvals
        xindex = {v: i for i, v in enumerate(xvals)}
        labels = [str(v) for v in xvals]

        def yval(q):
            return math.log2(q) if logy else q

        buckets = {}
        for f in files:
            q = per_file.get(f)
            if not q or q <= 0:
                continue
            xv = self._cat_value(f, xcol) if xcol else f
            buckets.setdefault(xv, []).append(q)

        # Decide whether the x labels need slanting: many ticks, or any long
        # label, would overlap horizontally. Estimate real widths and compare to
        # a rough per-tick budget so this adapts to the actual labels.
        fm = QFontMetrics(self.font())
        widest = max((fm.horizontalAdvance(s) for s in labels), default=0)
        angle = 45 if (len(labels) > 5 or widest > 60) else 0

        axis = RotatedAxisItem(orientation="bottom", angle=angle)
        plot = pg.PlotWidget(axisItems={"bottom": axis})
        style_plot(plot, pal)
        for name in ("left", "bottom"):
            plot.getAxis(name).enableAutoSIPrefix(False)
        plot.setLabel("left", "log2 quantity" if logy else "quantity")
        # Title names the split slice this sub-plot represents (the full path of
        # column=value choices that led here), so every panel is self-describing.
        if path:
            title = " / ".join((v if v != "" else "(blank)") for _, v in path)
            plot.setTitle(title, color=pal["fg"], size="10pt")

        pts = pal["points"]
        sx, sy = [], []
        for xv, qs in buckets.items():
            i = xindex.get(xv)
            if i is None:
                continue
            for q in qs:
                sx.append(i)
                sy.append(yval(q))

        plot.plot(sx, sy, pen=None, symbol="o", symbolSize=7,
                  symbolBrush=pg.mkBrush(pts[0], pts[1], pts[2], 220),
                  symbolPen=pg.mkPen(pal["fg"], width=0.5))

        # feed the shared-axis extents (see _link_facet_axes)
        if sy:
            self._y_lo = min(self._y_lo, min(sy))
            self._y_hi = max(self._y_hi, max(sy))
        self._x_max = max(self._x_max, len(labels) - 1)
        self._leaf_plots.append(plot)

        axis.setTicks([list(enumerate(labels))])
        if angle:
            # reserve vertical room for the slanted labels so they aren't clipped
            axis.setHeight(int(widest * 0.72) + 26)
        plot.setLabel("bottom", xcol if xcol else "file")
        return plot

    # ---- theming ---------------------------------------------------------

    def apply_theme(self, theme):
        self.theme = theme
        if not self._active:
            return
        if getattr(self, "theme_btn", None) is not None:
            self.theme_btn.setText("Light Mode" if theme != "light" else "Dark Mode")
        if getattr(self, "fold_status", None) is not None:
            # full-contrast (black on light, white on dark), not muted grey
            fg = "#101216" if theme == "light" else "#ffffff"
            self.fold_status.setStyleSheet(f"color: {fg}; font-weight: bold;")
        pal = palette(theme)
        if getattr(self, "fold_plot", None) is not None:
            style_plot(self.fold_plot, pal)
        self._refresh_fold_change()
        self._rebuild_facets()
        if getattr(self, "table", None) is not None:
            self._refresh_table()   # recolours the grouped header + averaged cells
