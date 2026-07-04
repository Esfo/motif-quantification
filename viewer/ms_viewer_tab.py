"""Tab 1 - MS viewing workspace, built from movable/resizable dock panels.

Layout (the default the spec drew):
    left   : list 1 proteins / list 2 peptides / list 3 PSMs, each with an
             "All" button; selecting in one filters the others.
    panel 1: m/z x intensity (2D) or m/z x time x intensity (3D); centroid or
             profile. The 2D/3D and centroid/profile toggles are fixed-size so
             they never jump when the label flips.
    panel 2: RT (x) x m/z (y) map of the exact window shown in panel 1. Panning
             panel 1's m/z (in 2D) and panel 2's m/z stay in sync.
    panel 3: MS1 theoretical-vs-experimental isotope overlay (when a peptide is
             identified) or the MS2 spectrum.
    table 1: per-line metrics for the selected distribution.

Panes are QDockWidgets, so they can be dragged, floated, stacked and resized.
The arrangement is persisted via QSettings between runs; "Reset layout"
restores the default captured at first build.

Staged (documented in ARCHITECTURE.md, not yet wired): distribution/line
selection colouring, charge-search navigation with history, profile<->centroid
peak linking, the colour-gradient pickers, and sequence coverage. The structure
here gives each of those a defined home.
"""

import re
import traceback

import numpy as np
from PySide6.QtCore import (
    QAbstractTableModel,
    QEvent,
    QModelIndex,
    QSortFilterProxyModel,
    Qt,
    QThread,
    QTimer,
    Signal,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QPushButton,
    QStackedWidget,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTabBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    import distinctipy
    HAVE_DISTINCTIPY = True
except Exception:
    distinctipy = None
    HAVE_DISTINCTIPY = False

try:
    from .mzml_store import MzmlStore, scan_arrays
    from .plots import plot_points, plot_spectrum
    from .region_view import HAVE_GL, gl
    from .session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
    from . import isotopes
    from . import coverage as seqcoverage
    from . import fragments as seqfragments_module
    from .theming import palette, style_plot, style_gl
    from .distributions_db import DistributionsDB
    from .points_store import PointsStore
except ImportError:
    from mzml_store import MzmlStore, scan_arrays
    from plots import plot_points, plot_spectrum
    from region_view import HAVE_GL, gl
    from session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
    import isotopes
    import coverage as seqcoverage
    import fragments as seqfragments_module
    from theming import palette, style_plot, style_gl
    from distributions_db import DistributionsDB
    from points_store import PointsStore


class NumericItem(QTableWidgetItem):
    """Table item that sorts by a numeric value when one is set, else by text.

    Used for Table 2's q / coverage-score columns so header-click sorting orders
    them by magnitude rather than lexically ('0.9' before '0.1' etc.)."""

    def __init__(self, text, value=None):
        super().__init__(text)
        self._value = value

    def __lt__(self, other):
        a = self._value
        b = getattr(other, "_value", None)
        if a is not None and b is not None:
            return a < b
        return super().__lt__(other)


def plain_seq(peptide):
    """Strip flanks and modifications to a bare uppercase residue sequence."""
    value = peptide or ""
    if len(value) >= 5 and value[1] == "." and value[-2] == ".":
        value = value[2:-2]
    value = re.sub(r"\[[^\]]*\]|\([^\)]*\)|\{[^\}]*\}", "", value)
    return re.sub(r"[^A-Za-z]", "", value).upper()


def peptide_mod_mass(peptide):
    """Total modification delta mass encoded in a peptide string, e.g. the
    +15.994915 in ``DVFLGM[+15.994915]FLYEYAR`` (oxidation) or the +57.021465
    of a fixed carbamidomethyl. Sums every signed number inside [], () or {}.

    The theoretical MS1 overlay is computed from the bare (modification-stripped)
    sequence, so its isotope m/z must be shifted by mod_mass / charge or it
    lands ~mod_mass/z Th away from the real, modified precursor."""
    total = 0.0
    for token in re.findall(r"[\[\(\{]([^\]\)\}]*)[\]\)\}]", peptide or ""):
        m = re.search(r"[-+]?\d+(?:\.\d+)?", token)
        if m:
            try:
                total += float(m.group(0))
            except ValueError:
                pass
    return total


class EvidenceWorker(QThread):
    """Reads everything a selection needs from the mzML off the UI thread.

    All mzML access for one selection happens here (metadata parse, the MS1 scan
    window, and the RT x m/z region), so the (multi-GB) file is never touched on
    the UI thread and the two reads never race on the same reader.
    """

    done = Signal(object)

    def __init__(self, centroid, store, scan, rt, mz_min, mz_max,
                 rt_start, rt_end, mz_bins, mode):
        super().__init__()
        self.centroid = centroid
        self.store = store
        self.scan = scan
        self.rt = rt
        self.mz_min = mz_min
        self.mz_max = mz_max
        self.rt_start = rt_start
        self.rt_end = rt_end
        self.mz_bins = mz_bins
        self.mode = mode

    def run(self):
        try:
            self.centroid.load_metadata()
            if self.store is not self.centroid:
                self.store.load_metadata()

            if self.rt is not None:
                ms1 = self.centroid.nearest_ms1_by_rt(self.rt)
            else:
                ms1 = self.centroid.preceding_ms1_for_scan(self.scan)

            scan_mz = scan_int = None
            ms1_number = ms1.number if ms1 is not None else None
            if self.rt is not None:
                # Overlay scan must come from the SAME store as the points (its
                # scan numbering matches); the points store keys by ms1_index, not
                # the centroid's scan numbers.
                src_ms1 = self.store.nearest_ms1_by_rt(self.rt)
                if src_ms1 is not None:
                    scan_mz, scan_int = self.store.scan_window_by_number(
                        src_ms1.number, self.mz_min, self.mz_max)
            elif ms1 is not None:
                scan_mz, scan_int = self.centroid.scan_window_by_number(
                    ms1.number, self.mz_min, self.mz_max)

            points = None
            region = None
            ms2 = []
            if self.rt is not None:
                points = self.store.extract_points(self.mz_min, self.mz_max, self.rt_start, self.rt_end)
                region = self.store.extract_region(
                    self.mz_min, self.mz_max, self.rt_start, self.rt_end,
                    mz_bins=self.mz_bins, mode=self.mode)
                # MS2 scans always come from the centroid run (where they live).
                ms2 = [{"rt": s.rt, "mz": s.precursor_mz, "number": s.number, "id": s.spectrum_id,
                        "iso_low": s.iso_low, "iso_high": s.iso_high}
                       for s in self.centroid.ms2_in_rt(self.rt_start, self.rt_end)]

            self.done.emit({"ms1_number": ms1_number, "scan_mz": scan_mz,
                            "scan_int": scan_int, "points": points, "region": region, "ms2": ms2})
        except Exception as exc:
            self.done.emit({"error": f"{exc}\n{traceback.format_exc()}"})


class Table1Worker(QThread):
    """Reads Table 1's line metrics off the UI thread.

    The distribution lookup + member query hit the (read-only) distributions
    sqlite, which can be slow; doing it here keeps selection responsive. SQLite
    connections are thread-bound, so the worker opens its OWN ``DistributionsDB``
    on ``db_path`` rather than touching the UI thread's connection. ``token``
    lets the tab discard results from a superseded selection (latest-wins).
    """

    done = Signal(object)

    def __init__(self, db_path, token, distribution_id=None, window=None):
        super().__init__()
        self.db_path = db_path
        self.token = token
        self.distribution_id = distribution_id
        self.window = window  # dict(mz_min, mz_max, rt_start, rt_end, charge)

    def run(self):
        result = {"token": self.token, "distribution_id": None, "rows": []}
        try:
            db = DistributionsDB(self.db_path)
            try:
                did = self.distribution_id
                if did is None and self.window is not None:
                    dists = db.distributions_in_window(limit=1, **self.window)
                    if dists:
                        did = dists[0]["distribution_id"]
                if did is not None:
                    result["distribution_id"] = did
                    result["rows"] = db.distribution_members(did)
            finally:
                db.close()
        except Exception as exc:
            result["error"] = f"{exc}\n{traceback.format_exc()}"
        self.done.emit(result)


class FragmentWorker(QThread):
    """Computes a peptide's theoretical b/y fragment ions and matches them to the
    MS2 spectrum off the UI thread.

    The fragment-isotope walk (``fragments.peptide_fragment_ions``) is heavy, so
    it runs here; the result drives the green/red overlay in panel 3. ``token``
    discards a superseded peptide selection (latest-wins).
    """

    done = Signal(object)

    def __init__(self, seq, charge, iso_low, iso_high, peak_mz, peak_int, ppm,
                 token, dividing_threshold, subiso_depth):
        super().__init__()
        self.seq = seq
        self.charge = charge
        self.iso_low = iso_low
        self.iso_high = iso_high
        self.peak_mz = peak_mz
        self.peak_int = peak_int
        self.ppm = ppm
        self.token = token
        self.dividing_threshold = dividing_threshold
        self.subiso_depth = subiso_depth

    def run(self):
        result = {"token": self.token, "seq": self.seq,
                  "matched": [], "unmatched": (np.array([]), np.array([]))}
        try:
            entries = seqfragments_module.peptide_fragment_ions(
                self.seq, self.charge, self.iso_low, self.iso_high,
                dividingthreshold=self.dividing_threshold,
                subisotopomericdepth=self.subiso_depth)
            matched, unmatched = seqfragments_module.annotate_spectrum(
                entries, self.peak_mz, self.peak_int, ppm=self.ppm)
            result["matched"] = matched
            result["unmatched"] = unmatched
        except Exception as exc:
            result["error"] = f"{exc}\n{traceback.format_exc()}"
        self.done.emit(result)


class DbTableWorker(QThread):
    """Loads a whole Table-1 tab (distributions / charge distributions) off the
    UI thread so the tabs can populate by default without freezing on open.
    Opens its own read-only DistributionsDB (sqlite connections are thread-bound)."""

    done = Signal(object)

    def __init__(self, db_path, kind, token):
        super().__init__()
        self.db_path = db_path
        self.kind = kind          # "distributions" | "charge"
        self.token = token

    def run(self):
        result = {"kind": self.kind, "token": self.token, "rows": []}
        try:
            db = DistributionsDB(self.db_path)
            try:
                if self.kind == "distributions":
                    result["rows"] = db.all_distributions()
                elif self.kind == "charge":
                    result["rows"] = db.all_analytes_multicharge()
            finally:
                db.close()
        except Exception as exc:
            result["error"] = f"{exc}\n{traceback.format_exc()}"
        self.done.emit(result)


# Horizontal alignment of the m/z (x) axis across panel 1 (2D & 3D) and panel 2.
# Panel 2's plot starts at the MS2 strip (pure bands) + its OWN real RT axis;
# panel 1's 2D left axis and the 3D's left gutter are set to that same total
# (PLOT_LEFT) so a given m/z sits at the same screen x in all of them. Using a
# real (fixed-width) axis on BOTH panels is what makes it deterministic.
MS2_STRIP_W = 34       # MS2 strip: pure bands, no axis
P2_AXIS_W = 60         # panel 2's RT axis width (and panel 1's intensity axis)
PLOT_LEFT = MS2_STRIP_W + P2_AXIS_W   # = where every plot's m/z axis begins
# Panel 1's 2D plot gets the SAME left structure as panel 2 -- a strip-width
# spacer + an axis of P2_AXIS_W -- so the two plot areas are pixel-identical and
# their m/z axes line up exactly (a single 80px axis on panel 1 rendered a hair
# wider than strip+axis on panel 2, shifting panel 1 right).


DIST_PALETTE = [
    (76, 114, 176), (221, 132, 82), (85, 168, 104), (196, 78, 82),
    (129, 114, 179), (147, 120, 96), (218, 139, 195), (140, 140, 140),
    (204, 185, 116), (100, 181, 205),
]


# (field, header, kind) where kind: "f"=float 4g, "t"=time 2f, "i"=int, "s"=str.
LINE_METRIC_COLUMNS = [
    ("isotope_index", "iso", "i"),
    ("mz_mean", "mean m/z", "f"),
    ("mz_min", "min m/z", "f"),
    ("mz_max", "max m/z", "f"),
    ("rt_apex", "RT", "t"),
    ("rt_start", "min t", "t"),
    ("rt_end", "max t", "t"),
    ("n_points", "n pts", "i"),
]

# 'lines' tab: every line (feature) + the charge of its distribution (0 if none).
LINES_TAB_COLUMNS = [
    ("mz_mean", "mean m/z", "f"),
    ("mz_min", "min m/z", "f"),
    ("mz_max", "max m/z", "f"),
    ("rt_apex", "RT", "t"),
    ("rt_start", "min t", "t"),
    ("rt_end", "max t", "t"),
    ("n_points", "n pts", "i"),
    ("charge", "charge", "i"),
]

# 'distributions' tab (distributionregions): one row per distribution.
DIST_TAB_COLUMNS = [
    ("charge", "z", "i"),
    ("neutral_mass", "neutral mass", "f"),
    ("mono_mz", "mono m/z", "f"),
    ("rt_apex", "RT", "t"),
    ("rt_start", "min t", "t"),
    ("rt_end", "max t", "t"),
    ("n_members", "n iso", "i"),
    ("auc", "AUC", "f"),
    ("iso_score", "iso", "f"),
    ("status", "status", "s"),
]

# 'charge distributions' tab (chargeregions): one row per analyte (expandable).
CHARGE_TAB_COLUMNS = [
    ("neutral_mass", "neutral mass", "f"),
    ("charge_min", "z min", "i"),
    ("charge_max", "z max", "i"),
    ("n_distributions", "n dist", "i"),
    ("rt_apex", "RT", "t"),
    ("rt_start", "min t", "t"),
    ("rt_end", "max t", "t"),
    ("auc", "AUC", "f"),
]


def _gtick(value):
    """Compact y-tick label so values stay narrow (no overlap with the row
    label): drop the exponent's '+'/leading zero, e.g. 1.98e+05 -> 1.98e5."""
    if value == 0:
        return "0"
    a = abs(value)
    if a >= 1e4 or a < 1e-2:
        s = f"{value:.2e}"
        return s.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    return f"{value:.3g}"


def _fmt(value, kind):
    if value is None:
        return ""
    if kind == "f":
        return f"{value:.4f}"
    if kind == "t":
        return f"{value:.2f}"
    if kind == "i":
        return str(int(value))
    return str(value)


class SimpleTableModel(QAbstractTableModel):
    """Read-only table over a list of dicts. DisplayRole is formatted text;
    EditRole returns the raw value so a QSortFilterProxyModel sorts numerically."""

    def __init__(self, columns, rows=None, parent=None):
        super().__init__(parent)
        self._columns = columns
        self._rows = rows or []
        self._loading = False

    def set_rows(self, rows):
        self.beginResetModel()
        self._loading = False
        self._rows = rows or []
        self.endResetModel()

    def set_loading(self, loading=True):
        """Show a single 'loading…' placeholder row while the data loads."""
        self.beginResetModel()
        self._loading = loading
        if loading:
            self._rows = []
        self.endResetModel()

    def row_dict(self, source_row):
        if 0 <= source_row < len(self._rows):
            return self._rows[source_row]
        return None

    def rowCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        return 1 if self._loading else len(self._rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._columns)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        if self._loading:
            if role == Qt.DisplayRole and index.column() == 0:
                return "loading…"
            return None
        field, _, kind = self._columns[index.column()]
        value = self._rows[index.row()].get(field)
        if role == Qt.DisplayRole:
            return _fmt(value, kind)
        if role == Qt.EditRole:
            return value
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return self._columns[section][1]
        return None


def fixed_toggle(off_text, on_text, width=70):
    """A two-state button whose label flips but whose size never changes."""
    button = QPushButton(off_text)
    button.setCheckable(True)
    button.setFixedWidth(width)
    button._off = off_text
    button._on = on_text
    return button


class FlipButton(QPushButton):
    """A mode button that cycles through labelled states on each press.

    Unlike a checkable toggle (which reads as 'pressed' / 'not pressed'), this is
    a plain push button: it always shows the label of the CURRENT mode, and each
    click advances to the next mode and fires ``on_change``. Two states make it a
    binary switch; more make it a cycle (e.g. the four noise levels). The width is
    fixed to the widest label so cycling never resizes the button.
    """

    def __init__(self, states, on_change=None, extra_pad=24, parent=None):
        super().__init__(parent)
        self._states = list(states)          # [(key, label), ...]
        self._idx = 0
        self._on_change = on_change
        fm = self.fontMetrics()
        width = max(fm.horizontalAdvance(label) for _, label in self._states) + extra_pad
        self.setFixedWidth(width)
        self.setText(self._states[0][1])
        self.clicked.connect(self._advance)

    def _advance(self):
        self._idx = (self._idx + 1) % len(self._states)
        self.setText(self._states[self._idx][1])
        if self._on_change is not None:
            self._on_change()

    def key(self):
        return self._states[self._idx][0]

    def index(self):
        return self._idx

    def set_index(self, i):
        self._idx = i % len(self._states)
        self.setText(self._states[self._idx][1])


class ArrowCycle(QWidget):
    """A ◀ label ▶ cycler: the same state machine as FlipButton, but driven by
    left/right arrows with a plain (non-button) label showing the current state.
    Exposes key()/index()/set_index() so it's a drop-in for FlipButton."""

    def __init__(self, states, on_change=None, parent=None):
        super().__init__(parent)
        self._states = list(states)
        self._idx = 0
        self._on_change = on_change
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        self.prev_btn = QPushButton("◀")
        self.prev_btn.setFixedWidth(22)
        self.next_btn = QPushButton("▶")
        self.next_btn.setFixedWidth(22)
        self.label = QLabel(self._states[0][1])
        self.label.setAlignment(Qt.AlignCenter)
        fm = self.label.fontMetrics()
        self.label.setFixedWidth(
            max(fm.horizontalAdvance(lbl) for _, lbl in self._states) + 10)
        self.prev_btn.clicked.connect(lambda: self._step(-1))
        self.next_btn.clicked.connect(lambda: self._step(1))
        lay.addWidget(self.prev_btn)
        lay.addWidget(self.label)
        lay.addWidget(self.next_btn)

    def _step(self, delta):
        self._idx = (self._idx + delta) % len(self._states)
        self.label.setText(self._states[self._idx][1])
        if self._on_change is not None:
            self._on_change()

    def key(self):
        return self._states[self._idx][0]

    def index(self):
        return self._idx

    def set_index(self, i):
        self._idx = i % len(self._states)
        self.label.setText(self._states[self._idx][1])


# Per-class noise rendering: brighter/larger for substantial noise lines, fainter
# and smaller down to lone stray points. Keys match the _noise_class values.
#   1 = noise line (a feature/line with >=5 points, in no distribution)
#   2 = small line (a feature with 2..4 points, in no distribution)
#   3 = single point (a 1-point feature, or a stray datapoint in no feature)
NOISE_STYLE = {
    1: ((205, 205, 205, 215), 2.5),
    2: ((170, 180, 205, 180), 2.0),
    3: ((150, 150, 150, 140), 1.5),
}


class MSViewerTab(QMainWindow):
    def __init__(self, session, distributions_db=None, xics_ppm=10.0, xics_rt_window=0.8, theme="dark"):
        super().__init__()
        self.session = session
        self.db = distributions_db
        self.xics_ppm = float(xics_ppm)
        # Fragment-match tolerance for MS2 ion annotation = the search fragment
        # tolerance (Sage fragment_tol, ±20 ppm in execution.xsh).
        self.frag_ppm = 20.0
        # Precursor tolerance for picking Table-2 candidates = the search
        # precursor tolerance (Sage precursor_tol, ±10 ppm) together with the
        # search isotope-error offsets (Sage isotope_errors, -1..+2). These are
        # the search parameters from execution.xsh, not a hard-coded Da window.
        self.precursor_ppm = 10.0
        self.precursor_isotope_errors = (-1, 2)
        # MS2-strip identification acceptance criteria: a scan's line is green if
        # it has a PSM with q-value <= this FDR (0.1% default), else red. Editable
        # via the "acceptance criteria" field on the panel-1 bar (value is a %).
        self._fdr_threshold = 0.001
        self.rt_half = float(xics_rt_window)
        self.mz_half = 2.5
        self.theme = theme

        self._centroid = {}
        self._profile = {}
        self.current = None
        self.psm_rows = []
        self.worker = None
        self._pending = None
        self._win = None
        self.window = None        # [mz_min, mz_max, rt_start, rt_end] view (source of truth)
        self._loaded_window = None # the padded region actually extracted/cached
        self._p2_scatters = []    # (ScatterPlotItem, base_size) for zoom rescaling
        self._p1_scatter_item = None  # panel-1 2D centroid scatter (for rescaling)
        self._assigned = None     # bool mask: raw points in a validated/ambiguous distribution
        self._groups = []         # [(distribution_id, colour, point_mask)]
        self._rec_assigned = None # bool mask: raw points in a recovered (less-confident) distribution
        self._rec_groups = []     # recovered distribution groups (shown at noise mode >= 1)
        self._noise_class = None  # int8 per point: 0 assigned, 1 line, 2 small, 3 single
        self._last_region = None
        self.center = None        # (mz_center, rt_center) for the ± controls
        self._guard = False       # suppress range-change handling during programmatic set
        self._dist_colors = {}    # distribution_id -> stable RGB colour
        # A pool of visually-distinct colours (distinctipy) assigned per
        # distribution in first-seen order; excludes black/white (white is the
        # 3D peak-tip colour). Falls back to the fixed palette if unavailable.
        if HAVE_DISTINCTIPY:
            try:
                pool = distinctipy.get_colors(48, exclude_colors=[(0, 0, 0), (1, 1, 1)], rng=0)
                self._color_pool = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in pool]
            except Exception:
                self._color_pool = list(DIST_PALETTE)
        else:
            self._color_pool = list(DIST_PALETTE)
        self.assumed_charge = None
        # Selected-distribution state: the dotted border + charge-search lock-on.
        self._selected_dist_id = None
        self._selected_charge_group = None   # {charge: {distribution, features}} for the analyte
        self._selected_bbox = None           # (mz_min, mz_max, rt_start, rt_end) of the selection
        self._sel_border_item = None         # the dotted-rectangle item on panel 2
        self._border_color = None            # None = normal (fg); "red" = hypothetical box
        self._selected_rt_band = None        # (rt_start, rt_end) of the analyte, for hypotheticals
        self._selected_n_members = None      # isotope count of the selection, for hypothetical width
        self._panel3_mode = "ms1"  # "ms1" (envelope/charge grid) or "ms2" (spectrum)
        self._ms2_scan = None      # the MS2 scan currently shown in panel 3
        self._ms2_band = None      # (low, high, rt) isolation band kept on panel 2
        self._ms1_theo = None      # {seq, charge} of the Table-2 peptide overlaid on panel 1
        self._ms1_theo_item = None # the BarGraphItem for that overlay
        self._nav = []            # navigation history of (window, charge)
        self._nav_i = -1
        self._suppress_record = False
        try:
            self._cmap = pg.colormap.get("viridis")
        except Exception:
            self._cmap = None

        from PySide6.QtCore import QTimer
        self._reload_timer = QTimer(self)
        self._reload_timer.setSingleShot(True)
        self._reload_timer.setInterval(140)
        self._reload_timer.timeout.connect(self.do_extract)

        # Debounced panel-1 redraw: panel 1 (and the 3D) must be re-filtered to
        # panel 2's visible window as it's zoomed/panned, even when no reload is
        # needed (zoom within the cached region).
        self._p1_view_timer = QTimer(self)
        self._p1_view_timer.setSingleShot(True)
        self._p1_view_timer.setInterval(60)
        self._p1_view_timer.timeout.connect(self._redraw_panel1_view)

        self.setDockNestingEnabled(True)
        self.build_lists_dock()
        self.build_panel1_dock()
        self.build_panel2_dock()
        self.build_panel3_dock()
        self.build_table1_dock()
        self.build_table2_dock()
        self.arrange_default()
        # Panes are resize-only: not movable, floatable, or closable. The user
        # can still drag the splitter edges between them to resize, but can't
        # rearrange/float them out of the layout.
        for dock in (self.dock_lists, self.dock_panel1, self.dock_panel2,
                     self.dock_panel3, self.dock_table1, self.dock_table2):
            dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        self._default_state = self.saveState()
        self.apply_theme(theme)
        # Table 2 is MS2-only: hidden until the user clicks an MS2 point, so the
        # MS1 panel 3 occupies the full panel-3 + table-2 space.
        self.dock_table2.hide()
        # Load the distributions / charge-distributions tabs up front (in the
        # background, with a 'loading…' placeholder) instead of lazily on click.
        self._eager_load_table1_tabs()

    # ---- docks -----------------------------------------------------------

    def build_lists_dock(self):
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(4, 4, 4, 4)

        # Tab 1 is a single-file view: a file must be chosen first. Default to
        # the first file so the lists are never empty/ambiguous on open.
        # Only files that have their OWN distributions sqlite are openable; a file
        # without one is omitted rather than shown with another file's overlay.
        self.file_combo = QComboBox()
        for row in self.session.files():
            name = row.get("filename", "")
            if name and self.session.distributions_db_for(name) is not None:
                self.file_combo.addItem(name, name)
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)
        layout.addWidget(self.file_combo)

        self.protein_list = self._titled_list(layout, "proteins", self.on_protein_selected, self.show_all_proteins)
        self.peptide_list = self._titled_list(layout, "peptides", self.on_peptide_selected, self.show_all_peptides)
        self.psm_list = self._titled_list(layout, "PSMs", self.on_psm_selected, self.show_all_psms)

        dock = QDockWidget("", self)
        dock.setObjectName("dock_lists")
        dock.setTitleBarWidget(QWidget())   # drop the "Lists" title bar (dead space)
        dock.setWidget(container)
        self.dock_lists = dock

        self.current_file = self.file_combo.currentData()
        self._set_db_for_current_file()
        self.repopulate_active_list()

    def _set_db_for_current_file(self):
        """Bind self.db to the sqlite for the currently selected file (or None)."""
        if getattr(self, "db", None) is not None:
            try:
                self.db.close()
            except Exception:
                pass
        db_path = self.session.distributions_db_for(self.current_file or "")
        self.db = DistributionsDB(db_path) if db_path is not None else None
        if hasattr(self, "table1_tabs"):
            self.reset_table1_tabs()

    def on_file_changed(self):
        self.current_file = self.file_combo.currentData()
        self._set_db_for_current_file()
        self.repopulate_active_list()

    def _titled_list(self, layout, title, on_select, on_all):
        header = QHBoxLayout()
        header.addWidget(QLabel(title))
        header.addStretch(1)
        all_button = QPushButton("All")
        all_button.setFixedWidth(40)
        all_button.clicked.connect(lambda: on_all())
        header.addWidget(all_button)
        layout.addLayout(header)

        listw = QListWidget()
        listw.setSelectionMode(QAbstractItemView.SingleSelection)
        # Keep the selection highlighted even when the list loses focus (Qt
        # otherwise renders the inactive selection in a muted grey, so the blue
        # "selected" cue disappeared when you clicked a panel/table).
        listw.setStyleSheet(
            "QListWidget::item:selected,"
            "QListWidget::item:selected:!active"
            " { background-color: #2f6fb3; color: white; }"
        )
        listw.itemSelectionChanged.connect(on_select)
        layout.addWidget(listw, stretch=1)
        return listw

    def build_panel1_dock(self):
        self.p1_2d = pg.PlotWidget()
        self.p1_load_overlay = self._make_loading_overlay(self.p1_2d)
        self.p1_2d.setLabel("bottom", "m/z")
        self.p1_2d.setLabel("left", "intensity")
        # No "x0.000#" SI-prefix multiplier on the m/z axis: with no data loaded
        # the default range makes pyqtgraph pick a tiny prefix and append it to
        # the label. m/z is an absolute value, never scaled.
        self.p1_2d.getAxis("bottom").enableAutoSIPrefix(False)
        # NOTE: clip-to-view + 'peak' auto-downsampling was culling scatter points
        # as you zoomed in (they "disappeared"). Disabled so every datapoint stays
        # drawn at any zoom; the window is bounded so the point count stays sane.
        self.p1_2d.setClipToView(False)
        self.p1_2d.setDownsampling(auto=False)
        # 2D panel 1: only the m/z (x) axis is interactive; y stays auto-scaled.
        self.p1_2d.setMouseEnabled(x=True, y=False)
        # Wheel over the y-axis label strip (left of the plot) scrolls intensity.
        self.p1_2d.viewport().installEventFilter(self)
        # Double-click re-fits the intensity (y) axis to the data.
        self.p1_2d.scene().sigMouseClicked.connect(self._on_p1_clicked)

        # Spawn / "align 3D" view: a FRONT-ON, near-orthographic view that looks
        # exactly like the 2D plot -- m/z on the horizontal axis (aligned with
        # panel 2), intensity vertical with the 0 baseline at the BOTTOM. Mapping:
        # m/z -> GL-x, time -> GL-y, intensity -> GL-z (0 -> bottom). elevation 0
        # + azimuth -90 looks straight along +time, so screen-right is +m/z and
        # screen-up is +intensity. A small FOV makes it near-orthographic so the
        # m/z axis is linear and lines up with panel 2. From here the user can
        # orbit to reveal time; "align 3D" returns to exactly this.
        self._p1_3d_fov = 2.0
        # pyqtgraph's fov is the HORIZONTAL field of view, so with the m/z data
        # mapped to x in [-1, 1] this distance makes m/z fill the pane width
        # edge-to-edge (aligned with panel 2). Intensity is scaled per-draw to the
        # pane's height/width so it fills vertically with the baseline at the
        # bottom. elevation 0 + azimuth -90 = front-on (looks like the 2D plot).
        self._p1_3d_dist = 1.0 / float(np.tan(np.radians(self._p1_3d_fov / 2.0)))
        self._p1_3d_cam = dict(distance=self._p1_3d_dist, elevation=0.0, azimuth=-90.0)
        if HAVE_GL:
            self.p1_3d = gl.GLViewWidget()
            self.p1_3d.opts["fov"] = self._p1_3d_fov
            self.p1_3d.setCameraPosition(**self._p1_3d_cam)
            from PySide6.QtWidgets import QSizePolicy
            self.p1_3d.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.p1_scatter = gl.GLScatterPlotItem(pos=np.zeros((1, 3), dtype=np.float32), size=4.0)
            self.p1_scatter.setVisible(False)
            self.p1_3d.addItem(self.p1_scatter)
            # The GL widget fills the WHOLE panel (same area as the 2D plot widget);
            # the m/z data is offset within it (in draw_panel1_3d) so it starts at
            # the strip width, matching panel 2 -- the empty left gutter mirrors the
            # 2D plot's left axis area.
            self.p1_3d_widget = self.p1_3d
        else:
            self.p1_3d = QLabel("3D needs pyqtgraph OpenGL (pip install PyOpenGL)")
            self.p1_3d.setAlignment(Qt.AlignCenter)
            self.p1_3d_widget = self.p1_3d

        # Wrap the 2D plot with a strip-width spacer so its left structure matches
        # panel 2's (spacer + axis), making the m/z axes pixel-identical.
        p1_2d_row = QWidget()
        _hb = QHBoxLayout(p1_2d_row)
        _hb.setContentsMargins(0, 0, 0, 0)
        _hb.setSpacing(0)
        self._p1_2d_spacer = QWidget()
        self._p1_2d_spacer.setFixedWidth(MS2_STRIP_W)
        self._p1_2d_spacer.setAutoFillBackground(True)
        _hb.addWidget(self._p1_2d_spacer)
        _hb.addWidget(self.p1_2d, stretch=1)

        self.p1_stack = QStackedWidget()
        self.p1_stack.addWidget(p1_2d_row)           # index 0 = 2D
        self.p1_stack.addWidget(self.p1_3d_widget)   # index 1 = 3D

        # Plain mode buttons (FlipButton): each shows the current mode and cycles
        # on press -- no sunken 'pressed' look. dim/source/colour are binary; noise
        # cycles through four levels of progressively finer noise.
        self.dim_toggle = FlipButton([("2D", "2D"), ("3D", "3D")], on_change=self.toggle_dimension)
        self.source_toggle = FlipButton([("centroid", "Centroid"), ("profile", "Profile")],
                                        on_change=self.refresh)
        # Noise levels: none -> noise lines -> + small lines -> + single points.
        # Levels: none -> recovered (less confident) distributions -> noise lines
        # -> + small lines -> + single points. Each level is cumulative. Driven by
        # the ◀ ▶ arrows (not a button) -- see the bar layout below.
        self.noise_toggle = ArrowCycle(
            [("none", "No Noise"), ("recovered", "Less Confident Distributions"),
             ("line", "Line Noise"), ("small", "Small Line Noise"),
             ("single", "Single Point Noise")],
            on_change=self._on_noise_changed)
        self.logcolor_toggle = FlipButton([("lin", "Lin Color"), ("log", "Log Color")],
                                          on_change=self._on_logcolor)
        # Theoretical-MS1 overlay mode (only matters once a Table-2 peptide is
        # selected): 'raw' keeps every isotopomer separate; 'summed' merges all
        # isotopomers at the same M+N into one signal.
        self.ms1theo_toggle = FlipButton(
            [("raw", "Raw Isotopomers"), ("summed", "Summed M+N")],
            on_change=self._on_ms1theo_changed)
        # Snap the 3D view to the aligned top-down (2D-like) orientation.
        self.reset3d_button = QPushButton("Align 3D")
        self.reset3d_button.setFixedWidth(70)
        self.reset3d_button.setToolTip("align the 3D view to a top-down, panel-2-aligned 2D view")
        self.reset3d_button.clicked.connect(self.reset_3d_view)
        self.reset3d_button.setEnabled(False)

        # Theme toggle, to the right of "align 3D". main_window wires
        # on_theme_toggle; the label tracks the current theme.
        self.on_theme_toggle = None
        self.theme_btn = QPushButton("Light Mode" if self.theme == "dark" else "Dark Mode")
        self.theme_btn.setFixedWidth(80)
        self.theme_btn.clicked.connect(lambda: self.on_theme_toggle and self.on_theme_toggle())

        # The "loading" indicator now lives in the tab bar (main_window sets
        # self.tab_loading_label); nothing above panel 1 anymore.
        self.tab_loading_label = None

        # FDR acceptance-criteria spin box for the MS2-strip green/red gate:
        # "[0.10] % FDR". Scroll the cursor over it to step the value. Percentage.
        self.fdr_edit = QDoubleSpinBox()
        self.fdr_edit.setDecimals(2)
        self.fdr_edit.setRange(0.0, 100.0)
        self.fdr_edit.setSingleStep(0.1)
        self.fdr_edit.setValue(0.1)
        self.fdr_edit.setFixedWidth(60)
        self.fdr_edit.setToolTip("MS2 lines are green when a PSM passes this FDR "
                                 "(percent), red otherwise")
        self.fdr_edit.valueChanged.connect(self._on_fdr_changed)
        self.fdr_unit = QLabel("% FDR")

        # Navigation-history arrows: same design as the charge arrows, placed to
        # their left ("history: ◀ ▶"). Moved here from the (now removed) top
        # toolbar.
        self.hist_back_btn = QPushButton("◀")
        self.hist_back_btn.setFixedWidth(28)
        self.hist_back_btn.setToolTip("navigation history: back")
        self.hist_back_btn.clicked.connect(self.nav_back)
        self.hist_fwd_btn = QPushButton("▶")
        self.hist_fwd_btn.setFixedWidth(28)
        self.hist_fwd_btn.setToolTip("navigation history: forward")
        self.hist_fwd_btn.clicked.connect(self.nav_forward)

        # Charge-search arrows: right-aligned, to the right of "align 3D".
        self.charge_prev_btn = QPushButton("◀")
        self.charge_prev_btn.setFixedWidth(28)
        self.charge_prev_btn.setToolTip("charge search: one charge higher at the same RT")
        self.charge_prev_btn.clicked.connect(lambda: self.charge_step(1))
        self.charge_next_btn = QPushButton("▶")
        self.charge_next_btn.setFixedWidth(28)
        self.charge_next_btn.setToolTip("charge search: one charge lower at the same RT")
        self.charge_next_btn.clicked.connect(lambda: self.charge_step(-1))

        bar = QHBoxLayout()
        bar.addWidget(self.dim_toggle)
        bar.addWidget(self.source_toggle)
        bar.addWidget(self.logcolor_toggle)
        bar.addWidget(self.ms1theo_toggle)
        bar.addWidget(self.reset3d_button)
        bar.addWidget(self.theme_btn)
        bar.addSpacing(10)
        # Noise cycler (◀ label ▶): right of "Light Mode", left of the
        # acceptance-criteria field. Replaces the old 'loading' text slot.
        bar.addWidget(self.noise_toggle)
        bar.addStretch(1)
        bar.addWidget(self.fdr_edit)
        bar.addWidget(self.fdr_unit)
        bar.addSpacing(12)
        bar.addWidget(QLabel("History"))
        bar.addWidget(self.hist_back_btn)
        bar.addWidget(self.hist_fwd_btn)
        bar.addSpacing(12)
        bar.addWidget(QLabel("Charge"))
        bar.addWidget(self.charge_prev_btn)
        bar.addWidget(self.charge_next_btn)

        # Hide pyqtgraph's in-plot auto-range "A" button: fit-to-data makes no
        # sense for the window-driven panels and it kept overlapping the data.
        self.p1_2d.getPlotItem().hideButtons()

        container = QWidget()
        container.setObjectName("panel1_frame")
        layout = QVBoxLayout(container)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.addLayout(bar)
        layout.addWidget(self.p1_stack, stretch=1)

        dock = QDockWidget("", self)
        dock.setObjectName("dock_panel1")
        dock.setTitleBarWidget(QWidget())   # drop the empty title bar (dead space)
        dock.setWidget(container)
        self.dock_panel1 = dock

    # ---- 3D side labels --------------------------------------------------

    def _recolor_gl(self, pal):
        """No 3D labels anymore (removed per the user); nothing to recolour."""
        return

    def reset_3d_view(self):
        """Return the 3D view to the EXACT front-on, panel-2-aligned spawn view --
        including resetting the pan centre, which ctrl+drag moves (that's why it
        didn't realign after panning)."""
        if HAVE_GL:
            try:
                self.p1_3d.opts["center"] = pg.Vector(0.0, 0.0, 0.0)
            except Exception:
                pass
            self.p1_3d.opts["fov"] = self._p1_3d_fov
            self.p1_3d.setCameraPosition(**self._p1_3d_cam)
            # Re-fit the data to the (possibly resized) pane.
            if getattr(self, "_p1_3d_inputs", None) is not None:
                self.draw_panel1_3d(*self._p1_3d_inputs)

    def build_panel2_dock(self):
        # Thin MS2 strip to the left of panel 2: MS2 scans as clickable points,
        # RT-aligned with panel 2 (shared y). It fits in the space panel 1's wide
        # y-axis labels leave between panels 1 and 2.
        self.ms2_plot = pg.PlotWidget()
        # Pure band strip (no axis), fixed width. The RT scale lives on panel 2's
        # own left axis to its right; using real axes on both panels is what keeps
        # the m/z (x) axes deterministically aligned.
        self.ms2_plot.setFixedWidth(MS2_STRIP_W)
        self.ms2_plot.setMouseEnabled(x=False, y=True)
        self.ms2_plot.getPlotItem().hideAxis("bottom")
        self.ms2_plot.getPlotItem().hideAxis("left")
        # Fix the x range so the strip never collapses when RT (y) is zoomed.
        self.ms2_plot.getViewBox().setXRange(0, 1, padding=0)
        self.ms2_plot.getViewBox().setMouseEnabled(x=False, y=True)
        self.ms2_plot.getViewBox().disableAutoRange(axis=pg.ViewBox.XAxis)
        # MS2 scans = solid horizontal lines spanning the strip at each RT, drawn
        # as ONE curve (NaN-separated segments) with a fixed pixel width so they
        # are always clearly visible as lines (never dots, never vanish) and stay
        # a consistent thickness at any zoom -- zooming in just spreads them apart
        # so individual scans become distinguishable. Clicking the strip selects
        # the nearest line by RT (scene click; a single curve isn't per-segment
        # clickable).
        # Two curves: GREEN for scans with an identified peptide passing the FDR
        # acceptance criteria, RED for scans without one.
        self._ms2_scans = []    # sorted [(rt, scan_dict)] for click hit-testing
        self.ms2_curve = pg.PlotCurveItem([], [], connect="finite",
                                          pen=pg.mkPen(60, 200, 100, 255, width=3))
        self.ms2_curve_red = pg.PlotCurveItem([], [], connect="finite",
                                              pen=pg.mkPen(220, 70, 70, 255, width=3))
        self.ms2_plot.addItem(self.ms2_curve_red)
        self.ms2_plot.addItem(self.ms2_curve)
        self.ms2_plot.scene().sigMouseClicked.connect(self._on_ms2_strip_clicked)
        self._ms2_all = []      # all loaded MS2 scans; filtered to the view below

        # Panel 2: m/z on x (aligned with panel 1), RT on its own real left axis
        # (fixed width) sitting to the right of the band strip.
        self.p2 = pg.PlotWidget()
        self.p2_load_overlay = self._make_loading_overlay(self.p2)
        self.p2.setLabel("bottom", "m/z")
        self.p2.setLabel("left", "RT", units="min")
        self.p2.getAxis("bottom").enableAutoSIPrefix(False)   # no "x0.000#" on m/z
        self.p2.getAxis("left").setWidth(P2_AXIS_W)
        self.p2_image = pg.ImageItem()
        if self._cmap is not None:
            self.p2_image.setColorMap(self._cmap)
        self.p2.addItem(self.p2_image)

        # Panel 1's 2D plot is wrapped with a strip-width spacer (see
        # build_panel1_dock) and its left axis is P2_AXIS_W -- IDENTICAL to panel
        # 2's strip(spacer) + RT axis -- so both plot areas start and span exactly
        # the same screen x and the m/z axes line up to the pixel.
        self.p1_2d.getAxis("left").setWidth(P2_AXIS_W)

        # No in-plot auto-range "A" buttons (they overlapped the data).
        self.p2.getPlotItem().hideButtons()
        self.ms2_plot.getPlotItem().hideButtons()

        # Remove the GraphicsView frames so the plot areas have NO border inset.
        # Panel 2's plot sits behind two framed widgets (strip + p2) vs panel 1's
        # one, so leftover 1px frames are exactly the kind of asymmetry that
        # nudges the m/z axes out of alignment. Zero them everywhere.
        for _w in (self.p1_2d, self.p2, self.ms2_plot):
            _w.setFrameStyle(0)
            _w.setContentsMargins(0, 0, 0, 0)
            try:
                _w.getPlotItem().getViewBox().setDefaultPadding(0.0)
            except Exception:
                pass

        # Band on panel 2 showing a hovered MS2 scan's isolation window: a thick
        # horizontal line at the scan's RT spanning the exact m/z range that was
        # isolated, same colour as the MS2 strip lines at ~50% opacity (3 px, so
        # it reads as a band the same width as the left strip's RT bands).
        self.p2_ms2_band = pg.PlotCurveItem(
            pen=pg.mkPen(255, 170, 50, 128, width=3))
        self.p2_ms2_band.setZValue(30)
        self.p2.addItem(self.p2_ms2_band)
        # Hover over the MS2 strip -> show that scan's isolation band on panel 2.
        self.ms2_plot.scene().sigMouseMoved.connect(self._on_ms2_hover)

        # MS2 strip shares panel 2's RT (y) axis.
        self.ms2_plot.setYLink(self.p2)
        # panel 1 (2D) and panel 2 share the m/z (x) axis -> link them.
        self.p1_2d.setXLink(self.p2)
        self.p2.sigXRangeChanged.connect(self.on_view_range_changed)
        self.p2.sigYRangeChanged.connect(self.on_view_range_changed)
        # Clicking empty space in panel 2 deselects the distribution (and clears
        # the theoretical overlay). A click ON a distribution's dots is handled
        # by the scatter's sigClicked (on_panel2_dist_clicked); the two are
        # disambiguated in _on_p2_clicked via a deferred check.
        self.p2.scene().sigMouseClicked.connect(self._on_p2_clicked)
        # Plain wheel over panel 2 pans the window (m/z / time); Ctrl+wheel zooms
        # (handled in eventFilter).
        self.p2.viewport().installEventFilter(self)

        self.p2_loading = QLabel("")
        self.p2_loading.setStyleSheet("color: black;")

        row_widget = QWidget()
        row = QHBoxLayout(row_widget)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        row.addWidget(self.ms2_plot)
        row.addWidget(self.p2, stretch=1)

        container = QWidget()
        container.setObjectName("panel2_frame")
        layout = QVBoxLayout(container)
        layout.setContentsMargins(2, 0, 2, 2)
        # p2_loading is intentionally NOT added to the layout: its row was the
        # dead-space border between panel 1 and panel 2. The loading state still
        # shows in panel 1's bar (and panel 3).
        layout.addWidget(row_widget, stretch=1)

        dock = QDockWidget("", self)
        dock.setObjectName("dock_panel2")
        dock.setTitleBarWidget(QWidget())   # drop the empty title bar (dead space)
        dock.setWidget(container)
        self.dock_panel2 = dock

    def build_panel3_dock(self):
        # Single plot (isotope overlay / MS2 spectrum) and the multi-distribution
        # charge-comparison grid live in a stack; panel 3 shows whichever fits.
        self.p3 = pg.PlotWidget()
        self.p3.setLabel("bottom", "m/z")
        self.p3.setLabel("left", "intensity")
        # Panel 3 zooms like panel 1: dragging/scrolling inside the plot moves the
        # m/z (x) axis only; wheel over the y-axis strip scrolls intensity with the
        # baseline pinned at 0 (handled in eventFilter); double-click resets.
        self.p3.setMouseEnabled(x=True, y=False)
        self.p3.getAxis("left").setWidth(60)
        self.p3.viewport().installEventFilter(self)
        self.p3.scene().sigMouseClicked.connect(self._on_p3_clicked)
        self.p3.getPlotItem().hideButtons()
        self.p3_grid = pg.GraphicsLayoutWidget()
        self.p3_grid.scene().sigMouseClicked.connect(self._on_grid_clicked)
        # Wheel over a grid cell's y-axis zooms that cell's intensity (handled in
        # eventFilter); wheel over the plot area zooms m/z (cells are x-only).
        self.p3_grid.viewport().installEventFilter(self)
        self.p3_stack = QStackedWidget()
        self.p3_stack.addWidget(self.p3)        # 0 = single plot
        self.p3_stack.addWidget(self.p3_grid)   # 1 = charge grid

        # p3_title still exists (error/status sink for setText calls) but is NOT
        # shown -- the user wants no descriptive caption above panel 3.
        self.p3_title = QLabel("")
        self.p3_title.setVisible(False)
        # p3_loading is intentionally NOT added to the layout: its row was the
        # dead-space border above panel 3. The loading state still updates this
        # label (harmless) but it no longer occupies a row, so the plot fills the
        # space.
        self.p3_loading = QLabel("")
        self.p3_loading.setStyleSheet("color: black;")

        # MS1 / MS2 selector for the current distribution: MS1 = isotope
        # envelope / charge grid, MS2 = the identification's fragment spectrum.
        # Reflects _panel3_mode (kept in sync by _sync_panel3_tab) and drives it
        # on a click (_on_panel3_tab_changed).
        self.p3_tabs = QTabBar()
        self.p3_tabs.addTab("MS1")   # index 0
        self.p3_tabs.addTab("MS2")   # index 1
        self.p3_tabs.setExpanding(False)
        self.p3_tabs.setDrawBase(False)
        self._syncing_p3_tabs = False
        self.p3_tabs.currentChanged.connect(self._on_panel3_tab_changed)

        # Dropdown of every MS2 scan taken on the current MS1 distribution (like
        # the file combo). Only shown in MS2 mode; picking one renders it and
        # highlights it on panel 2. Populated by _populate_ms2_combo.
        self.p3_ms2_combo = QComboBox()
        self.p3_ms2_combo.setVisible(False)
        self._syncing_ms2_combo = False
        self.p3_ms2_combo.currentIndexChanged.connect(self._on_ms2_combo_changed)

        container = QWidget()
        container.setObjectName("panel3_frame")
        layout = QVBoxLayout(container)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)
        layout.addWidget(self.p3_tabs)
        layout.addWidget(self.p3_ms2_combo)
        layout.addWidget(self.p3_stack, stretch=1)

        dock = QDockWidget("", self)
        dock.setObjectName("dock_panel3")
        dock.setTitleBarWidget(QWidget())   # drop the empty title bar (dead space)
        dock.setWidget(container)
        self.dock_panel3 = dock

    def _make_table_view(self, columns, on_double_click):
        """A sortable QTableView over a SimpleTableModel (via a sort proxy)."""
        view = QTableView()
        model = SimpleTableModel(columns, parent=view)
        proxy = QSortFilterProxyModel(view)
        proxy.setSortRole(Qt.EditRole)
        proxy.setSourceModel(model)
        view.setModel(proxy)
        view.setEditTriggers(QAbstractItemView.NoEditTriggers)
        view.setSelectionBehavior(QAbstractItemView.SelectRows)
        view.verticalHeader().setVisible(False)
        view.setSortingEnabled(True)
        view.setSelectionBehavior(QAbstractItemView.SelectRows)
        view.doubleClicked.connect(
            lambda proxy_index, m=model, p=proxy, cb=on_double_click:
            cb(m.row_dict(p.mapToSource(proxy_index).row()))
        )
        self._install_row_copy(view)
        return view, model

    def _install_row_copy(self, view):
        """Ctrl+C copies the whole selected row(s) (tab-separated)."""
        from PySide6.QtGui import QShortcut, QKeySequence
        sc = QShortcut(QKeySequence.Copy, view)
        sc.activated.connect(lambda v=view: self._copy_rows(v))

    def _copy_rows(self, view):
        from PySide6.QtWidgets import QApplication
        sel = view.selectionModel()
        model = view.model()
        if sel is None or model is None:
            return
        rows = sorted({i.row() for i in sel.selectedIndexes()})
        lines = []
        for r in rows:
            cells = [str(model.index(r, c).data() or "") for c in range(model.columnCount())]
            lines.append("\t".join(cells))
        if lines:
            QApplication.clipboard().setText("\n".join(lines))

    def build_table1_dock(self):
        # 'current' tab: the selected distribution's lines.
        self.table1 = QTableWidget()
        self.table1.setColumnCount(len(LINE_METRIC_COLUMNS))
        self.table1.setHorizontalHeaderLabels([h for _, h, _ in LINE_METRIC_COLUMNS])
        self.table1.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table1.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table1.verticalHeader().setVisible(False)
        self._install_row_copy(self.table1)

        # 'distributions' tab.
        self.dists_view, self.dists_model = self._make_table_view(
            DIST_TAB_COLUMNS, self._on_distribution_activated)
        self.dists_view.clicked.connect(self._on_dist_clicked)

        # 'charge distributions' tab: flat multi-charge analytes (single rows).
        self.charge_view, self.charge_model = self._make_table_view(
            CHARGE_TAB_COLUMNS, self._on_charge_activated)
        self.charge_view.clicked.connect(self._on_charge_clicked)

        self.table1_tabs = QTabWidget()
        self.table1_tabs.addTab(self.table1, "current")
        self.table1_tabs.addTab(self.dists_view, "distributions")
        self.table1_tabs.addTab(self.charge_view, "charge distributions")
        self.table1_tabs.currentChanged.connect(self._on_table1_tab_changed)
        self._table1_loaded = set()

        # 'All' button on the tab row (reloads whichever tab is active).
        all_btn = QPushButton("All")
        all_btn.setFixedWidth(44)
        all_btn.clicked.connect(self._reload_current_tab)
        self.table1_tabs.setCornerWidget(all_btn, Qt.TopRightCorner)
        # 'loading…' shows as a placeholder row WITHIN each table (models /
        # QTableWidget), not as a tab-corner widget (which pushed the tabs over).
        self.table1_loading = None

        dock = QDockWidget("", self)
        dock.setObjectName("dock_table1")
        dock.setTitleBarWidget(QWidget())   # reclaim the title-bar head space
        dock.setWidget(self.table1_tabs)
        self.dock_table1 = dock

    def build_table2_dock(self):
        # Shown for MS2: candidate PSMs for the sampled precursor + the b/y
        # fragment sequence coverage each gets against the MS2 spectrum (the
        # coverage concept from sequencecoverageconcept.py). Lives under panel 3.
        self.table2 = QTableWidget()
        # Table 2 compares the candidate PEPTIDES for the sampled precursor (not
        # proteins): each peptide's q-value and its coverage score (the b/y
        # fragment coverage metric from sequencecoverageconcept.py).
        cols = ["peptide", "q", "coverage"]
        self.table2.setColumnCount(len(cols))
        self.table2.setHorizontalHeaderLabels(cols)
        self.table2.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table2.verticalHeader().setVisible(False)   # no 1,2,3 index column
        # Click a header to sort by that column; click again to reverse.
        self.table2.setSortingEnabled(True)
        self.table2.horizontalHeader().setSortIndicatorShown(True)
        self.table2.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table2.setSelectionMode(QAbstractItemView.SingleSelection)
        # Selecting a candidate peptide overlays its fragment ions on panel 3.
        self.table2.itemSelectionChanged.connect(self._on_table2_peptide_selected)
        dock = QDockWidget("", self)
        dock.setObjectName("dock_table2")
        dock.setTitleBarWidget(QWidget())   # drop the title bar (dead space)
        dock.setWidget(self.table2)
        self.dock_table2 = dock

    def arrange_default(self):
        # Three columns: lists | [panel1 / panel2 / table1] | panel3
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_lists)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_panel1)
        self.splitDockWidget(self.dock_panel1, self.dock_panel3, Qt.Horizontal)
        self.splitDockWidget(self.dock_panel1, self.dock_panel2, Qt.Vertical)
        self.splitDockWidget(self.dock_panel2, self.dock_table1, Qt.Vertical)
        self.splitDockWidget(self.dock_panel3, self.dock_table2, Qt.Vertical)
        self.resizeDocks([self.dock_lists], [320], Qt.Horizontal)
        self.resizeDocks([self.dock_panel1, self.dock_panel3], [600, 460], Qt.Horizontal)
        # Give panel 2 the most vertical space (it was leaving dead space).
        self.resizeDocks([self.dock_panel1, self.dock_panel2, self.dock_table1],
                         [300, 460, 240], Qt.Vertical)

    def reset_layout(self):
        self.restoreState(self._default_state)

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)
        if getattr(self, "theme_btn", None) is not None:
            self.theme_btn.setText("Light Mode" if theme == "dark" else "Dark Mode")
        for plot in (self.p1_2d, self.p2, self.p3, self.ms2_plot):
            style_plot(plot, pal)
        try:
            self.p3_grid.setBackground(pal["bg"])
        except Exception:
            pass
        # Force the charge grid to rebuild on a theme change (its axis pens, tick
        # text and titles are coloured at build time, and the rebuild is otherwise
        # skipped when the same distribution is showing).
        self._grid_dist_id = None
        # Fill the panel-1 2D left spacer with the same plot background so it
        # blends in (instead of an awkward blank gap next to panel 2's strip).
        if getattr(self, "_p1_2d_spacer", None) is not None:
            self._p1_2d_spacer.setStyleSheet(f"background-color: {pal['bg']};")
        # Outline each panel (1/2/3) so their boundaries are visible -- black on
        # the light theme, light grey on the dark theme. Scoped by objectName so
        # the border doesn't cascade onto the child plot widgets.
        border = "#000000" if theme == "light" else "#9a9a9a"
        for name in ("panel1_frame", "panel2_frame", "panel3_frame"):
            w = self.findChild(QWidget, name)
            if w is not None:
                w.setStyleSheet(f"QWidget#{name} {{ border: 1px solid {border}; }}")
        # Keep the loading indicators legible against the tab-bar background.
        for lbl in (getattr(self, "tab_loading_label", None),
                    getattr(self, "table1_loading", None)):
            if lbl is not None:
                lbl.setStyleSheet(f"color: {pal['fg']}; padding: 0 8px;")
        self._recolor_gl(pal)   # 3D side labels (exist with or without GL)
        if HAVE_GL:
            style_gl(self.p1_3d, pal)
        # Re-render the data-bearing panels so theme-coloured bits (panel 1/3
        # titles, the MS2 spectrum data, the charge-grid axes) follow the theme.
        # The panel-3 rebuild routes through populate_table2, which clears the
        # theoretical-MS1 selection -- save/restore it across the redraw.
        saved_theo = self._ms1_theo
        if self.current is not None:
            if isinstance(getattr(self, "_last_points", None), dict):
                self._redraw_from_cache()
            self._redraw_panel3()   # rebuilds the charge grid with the new theme
        self._ms1_theo = saved_theo
        self._draw_ms1_theo_overlay()   # recolor the theoretical MS1 overlay
        # The panel-3 redraw wiped the MS2 fragment green/red annotation; restore
        # it from the last result (no need to re-run the worker).
        if (self._panel3_mode == "ms2"
                and getattr(self, "_last_frag", None) is not None):
            self._draw_fragment_overlay(*self._last_frag)

    def _redraw_panel3(self):
        """Re-draw whatever panel 3 is currently showing (MS2 spectrum / MS1
        envelope / charge grid) so its theme-dependent colours update."""
        if self.current is None:
            return
        if self._panel3_mode == "ms2" and self._ms2_scan is not None:
            self.render_ms2(self._ms2_scan)
        else:
            self.draw_panel3_ms1(self.current,
                                 getattr(self, "_last_scan_mz", None),
                                 getattr(self, "_last_scan_int", None))

    def _ms2_scan_for_current(self):
        """The MS2 scan for the current identification: the PSM's own scan if
        it's in the loaded window, else the nearest MS2 by RT. Mirrors the
        lookup in _apply_ms2_focus. Returns None when no MS2 is loaded."""
        ms2 = getattr(self, "_ms2_all", None)
        if not ms2 or not isinstance(self.current, dict):
            return None
        scan_no = str(self.current.get("scan", "") or "")
        for m in ms2:
            if str(m.get("number", "")) == scan_no:
                return m
        rt = self.current.get("rt")
        if rt is None:
            return None
        return min(ms2, key=lambda m: abs((m.get("rt") or 0.0) - rt))

    def _ms2_scans_on_distribution(self):
        """Every loaded MS2 scan taken ON the selected MS1 distribution: its
        precursor isolation window overlaps the distribution's m/z span and its
        RT falls inside the distribution's RT band. Sorted by RT. Falls back to
        the whole loaded window's RT if no distribution bbox is known."""
        ms2 = getattr(self, "_ms2_all", None)
        if not ms2:
            return []
        bbox = getattr(self, "_selected_bbox", None)
        if bbox is not None:
            mz0, mz1, rt0, rt1 = bbox
        elif self.window is not None:
            mz0, mz1, rt0, rt1 = self.window
        else:
            return []
        out = []
        for m in ms2:
            rt = m.get("rt")
            if rt is None or not (rt0 <= rt <= rt1):
                continue
            lo, hi = m.get("iso_low"), m.get("iso_high")
            pmz = m.get("mz")
            if lo is not None and hi is not None:
                if hi < mz0 or lo > mz1:      # isolation window misses the dist
                    continue
            elif pmz is not None and not (mz0 <= pmz <= mz1):
                continue
            out.append(m)
        return sorted(out, key=lambda m: (m.get("rt") or 0.0))

    def _default_ms2_scan(self):
        """The MS2 scan to show first when entering the MS2 view: the
        identification's own scan if it's on the distribution, else the first
        scan on the distribution, else the nearest-by-RT fallback."""
        scans = self._ms2_scans_on_distribution()
        if isinstance(self.current, dict):
            scan_no = str(self.current.get("scan", "") or "")
            for m in scans:
                if str(m.get("number", "")) == scan_no:
                    return m
        if scans:
            return scans[0]
        return self._ms2_scan_for_current()

    def _populate_ms2_combo(self, current):
        """Fill the MS2-scan dropdown with the scans on the current distribution
        and select ``current`` (a scan dict). No-op-safe; guarded so filling it
        doesn't fire the change handler."""
        combo = getattr(self, "p3_ms2_combo", None)
        if combo is None:
            return
        scans = self._ms2_scans_on_distribution()
        if current is not None and not any(
                str(m.get("number", "")) == str(current.get("number", ""))
                for m in scans):
            scans = [current] + scans   # always include the shown scan
        self._syncing_ms2_combo = True
        combo.clear()
        sel = 0
        for i, m in enumerate(scans):
            rt = m.get("rt")
            combo.addItem(
                f"scan {m.get('number', '?')}   rt={rt:.2f}   m/z={m.get('mz')}"
                if rt is not None else f"scan {m.get('number', '?')}", m)
            if current is not None and str(m.get("number", "")) == str(
                    current.get("number", "")):
                sel = i
        if combo.count():
            combo.setCurrentIndex(sel)
        self._syncing_ms2_combo = False

    def _on_ms2_combo_changed(self, index):
        """User picked a different MS2 scan from the dropdown: render it (which
        also moves the panel-2 isolation band to that scan)."""
        if getattr(self, "_syncing_ms2_combo", False) or index < 0:
            return
        m = self.p3_ms2_combo.itemData(index)
        if m is not None:
            self.render_ms2(m)

    def _sync_panel3_tab(self):
        """Reflect _panel3_mode in the MS1/MS2 tab bar without re-triggering the
        change handler, and show the MS2-scan dropdown only in MS2 mode. MS2 is
        disabled when no MS2 scan is available for the distribution."""
        tabs = getattr(self, "p3_tabs", None)
        if tabs is None:
            return
        has_ms2 = (self._ms2_scan is not None
                   or bool(self._ms2_scans_on_distribution())
                   or self._ms2_scan_for_current() is not None)
        self._syncing_p3_tabs = True
        tabs.setTabEnabled(1, has_ms2)
        tabs.setCurrentIndex(1 if self._panel3_mode == "ms2" else 0)
        self._syncing_p3_tabs = False
        combo = getattr(self, "p3_ms2_combo", None)
        if combo is not None:
            combo.setVisible(self._panel3_mode == "ms2" and combo.count() > 0)

    def _on_panel3_tab_changed(self, index):
        """User clicked the MS1 / MS2 tab: switch panel 3's view for the current
        distribution."""
        if self._syncing_p3_tabs or self.current is None:
            return
        if index == 1:   # MS2
            scan = self._ms2_scan or self._default_ms2_scan()
            if scan is None:
                self._sync_panel3_tab()   # nothing to show -> snap back to MS1
                return
            self.render_ms2(scan)
        else:            # MS1
            self._panel3_mode = "ms1"
            self.draw_panel3_ms1(self.current,
                                 getattr(self, "_last_scan_mz", None),
                                 getattr(self, "_last_scan_int", None))
            self._sync_panel3_tab()

    # ---- list population + cross-linking ---------------------------------

    def _fill(self, listw, entries, preserve=False):
        # Optionally keep the current selection (by text) and scroll position,
        # so clicking "All" doesn't jump the list or lose the highlight.
        sel = listw.currentItem().text() if (preserve and listw.currentItem()) else None
        scroll = listw.verticalScrollBar().value() if preserve else 0

        listw.blockSignals(True)
        listw.clear()
        restore_row = -1
        for i, (text, data) in enumerate(entries):
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, data)
            listw.addItem(item)
            if sel is not None and text == sel:
                restore_row = i
        if restore_row >= 0:
            listw.setCurrentRow(restore_row)
        listw.blockSignals(False)

        if preserve:
            listw.verticalScrollBar().setValue(scroll)

    def _filter(self, text):
        return True  # filter box removed; lists show everything for the file

    # All list content is scoped to the selected file (single-file view).
    def file_psms(self):
        rows = []
        for row in self.session.load_psms(self.current_file or ""):
            row = dict(row)
            row["filename"] = self.current_file
            rows.append(row)
        return rows

    def _peptide_label(self, row):
        n = str(row.get("n_psms", "") or "").strip()
        pep = row.get("peptide", "")
        return f"{pep}   ({n})" if n not in ("", "0") else pep

    def _identified_plains(self):
        """Plain sequences of peptides actually identified (have a PSM) in the
        current file -- used to keep the lists to file-identified entries only
        (LFQ-only / transferred peptides are excluded)."""
        return {plain_seq(r.get("peptide", "")) for r in self.file_psms()}

    def identified_peptides(self):
        ident = self._identified_plains()
        return [r for r in self.session.file_peptides(self.current_file or "")
                if plain_seq(r.get("peptide", "")) in ident]

    def identified_proteins(self):
        ident = self._identified_plains()
        out = []
        for r in self.session.file_proteins(self.current_file or ""):
            peps = [p for p in str(r.get("peptides", "")).split(";") if p]
            if any(plain_seq(p) in ident for p in peps):
                out.append(r)
        return out

    def show_all_proteins(self, preserve=True):
        # Only proteins with a peptide identified in THIS file (no LFQ-only).
        rows = self.identified_proteins()
        self._fill(self.protein_list, [(r.get("protein_id", ""), r) for r in rows
                                       if r.get("protein_id") and self._filter(r["protein_id"])],
                   preserve=preserve)

    def show_all_peptides(self, preserve=True):
        # Only peptides identified (with a PSM) in THIS file (no LFQ-only).
        rows = self.identified_peptides()
        self._fill(self.peptide_list, [(self._peptide_label(r), r) for r in rows
                                       if r.get("peptide") and self._filter(r["peptide"])],
                   preserve=preserve)

    def show_all_psms(self, preserve=True):
        self.psm_rows = [r for r in self.file_psms() if self._filter(r.get("peptide", ""))]
        self._fill(self.psm_list, [(f"{r.get('scan','')}  {r.get('peptide','')}", r) for r in self.psm_rows],
                   preserve=preserve)

    def repopulate_active_list(self):
        self.show_all_proteins(preserve=False)
        self.show_all_peptides(preserve=False)
        self.show_all_psms(preserve=False)

    def on_protein_selected(self):
        items = self.protein_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        peptides = set(p for p in str(row.get("peptides", "")).split(";") if p)
        plains = {plain_seq(p) for p in peptides}
        # show this protein's peptides, but only those identified in this file.
        matched = [r for r in self.identified_peptides()
                   if r.get("peptide") in peptides or plain_seq(r.get("peptide", "")) in plains]
        self._fill(self.peptide_list, [(r.get("peptide", ""), r) for r in matched])
        # If the protein resolves to a single peptide, auto-select it (which in
        # turn may auto-load its single PSM).
        if len(matched) == 1:
            self.peptide_list.setCurrentRow(0)

    def on_peptide_selected(self):
        items = self.peptide_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        plain = plain_seq(row.get("peptide", ""))
        # cross-link proteins of this peptide
        proteins = [p for p in str(row.get("proteins", "")).split(";") if p]
        if proteins:
            self._fill(self.protein_list, [(p, {"protein_id": p, "peptides": row.get("peptide", "")}) for p in proteins])
        # scope the PSM list to this peptide (within the file); don't auto-load
        # the (potentially huge) mzML until a PSM is explicitly chosen.
        self.psm_rows = [r for r in self.file_psms() if plain_seq(r.get("peptide", "")) == plain]
        self._fill(self.psm_list, [(f"{r.get('scan','')}  {r.get('peptide','')}", r) for r in self.psm_rows])
        # If there's exactly one PSM, auto-load it instead of making the user pick.
        if len(self.psm_rows) == 1:
            self.psm_list.setCurrentRow(0)

    def on_psm_selected(self):
        items = self.psm_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        if "scan" not in row:
            return
        try:
            self.update_evidence(row)
        except Exception as exc:
            import traceback
            self.p3_title.setText(f"evidence error: {exc}")
            traceback.print_exc()

    def focus_identification(self, filename, protein_id, peptide_plain, show_ms2=True):
        """Drive the lists to a specific identification, as if the user had
        clicked file → protein → peptide → PSM. Used by the Proteins tab when a
        peptide is double-clicked.

        With ``show_ms2`` the full MS2 evidence view is reconstructed once the
        (async) window load finishes: panel 3 shows the PSM's MS2 spectrum, panel
        1 gets the peptide's theoretical MS1 overlay, and panel 2 selects the
        matching distribution and draws the MS2 isolation band (see
        ``_apply_ms2_focus``, run from ``on_evidence_done``)."""
        if filename:
            idx = self.file_combo.findData(filename)
            if idx >= 0 and idx != self.file_combo.currentIndex():
                self.file_combo.setCurrentIndex(idx)   # triggers repopulate
            else:
                self.repopulate_active_list()
        # protein
        if protein_id:
            for i in range(self.protein_list.count()):
                if self.protein_list.item(i).text() == protein_id:
                    self.protein_list.setCurrentRow(i)   # populates peptide list
                    break
        # peptide (match on the plain, modification-stripped sequence)
        if peptide_plain:
            for i in range(self.peptide_list.count()):
                row = self.peptide_list.item(i).data(Qt.UserRole) or {}
                if plain_seq(row.get("peptide", "")) == peptide_plain:
                    self.peptide_list.setCurrentRow(i)
                    self.peptide_list.scrollToItem(self.peptide_list.item(i))
                    break

        # Pick the PSM to load: the best-scoring one for this peptide (lowest q).
        target_row, target_i = None, -1
        for i in range(self.psm_list.count()):
            row = self.psm_list.item(i).data(Qt.UserRole) or {}
            if plain_seq(row.get("peptide", "")) != peptide_plain:
                continue
            q = safe_float(row.get("percolator_q"))
            if target_row is None or (q is not None and
                                      (safe_float(target_row.get("percolator_q")) is None
                                       or q < safe_float(target_row.get("percolator_q")))):
                target_row, target_i = row, i
        if target_i < 0 and self.psm_list.count():
            target_row, target_i = self.psm_list.item(0).data(Qt.UserRole), 0

        # Arm the MS2 reconstruction BEFORE loading so on_evidence_done applies it
        # when the window finishes extracting.
        if show_ms2 and target_row is not None:
            self._pending_focus = {
                "scan": str(target_row.get("scan", "") or ""),
                "seq": peptide_plain,
                "charge": peptide_charge(target_row),
                "row": target_row,
            }
        else:
            self._pending_focus = None

        if target_i >= 0:
            if self.psm_list.currentRow() == target_i:
                self.on_psm_selected()   # re-load even if already selected
            else:
                self.psm_list.setCurrentRow(target_i)   # triggers on_psm_selected

    def _apply_ms2_focus(self, pf):
        """Reconstruct the MS2 evidence view for a Proteins-tab jump: render the
        PSM's MS2 spectrum, overlay the peptide's theoretical MS1 on panel 1, and
        select the matching distribution + MS2 band on panel 2. Runs after the
        window load so ``self._ms2_all`` is populated."""
        scan_no = pf.get("scan")
        scan = None
        for m in getattr(self, "_ms2_all", []):
            if str(m.get("number", "")) == scan_no:
                scan = m
                break
        if scan is None:
            # nearest MS2 by RT as a fallback so something relevant shows
            if getattr(self, "_ms2_all", None) and isinstance(self.current, dict):
                rt = self.current.get("rt")
                if rt is not None:
                    scan = min(self._ms2_all,
                               key=lambda m: abs((m.get("rt") or 0.0) - rt))

        # The theoretical MS1 overlay on panel 1 needs only the peptide's
        # sequence + charge, so draw it FIRST -- independent of whether an MS2
        # scan can be located. A validated ID whose MS2 scan isn't in the loaded
        # window (or whose stored scan number doesn't match) must still show its
        # MS1 isotope distribution. The MS2-dependent reconstruction below is
        # gated on ``scan`` separately.
        seq = pf.get("seq") or ""
        charge = pf.get("charge") or (scan.get("charge") if scan else None) or 2
        if len(seq) >= 2:
            row = pf.get("row") or {}
            # Select the matching distribution FIRST so the theoretical overlay
            # normalizes to it, then set the overlay and force a panel-1 redraw
            # (draw_panel1 re-adds the overlay) so it reliably appears. With no
            # MS2 scan the distribution is matched off the PSM's own RT.
            try:
                self._select_distribution_for_candidate(row, charge, scan)
            except Exception:
                pass
            self._ms1_theo = {"seq": seq, "charge": charge,
                              "mod_mass": peptide_mod_mass(row.get("peptide", ""))}
            self._redraw_panel1_view()
            self._draw_ms1_theo_overlay()

        if scan is None:
            return
        self.render_ms2(scan)

    # ---- store access ----------------------------------------------------

    def centroid_store(self, filename):
        path = self.session.centroid_path(filename)
        if path is None:
            return None
        key = str(path)
        if key not in self._centroid:
            self._centroid[key] = MzmlStore(path)
        return self._centroid[key]

    def profile_store(self, filename):
        path = self.session.profile_path(filename)
        if path is None:
            return None
        key = str(path)
        if key not in self._profile:
            self._profile[key] = MzmlStore(path)
        return self._profile[key]

    def points_store(self, filename):
        """Option B: a store that serves raw centroids from the file's
        distributions sqlite (scan_points), avoiding mzML re-decode on zoom.
        Returns None when the sqlite has no scan_points (falls back to mzML)."""
        sqlite_path = self.session.distributions_db_for(filename) if self.session else None
        if sqlite_path is None:
            return None
        key = str(sqlite_path)
        cache = getattr(self, "_points", None)
        if cache is None:
            cache = self._points = {}
        if key not in cache:
            cache[key] = PointsStore(key) if PointsStore.has_points(key) else None
        return cache[key]

    # ---- the selected match -> panels ------------------------------------

    def update_evidence(self, row):
        filename = row.get("filename", "")
        charge = peptide_charge(row)
        neutral_mass = peptide_mass(row)
        rt = peptide_rt(row)
        targets = isotope_mzs(neutral_mass, charge, n=6) if (neutral_mass and charge) else []
        mz_center = sum(targets) / len(targets) if targets else (neutral_mass or 500.0)

        self.current = {
            "row": row, "filename": filename, "scan": row.get("scan", ""),
            "charge": charge, "neutral_mass": neutral_mass, "rt": rt,
            "targets": targets, "mz_center": mz_center,
            "centroid": self.centroid_store(filename), "profile": self.profile_store(filename),
            "points": self.points_store(filename),
        }
        self.center = (mz_center, rt)
        self.assumed_charge = charge or 1
        # A fresh selection returns panel 3 to its MS1 view (the MS2 spectrum is
        # only shown after the user clicks an MS2 point).
        self._panel3_mode = "ms1"
        self._ms2_scan = None
        self._ms2_band = None   # drop the isolation band
        self._apply_ms2_band()
        # Overlay THIS peptide's theoretical MS1 isotope distribution on panel 1
        # directly from the list selection (not only via a Table-2 / Proteins-tab
        # jump). draw_panel1 re-adds it after the async extract, and
        # _on_table1_loaded re-normalizes it once the distribution is known.
        plain = plain_seq(row.get("peptide", ""))
        self._ms1_theo = ({"seq": plain, "charge": charge or 1,
                           "mod_mass": peptide_mod_mass(row.get("peptide", ""))}
                          if len(plain) >= 2 else None)
        self.render_table1(self.current)
        # Initialize the window from the ± controls and snap the views to it.
        rt_start = max(0.0, rt - self.rt_half) if rt is not None else 0.0
        rt_end = rt + self.rt_half if rt is not None else 1.0
        self.set_window([mz_center - self.mz_half, mz_center + self.mz_half, rt_start, rt_end],
                        set_view=True)

    def use_profile(self):
        return self.source_toggle.key() == "profile"

    def noise_mode(self):
        """Current level: 0 none, 1 recovered distributions, 2 lines, 3 +small
        lines, 4 +single points (cumulative)."""
        return self.noise_toggle.index()

    def recovered_visible(self):
        """Recovered (less-confident) distributions show from level 1 onward."""
        return self.noise_mode() >= 1

    def noise_visible_mask(self, n):
        """Boolean over the loaded points: which raw noise points the current
        level reveals. Noise classes start at level 2 (level 1 is recovered
        distributions), so class c (1 line, 2 small, 3 single) shows at mode>=c+1."""
        mode = self.noise_mode()
        cls = getattr(self, "_noise_class", None)
        if mode <= 1 or cls is None or cls.size != n:
            return np.zeros(n, dtype=bool)
        return (cls >= 1) & (cls <= mode - 1)

    def _on_noise_changed(self):
        # Pure redraw from cached data -- no re-extraction needed.
        self._redraw_from_cache()

    def use_logcolor(self):
        return self.logcolor_toggle.key() == "log"

    def _on_logcolor(self):
        # Recolour the cached 3D points without re-extraction.
        if HAVE_GL and getattr(self, "_p1_3d_inputs", None) is not None:
            self.draw_panel1_3d(*self._p1_3d_inputs)

    def _redraw_from_cache(self):
        """Redraw panels 1 & 2 from the last extracted points (e.g. after the
        noise toggle) without hitting the mzML again."""
        if not isinstance(getattr(self, "_last_points", None), dict) or self._win is None:
            return
        mz_min, mz_max, rt_start, rt_end = self._win
        self.draw_panel1(self._last_points, self._last_region, mz_min, mz_max, rt_start, rt_end)
        self.draw_panel2(self._last_points, (mz_min, mz_max, rt_start, rt_end))

    # ---- window-driven extraction ----------------------------------------

    def set_window(self, window, set_view=False):
        """Set the m/z x RT window (source of truth) and reload it.

        ``set_view`` snaps the panel views to the window (used on selection /
        ± changes). User drag/zoom calls this with set_view=False so the view the
        user produced is what gets reloaded.
        """
        self.window = [float(window[0]), float(window[1]), float(window[2]), float(window[3])]
        if set_view:
            self._guard = True
            try:
                self.p2.setXRange(self.window[0], self.window[1], padding=0)
                self.p2.setYRange(self.window[2], self.window[3], padding=0)
            finally:
                self._guard = False
        self.refresh()

    def on_view_range_changed(self, *_):
        if self._guard or self.window is None:
            return
        # Read the current window from panel 2 (m/z = x, RT = y).
        (mz0, mz1) = self.p2.getViewBox().viewRange()[0]
        (rt0, rt1) = self.p2.getViewBox().viewRange()[1]
        self.window = [mz0, mz1, max(0.0, rt0), rt1]
        # Inverse-scale the point sizes so dots stay visible as you zoom in.
        self._rescale_points()
        # Keep the MS2 strip in sync with what's visible in panel 2.
        self._refresh_ms2_visible()
        # Re-filter panel 1 (and the 3D) to the new visible window (debounced).
        self._p1_view_timer.start()
        # The big "data disappears on zoom" fix: we extract a PADDED region and
        # cache it. Zooming/panning *within* that cached region is a pure view
        # operation -- no re-extraction (which is what used to come back empty
        # when the RT view got narrower than the MS1 scan spacing). Only reload
        # when the view leaves the cached region.
        if self._loaded_window is not None and self._within(self.window, self._loaded_window):
            return
        self._reload_timer.start()

    @staticmethod
    def _within(inner, outer, eps=1e-9):
        return (inner[0] >= outer[0] - eps and inner[1] <= outer[1] + eps and
                inner[2] >= outer[2] - eps and inner[3] <= outer[3] + eps)

    @staticmethod
    def _padded(window):
        """Grow a view window into the region we actually extract, so there is
        margin to zoom/pan into before another reload is needed."""
        mz0, mz1, rt0, rt1 = window
        # Modest margin: enough to pan/zoom a bit before reloading, without
        # extracting a huge region every time (keeps it snappy). Zoom-IN within
        # the cache never reloads regardless of this size.
        mz_pad = max((mz1 - mz0) * 0.5, 0.4)
        rt_pad = max((rt1 - rt0) * 0.6, 0.05)
        return [mz0 - mz_pad, mz1 + mz_pad, max(0.0, rt0 - rt_pad), rt1 + rt_pad]

    def _rescale_points(self):
        """Scale scatter point sizes up as the view zooms in (relative to the
        cached region span), so datapoints stay visible instead of shrinking."""
        lw = self._loaded_window
        if lw is None:
            return
        lmz = max(lw[1] - lw[0], 1e-9)
        lrt = max(lw[3] - lw[2], 1e-9)
        try:
            (vmz0, vmz1) = self.p2.getViewBox().viewRange()[0]
            (vrt0, vrt1) = self.p2.getViewBox().viewRange()[1]
        except Exception:
            return
        # Keep point sizes pinned to their base (the "starting condition" look the
        # user prefers). The previous zoom-upscaling keyed off _loaded_window,
        # which grows when a zoom-out triggers a reload, so it was not reversible:
        # zooming back in returned bigger dots than the pristine view.
        for scatter, base in self._p2_scatters:
            try:
                scatter.setSize(base)
            except Exception:
                pass
        for scatter, base in getattr(self, "_p1_scatters", []):
            try:
                scatter.setSize(base)
            except Exception:
                pass

    def _redraw_panel1_view(self):
        """Re-draw panel 1 (filtered to panel 2's visible window) from cached
        points -- keeps panel 1 + the 3D matched to panel 2 on zoom/pan."""
        if isinstance(getattr(self, "_last_points", None), dict) and self._win is not None:
            self.draw_panel1(self._last_points, self._last_region, *self._win)

    def refresh(self):
        self._reload_timer.start()

    def do_extract(self):
        cur = self.current
        if cur is None or self.window is None:
            return
        centroid = cur["centroid"]
        if centroid is None:
            self.p3_title.setText(f"no centroid mzML for {cur['filename']}")
            return
        self._record_nav()
        # Extract a padded region around the view and cache its extent so zoom/pan
        # within it never re-extracts (and never comes back empty).
        load = self._padded(self.window)
        self._loaded_window = load
        mz_min, mz_max, rt_start, rt_end = load
        if self.use_profile() and cur["profile"]:
            store = cur["profile"]
        else:
            # Option B: read centroids from the sqlite points store when present,
            # else fall back to decoding the centroid mzML.
            store = cur.get("points") or centroid
        self._win = (mz_min, mz_max, rt_start, rt_end)
        self._pending = dict(
            centroid=centroid, store=store, scan=cur["scan"], rt=cur["rt"],
            mz_min=mz_min, mz_max=mz_max, rt_start=rt_start, rt_end=rt_end,
            mz_bins=400, mode="profile" if (self.use_profile() and cur["profile"]) else "centroid")
        self._set_loading(True, cur["row"].get("peptide", "") or "region")
        self._start_evidence()

    def _make_loading_overlay(self, parent):
        """A 'loading' badge floated ON TOP of a plot, hidden until shown.

        Parent it to the plot's VIEWPORT, not the PlotWidget (a QGraphicsView):
        a widget parented to the view itself renders *behind* the viewport, so
        the badge would be hidden by the plotted scene and never appear."""
        host = parent.viewport() if hasattr(parent, "viewport") else parent
        lbl = QLabel("Loading…", host)
        lbl.setStyleSheet(
            "background: rgba(0,0,0,150); color: white; padding: 3px 10px;"
            " border-radius: 4px; font-weight: bold;")
        lbl.hide()
        return lbl

    def _set_loading(self, on, context=""):
        """Show/clear a 'loading' badge on top of panels 1 & 2 while the data
        worker runs."""
        for ov in (getattr(self, "p1_load_overlay", None),
                   getattr(self, "p2_load_overlay", None)):
            if ov is None:
                continue
            if on:
                ov.adjustSize()
                pw = ov.parent().width()
                ov.move(max(4, (pw - ov.width()) // 2), 6)
                ov.show()
                ov.raise_()
            else:
                ov.hide()

    def _start_evidence(self):
        if self.worker is not None and self.worker.isRunning():
            return
        if self._pending is None:
            return
        params = self._pending
        self._pending = None
        self.worker = EvidenceWorker(**params)
        self.worker.done.connect(self.on_evidence_done)
        self.worker.start()

    def on_evidence_done(self, result):
        if self._pending is not None:
            self._start_evidence()
            return
        cur = self.current
        if cur is None or not isinstance(result, dict):
            self._set_loading(False)
            return
        if "error" in result:
            self._set_loading(False)
            self.p3_title.setText(f"evidence error: {result['error'].splitlines()[0]}")
            return

        mz_min, mz_max, rt_start, rt_end = self._win
        points = result.get("points")
        region = result.get("region")
        scan_mz = result.get("scan_mz")
        scan_int = result.get("scan_int")
        self._last_points = points
        self._last_region = region
        self._last_scan_mz = scan_mz
        self._last_scan_int = scan_int
        # Compute distribution membership once (shared by panels 1 and 2, and the
        # noise toggle). In profile mode use peak-finding so profile points are
        # linked to the centroid they'd reduce to (roadmap 1.9); otherwise match
        # centroid points directly.
        win = (mz_min, mz_max, rt_start, rt_end)
        if self.use_profile() and cur.get("profile"):
            self._assigned, self._groups = self._assignment_profile(points, win, status="primary")
            self._rec_assigned, self._rec_groups = self._assignment_profile(points, win, status="recovered")
        else:
            self._assigned, self._groups = self._assignment(points, win, status="primary")
            self._rec_assigned, self._rec_groups = self._assignment(points, win, status="recovered")
        # Noise = points in NO distribution (validated OR recovered): pass the
        # combined assignment so recovered members are not double-counted as noise.
        combined = self._assigned
        if combined is not None and self._rec_assigned is not None:
            combined = combined | self._rec_assigned
        elif self._rec_assigned is not None:
            combined = self._rec_assigned
        self._noise_class = self._compute_noise_class(points, win, combined)

        self.draw_panel1(points, region, mz_min, mz_max, rt_start, rt_end)
        self.draw_panel2(points, (mz_min, mz_max, rt_start, rt_end))
        self.draw_ms2_strip(result.get("ms2", []))
        # Keep the MS2 spectrum up if the user is in MS2 mode; otherwise (re)draw
        # the MS1 envelope / charge grid.
        if self._panel3_mode == "ms2" and getattr(self, "_pending_charge_refocus", False):
            # A charge step moved the selection to a linked distribution: shift
            # the MS2 view to that distribution's scan (now that _ms2_all is
            # reloaded for the new m/z region), not the stale previous scan.
            self._pending_charge_refocus = False
            scan = self._default_ms2_scan()
            if scan is not None:
                self.render_ms2(scan)
            else:
                # No MS2 on the new distribution -> show its MS1 view instead.
                self._panel3_mode = "ms1"
                self.draw_panel3_ms1(cur, scan_mz, scan_int)
        elif self._panel3_mode == "ms2" and self._ms2_scan is not None:
            self.render_ms2(self._ms2_scan)
        else:
            self.draw_panel3_ms1(cur, scan_mz, scan_int)
        self._set_loading(False)

        # A Proteins-tab jump armed an MS2 reconstruction: apply it now that the
        # window (and self._ms2_all) is loaded.
        pf = getattr(self, "_pending_focus", None)
        if pf is not None:
            self._pending_focus = None
            try:
                self._apply_ms2_focus(pf)
            except Exception:
                import traceback
                traceback.print_exc()

    def draw_panel1(self, points, region, mz_min, mz_max, rt_start, rt_end):
        self.p1_2d.clear()
        self._p1_scatters = []   # (ScatterPlotItem, base_size) for zoom rescaling
        if isinstance(points, dict) and points["mz"].size:
            mz = points["mz"]; rt = points["rt"]; inten = points["intensity"]
            # Panel 1 collapses RT onto m/z-vs-intensity, so it MUST be filtered to
            # panel 2's currently-visible window -- otherwise it shows points from
            # the whole loaded RT range (lots of distributions/colours) while panel
            # 2 is zoomed to one. This is what caused the mismatch.
            if self.window is not None:
                vmz0, vmz1, vrt0, vrt1 = self.window
                vm = (mz >= vmz0) & (mz <= vmz1) & (rt >= vrt0) & (rt <= vrt1)
            else:
                vm = np.ones(mz.size, dtype=bool)
            self.p1_2d.showGrid(x=True, y=True, alpha=0.25)
            # Colour each datapoint by the panel-2 distribution it belongs to
            # (one scatter per distribution); points in no distribution are grey
            # and only shown when noise is on. m/z vs intensity, every datapoint.
            shown_max = 0.0
            for _did, color, feat_masks in (self._groups or []):
                idx = np.zeros(mz.size, dtype=bool)
                for fm in feat_masks:
                    idx |= fm
                idx &= vm
                if idx.any():
                    # data=did + sigClicked makes a distribution clickable in
                    # panel 1 too, opening it in panel 3 exactly like panel 2.
                    sc = pg.ScatterPlotItem(x=mz[idx], y=inten[idx], size=3, pen=None,
                                            brush=pg.mkBrush(*color),
                                            data=[_did] * int(idx.sum()))
                    sc.sigClicked.connect(self.on_panel2_dist_clicked)
                    self.p1_2d.addItem(sc)
                    self._p1_scatters.append((sc, 3))
                    shown_max = max(shown_max, float(inten[idx].max()))
            # Recovered (less-confident) distributions: same colours, fainter, only
            # from noise level 1 onward.
            if self.recovered_visible():
                for _did, color, feat_masks in (self._rec_groups or []):
                    idx = np.zeros(mz.size, dtype=bool)
                    for fm in feat_masks:
                        idx |= fm
                    idx &= vm
                    if idx.any():
                        sc = pg.ScatterPlotItem(x=mz[idx], y=inten[idx], size=3, pen=None,
                                                brush=pg.mkBrush(*color, 140),
                                                data=[_did] * int(idx.sum()))
                        sc.sigClicked.connect(self.on_panel2_dist_clicked)
                        self.p1_2d.addItem(sc)
                        self._p1_scatters.append((sc, 3))
                        shown_max = max(shown_max, float(inten[idx].max()))
            noise_vis = self.noise_visible_mask(mz.size) & vm
            if noise_vis.any():
                # One scatter per noise class so the levels are visually distinct.
                for cls, (color, size) in NOISE_STYLE.items():
                    nm = noise_vis & (self._noise_class == cls)
                    if nm.any():
                        sc = pg.ScatterPlotItem(x=mz[nm], y=inten[nm], size=size, pen=None,
                                                brush=pg.mkBrush(*color))
                        self.p1_2d.addItem(sc)
                        self._p1_scatters.append((sc, size))
                        shown_max = max(shown_max, float(inten[nm].max()))
            if self._p1_scatters:
                self.p1_2d.getViewBox().setYRange(0, (shown_max or 1.0) * 1.05, padding=0)
            # Experimental points sit ABOVE the theoretical MS1 overlay (z=0), so
            # give every scatter a positive z; the overlay stays visually behind.
            for _sc, _base in self._p1_scatters:
                _sc.setZValue(1)
        # No title on panel 1 (per the user).
        # Rescale dot sizes for the current zoom, then build the 3D if shown.
        self._rescale_points()
        # The 3D scatter is expensive; only build it when 3D is actually shown.
        # Cache the inputs so toggling to 3D can render without a reload.
        self._p1_3d_inputs = (points, region, mz_min, mz_max, rt_start, rt_end)
        if HAVE_GL and self.dim_toggle.key() == "3D":
            self.draw_panel1_3d(points, region, mz_min, mz_max, rt_start, rt_end)
        # p1_2d.clear() above dropped any theoretical-MS1 overlay; re-add it.
        self._draw_ms1_theo_overlay()

    def draw_panel1_3d(self, points, region, mz_min, mz_max, rt_start, rt_end):
        if not HAVE_GL:
            return
        # The 3D view is just the individual datapoints (no surface).
        if not (isinstance(points, dict) and points["mz"].size):
            self.p1_scatter.setVisible(False)
            return

        mz = points["mz"]; rt = points["rt"]; inten = points["intensity"]

        # Per-point base colour = the panel-2 distribution colour (grey if the
        # point isn't in any distribution), so the 3D dots match panel 2.
        base = np.full((mz.size, 3), 0.5, dtype=np.float32)
        for _did, color, feat_masks in (self._groups or []):
            c = np.array(color, dtype=np.float32) / 255.0
            for fm in feat_masks:
                base[fm] = c
        if self.recovered_visible():
            for _did, color, feat_masks in (self._rec_groups or []):
                c = np.array(color, dtype=np.float32) / 255.0
                for fm in feat_masks:
                    base[fm] = c

        # Show EXACTLY panel 2's visible window (not the padded load region), so
        # the 3D mass/time extent matches panel 2 instead of being "way off".
        if self.window is not None:
            mz_min, mz_max, rt_start, rt_end = self.window
        view = (mz >= mz_min) & (mz <= mz_max) & (rt >= rt_start) & (rt <= rt_end)
        # Honour the noise level: keep assigned points plus only the noise classes
        # the current mode reveals.
        if self._assigned is not None:
            keep = self._assigned | self.noise_visible_mask(mz.size)
            if self.recovered_visible() and self._rec_assigned is not None:
                keep = keep | self._rec_assigned
            view &= keep
        mz, rt, inten, base = mz[view], rt[view], inten[view], base[view]
        if mz.size == 0:
            self.p1_scatter.setVisible(False)
            return

        # Cap the scatter (intensity-priority so the peaks survive) to keep orbit
        # smooth.
        MAX_3D_POINTS = 5000
        if mz.size > MAX_3D_POINTS:
            keep = np.argpartition(inten, -MAX_3D_POINTS)[-MAX_3D_POINTS:]
            mz, rt, inten, base = mz[keep], rt[keep], inten[keep], base[keep]

        zmax = float(inten.max()) or 1.0
        tnorm = (inten / zmax).astype(np.float32)   # 0..1 intensity for colour
        # White-tip gradient: low intensity keeps the distribution colour, high
        # intensity blends to white, so peak tips are white. Log option supported.
        if self.use_logcolor():
            t = np.log1p(inten); t = (t / (t.max() or 1.0)).astype(np.float32)
        else:
            t = tnorm
        t = t[:, None]
        rgb = base * (1.0 - t) + np.float32(1.0) * t
        colors = np.column_stack([rgb, np.ones(rgb.shape[0], dtype=np.float32)]).astype(np.float32)

        # m/z -> GL-x in [-1, 1] (fills the pane width, aligned with panel 2),
        # time -> GL-y (depth), intensity -> GL-z. With fov being horizontal, the
        # visible vertical half-extent is (height/width), so scale intensity to
        # that so it fills the pane height with the 0 baseline at the bottom.
        mz_span = max(mz_max - mz_min, 1e-6)
        rt_span = max(rt_end - rt_start, 1e-6)
        try:
            gl_w = max(self.p1_3d.width(), 1)
            gl_h = max(self.p1_3d.height(), 1)
        except Exception:
            gl_w, gl_h = 1, 1
        vh = gl_h / gl_w
        # The GL widget fills the whole panel, but m/z must start at PLOT_LEFT to
        # line up with panel 2's plot. Screen x in [0, gl_w] maps to GL-x in
        # [-1, 1]; the data's left edge sits at the PLOT_LEFT fraction.
        left_frac = min(max(PLOT_LEFT / gl_w, 0.0), 0.9)
        x_left = -1.0 + 2.0 * left_frac
        x = (x_left + (mz - mz_min) / mz_span * (1.0 - x_left)).astype(np.float32)
        y = ((rt - rt_start) / rt_span * 2 - 1).astype(np.float32)
        # baseline (tnorm=0) -> -vh (bottom), peak (tnorm=1) -> ~+vh (top), with a
        # little headroom so the tallest peak isn't clipped at the very top.
        z = ((tnorm * 1.9 - 1.0) * vh).astype(np.float32)
        pos = np.column_stack([x, y, z]).astype(np.float32)
        try:
            self.p1_scatter.setData(pos=pos, color=colors, size=4.0)
            self.p1_scatter.setVisible(True)
        except Exception:
            pass

    def _features_points(self, feats, cur):
        """Raw points covering a set of features (one charge state's lines),
        read directly from the store so charge states outside the panel-1 window
        still plot. One read over the features' combined m/z x RT span, cached."""
        if not feats:
            return None
        store = (cur.get("profile") if (self.use_profile() and cur.get("profile"))
                 else (cur.get("points") or cur.get("centroid")))
        if store is None:
            return None
        mzlo = min(f["mz_min"] for f in feats); mzhi = max(f["mz_max"] for f in feats)
        rtlo = min(f["rt_start"] for f in feats); rthi = max(f["rt_end"] for f in feats)
        key = (id(store), round(mzlo, 4), round(mzhi, 4), round(rtlo, 4), round(rthi, 4))
        cache = getattr(self, "_feat_points_cache", None)
        if cache is None:
            cache = self._feat_points_cache = {}
        if key in cache:
            return cache[key]
        try:
            store.load_metadata()
            pts = store.extract_points(mzlo, mzhi, rtlo, rthi)
        except Exception:
            pts = None
        if len(cache) > 64:
            cache.clear()
        cache[key] = pts
        return pts

    def distribution_color(self, distribution_id):
        """Stable colour per distribution (assigned in first-seen order) from the
        distinctipy pool."""
        if distribution_id not in self._dist_colors:
            i = len(self._dist_colors)
            pool = self._color_pool or DIST_PALETTE
            self._dist_colors[distribution_id] = pool[i % len(pool)]
        return self._dist_colors[distribution_id]

    def _assignment(self, points, window, status=None):
        """Membership of each raw point in a sqlite distribution. Returns
        ``(assigned_mask, groups)`` where groups is a list of
        ``(distribution_id, colour, [feature_point_masks])`` -- one mask per
        *line* (feature), since a connect-the-dots line must follow a single
        line's trace, not jump across a whole distribution. ``status`` selects the
        confidence tier ('primary', 'recovered', or None=all). Computed once per
        reload and shared by panels 1 and 2 + the noise toggle."""
        if not isinstance(points, dict) or points["mz"].size == 0 or self.db is None:
            return None, []
        mz = points["mz"]; rt = points["rt"]
        mz_min, mz_max, rt_start, rt_end = window
        assigned = np.zeros(mz.size, dtype=bool)
        groups = []
        for dist in self.db.distributions_in_window(mz_min, mz_max, rt_start, rt_end, status=status):
            did = dist["distribution_id"]
            feat_masks = []
            for feat in self.db.distribution_members(did):
                m = ((mz >= feat["mz_min"]) & (mz <= feat["mz_max"]) &
                     (rt >= feat["rt_start"]) & (rt <= feat["rt_end"]) & (~assigned))
                if m.any():
                    assigned |= m
                    feat_masks.append(m)
            if feat_masks:
                groups.append((did, self.distribution_color(did), feat_masks))
        return assigned, groups

    def _compute_noise_class(self, points, window, assigned):
        """Classify every loaded raw point into a noise level.

        0 = assigned (belongs to a distribution; always shown);
        1 = noise line   (a feature/line with >=5 points, in no distribution);
        2 = small line   (a feature with 2..4 points, in no distribution);
        3 = single point (a 1-point feature, or a stray datapoint in no feature).

        Unassigned points default to 3 (lone/stray) and are promoted to a stronger
        (lower) class when they fall inside a larger noise feature's box; a point
        in overlapping noise features takes the strongest (largest-line) class.
        """
        if not isinstance(points, dict) or points["mz"].size == 0:
            return None
        mz = points["mz"]; rt = points["rt"]
        cls = np.where(assigned if assigned is not None else False, 0, 3).astype(np.int8)
        if self.db is None or assigned is None:
            return cls
        unassigned = ~assigned
        if not unassigned.any():
            return cls
        for feat in self.db.noise_features_in_window(*window):
            npoints = feat["n_points"]
            c = 1 if npoints >= 5 else (2 if npoints >= 2 else 3)
            if c == 3:
                continue  # already the default for unassigned points
            m = (unassigned
                 & (mz >= feat["mz_min"]) & (mz <= feat["mz_max"])
                 & (rt >= feat["rt_start"]) & (rt <= feat["rt_end"]))
            if m.any():
                cls[m] = np.minimum(cls[m], np.int8(c))
        return cls

    def _assignment_profile(self, points, window, status=None):
        """Profile-mode membership via peak finding (roadmap 1.9): run the same
        centroiding peak-detection (``axis_peaks``) per scan, take each peak's
        apex as the centroid, match THAT centroid to a sqlite feature, and assign
        EVERY profile datapoint under the peak (left..right) to that feature's
        distribution -- so profile points are linked back to the centroid they'd
        be reduced to and coloured to match, instead of all reading as noise."""
        try:
            from .peaks import axis_peaks
        except ImportError:
            from peaks import axis_peaks
        if not isinstance(points, dict) or points["mz"].size == 0 or self.db is None:
            return None, []
        mz = points["mz"]; rt = points["rt"]; inten = points["intensity"]
        n = mz.size
        mz_min, mz_max, rt_start, rt_end = window
        # Flatten the features into arrays for vectorised apex->feature matching.
        flist = []
        for dist in self.db.distributions_in_window(mz_min, mz_max, rt_start, rt_end, status=status):
            did = dist["distribution_id"]
            color = self.distribution_color(did)
            for feat in self.db.distribution_members(did):
                flist.append((did, color, feat))
        assigned = np.zeros(n, dtype=bool)
        if not flist:
            return assigned, []
        fmzlo = np.array([f["mz_min"] for _, _, f in flist])
        fmzhi = np.array([f["mz_max"] for _, _, f in flist])
        frtlo = np.array([f["rt_start"] for _, _, f in flist])
        frthi = np.array([f["rt_end"] for _, _, f in flist])
        feat_pt = [np.zeros(n, dtype=bool) for _ in flist]
        for r in np.unique(rt):
            si = np.where(rt == r)[0]
            si = si[np.argsort(mz[si])]
            smz = mz[si]; sint = inten[si]
            for l, m, rg in axis_peaks(sint):
                if m < 0 or m >= smz.size:
                    continue
                apex = smz[m]
                cand = np.where((fmzlo <= apex) & (apex <= fmzhi)
                                & (frtlo <= r) & (r <= frthi))[0]
                if cand.size == 0:
                    continue
                pts = si[l:rg + 1]
                feat_pt[cand[0]][pts] = True
                assigned[pts] = True
        by_did = {}
        for fi, (did, color, _f) in enumerate(flist):
            if feat_pt[fi].any():
                by_did.setdefault(did, (color, []))[1].append(feat_pt[fi])
        groups = [(did, color, masks) for did, (color, masks) in by_did.items()]
        return assigned, groups

    def draw_panel2(self, points, window):
        # Connect-the-dots: raw points (m/z x, RT y) as small dots + thin
        # connecting lines. Each *line* (feature) is connected along its own trace
        # (sorted by RT) -- NOT across the whole distribution. Colour is per
        # distribution. Consolidated into one curve + one scatter per distribution
        # (instead of per feature) to keep it fast. Clicking a dot selects the
        # distribution and opens its MS1 panel 3.
        self.p2.clear()
        self._sel_border_item = None   # cleared by p2.clear(); re-added below
        self._p2_scatters = []   # (ScatterPlotItem, base_size) for zoom rescaling
        self.p2.addItem(self.p2_ms2_band)   # survives clear
        self._apply_ms2_band()              # restore the current scan's band
        if not isinstance(points, dict) or points["mz"].size == 0:
            self._render_selection_border()
            return
        mz = points["mz"]
        rt = points["rt"]
        nan = np.array([np.nan])

        for did, color, feat_masks in (self._groups or []):
            # One polyline for the whole distribution, but each line's points are
            # contiguous and separated by NaN so the curve never jumps between
            # different lines.
            xs, ys, smz, srt = [], [], [], []
            for fm in feat_masks:
                o = np.argsort(rt[fm])
                xs.append(mz[fm][o]); xs.append(nan)
                ys.append(rt[fm][o]); ys.append(nan)
                smz.append(mz[fm]); srt.append(rt[fm])
            self.p2.addItem(pg.PlotCurveItem(
                np.concatenate(xs), np.concatenate(ys), connect="finite",
                pen=pg.mkPen(color=(*color, 200), width=0.5)))
            cmz = np.concatenate(smz); crt = np.concatenate(srt)
            scatter = pg.ScatterPlotItem(
                x=cmz, y=crt, size=3, pen=None,
                brush=pg.mkBrush(*color, 230), data=[did] * cmz.size)
            scatter.sigClicked.connect(self.on_panel2_dist_clicked)
            self.p2.addItem(scatter)
            self._p2_scatters.append((scatter, 3))

        # Recovered (less-confident) distributions: same colours, fainter dots +
        # dashed connectors, only from noise level 1 onward.
        if self.recovered_visible():
            for did, color, feat_masks in (self._rec_groups or []):
                xs, ys, smz, srt = [], [], [], []
                for fm in feat_masks:
                    o = np.argsort(rt[fm])
                    xs.append(mz[fm][o]); xs.append(nan)
                    ys.append(rt[fm][o]); ys.append(nan)
                    smz.append(mz[fm]); srt.append(rt[fm])
                self.p2.addItem(pg.PlotCurveItem(
                    np.concatenate(xs), np.concatenate(ys), connect="finite",
                    pen=pg.mkPen(color=(*color, 130), width=0.5, style=Qt.DashLine)))
                cmz = np.concatenate(smz); crt = np.concatenate(srt)
                sc = pg.ScatterPlotItem(x=cmz, y=crt, size=2.5, pen=None,
                                        brush=pg.mkBrush(*color, 150), data=[did] * cmz.size)
                sc.sigClicked.connect(self.on_panel2_dist_clicked)
                self.p2.addItem(sc)
                self._p2_scatters.append((sc, 2.5))

        # Noise points -> one scatter per visible noise class (level-dependent).
        noise_vis = self.noise_visible_mask(mz.size)
        if noise_vis.any():
            for cls, (color, size) in NOISE_STYLE.items():
                nm = noise_vis & (self._noise_class == cls)
                if nm.any():
                    sc = pg.ScatterPlotItem(x=mz[nm], y=rt[nm], size=size, pen=None,
                                            brush=pg.mkBrush(*color))
                    self.p2.addItem(sc)
                    self._p2_scatters.append((sc, size))

        self._render_selection_border()   # dotted rect around the selected distribution
        self._rescale_points()

    def on_panel2_dist_clicked(self, _scatter, points):
        """Clicking a distribution's dots in panel 2 selects it and brings up the
        MS1 panel 3 for that distribution."""
        if points is None or len(points) == 0 or self.current is None:
            return
        did = points[0].data()
        if did is None:
            return
        # Bump the distribution-click sequence so the deferred empty-space check
        # in _on_p2_clicked knows a distribution was hit during this dispatch
        # (robust to scatter/scene signal ordering, and to panel-1 dot clicks
        # which route here too but fire no panel-2 scene click).
        self._dist_click_seq = getattr(self, "_dist_click_seq", 0) + 1
        # A real distribution is now selected: any theoretical overlay from a
        # prior peptide / charge-scroll hypothetical goes away.
        self._ms1_theo = None
        self._draw_ms1_theo_overlay()
        # Select it: dotted border + charge-search anchor (no view move on click).
        self._set_selected(did)
        # Selecting an MS1 distribution returns panel 3 to its MS1 view.
        self._panel3_mode = "ms1"
        self._ms2_scan = None
        self._ms2_band = None
        self._apply_ms2_band()
        # Refresh table 1 to this distribution's members, then redraw panel 3.
        try:
            self.table1_for_distribution(did)
        except Exception:
            pass
        self.draw_panel3_ms1(self.current,
                             getattr(self, "_last_scan_mz", None),
                             getattr(self, "_last_scan_int", None))

    def _on_p2_clicked(self, event):
        """Panel-2 scene click. A click on a distribution's dots is handled by
        the scatter's sigClicked (on_panel2_dist_clicked), which bumps
        _dist_click_seq; anything else is empty space -> deselect. Deferred to
        the next event-loop tick so the scatter signal (fired in the same click
        dispatch, order unspecified) is accounted for before we decide."""
        try:
            if event.double() or event.button() != Qt.LeftButton:
                return
        except Exception:
            pass
        seq = getattr(self, "_dist_click_seq", 0)
        QTimer.singleShot(0, lambda: self._resolve_p2_click(seq))

    def _resolve_p2_click(self, seq):
        # A distribution was clicked during this dispatch iff the sequence moved.
        if getattr(self, "_dist_click_seq", 0) != seq:
            return
        self._deselect_distribution()

    def _deselect_distribution(self):
        """Drop the current distribution selection: remove the panel-2 border,
        clear the charge-search anchor, and hide the theoretical MS1 overlay."""
        self._selected_dist_id = None
        self._selected_charge_group = None
        self._selected_bbox = None
        self._border_color = None
        self._clear_selection()          # removes the dotted border
        self._ms1_theo = None            # theoretical distributions disappear
        self._draw_ms1_theo_overlay()

    def draw_ms2_strip(self, ms2):
        # Keep every loaded MS2 scan; the strip only shows the ones whose RT *and*
        # precursor m/z fall inside panel 2's current view (so it tracks zoom and
        # doesn't show "tons" of scans whose precursor isn't even on screen).
        self._ms2_all = [m for m in ms2 if m.get("rt") is not None]
        self._refresh_ms2_visible()

    def _refresh_ms2_visible(self):
        """Filter the loaded MS2 scans to panel 2's current m/z x RT view and
        redraw the strip lines."""
        if not hasattr(self, "ms2_curve"):
            return
        try:
            (mz0, mz1) = self.p2.getViewBox().viewRange()[0]
            (rt0, rt1) = self.p2.getViewBox().viewRange()[1]
        except Exception:
            mz0, mz1, rt0, rt1 = -np.inf, np.inf, -np.inf, np.inf
        visible = []
        for m in self._ms2_all:
            rt = m.get("rt"); pmz = m.get("mz")
            if rt is None or not (rt0 <= rt <= rt1):
                continue
            # If we know the precursor m/z, require it inside the view too.
            if pmz is not None and not (mz0 <= pmz <= mz1):
                continue
            visible.append(m)
        scans = sorted(((m["rt"], m) for m in visible), key=lambda t: t[0])
        self._ms2_scans = scans
        ident = self._identified_scans()
        gx, gy, rx, ry = [], [], [], []
        for rt, m in scans:
            if str(m.get("number", "")) in ident:
                gx += [0.0, 1.0, np.nan]; gy += [rt, rt, np.nan]
            else:
                rx += [0.0, 1.0, np.nan]; ry += [rt, rt, np.nan]
        self.ms2_curve.setData(np.array(gx, dtype=float), np.array(gy, dtype=float))
        self.ms2_curve_red.setData(np.array(rx, dtype=float), np.array(ry, dtype=float))
        self.ms2_plot.getViewBox().setXRange(0, 1, padding=0)
        # Inverse-scale the line thickness with RT zoom: the more you zoom in, the
        # THICKER each MS2 line gets (so it never fades to nothing). Reference is
        # the full loaded RT span.
        try:
            (rv0, rv1) = self.p2.getViewBox().viewRange()[1]
            view_span = max(rv1 - rv0, 1e-9)
            full_span = (self._loaded_window[3] - self._loaded_window[2]
                         if self._loaded_window else view_span)
            factor = min(max((full_span / view_span) ** 0.6, 1.0), 6.0)
            self.ms2_curve.setPen(pg.mkPen(60, 200, 100, 255, width=3.0 * factor))
            self.ms2_curve_red.setPen(pg.mkPen(220, 70, 70, 255, width=3.0 * factor))
        except Exception:
            pass

    def _on_fdr_changed(self):
        """FDR spin box changed: update the threshold and recolor both the MS2
        strip lines AND the on-panel-2 isolation band."""
        self._fdr_threshold = max(0.0, self.fdr_edit.value()) / 100.0
        self._ident_cache_key = None   # invalidate cache
        self._refresh_ms2_visible()
        self._apply_ms2_band()

    def _identified_scans(self):
        """Set of MS2 scan numbers that have a PSM passing the FDR acceptance
        criteria (q-value <= threshold). Cached per file + threshold."""
        thr = getattr(self, "_fdr_threshold", 0.001)
        cache_key = (self.current_file, thr)
        if getattr(self, "_ident_cache_key", None) == cache_key:
            return self._ident_cache
        ident = set()
        for r in self.file_psms():
            q = safe_float(r.get("percolator_q"))
            scan_no = str(r.get("scan", "") or "")
            # q is None -> treat as identified (no q column); else gate on FDR.
            if scan_no and (q is None or q <= thr):
                ident.add(scan_no)
        self._ident_cache_key = cache_key
        self._ident_cache = ident
        return ident

    def _band_color(self, scan):
        """Green if the scan has an identified peptide (passes FDR), else red;
        50% opacity to match the strip lines."""
        ident = self._identified_scans()
        if scan is not None and str(scan.get("number", "")) in ident:
            return (60, 200, 100, 150)
        return (220, 70, 70, 150)

    def _draw_band(self, low, high, rt, scan):
        self.p2_ms2_band.setPen(pg.mkPen(*self._band_color(scan), width=3))
        self.p2_ms2_band.setData([low, high], [rt, rt])

    def _apply_ms2_band(self):
        """(Re)draw the persistent isolation band on panel 2 for the current MS2
        scan (coloured green/red by identification), or clear it when none."""
        band = getattr(self, "_ms2_band", None)
        if band is not None and band[0] is not None and band[2] is not None:
            low, high, rt = band
            self._draw_band(low, high, rt, getattr(self, "_ms2_scan", None))
        else:
            self.p2_ms2_band.setData([], [])

    def _on_ms2_hover(self, scene_pos):
        """Hovering an MS2 line on the strip previews its isolation band on panel
        2; leaving the strip reverts to the band of the scan being viewed."""
        if not self._ms2_scans:
            self._apply_ms2_band()
            return
        try:
            vb = self.ms2_plot.getViewBox()
            if not vb.sceneBoundingRect().contains(scene_pos):
                self._apply_ms2_band()
                return
            y = vb.mapSceneToView(scene_pos).y()
            (y0, y1) = vb.viewRange()[1]
        except Exception:
            return
        rt_arr = np.array([rt for rt, _ in self._ms2_scans])
        i = int(np.argmin(np.abs(rt_arr - y)))
        scan = self._ms2_scans[i][1]
        # m/z span = the isolation window; fall back to a small default width.
        low = scan.get("iso_low")
        high = scan.get("iso_high")
        if low is None or high is None:
            mz = scan.get("mz")
            low, high = (mz - 0.5, mz + 0.5) if mz is not None else (None, None)
        if abs(rt_arr[i] - y) <= max(abs(y1 - y0) * 0.03, 1e-4) and low is not None:
            self._draw_band(low, high, scan["rt"], scan)
        else:
            self._apply_ms2_band()

    def _on_ms2_strip_clicked(self, event):
        """Click anywhere on the MS2 strip -> select the nearest MS2 line by RT."""
        if not self._ms2_scans:
            return
        try:
            vb = self.ms2_plot.getViewBox()
            pos = vb.mapSceneToView(event.scenePos())
        except Exception:
            return
        y = pos.y()
        rt_arr = np.array([rt for rt, _ in self._ms2_scans])
        i = int(np.argmin(np.abs(rt_arr - y)))
        # Only act if the click is near a line: within ~3% of the visible RT span.
        (y0, y1) = vb.viewRange()[1]
        tol = max(abs(y1 - y0) * 0.03, 1e-4)
        if abs(rt_arr[i] - y) <= tol:
            self.render_ms2(self._ms2_scans[i][1])

    def render_ms2(self, scan):
        """Switch panel 3 to the MS2 spectrum for ``scan`` and remember it, so a
        later background reload keeps showing MS2 instead of snapping back to the
        MS1 view."""
        cur = self.current
        if cur is None or scan is None or cur.get("centroid") is None:
            return
        self._panel3_mode = "ms2"
        self._ms2_scan = scan
        # Keep the isolation band on panel 2 for the scan being viewed (persists
        # through redraws + Table-2 peptide selection, not just on hover).
        low, high = scan.get("iso_low"), scan.get("iso_high")
        if low is None or high is None:
            mz = scan.get("mz")
            # No isolation window in the mzML -> a visible default width around the
            # precursor (a zero-width line would not render).
            low, high = (mz - 0.5, mz + 0.5) if mz is not None else (None, None)
        self._ms2_band = (low, high, scan.get("rt")) if low is not None else None
        self._apply_ms2_band()
        # Table 2 (other candidate PSMs) only appears for the MS2 view.
        self.dock_table2.show()
        try:
            spectrum = cur["centroid"].get_scan_by_id(scan["id"])
            if spectrum is None:
                self.p3_title.setText(f"MS2 scan {scan.get('number','')} not found")
                return
            mz, inten = scan_arrays(spectrum)
            # Remember the spectrum so a Table-2 peptide selection can overlay its
            # fragment ions without re-reading the scan.
            self._ms2_mz, self._ms2_int = mz, inten
            self._frag_overlay = []          # cleared by plot_spectrum's clear()
            self._frag_token = getattr(self, "_frag_token", 0) + 1
            self.p3_stack.setCurrentIndex(0)
            # Default title until a peptide is picked in Table 2 (then it becomes
            # just the peptide).
            plot_spectrum(self.p3, mz, inten, color=palette(self.theme)["fg"],
                          title=self._ms2_title(scan))
            # Ground the baseline at 0 (no gap below the sticks).
            top = float(np.max(inten)) * 1.05 if len(inten) else 1.0
            self.p3.getViewBox().setYRange(0.0, top, padding=0)
            self._fit_p3_ms2_xrange()   # fit m/z to the data (labels added later)
            self.p3_title.setText(
                f"MS2 scan {scan.get('number','')}  rt={scan['rt']:.3f}  precursor m/z={scan.get('mz')}")
            self.populate_table2(scan, mz)
        except Exception as exc:
            self.p3_title.setText(f"MS2 load error: {exc}")
        # Refresh the MS2-scan dropdown to the distribution's scans, selecting
        # the one now shown. (Guarded so it doesn't re-enter render_ms2.)
        self._populate_ms2_combo(scan)
        self._sync_panel3_tab()

    def _ms2_title(self, scan, peptide=None):
        """Panel-3 MS2 title. Always shows the precursor isolation window +
        RT (the MS1 scan parameters, so the plotted lines can be eyeballed);
        prefixes the peptide once one is picked in Table 2."""
        low, high = scan.get("iso_low"), scan.get("iso_high")
        if low is None or high is None:
            mz = scan.get("mz")
            low, high = (mz - 0.5, mz + 0.5) if mz is not None else (None, None)
        rt = scan.get("rt")
        if low is not None and high is not None and rt is not None:
            win = f"{low:.4f} - {high:.4f} at RT {rt:.3f}"
        else:
            win = ""
        if peptide:
            return f"{peptide}    {win}".strip()
        return win or "MS2"

    def populate_table2(self, scan, peak_mz=None):
        # Candidate PSMs near this precursor m/z in the current file. The coverage
        # column is computed from the SAME generate-and-match used to annotate the
        # MS2 spectrum in panel 3 (fragments.peptide_fragment_ions ->
        # annotate_spectrum -> coverage_metrics matchcounts), so the table value
        # and the green ions in panel 3 always agree.
        prec = scan.get("mz")
        peaks = np.asarray(peak_mz, dtype=float) if peak_mz is not None else None
        ints = getattr(self, "_ms2_int", None)
        low = scan.get("iso_low")
        high = scan.get("iso_high")
        if low is None or high is None:
            base = prec or 0.0
            low, high = base - 1.0, base + 1.0
        # Candidates = the peptides the SEARCH would have considered for this
        # precursor: any file PSM whose precursor m/z falls within the search
        # precursor tolerance of this scan's precursor, allowing the search's
        # isotope-error offsets (Sage precursor_tol + isotope_errors). No
        # hard-coded Da window.
        NEUTRON = 1.0033548
        iso_lo, iso_hi = self.precursor_isotope_errors
        scan_no = str(scan.get("number", "") or "")
        rows = []
        for r in self.file_psms():
            # A PSM identified FROM this exact scan is a candidate by definition:
            # include it unconditionally, before any precursor-m/z math. This is
            # what guarantees the identified peptide the user is looking at is in
            # the table even when its row lacks calc/exp mass or its computed m/z
            # drifts outside the precursor tolerance.
            if scan_no and str(r.get("scan", "") or "") == scan_no:
                rows.append(r)
                continue
            try:
                row_mz = peptide_mass(r)
                z = peptide_charge(r) or 1
                psm_mz = (row_mz / z + 1.007276) if row_mz else None
            except Exception:
                psm_mz = None
            if prec is None:
                rows.append(r)
                continue
            if psm_mz is None:
                continue
            for k in range(iso_lo, iso_hi + 1):
                target = prec + k * NEUTRON / z
                if abs(psm_mz - target) <= target * self.precursor_ppm / 1e6:
                    rows.append(r)
                    break
        # Compute each candidate's fragment coverage (matchcounts), then rank the
        # peptides best-first so the most distinguishable candidate is at the top.
        reports = []
        for r in rows:
            mc = 0
            if peaks is not None and peaks.size and ints is not None:
                seq = plain_seq(r.get("peptide", ""))
                z = peptide_charge(r) or scan.get("charge") or 2
                try:
                    ents = seqfragments_module.peptide_fragment_ions(
                        seq, z, low, high, dividingthreshold=0.1,
                        subisotopomericdepth=0.5)
                    matched, _ = seqfragments_module.annotate_spectrum(
                        ents, peaks, ints, ppm=self.frag_ppm)
                    labels = [row["ion"] for m in matched for row in m.get("rows", [])]
                    mc = seqcoverage.coverage_metrics(seq, labels)["matchcounts"]
                except Exception:
                    mc = 0
            reports.append((r, {"matchcounts": mc}))
        # Rank by matchcounts (the reference's secondfinalmetrics): higher = more
        # fragment coverage = more confident. Best-first by default.
        reports.sort(key=lambda rr: rr[1]["matchcounts"], reverse=True)

        # Disable sorting while filling so row indices stay put; the default
        # best-first order shows once re-enabled (until the user clicks a header
        # to re-sort).
        self.table2.setSortingEnabled(False)
        # Keep the current theoretical MS1 overlay on panel 1 when the MS2 view
        # (re)loads -- e.g. switching to the MS2 tab -- so it persists from the
        # selected peptide. It's replaced only when a new Table-2 candidate is
        # picked (_on_table2_peptide_selected) or a new peptide is selected in
        # the list (update_evidence), not dropped here.
        self.table2.clearSelection()
        self.table2.setRowCount(len(reports))
        for i, (r, rep) in enumerate(reports):
            pep_item = QTableWidgetItem(str(r.get("peptide", "")))
            pep_item.setData(Qt.UserRole, r)   # full PSM row for selection lookup
            self.table2.setItem(i, 0, pep_item)
            self.table2.setItem(i, 1, NumericItem(str(r.get("percolator_q", "")),
                                                  safe_float(r.get("percolator_q"))))
            if peaks is not None and peaks.size:
                cov_val = rep.get("matchcounts", 0)
                cov = str(cov_val)
            else:
                cov, cov_val = "", None
            self.table2.setItem(i, 2, NumericItem(cov, cov_val))
        self.table2.setSortingEnabled(True)

    def _on_table2_peptide_selected(self):
        """A candidate peptide was clicked in Table 2 -> overlay its theoretical
        b/y fragment ions on the MS2 spectrum (green = matched, red = absent),
        and set the panel-3 title to just that peptide."""
        items = self.table2.selectedItems()
        scan = getattr(self, "_ms2_scan", None)
        if not items or scan is None or self._panel3_mode != "ms2":
            return
        row = items[0].row()
        pep_item = self.table2.item(row, 0)
        r = pep_item.data(Qt.UserRole) if pep_item is not None else None
        if not r:
            return
        seq = plain_seq(r.get("peptide", ""))
        if len(seq) < 2:
            return
        charge = peptide_charge(r) or scan.get("charge") or 2
        # Auto-select the matching MS1 distribution in panel 2 (dotted border)
        # WITHOUT flipping panel 3 back to its MS1 view (we stay on the MS2
        # spectrum). The band + border survive panel-2 redraws.
        self._select_distribution_for_candidate(r, charge, scan)
        # Also overlay this peptide's theoretical MS1 distribution on panel 1.
        self._ms1_theo = {"seq": seq, "charge": charge,
                          "mod_mass": peptide_mod_mass(r.get("peptide", ""))}
        self._draw_ms1_theo_overlay()
        low = scan.get("iso_low")
        high = scan.get("iso_high")
        if low is None or high is None:
            # No isolation window in the mzML -> use a default ~1 m/z window so at
            # least the monoisotopic precursor is included.
            prec = scan.get("mz") or 0.0
            low, high = prec - 1.0, prec + 1.0
        self._frag_token = getattr(self, "_frag_token", 0) + 1
        # Title stays the MS1 isolation window / RT (no peptide), so the plotted
        # lines can be verified against the scan parameters.
        self.p3.setTitle(self._ms2_title(scan), color=palette(self.theme)["fg"])
        worker = FragmentWorker(
            seq, charge, low, high,
            getattr(self, "_ms2_mz", np.array([])),
            getattr(self, "_ms2_int", np.array([])),
            self.frag_ppm, self._frag_token,
            dividing_threshold=0.1, subiso_depth=0.5)
        worker.done.connect(self._on_fragments_ready)
        self._frag_workers = [w for w in getattr(self, "_frag_workers", [])
                              if w.isRunning()]
        self._frag_workers.append(worker)
        worker.start()

    def _select_distribution_for_candidate(self, r, charge, scan):
        """Select the MS1 distribution matching a Table-2 candidate (dotted border
        on panel 2) without redrawing panel 3 (keeps the MS2 view up).

        Targets the candidate's actual mono m/z (within the precursor tolerance)
        and RT, and picks the CLOSEST distribution -- not just the highest-quality
        one in a wide window, which otherwise stuck on the first selection."""
        if self.db is None:
            return
        try:
            neutral = peptide_mass(r)
            if not neutral:
                return
            z = max(1, int(charge or 1))
            mono_mz = neutral / z + 1.007276
            rt = scan.get("rt") if scan else None
            if rt is None and isinstance(self.current, dict):
                rt = self.current.get("rt")
            if rt is None:
                return
            # Tight m/z window around the candidate's mono m/z (search precursor
            # tolerance, with a small absolute floor for isotope/centroiding slack).
            tol = max(mono_mz * self.precursor_ppm / 1e6, 0.02)
            dists = self.db.distributions_in_window(
                mz_min=mono_mz - tol, mz_max=mono_mz + tol,
                rt_start=rt - self.rt_half, rt_end=rt + self.rt_half,
                charge=z, limit=50)
            # If this candidate maps to a distribution, select it. If NOT, keep
            # whatever was already selected (don't clear it) -- the selection
            # stays true while candidates for this precursor are compared; only a
            # candidate that genuinely maps to a different distribution moves it.
            if dists:
                best = min(dists, key=lambda d: (
                    abs((d.get("mono_mz") or mono_mz) - mono_mz),
                    abs((d.get("rt_apex") or rt) - rt)))
                # _set_selected only draws the border + anchors charge search; it
                # does NOT touch panel 3, so the MS2 spectrum stays up.
                self._set_selected(best["distribution_id"])
        except Exception:
            pass

    def _on_ms1theo_changed(self):
        """raw <-> summed switch for the theoretical MS1 overlay."""
        self._draw_ms1_theo_overlay()

    def _p1_experimental_max(self):
        """Max experimental intensity currently drawn in panel 1 (2D)."""
        ymax = 0.0
        for scatter, _base in getattr(self, "_p1_scatters", []):
            try:
                ys = scatter.getData()[1]
                if ys is not None and len(ys):
                    ymax = max(ymax, float(np.nanmax(ys)))
            except Exception:
                pass
        return ymax

    def _distribution_experimental_max(self, did):
        """Max experimental intensity of just ONE distribution's points (within
        the visible window), used to normalize its theoretical MS1 overlay."""
        points = getattr(self, "_last_points", None)
        if did is None or not isinstance(points, dict):
            return 0.0
        mz = points.get("mz")
        rt = points.get("rt")
        inten = points.get("intensity")
        if inten is None or mz is None or not len(inten):
            return 0.0
        if self.window is not None:
            vmz0, vmz1, vrt0, vrt1 = self.window
            vm = (mz >= vmz0) & (mz <= vmz1) & (rt >= vrt0) & (rt <= vrt1)
        else:
            vm = np.ones(mz.size, dtype=bool)
        idx = np.zeros(mz.size, dtype=bool)
        for g_did, _color, feat_masks in ((self._groups or []) + (self._rec_groups or [])):
            if g_did == did:
                for fm in feat_masks:
                    idx |= fm
                break
        idx &= vm
        return float(inten[idx].max()) if idx.any() else 0.0

    def _draw_ms1_theo_overlay(self):
        """Overlay the selected Table-2 peptide's theoretical MS1 isotope
        distribution on panel 1 (2D only) as a 50%-opacity bar chart, normalized
        so the tallest theoretical bar matches the tallest experimental peak of
        the DISTRIBUTION it matched (not all the data visible in the plot)."""
        # Remove EVERY theoretical bar on panel 1, not just the last-tracked one.
        # The overlay is the only BarGraphItem panel 1 ever holds (experimental
        # data are ScatterPlotItems), so sweeping all bars guarantees a single
        # overlay even if a prior draw left a stray (e.g. an item orphaned by a
        # plot clear() whose tracked reference was then overwritten).
        try:
            for it in [i for i in self.p1_2d.getPlotItem().items
                       if isinstance(i, pg.BarGraphItem)]:
                self.p1_2d.removeItem(it)
        except Exception:
            pass
        self._ms1_theo_item = None
        # 2D only.
        if (self._ms1_theo is None
                or getattr(self, "p1_stack", None) is None
                or self.p1_stack.currentIndex() != 0):
            return
        # Normalize to the matched distribution's max; fall back to the overall
        # visible max, then to the current y-range top, so the overlay ALWAYS
        # shows (even when a distribution is selected / data is sparse).
        exp_max = self._distribution_experimental_max(
            getattr(self, "_selected_dist_id", None))
        if exp_max <= 0:
            exp_max = self._p1_experimental_max()
        if exp_max <= 0:
            try:
                exp_max = float(self.p1_2d.getViewBox().viewRange()[1][1])
            except Exception:
                exp_max = 1.0
        if exp_max <= 0:
            exp_max = 1.0
        try:
            mode = self.ms1theo_toggle.key()
            mzs, abund = isotopes.peptide_isotope_bars(
                self._ms1_theo["seq"], self._ms1_theo["charge"], mode=mode,
                dividing_threshold=0.1)
        except Exception as exc:
            # Don't silently drop the overlay: a computation failure (e.g. a
            # residue with no atomic composition) otherwise looks identical to
            # "no distribution here". Surface it so it's diagnosable.
            import sys
            print(f"[ms1-theo] no overlay for "
                  f"{self._ms1_theo.get('seq')!r}: {exc}", file=sys.stderr)
            return
        if mzs is None or len(mzs) == 0:
            return
        # The bars come from the bare sequence; shift them by the peptide's
        # modification delta mass / charge so a modified peptide (e.g. +15.99
        # oxidation) lines up with its real precursor instead of sitting
        # ~mod_mass/z Th to the left.
        mod_mass = float(self._ms1_theo.get("mod_mass", 0.0) or 0.0)
        z = max(1, int(self._ms1_theo["charge"] or 1))
        mzs = np.asarray(mzs, dtype=float)
        abund = np.asarray(abund, dtype=float)
        if mod_mass:
            mzs = mzs + mod_mass / z
        # Merge isotopomer bars that fall within a bar width of each other. In
        # raw mode, isotopomers sharing a nominal M+N (e.g. a 13C vs a 15N
        # substitution) sit ~0.002 Th apart -- unresolvable and rendered as two
        # overlapping bars at the same spot, which reads as the distribution
        # being "plotted twice". Collapse each cluster into one abundance-summed
        # bar so a position is only ever drawn once. Summed mode has no such
        # near-duplicates, so this leaves it unchanged.
        mzs, abund = self._merge_close_bars(mzs, abund, tol=0.02)
        theo_max = float(np.max(abund)) or 1.0
        # Normalize: tallest theoretical bar == tallest experimental peak.
        heights = np.asarray(abund, dtype=float) * (exp_max / theo_max)
        # 50%-opacity theme-adaptive bars (white on the dark plot, near-black on
        # the light plot), drawn BENEATH the experimental data so the measured
        # points sit visually in front.
        shade = 255 if self.theme == "dark" else 20
        bars = pg.BarGraphItem(x=np.asarray(mzs, dtype=float), height=heights,
                               width=0.02, brush=pg.mkBrush(shade, shade, shade, 128),
                               pen=pg.mkPen(shade, shade, shade, 128))
        bars.setZValue(0)   # beneath the experimental scatters (z=1), above bg
        self.p1_2d.addItem(bars)
        self._ms1_theo_item = bars

    @staticmethod
    def _merge_close_bars(mzs, abund, tol=0.02):
        """Collapse bars closer than ``tol`` m/z into one abundance-summed bar at
        their abundance-weighted mean m/z (so visually-coincident isotopomers
        aren't drawn on top of each other)."""
        if mzs is None or len(mzs) == 0:
            return mzs, abund
        order = np.argsort(mzs)
        mzs = np.asarray(mzs, dtype=float)[order]
        abund = np.asarray(abund, dtype=float)[order]
        out_mz, out_ab = [], []
        cluster_mz = [mzs[0]]
        cluster_ab = [abund[0]]
        for m, a in zip(mzs[1:], abund[1:]):
            if m - cluster_mz[-1] <= tol:
                cluster_mz.append(m)
                cluster_ab.append(a)
            else:
                tot = sum(cluster_ab) or 1.0
                out_mz.append(sum(mm * aa for mm, aa in zip(cluster_mz, cluster_ab)) / tot)
                out_ab.append(sum(cluster_ab))
                cluster_mz, cluster_ab = [m], [a]
        tot = sum(cluster_ab) or 1.0
        out_mz.append(sum(mm * aa for mm, aa in zip(cluster_mz, cluster_ab)) / tot)
        out_ab.append(sum(cluster_ab))
        return np.array(out_mz), np.array(out_ab)

    def _on_fragments_ready(self, result):
        if result.get("token") != getattr(self, "_frag_token", 0):
            return
        if result.get("error"):
            return
        matched = result.get("matched", [])
        unmatched = result.get("unmatched", (np.array([]), np.array([])))
        # Remember it so a theme redraw (which rebuilds the spectrum) can restore
        # the green/red annotation without re-running the worker.
        self._last_frag = (matched, unmatched)
        self._draw_fragment_overlay(matched, unmatched)

    def _clear_fragment_overlay(self):
        for item in getattr(self, "_frag_overlay", []):
            try:
                self.p3.removeItem(item)
            except Exception:
                pass
        self._frag_overlay = []

    @staticmethod
    def _sticks(mzs, ints):
        """NaN-separated (x, y) for drawing baseline-grounded spectrum sticks."""
        mzs = np.asarray(mzs, dtype=float)
        ints = np.asarray(ints, dtype=float)
        if mzs.size == 0:
            return np.array([]), np.array([])
        x = np.empty(mzs.size * 3)
        y = np.empty(mzs.size * 3)
        x[0::3] = mzs
        x[1::3] = mzs
        x[2::3] = np.nan
        y[0::3] = 0.0
        y[1::3] = ints
        y[2::3] = np.nan
        return x, y

    def _draw_fragment_overlay(self, matched, unmatched):
        """Recolour the ACTUAL MS2 spectrum: peaks matched by a theoretical b/y
        fragment ion are drawn green (and labelled with the ion + the MS1
        isotope it came from); the remaining real peaks are drawn red.
        Theoretical ions with no experimental match are not drawn."""
        self._clear_fragment_overlay()
        pal = palette(self.theme)
        fg = pal["fg"]
        dark = self.theme == "dark"
        green = (90, 220, 120) if dark else (0, 150, 0)
        red = (235, 90, 90) if dark else (205, 0, 0)
        overlay = []

        um_mz, um_int = unmatched
        ux, uy = self._sticks(um_mz, um_int)
        if ux.size:
            red_curve = pg.PlotCurveItem(ux, uy, pen=pg.mkPen(*red, width=1),
                                         connect="finite")
            red_curve.setZValue(18)
            self.p3.addItem(red_curve)
            overlay.append(red_curve)

        if matched:
            mx, my = self._sticks([m["mz"] for m in matched],
                                  [m["intensity"] for m in matched])
            green_curve = pg.PlotCurveItem(mx, my, pen=pg.mkPen(*green, width=2),
                                           connect="finite")
            green_curve.setZValue(20)
            self.p3.addItem(green_curve)
            overlay.append(green_curve)
            for m in matched:
                # Every ion that matched this peak is its own row (already sorted
                # by decreasing ppm error), with the ppm error to the right of the
                # M+N label. Stacked as one multi-line label above the peak.
                lines = []
                for row in m.get("rows", []):
                    isos = row.get("isotopes", [])
                    iso_txt = ",".join(f"M+{i}" for i in isos) if isos else ""
                    z = row.get("charge", 1)
                    ztxt = "" if z == 1 else f"({z}+)"
                    lines.append(
                        f"{row['ion']}{ztxt} {iso_txt}  {row['ppm']:+.1f} ppm".strip())
                if not lines:
                    continue
                label = pg.TextItem("\n".join(lines), color=fg, anchor=(0.5, 1.0))
                label.setPos(m["mz"], m["intensity"])
                label.setZValue(21)
                self.p3.addItem(label)
                overlay.append(label)
        self._frag_overlay = overlay
        # Widen the m/z axis so edge labels aren't clipped (below).
        self._fit_p3_ms2_xrange()

    def _fit_p3_ms2_xrange(self):
        """Set the MS2 plot's m/z bounds to whichever reaches further -- the peak
        data or the fragment-ion text labels. Labels are screen-fixed size and
        anchored (centred) on their peak, so a label on an edge peak overhangs
        the data; convert its pixel half-width to m/z and extend the range to
        fit. The m/z-per-pixel is referenced to the DATA span (not the live view)
        so widening can't feed back into a wider label and diverge."""
        if self._panel3_mode != "ms2":
            return
        mz = getattr(self, "_ms2_mz", None)
        if mz is None or not len(mz):
            return
        data_lo, data_hi = float(np.min(mz)), float(np.max(mz))
        span = max(data_hi - data_lo, 1e-9)
        vb = self.p3.getViewBox()
        try:
            px_w = float(vb.width())
        except Exception:
            px_w = 0.0
        xscale = span / px_w if px_w > 1 else 0.0
        lo, hi = data_lo, data_hi
        for it in getattr(self, "_frag_overlay", []):
            if not isinstance(it, pg.TextItem):
                continue
            try:
                x = float(it.pos().x())
                half = it.boundingRect().width() * xscale / 2.0   # anchor 0.5
                lo = min(lo, x - half)
                hi = max(hi, x + half)
            except Exception:
                pass
        pad = span * 0.02
        vb.setXRange(lo - pad, hi + pad, padding=0)

    # Wheel over a plot's y-axis strip scrolls intensity (y zoom) with the
    # baseline pinned at 0, even though y-drag inside the plot is disabled. Used
    # by both panel 1 (2D) and panel 3.
    def eventFilter(self, obj, event):
        # Guard: events can arrive while the docks are still being built (some
        # of these widgets don't exist yet).
        grid = getattr(self, "p3_grid", None)
        if grid is not None and obj is grid.viewport() and event.type() == QEvent.Wheel:
            if self._grid_axis_wheel(event):
                return True
            return super().eventFilter(obj, event)
        p1 = getattr(self, "p1_2d", None)
        p2 = getattr(self, "p2", None)
        p3 = getattr(self, "p3", None)
        # Panel 1 / panel 3 intensity (y-axis strip) zoom: wheel over the left
        # axis scrolls intensity with the baseline pinned at 0. Panel 1's
        # intensity behaviour is intentionally left untouched.
        target = None
        if p1 is not None and obj is p1.viewport():
            target = self.p1_2d
        elif p3 is not None and obj is p3.viewport():
            target = self.p3
        if target is not None and event.type() == QEvent.Wheel:
            axis = target.getPlotItem().getAxis("left")
            if event.position().x() < axis.width():
                vb = target.getViewBox()
                (y0, y1) = vb.viewRange()[1]
                factor = 0.85 if event.angleDelta().y() > 0 else 1.0 / 0.85
                # Keep the baseline pinned at 0 (intensity is always >= 0); only
                # the top of the y range moves, so the data baseline never lifts.
                vb.setYRange(0.0, y1 * factor, padding=0)
                return True

        # Panels 1 & 2 plot area: plain wheel pans the window (m/z on panel 1's
        # x-axis; m/z / time on panel 2), Ctrl+wheel falls through to pyqtgraph's
        # default zoom. Panel 1's intensity (y) axis is excluded above.
        pan_target = None
        if p1 is not None and obj is p1.viewport():
            pan_target = self.p1_2d
        elif p2 is not None and obj is p2.viewport():
            pan_target = self.p2
        if pan_target is not None and event.type() == QEvent.Wheel:
            if event.modifiers() & Qt.ControlModifier:
                # Let pyqtgraph zoom (current behaviour).
                return super().eventFilter(obj, event)
            vb = pan_target.getViewBox()
            dy = event.angleDelta().y()
            dx = event.angleDelta().x()
            if pan_target is self.p1_2d:
                # Panel 1: vertical wheel pans m/z (x). Intensity y stays fixed.
                self._pan_axis(vb, 0, dy)
            else:
                # Panel 2: horizontal wheel, Shift+wheel, or a wheel over the
                # bottom m/z axis pans m/z (x); vertical wheel over the plot
                # pans time (y).
                bottom = pan_target.getPlotItem().getAxis("bottom")
                plot_h = pan_target.viewport().height()
                over_mz_axis = event.position().y() > (plot_h - bottom.height())
                if (event.modifiers() & Qt.ShiftModifier) or over_mz_axis:
                    self._pan_axis(vb, 0, dy or dx)
                elif abs(dx) > abs(dy):
                    self._pan_axis(vb, 0, dx)
                else:
                    self._pan_axis(vb, 1, dy)
            return True
        return super().eventFilter(obj, event)

    @staticmethod
    def _pan_axis(vb, axis, delta):
        """Pan a viewbox along one axis (0=x, 1=y) by a fraction of the visible
        range proportional to the wheel delta. Positive delta moves the window
        toward higher values."""
        if not delta:
            return
        lo, hi = vb.viewRange()[axis]
        span = hi - lo
        shift = span * 0.10 * (delta / 120.0)
        if axis == 0:
            vb.setXRange(lo + shift, hi + shift, padding=0)
        else:
            vb.setYRange(lo + shift, hi + shift, padding=0)

    def _on_p1_clicked(self, event):
        # Double-click re-fits panel 1's intensity (y) axis to the visible data,
        # grounded at 0.
        if not event.double():
            return
        ymax = 0.0
        for scatter, _base in getattr(self, "_p1_scatters", []):
            try:
                ys = scatter.getData()[1]
                if ys is not None and len(ys):
                    ymax = max(ymax, float(np.nanmax(ys)))
            except Exception:
                pass
        if ymax > 0:
            self.p1_2d.getViewBox().setYRange(0.0, ymax * 1.05, padding=0)

    def _on_p3_clicked(self, event):
        # Double-click resets panel 3's zoom: auto-range m/z, ground intensity at 0.
        if not event.double():
            return
        vb = self.p3.getViewBox()
        vb.enableAutoRange(axis=vb.XAxis)
        (x0, x1) = vb.viewRange()[0]
        # Re-ground the y baseline at 0 against the data's max in view.
        try:
            ymax = 0.0
            for item in self.p3.getPlotItem().listDataItems():
                yd = item.yData
                if yd is not None and len(yd):
                    ymax = max(ymax, float(np.nanmax(yd)))
            if ymax > 0:
                vb.setYRange(0.0, ymax * 1.05, padding=0)
        except Exception:
            vb.enableAutoRange(axis=vb.YAxis)

    def draw_panel3_ms1(self, cur, scan_mz, scan_int):
        # Table 2 (MS2 candidates) only belongs to the MS2 view; hide it so the
        # MS1 panel 3 takes the full panel-3 + table-2 space.
        self.dock_table2.hide()
        # If this match maps to a distribution with charge states, show the
        # charge-comparison grid; otherwise fall back to the isotope overlay.
        dist_id = cur.get("distribution_id")
        if dist_id is not None:
            group = self.db.charge_group(dist_id)
            if group:
                # Skip rebuilding the grid (and the flicker) when it's already
                # showing this same distribution -- e.g. on a profile<->centroid
                # toggle or an in-window reload, where the grid data is unchanged.
                if not (self.p3_stack.currentIndex() == 1
                        and getattr(self, "_grid_dist_id", None) == dist_id):
                    self.draw_charge_grid(group, cur)
                    self._grid_dist_id = dist_id
                self.p3_stack.setCurrentIndex(1)
                self._sync_panel3_tab()
                return
        self._grid_dist_id = None
        self.p3_stack.setCurrentIndex(0)
        self.p3.clear()
        self.p3.setLabel("bottom", "m/z")
        plain = plain_seq(cur["row"].get("peptide", ""))
        charge = cur["charge"] or 1
        title = f"MS1 isotope envelope - {plain} z={charge}"

        exp_peak = 1.0
        if scan_mz is not None and len(scan_int):
            exp_peak = float(np.max(scan_int)) or 1.0
            plot_spectrum(self.p3, scan_mz, scan_int, title=title,
                          color=palette(self.theme)["fg"])

        if plain and set(plain) <= set("ACDEFGHIJKLMNPQRSTVWYUO"):
            try:
                t_mz, t_norm = isotopes.peptide_isotope_mzs(plain, charge)
                # Shift by the peptide's modification mass so the theoretical
                # envelope matches the modified precursor (same fix as panel 1).
                mod_mass = peptide_mod_mass(cur["row"].get("peptide", ""))
                if mod_mass:
                    t_mz = np.asarray(t_mz, dtype=float) + mod_mass / max(1, int(charge))
                t_y = t_norm * exp_peak
                x = np.empty(t_mz.size * 3); y = np.empty(t_mz.size * 3)
                x[0::3] = t_mz; x[1::3] = t_mz; x[2::3] = np.nan
                y[0::3] = 0.0; y[1::3] = t_y; y[2::3] = np.nan
                self.p3.plot(x, y, pen=pg.mkPen("#e85d58", width=2))
                self.p3_title.setText(title + "  (red = theoretical)")
            except Exception as exc:
                self.p3_title.setText(f"{title}  (theory failed: {exc})")
        else:
            self.p3_title.setText(title)
        self._sync_panel3_tab()

    def render_table1(self, cur):
        cur["distribution_id"] = None
        if self.db is None or cur["rt"] is None or not cur["charge"]:
            self._fill_table1([])
            return
        window = {
            "mz_min": cur["mz_center"] - self.mz_half,
            "mz_max": cur["mz_center"] + self.mz_half,
            "rt_start": cur["rt"] - self.rt_half,
            "rt_end": cur["rt"] + self.rt_half,
            "charge": cur["charge"],
        }
        # Window mode discovers the distribution id, which panel 3's charge grid
        # also needs -> let the worker result redraw panel 3 if it arrives late.
        self._start_table1_worker(window=window, redraw_panel3=True)

    def table1_for_distribution(self, distribution_id):
        """Fill table 1 with a specific distribution's line metrics (used when a
        distribution is clicked directly in panel 2)."""
        if self.db is None:
            self._fill_table1([])
            return
        self._start_table1_worker(distribution_id=distribution_id)

    def _start_table1_worker(self, distribution_id=None, window=None,
                             redraw_panel3=False):
        """Load Table 1's rows on a background thread (latest-wins by token)."""
        self._table1_token = getattr(self, "_table1_token", 0) + 1
        self._table1_redraw_panel3 = redraw_panel3
        # Show a 'loading…' placeholder row while the query runs.
        self.table1.setRowCount(1)
        self.table1.setItem(0, 0, QTableWidgetItem("loading…"))
        for j in range(1, self.table1.columnCount()):
            self.table1.setItem(0, j, QTableWidgetItem(""))
        worker = Table1Worker(str(self.db.path), self._table1_token,
                              distribution_id=distribution_id, window=window)
        worker.done.connect(self._on_table1_loaded)
        # Keep a reference so the QThread isn't garbage-collected mid-run.
        self._table1_workers = [w for w in getattr(self, "_table1_workers", [])
                                if w.isRunning()]
        self._table1_workers.append(worker)
        self._table_loading(1)
        worker.start()

    def _table_loading(self, delta):
        """Ref-counted 'loading…' indicator on the Table-1 tab row."""
        self._table_loading_count = max(
            0, getattr(self, "_table_loading_count", 0) + delta)
        if getattr(self, "table1_loading", None) is not None:
            self.table1_loading.setText(
                "<b>loading…</b>" if self._table_loading_count > 0 else "")

    def _on_table1_loaded(self, result):
        self._table_loading(-1)
        # Ignore results from a superseded selection.
        if result.get("token") != getattr(self, "_table1_token", 0):
            return
        did = result.get("distribution_id")
        if self.current is not None:
            self.current["distribution_id"] = did
        self._fill_table1(result.get("rows", []))
        # Auto-select the distribution this peptide was identified as: draw its
        # selection border on panel 2 and re-normalize panel 1's theoretical MS1
        # overlay to it. (Previously the id was only recorded, never selected,
        # so a plain list click showed neither the highlight nor the overlay at
        # the distribution's scale.)
        if did is not None and self.db is not None:
            self._set_selected(did)
            self._draw_ms1_theo_overlay()
            # Now that the distribution (and its bbox) is known, refresh the
            # MS1/MS2 tab enable + scan dropdown for the scans on it.
            self._sync_panel3_tab()
        # If panel 3 drew its MS1 view before the id was known (the sqlite query
        # finished after the mzML read), rebuild it now that we have the id.
        if (self._table1_redraw_panel3 and self._panel3_mode == "ms1"
                and did != getattr(self, "_grid_dist_id", None)):
            self._redraw_panel3()

    def _fill_table1(self, rows):
        self.table1.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, (field, _, kind) in enumerate(LINE_METRIC_COLUMNS):
                self.table1.setItem(i, j, QTableWidgetItem(_fmt(row.get(field), kind)))

    # ---- table-1 tabs: lines / distributions / charge distributions -------

    def _eager_load_table1_tabs(self):
        """Populate the distributions / charge-distributions tabs in the
        background on open (not lazily), showing 'loading…' until they arrive."""
        if self.db is None or getattr(self, "dists_model", None) is None:
            return
        self._table1_tab_token = getattr(self, "_table1_tab_token", 0) + 1
        # Mark them loaded so a tab click doesn't ALSO trigger the synchronous
        # lazy path while the background load is in flight.
        self._table1_loaded.update({"distributions", "charge distributions"})
        self._table1_tab_workers = []
        for kind, model in (("distributions", self.dists_model),
                            ("charge", self.charge_model)):
            model.set_loading(True)
            worker = DbTableWorker(str(self.db.path), kind, self._table1_tab_token)
            worker.done.connect(self._on_table1_tab_loaded)
            self._table1_tab_workers.append(worker)
            self._table_loading(1)
            worker.start()

    def _on_table1_tab_loaded(self, result):
        self._table_loading(-1)
        if result.get("token") != getattr(self, "_table1_tab_token", 0):
            return
        rows = result.get("rows", [])
        if result.get("kind") == "distributions":
            self.dists_model.set_rows(rows)
        elif result.get("kind") == "charge":
            self.charge_model.set_rows(rows)

    def reset_table1_tabs(self):
        """Reload tab data in the background (call when the file/db changes)."""
        self._table1_loaded = set()
        self._eager_load_table1_tabs()

    def _reload_current_tab(self):
        """'All' button: reload the full list for whichever tab is active."""
        if self.db is None or getattr(self, "table1_tabs", None) is None:
            return
        name = self.table1_tabs.tabText(self.table1_tabs.currentIndex())
        if name == "distributions":
            self.dists_model.set_rows(self.db.all_distributions())
            self._table1_loaded.add("distributions")
        elif name == "charge distributions":
            self.charge_model.set_rows(self.db.all_analytes_multicharge())
            self._table1_loaded.add("charge distributions")

    def _on_table1_tab_changed(self, index):
        if getattr(self, "table1_tabs", None) is None or self.db is None:
            return
        name = self.table1_tabs.tabText(index)
        if name in self._table1_loaded:
            return
        if name == "distributions":
            self.dists_model.set_rows(self.db.all_distributions())
        elif name == "charge distributions":
            self.charge_model.set_rows(self.db.all_analytes_multicharge())
        else:
            return
        self._table1_loaded.add(name)

    def _jump_to(self, mz_center, rt, distribution_id=None):
        """Centre panels 1 & 2 on (m/z, RT) and optionally select a distribution."""
        if not self.current_file or mz_center is None or rt is None:
            return
        filename = self.current_file
        if not isinstance(self.current, dict) or self.current.get("filename") != filename:
            self.current = {
                "row": {}, "filename": filename, "scan": "",
                "charge": None, "neutral_mass": None, "rt": rt,
                "targets": [], "mz_center": mz_center,
                "centroid": self.centroid_store(filename),
                "profile": self.profile_store(filename),
                "points": self.points_store(filename),
            }
        self.current["mz_center"] = mz_center
        self.current["rt"] = rt
        self.current["distribution_id"] = distribution_id
        self._panel3_mode = "ms1"
        self._ms2_scan = None
        # Select it first so the dotted border + bbox are ready when panel 2
        # redraws after the window move.
        if distribution_id is not None:
            self._set_selected(distribution_id)
        else:
            self._clear_selection()
        rt_start = max(0.0, rt - self.rt_half)
        rt_end = rt + self.rt_half
        self.set_window([mz_center - self.mz_half, mz_center + self.mz_half, rt_start, rt_end],
                        set_view=True)
        if distribution_id is not None:
            try:
                self.table1_for_distribution(distribution_id)
            except Exception:
                pass
            self.draw_panel3_ms1(self.current, getattr(self, "_last_scan_mz", None),
                                 getattr(self, "_last_scan_int", None))

    def _tab_index(self, name):
        for i in range(self.table1_tabs.count()):
            if self.table1_tabs.tabText(i) == name:
                return i
        return -1

    def _on_distribution_activated(self, row):
        # double-click a distribution -> show its lines (current tab) and jump
        if not row:
            return
        self.table1_tabs.setCurrentIndex(self._tab_index("current"))
        self._jump_to(row.get("mono_mz"), row.get("rt_apex"), row.get("distribution_id"))

    def _on_charge_activated(self, row):
        # double-click a charge distribution -> open ALL its member distributions
        # in the distributions tab and jump panels 1/2 to the analyte.
        if not row or self.db is None:
            return
        members = self.db.analyte_distributions(row.get("analyte_id"))
        self.dists_model.set_rows(members)
        self.table1_tabs.setCurrentIndex(self._tab_index("distributions"))
        # jump to a representative member (highest n_members, else first)
        rep = max(members, key=lambda d: d.get("n_members", 0)) if members else None
        if rep is not None:
            self._jump_to(rep.get("mono_mz"), rep.get("rt_apex"), rep.get("distribution_id"))

    # single click in a list -> load that distribution into panel 3 (no window
    # move; double-click still jumps panels 1/2).
    def _ensure_current_for_file(self):
        if not self.current_file:
            return False
        if not isinstance(self.current, dict) or self.current.get("filename") != self.current_file:
            self.current = {
                "row": {}, "filename": self.current_file, "scan": "",
                "charge": None, "neutral_mass": None, "rt": None,
                "targets": [], "mz_center": 500.0,
                "centroid": self.centroid_store(self.current_file),
                "profile": self.profile_store(self.current_file),
                "points": self.points_store(self.current_file),
            }
        return True

    def _show_distribution_in_panel3(self, distribution_id):
        if distribution_id is None or self.db is None or not self._ensure_current_for_file():
            return
        self._set_selected(distribution_id)
        self._panel3_mode = "ms1"
        self._ms2_scan = None
        try:
            self.table1_for_distribution(distribution_id)
        except Exception:
            pass
        self.draw_panel3_ms1(self.current, getattr(self, "_last_scan_mz", None),
                             getattr(self, "_last_scan_int", None))

    def _on_dist_clicked(self, index):
        row = self.dists_model.row_dict(self.dists_view.model().mapToSource(index).row())
        if row:
            # loads panel 3 + fills the 'current' tab with this distribution's lines
            self._show_distribution_in_panel3(row.get("distribution_id"))

    def _on_charge_clicked(self, index):
        row = self.charge_model.row_dict(self.charge_view.model().mapToSource(index).row())
        if row:
            self._show_distribution_in_panel3(row.get("rep_distribution_id"))

    # ---- panel 3 charge-comparison grid ----------------------------------

    # 8 rows, one column per charge state (port of the charge-state comparison
    # plot). Labels match the original axis labels.
    CHARGE_ROW_LABELS = ["retention time", "peak area", "charge distances",
                         "cross-charge", "intensity sum %", "adjacency",
                         "ppm to mean", "ppm error"]

    @staticmethod
    def _align(a, b):
        """Align two base-mass arrays by nearest match; returns (ai, bi, size)
        index offsets into a and b and the overlapping length (per the original
        argwhere-min alignment)."""
        if a.size == 0 or b.size == 0:
            return 0, 0, 0
        basediffs = np.abs(b - a[:, None])
        loc = np.argwhere(basediffs == basediffs.min())[0].tolist()
        minindex = min(loc)
        loc = [i - minindex for i in loc]
        amax = a.size - loc[0]
        bmax = b.size - loc[1]
        return loc[0], loc[1], max(0, min(amax, bmax))

    def draw_charge_grid(self, group, cur):
        """Charge-comparison grid: columns = the analyte's charge states, rows =
        the eight metrics from the charge-state-determination plot. Built on the
        sqlite features (mz_mean / height / area) plus the window raw points for
        the retention-time row."""
        self.p3_grid.clear()
        charges = sorted(group)
        self.p3_title.setText(
            f"charge comparison - {plain_seq(cur['row'].get('peptide',''))}  "
            f"analyte charges: {', '.join(map(str, charges))}")

        proton = isotopes.proton
        cols = {c: i for i, c in enumerate(charges)}
        # One opaque base colour per charge (its distribution colour). Rows that
        # compare against other charges (cross-charge, ppm error) colour their
        # bars by the OTHER charge's colour, like the reference `cols[nc]`.
        charge_color = {
            c: (self.distribution_color(group[c]["distribution"]["distribution_id"])
                if group[c].get("distribution") else DIST_PALETTE[cols[c] % len(DIST_PALETTE)])
            for c in charges
        }
        BAR_W = 0.03   # uniform bar width (m/z units) for the single-series rows

        # Per-charge arrays (sorted by isotope index).
        masses, intens, areas, bases = {}, {}, {}, {}
        for c in charges:
            feats = sorted(group[c]["features"], key=lambda f: f.get("isotope_index", 0))
            m = np.array([f.get("mz_mean", 0.0) for f in feats], dtype=float)
            masses[c] = m
            intens[c] = np.array([f.get("height", 0.0) for f in feats], dtype=float)
            areas[c] = np.array([f.get("area", 0.0) for f in feats], dtype=float)
            bases[c] = m * c - proton * c

        # Cross-charge mean base masses + summed intensities (size-ordered align).
        order = sorted(charges, key=lambda c: -bases[c].size)
        moving = bases[order[0]].copy()
        size0 = moving.size
        arraysums = np.zeros(size0)
        intsums = np.zeros(size0)
        divs = np.zeros(size0)
        for c in order:
            ai, bi, size = self._align(bases[c], moving)
            if size <= 0:
                continue
            arraysums[bi:bi + size] += bases[c][ai:ai + size]
            intsums[bi:bi + size] += intens[c][ai:ai + size]
            divs[bi:bi + size] += 1
        divs[divs == 0] = 1
        arraymeans = arraysums / divs

        points = getattr(self, "_last_points", None)
        self._grid_cells = {}
        row_first = {}   # row index -> leftmost cell (owns the shared y-axis labels)
        fg = palette(self.theme)["fg"]

        last_row = len(self.CHARGE_ROW_LABELS) - 1

        def cell(ri, ci):
            p = self.p3_grid.addPlot(row=ri, col=ci)
            p.showGrid(x=True, y=True, alpha=0.2)
            p.hideButtons()
            vb = p.getViewBox()
            # In-plot drag/scroll zooms the m/z (x) axis only; intensity (y) is
            # zoomed by scrolling over the y-axis (handled in eventFilter).
            vb.setMouseEnabled(x=True, y=False)
            if ri in self.GRID_ZERO_ROWS:
                vb.setLimits(yMin=0)   # baseline can never separate from 0
            left = p.getAxis("left"); bottom = p.getAxis("bottom")
            for ax in (left, bottom):
                ax.setPen(pg.mkPen(fg)); ax.setTextPen(pg.mkPen(fg))
                ax.enableAutoSIPrefix(False)   # no "1e6"-style SI prefixes
            if ri == 0:
                # Just the charge per column (the iso-score / ambiguous badge were
                # dropped per the user).
                p.setTitle(f"z={charges[ci]}", color=fg)

            # Y axis: only the leftmost column carries the axis (labels + values);
            # the other columns get a ZERO-width axis so there's no empty gap
            # between plots. Equal stretch over the remaining space keeps the
            # plots themselves the same width, and the left labels stay put.
            if ci == 0:
                left.setWidth(54)
                row_first[ri] = p
                p.setLabel("left", self.CHARGE_ROW_LABELS[ri], color=fg)
                # Ticks are inset from the view edges (pyqtgraph clips edge labels,
                # which is why only one value showed). 2 values (top, middle) for a
                # zero-baseline bar row; 3 (top, middle, bottom) otherwise.
                zero = ri in self.GRID_ZERO_ROWS

                def _ticks(*_args, ax=left, vb=vb, zero=zero):
                    (y0, y1) = vb.viewRange()[1]
                    span = y1 - y0
                    if span <= 0:
                        return
                    top = y0 + 0.92 * span
                    mid = y0 + 0.50 * span
                    ticks = [(top, _gtick(top)), (mid, _gtick(mid))]
                    if not zero:
                        ticks.append((y0 + 0.08 * span, _gtick(y0 + 0.08 * span)))
                    ax.setTicks([ticks])

                vb.sigYRangeChanged.connect(_ticks)
                p._tick_updater = _ticks
            else:
                left.setStyle(showValues=False)
                left.setWidth(0)                  # no empty axis -> no inter-column gap
                if ri in row_first:
                    p.setYLink(row_first[ri])

            # X axis: only the bottom row shows m/z values (3 ticks, no squish);
            # every other row hides its x labels to reclaim vertical space.
            if ri == last_row:
                def _xticks(*_args, ax=bottom, vb=vb):
                    (x0, x1) = vb.viewRange()[0]
                    xm = (x0 + x1) / 2.0
                    ax.setTicks([[(x0, f"{x0:.2f}"), (xm, f"{xm:.2f}"), (x1, f"{x1:.2f}")]])

                vb.sigXRangeChanged.connect(_xticks)
                p._xtick_updater = _xticks
            else:
                bottom.setStyle(showValues=False)
                bottom.setHeight(6)         # reclaim vertical space
            self._grid_cells[(ri, ci)] = p
            return p

        for c in charges:
            ci = cols[c]
            cmz = masses[c]
            cint = intens[c]
            carea = areas[c]
            cbase = bases[c]
            color = charge_color[c]

            try:
                # Row 0: retention time -- raw points (m/z vs rt) per line. Other
                # charge states sit at a different m/z than the loaded window, so
                # fetch each charge's own points directly (this is why the left-
                # most charge's RT plot used to come up empty).
                p0 = cell(0, ci)
                cpoints = self._features_points(group[c]["features"], cur)
                if isinstance(cpoints, dict) and cpoints["mz"].size:
                    for feat in group[c]["features"]:
                        m = ((cpoints["mz"] >= feat["mz_min"]) & (cpoints["mz"] <= feat["mz_max"]) &
                             (cpoints["rt"] >= feat["rt_start"]) & (cpoints["rt"] <= feat["rt_end"]))
                        if m.any():
                            o = np.argsort(cpoints["rt"][m])
                            p0.plot(cpoints["mz"][m][o], cpoints["rt"][m][o],
                                    pen=pg.mkPen(*color, width=1.5))
                            p0.addItem(pg.ScatterPlotItem(x=cpoints["mz"][m], y=cpoints["rt"][m],
                                                          size=4, pen=None, brush=pg.mkBrush(*color)))

                # Row 1: peak area -- linear bars rising from the 0 baseline.
                p1 = cell(1, ci)
                if cmz.size:
                    p1.addItem(pg.BarGraphItem(x=cmz, height=carea, y0=0.0,
                                               width=BAR_W, brush=pg.mkBrush(*color)))

                # Row 2: charge distances = diff(mz) * charge.
                p2 = cell(2, ci)
                if cmz.size > 1:
                    mids = cmz[:-1] + np.diff(cmz) / 2
                    p2.addItem(pg.BarGraphItem(x=mids, height=np.diff(cmz) * c,
                                               width=BAR_W, brush=pg.mkBrush(*color)))

                # Row 3: cross-charge intensity ratios vs every other charge.
                # Linear bars from the 0 baseline; each bar coloured by the OTHER
                # charge it compares to (ref cols[nc]), sub-divided side-by-side.
                p3 = cell(3, ci)
                nch = len(charges)
                subw = BAR_W / max(nch, 1)
                for nc in charges:
                    ai, bi, size = self._align(cbase, bases[nc])
                    if size <= 0:
                        continue
                    denom = cint[ai:ai + size]
                    ratio = np.divide(intens[nc][bi:bi + size], denom,
                                      out=np.ones(size), where=denom > 0)
                    bx = cmz[ai:ai + size] + (cols[nc] - (nch - 1) / 2.0) * subw
                    p3.addItem(pg.BarGraphItem(x=bx, height=ratio, y0=0.0,
                                               width=subw, brush=pg.mkBrush(*charge_color[nc])))

                # Row 4: intensity sum % -- linear bars from the 0 baseline.
                p4 = cell(4, ci)
                n = min(cint.size, intsums.size)
                if n:
                    denom = np.where(intsums[:n] > 0, intsums[:n], 1.0)
                    p4.addItem(pg.BarGraphItem(x=cmz[:n], height=cint[:n] / denom, y0=0.0,
                                               width=BAR_W, brush=pg.mkBrush(*color)))

                # Row 5: adjacency = signed adjacent intensity ratios (symlog-ish).
                p5 = cell(5, ci)
                if cint.size > 1:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        cr = cint[:-1] / cint[1:]
                    cr[~np.isfinite(cr)] = 0
                    cr[cr < 1] = -1.0 / np.where(cr[cr < 1] == 0, 1, cr[cr < 1])
                    mids = cmz[:-1] + np.diff(cmz) / 2
                    p5.addItem(pg.BarGraphItem(x=mids, height=cr, width=BAR_W, brush=pg.mkBrush(*color)))

                # Row 6: ppm of base mass to cross-charge mean (own colour).
                p6 = cell(6, ci)
                ai, bi, size = self._align(cbase, arraymeans)
                if size > 0:
                    base = cbase[ai:ai + size]
                    mean = arraymeans[bi:bi + size]
                    ppm = (mean - base) / np.where(base == 0, 1, base) * 1e6
                    p6.addItem(pg.BarGraphItem(x=cmz[ai:ai + size], height=ppm,
                                               width=BAR_W, brush=pg.mkBrush(*color)))

                # Row 7: ppm error vs each other charge's base mass, coloured by
                # the OTHER charge (ref cols[nc]); bars sub-divided side-by-side.
                p7 = cell(7, ci)
                others = [nc for nc in charges if nc != c]
                subw7 = BAR_W / max(len(others), 1)
                for k, nc in enumerate(others):
                    ai, bi, size = self._align(cbase, bases[nc])
                    if size <= 0:
                        continue
                    base = cbase[ai:ai + size]
                    other = bases[nc][bi:bi + size]
                    ppm = (other - base) / np.where(other == 0, 1, other) * 1e6
                    bx = cmz[ai:ai + size] + (k - (len(others) - 1) / 2.0) * subw7
                    p7.addItem(pg.BarGraphItem(x=bx, height=ppm,
                                               width=subw7, brush=pg.mkBrush(*charge_color[nc])))
            except Exception as exc:
                # one bad column shouldn't blank the whole grid
                import traceback
                traceback.print_exc()

        # Each column shares one m/z (x) axis -> link all rows in a column to the
        # top cell so panning/zooming mass moves them together (per-column, since
        # different charges have different m/z scales).
        for ci in range(len(charges)):
            top = self._grid_cells.get((0, ci))
            if top is None:
                continue
            for ri in range(1, len(self.CHARGE_ROW_LABELS)):
                c = self._grid_cells.get((ri, ci))
                if c is not None:
                    c.setXLink(top)
        # Each ROW shares one y axis (matplotlib sharey='row'): link every column
        # to the left-most cell, which is the only one showing values, so the left
        # y-axis represents the whole row.
        for ri in range(len(self.CHARGE_ROW_LABELS)):
            left = self._grid_cells.get((ri, 0))
            if left is None:
                continue
            for ci in range(1, len(charges)):
                c = self._grid_cells.get((ri, ci))
                if c is not None:
                    c.setYLink(left)
        self._grid_ncols = len(charges)
        # Force every charge column to the SAME width (equal stretch), so no
        # column's plots are wider than another's.
        try:
            glayout = self.p3_grid.ci.layout
            glayout.setHorizontalSpacing(0)   # tight gaps -> reclaim blank space
            glayout.setVerticalSpacing(1)
            glayout.setContentsMargins(2, 2, 2, 2)
            for ci in range(len(charges)):
                glayout.setColumnStretchFactor(ci, 1)
                glayout.setColumnPreferredWidth(ci, 0)
        except Exception:
            pass
        # Fit synchronously from the items' data (no deferred pass), so the grid
        # opens directly in the right state -- the deferred reset was causing a
        # one-frame flicker through the wrong auto-fit.
        # Remove the outline border on every bar and force full opacity so the
        # true colour reads cleanly at any zoom.
        for p in self._grid_cells.values():
            for it in p.items:
                if isinstance(it, pg.BarGraphItem):
                    it.setOpts(pen=None)
        self._fit_grid_cols()
        self._fit_grid_rows()
        for p in self._grid_cells.values():
            if hasattr(p, "_tick_updater"):
                p._tick_updater()
            if hasattr(p, "_xtick_updater"):
                p._xtick_updater()

    # Non-negative bar rows glued to a 0 baseline (peak area, cross-charge,
    # intensity sum %): y starts at 0 at the bottom, like the MS2 spectrum.
    GRID_ZERO_ROWS = {1, 3, 4}

    @staticmethod
    def _grid_item_yvals(plotitem):
        vals = []
        for it in plotitem.items:
            if isinstance(it, pg.BarGraphItem):
                h = it.opts.get("height")
                if h is not None:
                    vals.append(np.asarray(h, dtype=float).ravel())
            else:
                yd = getattr(it, "yData", None)
                if yd is not None and len(yd):
                    vals.append(np.asarray(yd, dtype=float).ravel())
        return vals

    @staticmethod
    def _grid_item_xvals(plotitem):
        vals = []
        for it in plotitem.items:
            if isinstance(it, pg.BarGraphItem):
                x = it.opts.get("x")
                if x is not None:
                    vals.append(np.asarray(x, dtype=float).ravel())
            else:
                xd = getattr(it, "xData", None)
                if xd is not None and len(xd):
                    vals.append(np.asarray(xd, dtype=float).ravel())
        return vals

    def _fit_grid_cols(self):
        """Fit each column's x-axis (m/z) to its own data, on the column's top
        (master) cell the rows are x-linked to. Deterministic = no flicker."""
        for ci in range(getattr(self, "_grid_ncols", 0)):
            top = self._grid_cells.get((0, ci))
            if top is None:
                continue
            vals = []
            for ri in range(len(self.CHARGE_ROW_LABELS)):
                p = self._grid_cells.get((ri, ci))
                if p is not None:
                    vals += self._grid_item_xvals(p)
            if not vals:
                continue
            allv = np.concatenate(vals)
            allv = allv[np.isfinite(allv)]
            if allv.size == 0:
                continue
            lo, hi = float(allv.min()), float(allv.max())
            if hi - lo < 1e-9:
                lo -= 0.5; hi += 0.5
            pad = (hi - lo) * 0.05
            top.getViewBox().setXRange(lo - pad, hi + pad, padding=0)

    def _fit_grid_rows(self):
        """Fit each row's y-axis to the UNION of all its columns' data, applied to
        the left (master) cell that the others are y-linked to. This is why both
        ppm rows (positive AND negative values) now fit on screen, and it's
        deterministic so double-click doesn't oscillate between two auto-fits."""
        for ri in range(len(self.CHARGE_ROW_LABELS)):
            master = self._grid_cells.get((ri, 0))
            if master is None:
                continue
            vals = []
            for ci in range(getattr(self, "_grid_ncols", 0)):
                p = self._grid_cells.get((ri, ci))
                if p is not None:
                    vals += self._grid_item_yvals(p)
            if not vals:
                continue
            allv = np.concatenate(vals)
            allv = allv[np.isfinite(allv)]
            if allv.size == 0:
                continue
            hi = float(allv.max())
            if ri in self.GRID_ZERO_ROWS:
                # Glue the baseline to 0 at the bottom; only headroom on top.
                lo = 0.0
                if hi <= 0:
                    hi = 1.0
                master.getViewBox().setYRange(lo, hi * 1.08, padding=0)
            else:
                lo = float(allv.min())
                if hi - lo < 1e-9:
                    lo -= 1.0; hi += 1.0
                pad = (hi - lo) * 0.08
                master.getViewBox().setYRange(lo - pad, hi + pad, padding=0)

    def _grid_axis_wheel(self, event):
        """Wheel over a grid cell's left (y) axis zooms that cell's intensity.
        For the 0-baseline rows the bottom stays pinned at 0 (only the top moves,
        so it 'goes down' rather than zooming symmetrically)."""
        try:
            scene_pos = self.p3_grid.mapToScene(event.position().toPoint())
        except Exception:
            return False
        factor = 0.85 if event.angleDelta().y() > 0 else 1.0 / 0.85
        for (ri, ci), p in getattr(self, "_grid_cells", {}).items():
            ax = p.getAxis("left")
            if not ax.isVisible() or not ax.sceneBoundingRect().contains(scene_pos):
                continue
            vb = p.getViewBox()
            (y0, y1) = vb.viewRange()[1]
            if ri in self.GRID_ZERO_ROWS:
                vb.setYRange(0.0, max(y1 * factor, 1e-12), padding=0)
            else:
                c = (y0 + y1) / 2.0
                h = (y1 - y0) / 2.0 * factor
                vb.setYRange(c - h, c + h, padding=0)
            return True
        return False

    def reset_charge_grid_zoom(self):
        # Deterministic reset (same as the initial fit): no flicker/oscillation.
        self._fit_grid_cols()
        self._fit_grid_rows()

    def _on_grid_clicked(self, event):
        if event.double():
            self.reset_charge_grid_zoom()

    # ---- toggles + sync --------------------------------------------------

    def toggle_dimension(self):
        to_3d = self.dim_toggle.key() == "3D"
        self.p1_stack.setCurrentIndex(1 if to_3d else 0)
        self.reset3d_button.setEnabled(to_3d and HAVE_GL)
        # Build the 3D view lazily on first switch (it's skipped while 2D shows).
        if to_3d and HAVE_GL and getattr(self, "_p1_3d_inputs", None) is not None:
            self.draw_panel1_3d(*self._p1_3d_inputs)
        # The theoretical MS1 overlay is 2D-only.
        self._draw_ms1_theo_overlay()

    def _recenter_window(self):
        if self.center is None:
            return
        mz_c, rt_c = self.center
        rt_start = max(0.0, rt_c - self.rt_half) if rt_c is not None else 0.0
        rt_end = rt_c + self.rt_half if rt_c is not None else 1.0
        self.set_window([mz_c - self.mz_half, mz_c + self.mz_half, rt_start, rt_end], set_view=True)

    def set_mz_half(self, value):
        self.mz_half = float(value)
        self._recenter_window()

    def set_rt_half(self, value):
        self.rt_half = float(value)
        self._recenter_window()

    # ---- distribution selection + dotted border --------------------------

    def _set_selected(self, distribution_id, keep_group=False):
        """Mark a distribution as selected: update the charge-search anchor
        (current charge/mass/RT) to it, cache its analyte's charge group, and
        compute its bounding box for the dotted border. Does not move the view."""
        if self.db is None or distribution_id is None:
            return
        dist = self.db.distribution(distribution_id)
        if dist is None:
            return
        members = self.db.distribution_members(distribution_id)
        if isinstance(self.current, dict):
            self.current["distribution_id"] = distribution_id
            self.current["charge"] = dist.get("charge")
            self.current["neutral_mass"] = dist.get("neutral_mass")
            self.current["rt"] = dist.get("rt_apex")
        self.assumed_charge = dist.get("charge")
        self._border_color = None   # a real distribution -> normal (fg) border
        if members:
            self._selected_bbox = (
                min(m["mz_min"] for m in members),
                max(m["mz_max"] for m in members),
                min(m["rt_start"] for m in members),
                max(m["rt_end"] for m in members),
            )
            self._selected_rt_band = (self._selected_bbox[2], self._selected_bbox[3])
            self._selected_n_members = len(members)
        else:
            self._selected_bbox = None
        if not keep_group:
            self._selected_charge_group = self.db.charge_group(distribution_id)
        self._selected_dist_id = distribution_id
        self._render_selection_border()

    def _set_hypothetical(self, charge):
        """Show a RED dotted box where the distribution at `charge` *would* sit if
        it existed: the theoretical isotope m/z span at the selection's neutral
        mass, over the analyte's RT band. Auto-fits to it like a real lock-on."""
        nm = self.current.get("neutral_mass") if isinstance(self.current, dict) else None
        if nm is None:
            self._clear_selection()
            self._recenter_for_charge()
            return
        n = self._selected_n_members or 6
        mzs = isotope_mzs(nm, max(1, charge), n=n)
        if not mzs:
            self._clear_selection()
            self._recenter_for_charge()
            return
        if self._selected_rt_band is not None:
            y0, y1 = self._selected_rt_band
        elif self.window is not None:
            y0, y1 = self.window[2], self.window[3]
        else:
            y0, y1 = 0.0, 1.0
        self._selected_bbox = (min(mzs), max(mzs), y0, y1)
        self._border_color = "red"
        self._selected_dist_id = None
        self._fit_to_selected()

    def _clear_selection(self):
        self._selected_bbox = None
        self._render_selection_border()

    def _render_selection_border(self):
        """Draw (or refresh) the dotted rectangle around the selected
        distribution on panel 2. A standalone item so it survives without a full
        redraw; draw_panel2 re-adds it after its clear()."""
        if getattr(self, "_sel_border_item", None) is not None:
            try:
                self.p2.removeItem(self._sel_border_item)
            except Exception:
                pass
            self._sel_border_item = None
        bbox = getattr(self, "_selected_bbox", None)
        if bbox is None:
            return
        x0, x1, y0, y1 = bbox
        px = (x1 - x0) * 0.04 + 1e-4
        py = (y1 - y0) * 0.04 + 1e-4
        xs = np.array([x0 - px, x1 + px, x1 + px, x0 - px, x0 - px])
        ys = np.array([y0 - py, y0 - py, y1 + py, y1 + py, y0 - py])
        color = (230, 70, 70) if self._border_color == "red" else palette(self.theme)["fg"]
        pen = pg.mkPen(color=color, width=1.4, style=Qt.DashLine)
        item = pg.PlotCurveItem(xs, ys, pen=pen)
        item.setZValue(50)
        self.p2.addItem(item)
        self._sel_border_item = item

    def _fit_to_selected(self):
        """Auto-fit panels 1 & 2 to the selected distribution's bounding box so it
        sits centred and fully framed (the charge-search lock-on)."""
        bbox = getattr(self, "_selected_bbox", None)
        if bbox is None:
            return
        x0, x1, y0, y1 = bbox
        mzpad = (x1 - x0) * 0.18 + 0.05
        rtpad = (y1 - y0) * 0.30 + 0.02
        self.set_window([x0 - mzpad, x1 + mzpad, max(0.0, y0 - rtpad), y1 + rtpad], set_view=True)

    # ---- charge search + navigation history ------------------------------

    def _recenter_for_charge(self):
        """Recentre the m/z window on the same neutral mass at the assumed charge,
        keeping the RT window -- this is the charge-search 'lock on'."""
        cur = self.current
        if cur is None or cur.get("neutral_mass") is None:
            return
        z = max(1, self.assumed_charge or 1)
        mzs = isotope_mzs(cur["neutral_mass"], z, n=6)
        center = sum(mzs) / len(mzs)
        if self.window is not None:
            rt_start, rt_end = self.window[2], self.window[3]
        elif cur["rt"] is not None:
            rt_start, rt_end = max(0.0, cur["rt"] - self.rt_half), cur["rt"] + self.rt_half
        else:
            rt_start, rt_end = 0.0, 1.0
        self.set_window([center - self.mz_half, center + self.mz_half, rt_start, rt_end], set_view=True)

    def charge_step(self, delta):
        if self.current is None:
            return
        step = 1 if delta > 0 else -1
        base = self.current.get("charge") or self.assumed_charge or 1
        target = max(1, base + step)
        # Follow the new charge on panel 1's theoretical overlay in BOTH cases
        # (a linked distribution or a hypothetical one): the theoretical MS1 is
        # the current peptide's isotope pattern at the target charge, i.e. it
        # sits at the presumed m/z for that charge and scales its abundances the
        # same way as any other theoretical overlay. Rebuild it from the peptide
        # if a prior action cleared it.
        self._set_theo_for_charge(target)
        grp = getattr(self, "_selected_charge_group", None)
        # If the selected distribution's analyte has a distribution at the target
        # charge, lock onto it: transfer the dotted border and auto-fit to it.
        if grp and target in grp and grp[target].get("distribution"):
            did = grp[target]["distribution"]["distribution_id"]
            self._set_selected(did, keep_group=True)
            # Fitting reloads the region (async); if the MS2 view is up, tell
            # on_evidence_done to shift panel 3 to THIS distribution's MS2
            # instead of re-rendering the old scan. (The MS1 view already
            # follows the new distribution id.)
            if self._panel3_mode == "ms2":
                self._pending_charge_refocus = True
            self._fit_to_selected()
            self._draw_ms1_theo_overlay()   # show the new-charge overlay at once
            return
        # Otherwise no distribution exists at this charge: show a RED box where the
        # hypothetical distribution would be, and fit to it.
        self.assumed_charge = target
        if isinstance(self.current, dict):
            self.current["charge"] = target
        self._set_hypothetical(target)
        # No distribution here, so no async reload is guaranteed to redraw the
        # overlay -- draw it now at the hypothetical charge's presumed m/z.
        self._draw_ms1_theo_overlay()

    def _set_theo_for_charge(self, charge):
        """Point the panel-1 theoretical MS1 overlay at the current peptide at
        ``charge`` (used by charge stepping). Rebuilds _ms1_theo from the current
        PSM if it was cleared, so scrolling into a charge with no distribution
        still shows a theoretical distribution at that charge's presumed m/z."""
        if isinstance(self._ms1_theo, dict):
            self._ms1_theo = {**self._ms1_theo, "charge": charge}
            return
        row = self.current.get("row", {}) if isinstance(self.current, dict) else {}
        plain = plain_seq(row.get("peptide", ""))
        if len(plain) >= 2:
            self._ms1_theo = {"seq": plain, "charge": charge,
                              "mod_mass": peptide_mod_mass(row.get("peptide", ""))}

    def set_charge(self, z):
        if self.current is None:
            return
        self.assumed_charge = max(1, int(z))
        self._recenter_for_charge()

    def _record_nav(self):
        if self.window is None:
            return
        entry = (list(self.window), self.assumed_charge)
        if self._suppress_record:
            self._suppress_record = False
            return
        if self._nav and self._nav[self._nav_i][0] == entry[0]:
            return
        self._nav = self._nav[:self._nav_i + 1]
        self._nav.append(entry)
        self._nav_i = len(self._nav) - 1

    def _apply_nav(self):
        window, charge = self._nav[self._nav_i]
        self.assumed_charge = charge
        self._suppress_record = True
        self.set_window(list(window), set_view=True)

    def nav_back(self):
        if self._nav_i > 0:
            self._nav_i -= 1
            self._apply_nav()

    def nav_forward(self):
        if self._nav_i < len(self._nav) - 1:
            self._nav_i += 1
            self._apply_nav()
