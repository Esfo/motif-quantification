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
from PySide6.QtCore import QEvent, Qt, QThread, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDockWidget,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QPushButton,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .mzml_store import MzmlStore, scan_arrays
    from .plots import plot_points, plot_spectrum
    from .region_view import HAVE_GL, gl
    from .session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
    from . import isotopes
    from .theming import palette, style_plot, style_gl
except ImportError:
    from mzml_store import MzmlStore, scan_arrays
    from plots import plot_points, plot_spectrum
    from region_view import HAVE_GL, gl
    from session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
    import isotopes
    from theming import palette, style_plot, style_gl


def plain_seq(peptide):
    """Strip flanks and modifications to a bare uppercase residue sequence."""
    value = peptide or ""
    if len(value) >= 5 and value[1] == "." and value[-2] == ".":
        value = value[2:-2]
    value = re.sub(r"\[[^\]]*\]|\([^\)]*\)|\{[^\}]*\}", "", value)
    return re.sub(r"[^A-Za-z]", "", value).upper()


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
            ms1_number = None
            if ms1 is not None:
                ms1_number = ms1.number
                scan_mz, scan_int = self.store.scan_window_by_number(ms1.number, self.mz_min, self.mz_max)

            points = None
            region = None
            ms2 = []
            if self.rt is not None:
                points = self.store.extract_points(self.mz_min, self.mz_max, self.rt_start, self.rt_end)
                region = self.store.extract_region(
                    self.mz_min, self.mz_max, self.rt_start, self.rt_end,
                    mz_bins=self.mz_bins, mode=self.mode)
                # MS2 scans always come from the centroid run (where they live).
                ms2 = [{"rt": s.rt, "mz": s.precursor_mz, "number": s.number, "id": s.spectrum_id}
                       for s in self.centroid.ms2_in_rt(self.rt_start, self.rt_end)]

            self.done.emit({"ms1_number": ms1_number, "scan_mz": scan_mz,
                            "scan_int": scan_int, "points": points, "region": region, "ms2": ms2})
        except Exception as exc:
            self.done.emit({"error": f"{exc}\n{traceback.format_exc()}"})


LINE_METRIC_COLUMNS = [
    ("isotope_index", "iso"),
    ("area", "AUC"),
    ("height", "max I"),
    ("n_points", "n pts"),
    ("rt_start", "min t"),
    ("rt_end", "max t"),
    ("mz_min", "min m/z"),
    ("mz_max", "max m/z"),
    ("mz_mean", "mean m/z"),
    ("rt_apex", "RT"),
]


def fixed_toggle(off_text, on_text, width=70):
    """A two-state button whose label flips but whose size never changes."""
    button = QPushButton(off_text)
    button.setCheckable(True)
    button.setFixedWidth(width)
    button._off = off_text
    button._on = on_text
    return button


class MSViewerTab(QMainWindow):
    def __init__(self, session, distributions_db=None, xics_ppm=10.0, xics_rt_window=0.8, theme="dark"):
        super().__init__()
        self.session = session
        self.db = distributions_db
        self.xics_ppm = float(xics_ppm)
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
        self.window = None        # [mz_min, mz_max, rt_start, rt_end] source of truth
        self.center = None        # (mz_center, rt_center) for the ± controls
        self._guard = False       # suppress range-change handling during programmatic set
        try:
            self._cmap = pg.colormap.get("viridis")
        except Exception:
            self._cmap = None

        from PySide6.QtCore import QTimer
        self._reload_timer = QTimer(self)
        self._reload_timer.setSingleShot(True)
        self._reload_timer.setInterval(140)
        self._reload_timer.timeout.connect(self.do_extract)

        self.setDockNestingEnabled(True)
        self.build_lists_dock()
        self.build_panel1_dock()
        self.build_panel2_dock()
        self.build_panel3_dock()
        self.build_table1_dock()
        self.arrange_default()
        self._default_state = self.saveState()
        self.apply_theme(theme)

    # ---- docks -----------------------------------------------------------

    def build_lists_dock(self):
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(4, 4, 4, 4)

        # Tab 1 is a single-file view: a file must be chosen first. Default to
        # the first file so the lists are never empty/ambiguous on open.
        self.file_combo = QComboBox()
        for row in self.session.files():
            name = row.get("filename", "")
            if name:
                self.file_combo.addItem(name, name)
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)
        layout.addWidget(QLabel("file"))
        layout.addWidget(self.file_combo)

        self.search = QLineEdit()
        self.search.setPlaceholderText("filter…")
        self.search.textChanged.connect(self.repopulate_active_list)
        layout.addWidget(self.search)

        self.protein_list = self._titled_list(layout, "proteins", self.on_protein_selected, self.show_all_proteins)
        self.peptide_list = self._titled_list(layout, "peptides", self.on_peptide_selected, self.show_all_peptides)
        self.psm_list = self._titled_list(layout, "PSMs", self.on_psm_selected, self.show_all_psms)

        dock = QDockWidget("Lists", self)
        dock.setObjectName("dock_lists")
        dock.setWidget(container)
        self.dock_lists = dock

        self.current_file = self.file_combo.currentData()
        self.repopulate_active_list()

    def on_file_changed(self):
        self.current_file = self.file_combo.currentData()
        self.search.clear()
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
        listw.itemSelectionChanged.connect(on_select)
        layout.addWidget(listw, stretch=1)
        return listw

    def build_panel1_dock(self):
        self.p1_2d = pg.PlotWidget()
        self.p1_2d.setLabel("bottom", "m/z")
        self.p1_2d.setLabel("left", "intensity")
        self.p1_2d.setClipToView(True)
        self.p1_2d.setDownsampling(auto=True, mode="peak")
        # 2D panel 1: only the m/z (x) axis is interactive; y stays auto-scaled.
        self.p1_2d.setMouseEnabled(x=True, y=False)
        # Wheel over the y-axis label strip (left of the plot) scrolls intensity.
        self.p1_2d.viewport().installEventFilter(self)

        if HAVE_GL:
            self.p1_3d = gl.GLViewWidget()
            # Deterministic orientation: time runs left->right, m/z front->back,
            # so the m/z axis reads the same direction as panel 2's x.
            self.p1_3d.setCameraPosition(distance=3.4, elevation=22, azimuth=-90)
            self.p1_surface = gl.GLSurfacePlotItem(
                x=np.array([0.0, 1.0], dtype=np.float32),
                y=np.array([0.0, 1.0], dtype=np.float32),
                z=np.zeros((2, 2), dtype=np.float32),
                shader="shaded", smooth=True,
            )
            self.p1_surface.setVisible(False)
            self.p1_scatter = gl.GLScatterPlotItem(pos=np.zeros((1, 3), dtype=np.float32), size=4.0)
            self.p1_scatter.setVisible(False)
            axis = gl.GLAxisItem()
            axis.setSize(2.2, 2.2, 1.2)
            self.p1_3d.addItem(axis)
            self.p1_3d.addItem(self.p1_surface)
            self.p1_3d.addItem(self.p1_scatter)
            # Axis name labels on the outside (x = time, y = m/z).
            for text, pos in (("time (min)", (1.3, -1.1, -1.0)), ("m/z", (-1.1, 1.3, -1.0))):
                try:
                    self.p1_3d.addItem(gl.GLTextItem(pos=pos, text=text, color=(220, 220, 220, 255)))
                except Exception:
                    pass
        else:
            self.p1_3d = QLabel("3D needs pyqtgraph OpenGL (pip install PyOpenGL)")
            self.p1_3d.setAlignment(Qt.AlignCenter)

        self.p1_stack = QStackedWidget()
        self.p1_stack.addWidget(self.p1_2d)   # index 0 = 2D
        self.p1_stack.addWidget(self.p1_3d)   # index 1 = 3D

        self.dim_toggle = fixed_toggle("3D", "2D")     # shows what it will switch TO
        self.dim_toggle.clicked.connect(self.toggle_dimension)
        self.source_toggle = fixed_toggle("profile", "centroid")
        self.source_toggle.clicked.connect(self.refresh)

        bar = QHBoxLayout()
        bar.addWidget(self.dim_toggle)
        bar.addWidget(self.source_toggle)
        bar.addStretch(1)

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.addLayout(bar)
        layout.addWidget(self.p1_stack, stretch=1)

        dock = QDockWidget("Panel 1 - spectrum", self)
        dock.setObjectName("dock_panel1")
        dock.setWidget(container)
        self.dock_panel1 = dock

    def build_panel2_dock(self):
        # Thin MS2 strip to the left of panel 2: MS2 scans as clickable points,
        # RT-aligned with panel 2 (shared y). It fits in the space panel 1's wide
        # y-axis labels leave between panels 1 and 2.
        self.ms2_plot = pg.PlotWidget()
        self.ms2_plot.setFixedWidth(64)
        self.ms2_plot.setMouseEnabled(x=False, y=True)
        self.ms2_plot.getPlotItem().hideAxis("bottom")
        self.ms2_plot.setLabel("left", "MS2 RT")
        self.ms2_plot.setXRange(0, 1, padding=0)
        self.ms2_scatter = pg.ScatterPlotItem(size=7, brush=pg.mkBrush(255, 180, 60, 220),
                                              pen=pg.mkPen("#222"))
        self.ms2_scatter.sigClicked.connect(self.on_ms2_clicked)
        self.ms2_plot.addItem(self.ms2_scatter)

        # Panel 2: m/z on x (aligned with panel 1), RT on y.
        self.p2 = pg.PlotWidget()
        self.p2.setLabel("bottom", "m/z")
        self.p2.setLabel("left", "RT", units="min")
        self.p2_image = pg.ImageItem()
        if self._cmap is not None:
            self.p2_image.setColorMap(self._cmap)
        self.p2.addItem(self.p2_image)

        # MS2 strip shares panel 2's RT (y) axis.
        self.ms2_plot.setYLink(self.p2)
        # panel 1 (2D) and panel 2 share the m/z (x) axis -> link them.
        self.p1_2d.setXLink(self.p2)
        self.p2.sigXRangeChanged.connect(self.on_view_range_changed)
        self.p2.sigYRangeChanged.connect(self.on_view_range_changed)

        container = QWidget()
        row = QHBoxLayout(container)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        row.addWidget(self.ms2_plot)
        row.addWidget(self.p2, stretch=1)

        dock = QDockWidget("Panel 2 - m/z x RT map (MS2 strip left)", self)
        dock.setObjectName("dock_panel2")
        dock.setWidget(container)
        self.dock_panel2 = dock

    def build_panel3_dock(self):
        self.p3 = pg.PlotWidget()
        self.p3.setLabel("bottom", "m/z")
        self.p3.setLabel("left", "intensity")
        self.p3_title = QLabel("Panel 3")
        self.p3_title.setStyleSheet("font-weight: bold;")

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.addWidget(self.p3_title)
        layout.addWidget(self.p3, stretch=1)

        dock = QDockWidget("Panel 3 - MS1 / MS2", self)
        dock.setObjectName("dock_panel3")
        dock.setWidget(container)
        self.dock_panel3 = dock

    def build_table1_dock(self):
        self.table1 = QTableWidget()
        self.table1.setColumnCount(len(LINE_METRIC_COLUMNS))
        self.table1.setHorizontalHeaderLabels([h for _, h in LINE_METRIC_COLUMNS])
        self.table1.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table1.setSelectionBehavior(QAbstractItemView.SelectRows)

        dock = QDockWidget("Table 1 - distribution lines", self)
        dock.setObjectName("dock_table1")
        dock.setWidget(self.table1)
        self.dock_table1 = dock

    def arrange_default(self):
        # Three columns: lists | [panel1 / panel2 / table1] | panel3
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_lists)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_panel1)
        self.splitDockWidget(self.dock_panel1, self.dock_panel3, Qt.Horizontal)
        self.splitDockWidget(self.dock_panel1, self.dock_panel2, Qt.Vertical)
        self.splitDockWidget(self.dock_panel2, self.dock_table1, Qt.Vertical)
        self.resizeDocks([self.dock_lists], [320], Qt.Horizontal)
        self.resizeDocks([self.dock_panel1, self.dock_panel3], [600, 460], Qt.Horizontal)

    def reset_layout(self):
        self.restoreState(self._default_state)

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)
        for plot in (self.p1_2d, self.p2, self.p3):
            style_plot(plot, pal)
        if HAVE_GL:
            style_gl(self.p1_3d, pal)

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
        t = self.search.text().strip().lower()
        return (not t) or (t in text.lower())

    # All list content is scoped to the selected file (single-file view).
    def file_psms(self):
        rows = []
        for row in self.session.load_psms(self.current_file or ""):
            row = dict(row)
            row["filename"] = self.current_file
            rows.append(row)
        return rows

    def _peptide_label(self, row):
        # Annotate LFQ-only peptides (quantified but with no PSM in this file).
        n = str(row.get("n_psms", "") or "").strip()
        pep = row.get("peptide", "")
        if n in ("", "0"):
            return f"{pep}   · LFQ-only"
        return f"{pep}   ({n})"

    def show_all_proteins(self, preserve=True):
        rows = self.session.file_proteins(self.current_file or "")
        self._fill(self.protein_list, [(r.get("protein_id", ""), r) for r in rows
                                       if r.get("protein_id") and self._filter(r["protein_id"])],
                   preserve=preserve)

    def show_all_peptides(self, preserve=True):
        rows = self.session.file_peptides(self.current_file or "")
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
        # show this protein's peptides (matched within the file's peptide table)
        matched = [r for r in self.session.file_peptides(self.current_file or "")
                   if r.get("peptide") in peptides or plain_seq(r.get("peptide", "")) in plains]
        self._fill(self.peptide_list, [(r.get("peptide", ""), r) for r in matched])

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
        }
        self.center = (mz_center, rt)
        self.render_table1(self.current)
        # Initialize the window from the ± controls and snap the views to it.
        rt_start = max(0.0, rt - self.rt_half) if rt is not None else 0.0
        rt_end = rt + self.rt_half if rt is not None else 1.0
        self.set_window([mz_center - self.mz_half, mz_center + self.mz_half, rt_start, rt_end],
                        set_view=True)

    def use_profile(self):
        return self.source_toggle.isChecked()

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
        # Read the current window from panel 2 (m/z = x, RT = y) and debounce.
        (mz0, mz1) = self.p2.getViewBox().viewRange()[0]
        (rt0, rt1) = self.p2.getViewBox().viewRange()[1]
        self.window = [mz0, mz1, max(0.0, rt0), rt1]
        self._reload_timer.start()

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
        mz_min, mz_max, rt_start, rt_end = self.window
        store = cur["profile"] if (self.use_profile() and cur["profile"]) else centroid
        self._win = (mz_min, mz_max, rt_start, rt_end)
        self._pending = dict(
            centroid=centroid, store=store, scan=cur["scan"], rt=cur["rt"],
            mz_min=mz_min, mz_max=mz_max, rt_start=rt_start, rt_end=rt_end,
            mz_bins=400, mode="profile" if (self.use_profile() and cur["profile"]) else "centroid")
        self.p3_title.setText(f"loading {cur['row'].get('peptide','')}…")
        self._start_evidence()

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
            return
        if "error" in result:
            self.p3_title.setText(f"evidence error: {result['error'].splitlines()[0]}")
            return

        mz_min, mz_max, rt_start, rt_end = self._win
        points = result.get("points")
        region = result.get("region")
        scan_mz = result.get("scan_mz")
        scan_int = result.get("scan_int")
        profile = bool(self.use_profile() and cur["profile"])

        self.draw_panel1(points, region, mz_min, mz_max, rt_start, rt_end, profile)
        self.draw_panel2(region, mz_min, mz_max, rt_start, rt_end)
        self.draw_ms2_strip(result.get("ms2", []))
        self.draw_panel3_ms1(cur, scan_mz, scan_int)

    def draw_panel1(self, points, region, mz_min, mz_max, rt_start, rt_end, profile):
        self.p1_2d.clear()
        if isinstance(points, dict) and points["mz"].size:
            self.p1_2d.showGrid(x=True, y=True, alpha=0.25)
            if profile:
                # Profile data is continuous: draw one curve per scan so the
                # envelope stays visible (and crisp) at any zoom, unlike fixed dots.
                rts = points["rt"]
                for r in np.unique(rts):
                    m = rts == r
                    order = np.argsort(points["mz"][m])
                    self.p1_2d.plot(points["mz"][m][order], points["intensity"][m][order],
                                    pen=pg.mkPen(120, 170, 255, 90))
            else:
                self.p1_2d.plot(points["mz"], points["intensity"], pen=None,
                                symbol="o", symbolSize=4, symbolPen=None,
                                symbolBrush=pg.mkBrush(120, 170, 255, 170))
            self.p1_2d.setTitle(f"{self.current['row'].get('peptide','')} z={self.current['charge']}  "
                                f"({points['mz'].size} pts)")
            self.p1_2d.getViewBox().setYRange(0, float(points["intensity"].max()) * 1.05, padding=0)
        else:
            self.p1_2d.setTitle("no MS1 points in window")
        self.draw_panel1_3d(points, region, mz_min, mz_max, rt_start, rt_end)

    def draw_panel1_3d(self, points, region, mz_min, mz_max, rt_start, rt_end):
        if not HAVE_GL:
            return
        mz_span = max(mz_max - mz_min, 1e-6)
        rt_span = max(rt_end - rt_start, 1e-6)

        if isinstance(region, dict) and region.get("z") is not None and region["z"].size:
            z = region["z"]
            rts = region["rts"]
            mz_grid = region["mz_grid"]
            zmax = float(z.max()) or 1.0
            zr = (z / zmax).astype(np.float32)            # (n_rt, n_mz)
            # Map surface to ACTUAL rt/mz (not even index spacing) so it lines up
            # with the scatter points -- fixes the profile points/surface mismatch.
            xs = ((rts - rt_start) / rt_span * 2 - 1).astype(np.float32)
            ys = ((mz_grid - mz_min) / mz_span * 2 - 1).astype(np.float32)
            try:
                self.p1_surface.setData(x=xs, y=ys, z=zr)
                self.p1_surface.setColor((0.30, 0.45, 0.70, 0.45))
                self.p1_surface.setVisible(True)
            except Exception:
                pass

        if isinstance(points, dict) and points["mz"].size:
            x = ((points["rt"] - rt_start) / rt_span * 2 - 1).astype(np.float32)
            y = ((points["mz"] - mz_min) / mz_span * 2 - 1).astype(np.float32)
            inten = points["intensity"]
            zmax = float(inten.max()) or 1.0
            z = (inten / zmax).astype(np.float32)
            pos = np.column_stack([x, y, z])
            if self._cmap is not None:
                colors = self._cmap.map(z, mode="float")
            else:
                colors = np.tile(np.array([1, 1, 1, 1], dtype=np.float32), (z.size, 1))
            try:
                self.p1_scatter.setData(pos=pos, color=colors, size=4.0)
                self.p1_scatter.setVisible(True)
            except Exception:
                pass

    def draw_panel2(self, region, mz_min, mz_max, rt_start, rt_end):
        if not isinstance(region, dict):
            return
        z = region.get("z")
        rts = region.get("rts")
        if z is None or z.size == 0 or rts.size == 0:
            self.p2_image.clear()
            return
        # x = m/z, y = RT -> image indexed [x=mz][y=rt] = z transposed.
        self.p2_image.setImage(np.log1p(z).T, autoLevels=True)
        rt_span = max(rts[-1] - rts[0], 1e-6) if rts.size > 1 else max(rt_end - rt_start, 1e-6)
        self.p2_image.setRect(pg.QtCore.QRectF(mz_min, rts[0], mz_max - mz_min, rt_span))

    def draw_ms2_strip(self, ms2):
        # Clickable MS2 scans as points at fixed x, RT on y (linked to panel 2).
        spots = [{"pos": (0.5, m["rt"]), "data": m} for m in ms2 if m.get("rt") is not None]
        self.ms2_scatter.setData(spots)

    def on_ms2_clicked(self, _scatter, points):
        if not points:
            return
        scan = points[0].data()
        cur = self.current
        if cur is None or cur["centroid"] is None:
            return
        try:
            spectrum = cur["centroid"].get_scan_by_id(scan["id"])
            if spectrum is None:
                return
            mz, inten = scan_arrays(spectrum)
            plot_spectrum(self.p3, mz, inten,
                          title=f"MS2  rt={scan['rt']:.3f}  precursor m/z={scan.get('mz')}")
            self.p3_title.setText(
                f"MS2 scan {scan.get('number','')}  rt={scan['rt']:.3f}  precursor m/z={scan.get('mz')}")
        except Exception as exc:
            self.p3_title.setText(f"MS2 load error: {exc}")

    # Wheel over panel 1's y-axis strip scrolls intensity (y zoom), even though
    # y-drag inside the plot is disabled.
    def eventFilter(self, obj, event):
        if obj is self.p1_2d.viewport() and event.type() == QEvent.Wheel:
            axis = self.p1_2d.getPlotItem().getAxis("left")
            if event.position().x() < axis.width():
                vb = self.p1_2d.getViewBox()
                (y0, y1) = vb.viewRange()[1]
                factor = 0.85 if event.angleDelta().y() > 0 else 1.0 / 0.85
                center = (y0 + y1) / 2
                half = (y1 - y0) / 2 * factor
                vb.setYRange(center - half, center + half, padding=0)
                return True
        return super().eventFilter(obj, event)

    def draw_panel3_ms1(self, cur, scan_mz, scan_int):
        self.p3.clear()
        self.p3.setLabel("bottom", "m/z")
        plain = plain_seq(cur["row"].get("peptide", ""))
        charge = cur["charge"] or 1
        title = f"MS1 isotope envelope - {plain} z={charge}"

        exp_peak = 1.0
        if scan_mz is not None and len(scan_int):
            exp_peak = float(np.max(scan_int)) or 1.0
            plot_spectrum(self.p3, scan_mz, scan_int, title=title)

        if plain and set(plain) <= set("ACDEFGHIKLMNPQRSTVWYUO"):
            try:
                t_mz, t_norm = isotopes.peptide_isotope_mzs(plain, charge)
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

    def render_table1(self, cur):
        rows = []
        if self.db is not None and cur["rt"] is not None and cur["charge"]:
            mz_min, mz_max = cur["mz_center"] - self.mz_half, cur["mz_center"] + self.mz_half
            dists = self.db.distributions_in_window(
                mz_min=mz_min, mz_max=mz_max,
                rt_start=cur["rt"] - self.rt_half, rt_end=cur["rt"] + self.rt_half,
                charge=cur["charge"], limit=1,
            )
            if dists:
                rows = self.db.distribution_members(dists[0]["distribution_id"])

        self.table1.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, (field, _) in enumerate(LINE_METRIC_COLUMNS):
                value = row.get(field, "")
                if isinstance(value, float):
                    value = f"{value:.4g}"
                self.table1.setItem(i, j, QTableWidgetItem(str(value)))

    # ---- toggles + sync --------------------------------------------------

    def toggle_dimension(self):
        to_3d = self.dim_toggle.isChecked()
        self.dim_toggle.setText("2D" if to_3d else "3D")
        self.p1_stack.setCurrentIndex(1 if to_3d else 0)

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
