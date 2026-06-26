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

import numpy as np
from PySide6.QtCore import Qt
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
    from .region_view import HAVE_GL, RegionWorker, gl
    from .session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
    from . import isotopes
    from .theming import palette, style_plot, style_gl
except ImportError:
    from mzml_store import MzmlStore, scan_arrays
    from plots import plot_points, plot_spectrum
    from region_view import HAVE_GL, RegionWorker, gl
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

        self.show_all_proteins()
        self.show_all_peptides()
        self.show_all_psms()

    def _titled_list(self, layout, title, on_select, on_all):
        header = QHBoxLayout()
        header.addWidget(QLabel(title))
        header.addStretch(1)
        all_button = QPushButton("All")
        all_button.setFixedWidth(40)
        all_button.clicked.connect(on_all)
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
        self.p1_2d.sigXRangeChanged.connect(self.on_panel1_mz_changed)

        if HAVE_GL:
            self.p1_3d = gl.GLViewWidget()
            self.p1_3d.setCameraPosition(distance=3.2, elevation=30, azimuth=-60)
            self.p1_surface = gl.GLSurfacePlotItem(
                x=np.array([0.0, 1.0], dtype=np.float32),
                y=np.array([0.0, 1.0], dtype=np.float32),
                z=np.zeros((2, 2), dtype=np.float32),
                shader="shaded", smooth=True,
            )
            self.p1_surface.setVisible(False)
            self.p1_scatter = gl.GLScatterPlotItem(pos=np.zeros((1, 3), dtype=np.float32), size=4.0)
            self.p1_scatter.setVisible(False)
            self.p1_3d.addItem(self.p1_surface)
            self.p1_3d.addItem(self.p1_scatter)
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
        self.p2 = pg.PlotWidget()
        self.p2.setLabel("bottom", "RT", units="min")
        self.p2.setLabel("left", "m/z")
        self.p2_image = pg.ImageItem()
        try:
            self.p2_image.setColorMap(pg.colormap.get("viridis"))
        except Exception:
            pass
        self.p2.addItem(self.p2_image)
        self.p2.sigYRangeChanged.connect(self.on_panel2_mz_changed)

        dock = QDockWidget("Panel 2 - RT x m/z map", self)
        dock.setObjectName("dock_panel2")
        dock.setWidget(self.p2)
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
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_lists)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_panel1)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_panel3)
        self.splitDockWidget(self.dock_panel1, self.dock_panel2, Qt.Vertical)
        self.splitDockWidget(self.dock_panel2, self.dock_table1, Qt.Vertical)
        self.resizeDocks([self.dock_lists], [320], Qt.Horizontal)

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

    def _fill(self, listw, entries):
        listw.blockSignals(True)
        listw.clear()
        for text, data in entries:
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, data)
            listw.addItem(item)
        listw.blockSignals(False)

    def _filter(self, text):
        t = self.search.text().strip().lower()
        return (not t) or (t in text.lower())

    def show_all_proteins(self):
        self.active_list = "protein"
        self._fill(self.protein_list, [(r.get("protein_id", ""), r) for r in self.session.global_proteins()
                                       if r.get("protein_id") and self._filter(r["protein_id"])])

    def show_all_peptides(self):
        self._fill(self.peptide_list, [(r.get("peptide", ""), r) for r in self.session.global_peptides()
                                       if r.get("peptide") and self._filter(r["peptide"])])

    def show_all_psms(self):
        rows = []
        for gp in self.session.global_peptides():
            rows.append((gp.get("peptide", ""), gp))
        self._fill(self.psm_list, [r for r in rows if self._filter(r[0])])

    def repopulate_active_list(self):
        self.show_all_proteins()
        self.show_all_peptides()
        self.show_all_psms()

    def on_protein_selected(self):
        items = self.protein_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        peptides = [p for p in str(row.get("peptides", "")).split(";") if p]
        self._fill(self.peptide_list, [(p, {"peptide": p}) for p in peptides])

    def on_peptide_selected(self):
        items = self.peptide_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        plain = plain_seq(row.get("peptide", ""))
        # cross-link proteins
        proteins = [p for p in str(row.get("proteins", "")).split(";") if p]
        if proteins:
            self._fill(self.protein_list, [(p, {"protein_id": p, "peptides": ""}) for p in proteins])
        # list the PSMs of this peptide
        self.psm_rows = self.psms_for_peptide(plain)
        self._fill(self.psm_list, [(f"{r.get('scan','')}  {r.get('peptide','')}", r) for r in self.psm_rows])
        if self.psm_rows:
            self.update_evidence(self.psm_rows[0])

    def on_psm_selected(self):
        items = self.psm_list.selectedItems()
        if not items:
            return
        row = items[0].data(Qt.UserRole)
        if "scan" in row:
            self.update_evidence(row)
        else:
            # a peptide row from the all-PSMs list
            plain = plain_seq(row.get("peptide", ""))
            self.psm_rows = self.psms_for_peptide(plain)
            self._fill(self.psm_list, [(f"{r.get('scan','')}  {r.get('peptide','')}", r) for r in self.psm_rows])
            if self.psm_rows:
                self.update_evidence(self.psm_rows[0])

    def psms_for_peptide(self, plain):
        rows = []
        files = set()
        for gp in self.session.global_peptides():
            if plain_seq(gp.get("peptide", "")) == plain:
                for f in str(gp.get("files", "")).split(";"):
                    if f:
                        files.add(f)
        if not files:
            files = {r.get("filename", "") for r in self.session.files()}
        for filename in sorted(files):
            for row in self.session.load_psms(filename):
                if plain_seq(row.get("peptide", "")) == plain:
                    row = dict(row)
                    row["filename"] = filename
                    rows.append(row)
        return rows

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
        self.refresh()

    def use_profile(self):
        return self.source_toggle.isChecked()

    def refresh(self):
        cur = self.current
        if cur is None:
            return
        centroid = cur["centroid"]
        if centroid is None:
            self.p3_title.setText(f"no centroid mzML for {cur['filename']}")
            return
        centroid.load_metadata()

        mz_center = cur["mz_center"]
        rt = cur["rt"]
        mz_min, mz_max = mz_center - self.mz_half, mz_center + self.mz_half

        # Panel 1 - 2D spectrum (m/z x intensity) of nearest MS1 scan in window
        store = cur["profile"] if (self.use_profile() and cur["profile"]) else centroid
        ms1 = centroid.nearest_ms1_by_rt(rt) if rt is not None else centroid.preceding_ms1_for_scan(cur["scan"])
        if ms1 is not None:
            mz, inten = store.scan_window_by_number(ms1.number, mz_min, mz_max)
            plot_points(self.p1_2d, mz, inten, title=f"{cur['row'].get('peptide','')} z={cur['charge']}")
        else:
            plot_points(self.p1_2d, [], [], title="MS1 scan unavailable")

        # Panel 2 - RT x m/z map (threaded extract)
        self.render_panel2(store, mz_min, mz_max, rt)

        # Panel 3 - MS1 isotope overlay (theoretical vs experimental)
        self.render_panel3_ms1(cur, ms1, store, mz_min, mz_max)

        # Table 1 - distribution lines (from sqlite if present)
        self.render_table1(cur)

    def render_panel2(self, store, mz_min, mz_max, rt):
        if rt is None:
            return
        if self.worker is not None and self.worker.isRunning():
            return
        params = dict(mz_min=mz_min, mz_max=mz_max,
                      rt_start=max(0.0, rt - self.rt_half), rt_end=rt + self.rt_half,
                      mz_bins=400, mode="profile" if self.use_profile() else "centroid")
        self._p2_win = (mz_min, mz_max, params["rt_start"], params["rt_end"])
        self.worker = RegionWorker(store, params)
        self.worker.done.connect(self.on_panel2_done)
        self.worker.start()

    def on_panel2_done(self, region):
        if not isinstance(region, dict) or "error" in region:
            return
        z = region.get("z")
        rts = region.get("rts")
        if z is None or z.size == 0 or rts.size == 0:
            self.p2_image.clear()
            return
        mz_min, mz_max, rt_start, rt_end = self._p2_win
        self.p2_image.setImage(np.log1p(z), autoLevels=True)
        rt_span = max(rts[-1] - rts[0], 1e-6) if rts.size > 1 else max(rt_end - rt_start, 1e-6)
        self.p2_image.setRect(pg.QtCore.QRectF(rts[0], mz_min, rt_span, mz_max - mz_min))

    def render_panel3_ms1(self, cur, ms1, store, mz_min, mz_max):
        self.p3.clear()
        self.p3.setLabel("bottom", "m/z")
        plain = plain_seq(cur["row"].get("peptide", ""))
        charge = cur["charge"] or 1
        title = f"MS1 isotope envelope - {plain} z={charge}"

        # experimental sticks
        exp_peak = 1.0
        if ms1 is not None:
            mz, inten = store.scan_window_by_number(ms1.number, mz_min, mz_max)
            if len(inten):
                exp_peak = float(np.max(inten)) or 1.0
                plot_spectrum(self.p3, mz, inten, title=title)

        # theoretical overlay scaled to the experimental peak height
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

    def on_panel1_mz_changed(self, _vb, rng):
        if getattr(self, "_syncing", False) or self.p1_stack.currentIndex() != 0:
            return
        self._syncing = True
        try:
            self.p2.setYRange(rng[0], rng[1], padding=0)
        finally:
            self._syncing = False

    def on_panel2_mz_changed(self, _vb, rng):
        if getattr(self, "_syncing", False):
            return
        self._syncing = True
        try:
            if self.p1_stack.currentIndex() == 0:
                self.p1_2d.setXRange(rng[0], rng[1], padding=0)
        finally:
            self._syncing = False

    def set_mz_half(self, value):
        self.mz_half = float(value)
        self.refresh()

    def set_rt_half(self, value):
        self.rt_half = float(value)
        self.refresh()
