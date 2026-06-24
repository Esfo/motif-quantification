"""Unified spectrum workspace.

One area that ties together everything for a single Sage-validated match:

  left   : browse by file / protein / peptide, then a PSM list
  right   : MS2 spectrum, MS1 scan, MS1 XIC, a 2D RT x m/z heatmap, and the
            bounded 3D profile -- all centered on the selected match.

The plots that share an axis can be zoom/pan synchronized (m/z plots together,
RT plots together); sync can be toggled off, and a reset re-centers everything.
The m/z and RT half-windows are driven from the top toolbar and update live.
Every panel is in a splitter, so panes resize freely.
"""

import re
import traceback

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .mzml_store import MzmlStore, scan_arrays
    from .plots import add_profile_line, plot_points, plot_spectrum, plot_traces, short_file_label
    from .region_view import RegionView
    from .session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float
except ImportError:
    from mzml_store import MzmlStore, scan_arrays
    from plots import add_profile_line, plot_points, plot_spectrum, plot_traces, short_file_label
    from region_view import RegionView
    from session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float


PSM_COLUMNS = [
    ("file", "file"),
    ("scan", "scan"),
    ("peptide", "peptide"),
    ("charge", "z"),
    ("percolator_q", "perc q"),
    ("hyperscore", "hyperscore"),
    ("rt", "rt"),
]


def plain_seq(peptide):
    value = peptide or ""

    if len(value) >= 5 and value[1] == "." and value[-2] == ".":
        value = value[2:-2]

    value = re.sub(r"\[[^\]]*\]|\([^\)]*\)|\{[^\}]*\}", "", value)
    return re.sub(r"[^A-Za-z]", "", value).upper()


class Workspace(QWidget):
    def __init__(self, session, xics_ppm=10.0, xics_rt_window=0.8, theme="dark"):
        super().__init__()

        self.session = session
        self.xics_ppm = float(xics_ppm)
        self.theme = theme

        self.mz_half = 2.5
        self.rt_half = float(xics_rt_window)
        self.sync = True

        self._centroid = {}
        self._profile = {}
        self.current = None
        self.psm_rows = []

        self.build_ui()
        self.populate_primary()
        self.set_sync(True)

    # ---- UI --------------------------------------------------------------

    def build_ui(self):
        self.browse_combo = QComboBox()
        self.browse_combo.addItems(["File", "Protein", "Peptide"])
        self.browse_combo.currentIndexChanged.connect(self.on_browse_changed)

        self.search = QLineEdit()
        self.search.setPlaceholderText("filter…")
        self.search.textChanged.connect(self.populate_primary)

        self.primary = QListWidget()
        self.primary.itemSelectionChanged.connect(self.on_primary_changed)

        self.secondary = QListWidget()
        self.secondary.itemSelectionChanged.connect(self.on_secondary_changed)
        self.secondary_label = QLabel("peptides")

        self.psm_table = QTableWidget()
        self.psm_table.setColumnCount(len(PSM_COLUMNS))
        self.psm_table.setHorizontalHeaderLabels([h for _, h in PSM_COLUMNS])
        self.psm_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.psm_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.psm_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.psm_table.itemSelectionChanged.connect(self.on_psm_changed)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(4, 4, 4, 4)
        left_layout.addWidget(self.browse_combo)
        left_layout.addWidget(self.search)
        left_layout.addWidget(self.primary, stretch=2)
        left_layout.addWidget(self.secondary_label)
        left_layout.addWidget(self.secondary, stretch=2)
        left_layout.addWidget(QLabel("PSMs"))
        left_layout.addWidget(self.psm_table, stretch=3)

        self.ms2_plot = pg.PlotWidget()
        self.ms2_plot.setLabel("bottom", "m/z")
        self.ms1_plot = pg.PlotWidget()
        self.ms1_plot.setLabel("bottom", "m/z")
        self.xic_plot = pg.PlotWidget()
        self.xic_plot.setLabel("bottom", "RT", units="min")

        for plot in (self.ms2_plot, self.ms1_plot, self.xic_plot):
            plot.setClipToView(True)
            plot.setDownsampling(auto=True, mode="peak")

        self.region = RegionView()
        self.region.set_controls_visible(False)

        right = QSplitter(Qt.Vertical)
        right.addWidget(self._titled("MS2 spectrum", self.ms2_plot))
        right.addWidget(self._titled("MS1 scan", self.ms1_plot))
        right.addWidget(self._titled("MS1 XIC", self.xic_plot))
        right.addWidget(self._titled("MS1 profile region (2D + 3D)", self.region))
        right.setStretchFactor(0, 2)
        right.setStretchFactor(1, 2)
        right.setStretchFactor(2, 2)
        right.setStretchFactor(3, 5)

        main = QSplitter(Qt.Horizontal)
        main.addWidget(left)
        main.addWidget(right)
        main.setStretchFactor(0, 2)
        main.setStretchFactor(1, 6)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(main)

        self.secondary_label.hide()
        self.secondary.hide()

    def _titled(self, title, widget):
        box = QWidget()
        v = QVBoxLayout(box)
        v.setContentsMargins(2, 2, 2, 2)
        v.setSpacing(1)
        label = QLabel(title)
        label.setStyleSheet("font-weight: bold;")
        v.addWidget(label)
        v.addWidget(widget, stretch=1)
        return box

    # ---- selectors -------------------------------------------------------

    def browse_mode(self):
        return self.browse_combo.currentText()

    def on_browse_changed(self):
        protein = self.browse_mode() == "Protein"
        self.secondary_label.setVisible(protein)
        self.secondary.setVisible(protein)
        self.populate_primary()

    def populate_primary(self):
        self.primary.blockSignals(True)
        self.primary.clear()
        text = self.search.text().strip().lower()
        mode = self.browse_mode()

        if mode == "File":
            for row in self.session.files():
                name = row.get("filename", "")
                if name and (not text or text in name.lower()):
                    self._add(self.primary, name, name)
        elif mode == "Protein":
            for row in self.session.global_proteins():
                pid = row.get("protein_id", "")
                if pid and (not text or text in pid.lower()):
                    self._add(self.primary, pid, row)
        else:
            for row in self.session.global_peptides():
                pep = row.get("peptide", "")
                if pep and (not text or text in pep.lower()):
                    self._add(self.primary, pep, row)

        self.primary.blockSignals(False)

    def _add(self, listw, text, data):
        item = QListWidgetItem(text)
        item.setData(Qt.UserRole, data)
        listw.addItem(item)

    def on_primary_changed(self):
        items = self.primary.selectedItems()
        if not items:
            return

        data = items[0].data(Qt.UserRole)
        mode = self.browse_mode()

        if mode == "File":
            self.show_psms(self.psms_for_file(data))
        elif mode == "Protein":
            self.populate_secondary(data)
        else:
            self.show_psms(self.psms_for_peptide(plain_seq(data.get("peptide", ""))))

    def populate_secondary(self, protein_row):
        self.secondary.blockSignals(True)
        self.secondary.clear()
        peptides = [p for p in str(protein_row.get("peptides", "")).split(";") if p]
        for pep in peptides:
            self._add(self.secondary, pep, pep)
        self.secondary.blockSignals(False)
        self.show_psms([])

    def on_secondary_changed(self):
        items = self.secondary.selectedItems()
        if not items:
            return
        self.show_psms(self.psms_for_peptide(plain_seq(items[0].data(Qt.UserRole))))

    # ---- PSM gathering ---------------------------------------------------

    def psms_for_file(self, filename):
        rows = []
        for row in self.session.load_psms(filename):
            row = dict(row)
            row["filename"] = filename
            rows.append(row)
        return rows

    def psms_for_peptide(self, plain):
        rows = []
        seen_files = set()

        for gp in self.session.global_peptides():
            if plain_seq(gp.get("peptide", "")) != plain:
                continue
            for f in str(gp.get("files", "")).split(";"):
                if f:
                    seen_files.add(f)

        if not seen_files:
            seen_files = {r.get("filename", "") for r in self.session.files()}

        for filename in sorted(seen_files):
            for row in self.session.load_psms(filename):
                if plain_seq(row.get("peptide", "")) == plain:
                    row = dict(row)
                    row["filename"] = filename
                    rows.append(row)
        return rows

    def show_psms(self, rows):
        self.psm_rows = rows
        self.psm_table.setUpdatesEnabled(False)
        self.psm_table.setRowCount(len(rows))

        for i, row in enumerate(rows):
            for j, (field, _) in enumerate(PSM_COLUMNS):
                if field == "file":
                    text = short_file_label(row.get("filename", ""))
                else:
                    text = str(row.get(field, ""))
                item = QTableWidgetItem(text)
                item.setData(Qt.UserRole, i)
                self.psm_table.setItem(i, j, item)

        self.psm_table.setUpdatesEnabled(True)

    # ---- the match -> plots ---------------------------------------------

    def on_psm_changed(self):
        items = self.psm_table.selectedItems()
        if not items:
            return
        i = items[0].data(Qt.UserRole)
        if i is None or i >= len(self.psm_rows):
            return
        try:
            self.update_evidence(self.psm_rows[i])
        except Exception as exc:
            self.region.status.setText(f"error: {exc}\n{traceback.format_exc()}")

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

    def update_evidence(self, row):
        filename = row.get("filename", "")
        scan = row.get("scan", "")
        charge = peptide_charge(row)
        neutral_mass = peptide_mass(row)
        rt = peptide_rt(row)

        centroid = self.centroid_store(filename)
        profile = self.profile_store(filename)

        targets = []
        if neutral_mass is not None and charge:
            targets = isotope_mzs(neutral_mass, charge, n=6)

        mz_center = sum(targets) / len(targets) if targets else (neutral_mass or 500.0)

        self.current = {
            "row": row,
            "filename": filename,
            "scan": scan,
            "charge": charge,
            "neutral_mass": neutral_mass,
            "rt": rt,
            "targets": targets,
            "mz_center": mz_center,
            "centroid": centroid,
            "profile": profile,
        }

        self.refresh_plots()

    def refresh_plots(self):
        cur = self.current
        if cur is None:
            return

        centroid = cur["centroid"]
        targets = cur["targets"]
        rt = cur["rt"]
        scan = cur["scan"]
        charge = cur["charge"]
        mz_center = cur["mz_center"]
        label = f"{cur['row'].get('peptide', '')} z={charge}"

        if centroid is None:
            self.region.status.setText(f"no centroid mzML found for {cur['filename']}")
            return

        centroid.load_metadata()

        # MS2 spectrum (mass x intensity)
        ms2 = centroid.get_scan_by_number(scan)
        if ms2 is not None:
            mz, inten = scan_arrays(ms2)
            plot_spectrum(self.ms2_plot, mz, inten, title=f"MS2 scan {scan}: {label}")
        else:
            plot_spectrum(self.ms2_plot, [], [], title=f"MS2 scan {scan} not found")

        # nearest MS1
        if rt is None:
            ms1 = centroid.preceding_ms1_for_scan(scan)
        else:
            ms1 = centroid.nearest_ms1_by_rt(rt)

        mz_min = mz_center - self.mz_half
        mz_max = mz_center + self.mz_half

        if ms1 is not None:
            cmz, cinten = centroid.scan_window_by_number(ms1.number, mz_min, mz_max)
            plot_points(self.ms1_plot, cmz, cinten, title=f"MS1 scan {ms1.number}")
            profile = cur["profile"]
            if profile is not None:
                try:
                    pmz, pint = profile.scan_window_by_number(ms1.number, mz_min, mz_max)
                    add_profile_line(self.ms1_plot, pmz, pint)
                except Exception:
                    pass
        else:
            plot_points(self.ms1_plot, [], [], title="MS1 scan unavailable")

        # XIC over isotopes
        if rt is not None and targets:
            xic = centroid.extract_xics(
                targets=targets,
                rt_start=max(0.0, rt - self.rt_half),
                rt_end=rt + self.rt_half,
                ppm=self.xics_ppm,
            )
            plot_traces(self.xic_plot, xic["rts"], xic["traces"], xic["targets"], title="MS1 XIC")
        else:
            plot_traces(self.xic_plot, [], [], [], title="MS1 XIC unavailable")

        # 2D + 3D bounded region around this match
        if rt is not None:
            profile = cur["profile"]
            self.region.configure_window(mz_center, self.mz_half, rt, self.rt_half)
            self.region.set_target(centroid, profile, mz_min, mz_max,
                                   max(0.0, rt - self.rt_half), rt + self.rt_half, label=label)
            self.region.source_combo.setCurrentIndex(1 if profile is not None else 0)
            self.region.pending = False
            self.region.render_region()

        self.reset_zoom()

    # ---- toolbar hooks ---------------------------------------------------

    def set_mz_half(self, value):
        self.mz_half = float(value)
        self.refresh_plots()

    def set_rt_half(self, value):
        self.rt_half = float(value)
        self.refresh_plots()

    def reset_zoom(self):
        for plot in (self.ms2_plot, self.ms1_plot, self.xic_plot, self.region.heatmap):
            plot.getViewBox().autoRange()

    def set_sync(self, on):
        self.sync = on
        if on:
            self.ms1_plot.setXLink(self.ms2_plot)
            self.xic_plot.setXLink(self.region.heatmap)
            self.reset_zoom()
        else:
            self.ms1_plot.setXLink(None)
            self.xic_plot.setXLink(None)

    def apply_theme(self, theme):
        self.theme = theme
        self.region.apply_theme(theme)
