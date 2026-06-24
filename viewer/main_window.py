from pathlib import Path
import traceback

from PySide6.QtCore import QEvent, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .mzml_store import MzmlStore, scan_arrays
    from .plots import add_profile_line, plot_points, plot_spectrum, plot_traces
    from .session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float, safe_int, ViewerSession
    from .views import OverviewView, PeptidesView, ProteinsView
    from .region_view import RegionView
except ImportError:
    from mzml_store import MzmlStore, scan_arrays
    from plots import add_profile_line, plot_points, plot_spectrum, plot_traces
    from session import isotope_mzs, peptide_charge, peptide_mass, peptide_rt, safe_float, safe_int, ViewerSession
    from views import OverviewView, PeptidesView, ProteinsView
    from region_view import RegionView


PSM_COLUMNS = [
    ("scan", "scan"),
    ("peptide", "peptide"),
    ("charge", "z"),
    ("percolator_q", "perc q"),
    ("percolator_score", "perc score"),
    ("sage_spectrum_q", "sage q"),
    ("hyperscore", "hyperscore"),
    ("exp_mass", "exp mass"),
    ("calc_mass", "calc mass"),
    ("rt", "rt"),
    ("matched_peaks", "matched"),
    ("matched_intensity_pct", "ms2 %"),
]

DISTRIBUTION_COLUMNS = [
    ("distribution_id", "distribution"),
    ("charge", "z"),
    ("neutral_mass", "mass"),
    ("mono_mz", "mono m/z"),
    ("rt_apex", "rt apex"),
    ("n_members", "members"),
    ("score", "score"),
    ("quality", "quality"),
    ("mass_error", "mass err"),
    ("rt_error", "rt err"),
]


class MainWindow(QMainWindow):
    def __init__(
        self,
        reorganized=None,
        distribution_db=None,
        centroid_dir=None,
        profile_dir=None,
        xics_ppm=10.0,
        xics_rt_window=0.8,
    ):
        super().__init__()

        self.setWindowTitle("Motif Quantification Viewer")

        self.distribution_db = distribution_db
        self.centroid_dir = centroid_dir
        self.profile_dir = profile_dir
        self.xics_ppm = float(xics_ppm)
        self.xics_rt_window = float(xics_rt_window)

        self.session = None
        self.current_filename = None
        self.current_psms = []
        self.centroid_stores = {}
        self.profile_stores = {}

        self.build_menu()

        # Double-clicking anywhere while no folder is loaded opens the picker.
        QApplication.instance().installEventFilter(self)

        # Always show the full (empty) GUI; load data if a folder was given.
        self.load_session(reorganized)

    def eventFilter(self, obj, event):
        if (
            event.type() == QEvent.MouseButtonDblClick
            and (self.session is None or self.session.is_empty)
        ):
            self.choose_reorganized()
            return True

        return super().eventFilter(obj, event)

    def build_menu(self):
        file_menu = self.menuBar().addMenu("&File")

        open_action = file_menu.addAction("&Open reorganized folder…")
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.choose_reorganized)

        self.reload_action = file_menu.addAction("&Reload")
        self.reload_action.setShortcut("Ctrl+R")
        self.reload_action.triggered.connect(self.reload_current)
        self.reload_action.setEnabled(False)

        file_menu.addSeparator()
        quit_action = file_menu.addAction("&Quit")
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)

    def choose_reorganized(self):
        start = str(self.session.reorganized) if self.session is not None else ""
        path = QFileDialog.getExistingDirectory(
            self, "Open reorganized search folder", start
        )

        if path:
            self.open_reorganized(path)

    def reload_current(self):
        if self.session is not None:
            self.open_reorganized(self.session.reorganized)

    def open_reorganized(self, reorganized):
        reorganized = Path(reorganized)

        if not (reorganized / "files.tsv").exists():
            QMessageBox.warning(
                self,
                "Not a reorganized folder",
                f"{reorganized}\n\nNo files.tsv found here. Pick the 'reorganized' "
                "directory that reorganize-results.py wrote.",
            )
            return

        self.load_session(reorganized)

    def load_session(self, reorganized):
        try:
            self.session = ViewerSession(
                reorganized=reorganized,
                distribution_db=self.distribution_db,
                centroid_dir=self.centroid_dir,
                profile_dir=self.profile_dir,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Failed to open folder", str(exc))
            return

        self.current_filename = None
        self.current_psms = []
        self.centroid_stores = {}
        self.profile_stores = {}

        self.create_spectra_widgets()
        self.build_layout()
        self.load_files()

        self.reload_action.setEnabled(not self.session.is_empty)

        if self.session.is_empty:
            self.setWindowTitle("Motif Quantification Viewer — no folder open (double-click to open)")
        else:
            self.setWindowTitle(f"Motif Quantification Viewer — {reorganized}")

    def create_spectra_widgets(self):
        self.file_combo = QComboBox()
        self.profile_checkbox = QCheckBox("profile overlay")
        self.profile_checkbox.setChecked(False)
        self.profile_checkbox.stateChanged.connect(self.on_current_psm_changed)

        self.psm_table = QTableWidget()
        self.psm_table.setColumnCount(len(PSM_COLUMNS))
        self.psm_table.setHorizontalHeaderLabels([label for _, label in PSM_COLUMNS])
        self.psm_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.psm_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.psm_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.psm_table.itemSelectionChanged.connect(self.on_current_psm_changed)

        self.distribution_table = QTableWidget()
        self.distribution_table.setColumnCount(len(DISTRIBUTION_COLUMNS))
        self.distribution_table.setHorizontalHeaderLabels([label for _, label in DISTRIBUTION_COLUMNS])
        self.distribution_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.distribution_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.distribution_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.detail_text = QPlainTextEdit()
        self.detail_text.setReadOnly(True)

        self.ms2_plot = pg.PlotWidget()
        self.ms1_trace_plot = pg.PlotWidget()
        self.ms1_scan_plot = pg.PlotWidget()

    def build_layout(self):
        self.region_view = RegionView()

        self.tabs = QTabWidget()
        self.tabs.addTab(OverviewView(self.session), "Overview")
        self.tabs.addTab(self.build_spectra_tab(), "Spectra")
        self.tabs.addTab(self.region_view, "3D Profile")
        self.tabs.addTab(PeptidesView(self.session), "Peptides")
        self.tabs.addTab(ProteinsView(self.session), "Proteins")
        self.tabs.setCurrentIndex(1)
        self.tabs.currentChanged.connect(self.on_tab_changed)
        self.setCentralWidget(self.tabs)

    def on_tab_changed(self, index):
        if self.tabs.widget(index) is self.region_view:
            self.region_view.render_if_pending()

    def build_spectra_tab(self):
        top_bar = QWidget()
        top_layout = QHBoxLayout(top_bar)
        top_layout.setContentsMargins(4, 4, 4, 4)
        top_layout.addWidget(QLabel("file"))
        top_layout.addWidget(self.file_combo, stretch=1)
        top_layout.addWidget(self.profile_checkbox)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(4, 4, 4, 4)
        left_layout.addWidget(top_bar)
        left_layout.addWidget(self.psm_table, stretch=1)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(4, 4, 4, 4)

        upper_right = QSplitter(Qt.Horizontal)
        upper_right.addWidget(self.detail_text)
        upper_right.addWidget(self.distribution_table)
        upper_right.setStretchFactor(0, 1)
        upper_right.setStretchFactor(1, 2)

        plot_splitter = QSplitter(Qt.Vertical)
        plot_splitter.addWidget(self.ms2_plot)
        plot_splitter.addWidget(self.ms1_trace_plot)
        plot_splitter.addWidget(self.ms1_scan_plot)
        plot_splitter.setStretchFactor(0, 1)
        plot_splitter.setStretchFactor(1, 1)
        plot_splitter.setStretchFactor(2, 1)

        right_layout.addWidget(upper_right, stretch=1)
        right_layout.addWidget(plot_splitter, stretch=4)

        main = QSplitter(Qt.Horizontal)
        main.addWidget(left)
        main.addWidget(right)
        main.setStretchFactor(0, 2)
        main.setStretchFactor(1, 5)

        return main

    def load_files(self):
        self.file_combo.blockSignals(True)
        self.file_combo.clear()

        rows = self.session.files()

        for row in rows:
            filename = row.get("filename", "")

            if filename:
                label = filename
                n_psms = row.get("n_psms", "")

                if n_psms:
                    label = f"{filename} ({n_psms} PSMs)"

                self.file_combo.addItem(label, filename)

        self.file_combo.blockSignals(False)
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)

        if self.file_combo.count():
            self.file_combo.setCurrentIndex(0)
            self.on_file_changed(0)

    def on_file_changed(self, index):
        filename = self.file_combo.itemData(index)

        if not filename:
            return

        self.current_filename = filename
        self.current_psms = self.session.load_psms(filename)
        self.populate_psm_table(self.current_psms)

        centroid_path = self.session.centroid_path(filename)
        profile_path = self.session.profile_path(filename)

        status = [
            f"file: {filename}",
            f"psms: {len(self.current_psms)}",
            f"centroid mzML: {centroid_path or '(not found)'}",
            f"profile mzML: {profile_path or '(not found)'}",
            f"reorganized: {self.session.reorganized}",
        ]

        if self.session.distribution_db:
            status.append(f"distribution db: {self.session.distribution_db}")

        self.detail_text.setPlainText("\n".join(status))

    def populate_psm_table(self, rows):
        self.psm_table.setRowCount(len(rows))

        for row_i, row in enumerate(rows):
            for col_i, (field, _) in enumerate(PSM_COLUMNS):
                value = row.get(field, "")
                item = QTableWidgetItem(str(value))
                item.setData(Qt.UserRole, row_i)

                if field in {"scan", "charge", "matched_peaks"}:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

                self.psm_table.setItem(row_i, col_i, item)

        self.psm_table.resizeColumnsToContents()

    def current_psm(self):
        selected = self.psm_table.selectedItems()

        if not selected:
            return None

        row_i = selected[0].row()

        if row_i < 0 or row_i >= len(self.current_psms):
            return None

        return self.current_psms[row_i]

    def on_current_psm_changed(self):
        row = self.current_psm()

        if row is None:
            return

        try:
            self.update_evidence(row)
        except Exception as exc:
            self.detail_text.setPlainText(
                "viewer error\n"
                f"{exc}\n\n"
                f"{traceback.format_exc()}"
            )

    def get_centroid_store(self):
        if not self.current_filename:
            return None

        path = self.session.centroid_path(self.current_filename)

        if path is None:
            return None

        key = str(path)

        if key not in self.centroid_stores:
            self.centroid_stores[key] = MzmlStore(path)

        return self.centroid_stores[key]

    def get_profile_store(self):
        if not self.current_filename:
            return None

        path = self.session.profile_path(self.current_filename)

        if path is None:
            return None

        key = str(path)

        if key not in self.profile_stores:
            self.profile_stores[key] = MzmlStore(path)

        return self.profile_stores[key]

    def update_evidence(self, row):
        scan_number = row.get("scan", "")
        charge = peptide_charge(row)
        neutral_mass = peptide_mass(row)
        rt = peptide_rt(row)

        detail_lines = []
        detail_lines.append("Selected ID")
        detail_lines.append("===========")

        for key in [
            "scan",
            "scan_native",
            "psm_id",
            "peptide",
            "peptide_plain",
            "peptide_flanked",
            "proteins",
            "charge",
            "exp_mass",
            "calc_mass",
            "precursor_ppm",
            "rt",
            "aligned_rt",
            "percolator_q",
            "percolator_score",
            "sage_spectrum_q",
            "hyperscore",
            "matched_peaks",
            "matched_intensity_pct",
        ]:
            detail_lines.append(f"{key}: {row.get(key, '')}")

        centroid_store = self.get_centroid_store()

        if centroid_store is None:
            self.detail_text.setPlainText("\n".join(detail_lines + ["", "No centroid mzML found."]))
            return

        centroid_store.load_metadata()

        ms2_scan = centroid_store.get_scan_by_number(scan_number)

        if ms2_scan is not None:
            ms2_mz, ms2_intensity = scan_arrays(ms2_scan)
            plot_spectrum(
                self.ms2_plot,
                ms2_mz,
                ms2_intensity,
                title=f"MS2 scan {scan_number}: {row.get('peptide', '')}",
            )
        else:
            plot_spectrum(
                self.ms2_plot,
                [],
                [],
                title=f"MS2 scan {scan_number} not found",
            )

        if rt is None:
            ms1_summary = centroid_store.preceding_ms1_for_scan(scan_number)
        else:
            ms1_summary = centroid_store.nearest_ms1_by_rt(rt)

        if neutral_mass is not None and charge is not None:
            targets = isotope_mzs(neutral_mass, charge, n=6)
        else:
            targets = []

        if rt is not None and targets:
            rt_start = max(0.0, rt - self.xics_rt_window)
            rt_end = rt + self.xics_rt_window

            xics = centroid_store.extract_xics(
                targets=targets,
                rt_start=rt_start,
                rt_end=rt_end,
                ppm=self.xics_ppm,
                abs_tol=0.01,
            )

            plot_traces(
                self.ms1_trace_plot,
                xics["rts"],
                xics["traces"],
                xics["targets"],
                title=f"MS1 centroid isotope traces, z={charge}, mass={neutral_mass:.4f}",
            )
        else:
            plot_traces(
                self.ms1_trace_plot,
                [],
                [],
                [],
                title="MS1 centroid isotope traces unavailable",
            )

        if ms1_summary is not None and targets:
            mz_min = min(targets) - 0.75
            mz_max = max(targets) + 0.75

            centroid_mz, centroid_intensity = centroid_store.scan_window_by_number(
                ms1_summary.number,
                mz_min,
                mz_max,
            )

            plot_points(
                self.ms1_scan_plot,
                centroid_mz,
                centroid_intensity,
                title=f"MS1 centroid scan {ms1_summary.number}, RT={ms1_summary.rt:.4f}",
            )

            if self.profile_checkbox.isChecked():
                profile_store = self.get_profile_store()

                if profile_store is not None:
                    profile_mz, profile_intensity = profile_store.scan_window_by_number(
                        ms1_summary.number,
                        mz_min,
                        mz_max,
                    )

                    add_profile_line(self.ms1_scan_plot, profile_mz, profile_intensity)
        else:
            plot_points(
                self.ms1_scan_plot,
                [],
                [],
                title="MS1 scan unavailable",
            )

        candidates = self.session.distribution_candidates(
            neutral_mass=neutral_mass,
            charge=charge,
            rt=rt,
            ppm=20.0,
            rt_window=max(1.0, self.xics_rt_window),
            limit=50,
        )

        self.populate_distribution_table(candidates)

        detail_lines.append("")
        detail_lines.append("Derived evidence target")
        detail_lines.append("=======================")
        detail_lines.append(f"neutral_mass: {neutral_mass if neutral_mass is not None else ''}")
        detail_lines.append(f"charge: {charge if charge is not None else ''}")
        detail_lines.append(f"rt: {rt if rt is not None else ''}")

        if targets:
            detail_lines.append("expected isotope m/z:")
            for isotope_i, mz_value in enumerate(targets):
                detail_lines.append(f"  M+{isotope_i}: {mz_value:.6f}")

        if ms1_summary is not None:
            detail_lines.append("")
            detail_lines.append("Nearest MS1")
            detail_lines.append("===========")
            detail_lines.append(f"scan: {ms1_summary.number}")
            detail_lines.append(f"rt: {ms1_summary.rt}")
            detail_lines.append(f"id: {ms1_summary.spectrum_id}")

        detail_lines.append("")
        detail_lines.append(f"distribution candidates: {len(candidates)}")

        self.detail_text.setPlainText("\n".join(detail_lines))

        self.arm_region_view(row, targets, neutral_mass, charge, rt, centroid_store)

    def arm_region_view(self, row, targets, neutral_mass, charge, rt, centroid_store):
        if rt is None:
            return

        if targets:
            mz_min = min(targets) - 0.75
            mz_max = max(targets) + 0.75
        elif neutral_mass is not None and charge:
            mono = neutral_mass / charge
            mz_min = mono - 1.0
            mz_max = mono + 4.0
        else:
            return

        rt_start = max(0.0, rt - self.xics_rt_window)
        rt_end = rt + self.xics_rt_window

        self.region_view.set_target(
            centroid_store=centroid_store,
            profile_store=self.get_profile_store(),
            mz_min=mz_min,
            mz_max=mz_max,
            rt_start=rt_start,
            rt_end=rt_end,
            label=f"scan {row.get('scan', '')} {row.get('peptide', '')} z={charge}",
        )

    def populate_distribution_table(self, rows):
        self.distribution_table.setRowCount(len(rows))

        for row_i, row in enumerate(rows):
            for col_i, (field, _) in enumerate(DISTRIBUTION_COLUMNS):
                value = row.get(field, "")

                if isinstance(value, float):
                    if field in {"score", "quality", "mass_error", "rt_error"}:
                        text = f"{value:.6g}"
                    else:
                        text = f"{value:.6f}"
                else:
                    text = str(value)

                item = QTableWidgetItem(text)
                item.setData(Qt.UserRole, row)
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.distribution_table.setItem(row_i, col_i, item)

        self.distribution_table.resizeColumnsToContents()
