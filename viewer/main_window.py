from pathlib import Path

from PySide6.QtCore import QEvent, QSettings, Qt
from PySide6.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .session import ViewerSession
    from .theming import palette, style_plot
    from .ms_viewer_tab import MSViewerTab
    from .distributions_db import DistributionsDB
    from .experimental import ExperimentalSetup
except ImportError:
    from session import ViewerSession
    from theming import palette, style_plot
    from ms_viewer_tab import MSViewerTab
    from distributions_db import DistributionsDB
    from experimental import ExperimentalSetup


def find_distributions_db(reorganized):
    """Locate a distributions sqlite at <project>/distributions/*.sqlite."""
    if reorganized is None:
        return None
    # reorganized = <project>/searches/reorganized
    project = Path(reorganized).resolve().parent.parent
    dist_dir = project / "distributions"
    if dist_dir.is_dir():
        for path in sorted(dist_dir.glob("*.sqlite")):
            return path
    return None


def find_experimental_setup(reorganized):
    if reorganized is None:
        return None
    project = Path(reorganized).resolve().parent.parent
    path = project / "experimental-setup"
    return path if path.exists() else None


class MainWindow(QMainWindow):
    def __init__(self, reorganized=None, distribution_db=None, centroid_dir=None,
                 profile_dir=None, xics_ppm=10.0, xics_rt_window=0.8):
        super().__init__()
        self.setWindowTitle("Motif Quantification Viewer")

        self.distribution_db = distribution_db
        self.centroid_dir = centroid_dir
        self.profile_dir = profile_dir
        self.xics_ppm = float(xics_ppm)
        self.xics_rt_window = float(xics_rt_window)

        self.session = None
        self.ms_tab = None
        self._opening = False
        self.theme = "dark"
        self.settings = QSettings("motif-quantification", "viewer")

        self.build_menu()
        self.build_toolbar()
        QApplication.instance().installEventFilter(self)
        self.restore_geometry()
        self.load_session(reorganized)

        # Autosave the dock layout periodically so it persists even if the app
        # is force-quit / Ctrl+C'd (closeEvent wouldn't run then).
        from PySide6.QtCore import QTimer
        self._layout_timer = QTimer(self)
        self._layout_timer.setInterval(4000)
        self._layout_timer.timeout.connect(self.save_layout)
        self._layout_timer.start()

    # ---- chrome ----------------------------------------------------------

    def build_menu(self):
        file_menu = self.menuBar().addMenu("&File")
        open_action = file_menu.addAction("&Open project / reorganized folder…")
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

        view_menu = self.menuBar().addMenu("&View")
        self.theme_action = view_menu.addAction("Switch to &light mode")
        self.theme_action.setShortcut("Ctrl+T")
        self.theme_action.triggered.connect(self.toggle_theme)
        reset_layout = view_menu.addAction("Reset &panel layout")
        reset_layout.triggered.connect(self.reset_layout)

    def build_toolbar(self):
        bar = self.addToolBar("controls")
        bar.setMovable(False)
        bar.addWidget(QLabel(" ± m/z "))
        self.mz_spin = self._spin(0.1, 25.0, 2, 2.5, 0.5, self.on_mz_changed)
        bar.addWidget(self.mz_spin)
        bar.addWidget(QLabel("  ± RT "))
        self.rt_spin = self._spin(0.02, 10.0, 2, self.xics_rt_window, 0.1, self.on_rt_changed)
        bar.addWidget(self.rt_spin)
        bar.addSeparator()
        bar.addWidget(QLabel(" charge "))
        back = bar.addAction("◀")
        back.setToolTip("charge search: one charge lower at the same RT")
        back.triggered.connect(lambda: self.ms_tab and self.ms_tab.charge_step(-1))
        fwd = bar.addAction("▶")
        fwd.setToolTip("charge search: one charge higher at the same RT")
        fwd.triggered.connect(lambda: self.ms_tab and self.ms_tab.charge_step(1))
        hist_back = bar.addAction("⟲")
        hist_back.setToolTip("navigation history: back")
        hist_back.triggered.connect(lambda: self.ms_tab and self.ms_tab.nav_back())
        hist_fwd = bar.addAction("⟳")
        hist_fwd.setToolTip("navigation history: forward")
        hist_fwd.triggered.connect(lambda: self.ms_tab and self.ms_tab.nav_forward())

        bar.addSeparator()
        reset_zoom = bar.addAction("Reset zoom")
        reset_zoom.setShortcut("Ctrl+0")
        reset_zoom.triggered.connect(self.on_reset_zoom)
        reset_layout = bar.addAction("Reset layout")
        reset_layout.triggered.connect(self.reset_layout)
        align3d = bar.addAction("Align 3D")
        align3d.setToolTip("reset the Panel 1 3D view to its default orientation")
        align3d.triggered.connect(lambda: self.ms_tab and self.ms_tab.reset_3d_view())
        bar.addSeparator()
        self.theme_action_tb = bar.addAction("Light mode")
        self.theme_action_tb.triggered.connect(self.toggle_theme)

    def _spin(self, lo, hi, decimals, value, step, handler):
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(value)
        spin.valueChanged.connect(handler)
        return spin

    # ---- toolbar handlers ------------------------------------------------

    def on_mz_changed(self, value):
        if self.ms_tab is not None:
            self.ms_tab.set_mz_half(value)

    def on_rt_changed(self, value):
        if self.ms_tab is not None:
            self.ms_tab.set_rt_half(value)

    def on_reset_zoom(self):
        if self.ms_tab is not None:
            for plot in (self.ms_tab.p1_2d, self.ms_tab.p2, self.ms_tab.p3):
                plot.getViewBox().autoRange()

    def reset_layout(self):
        if self.ms_tab is not None:
            self.ms_tab.reset_layout()

    def toggle_theme(self):
        self.apply_theme("light" if self.theme == "dark" else "dark")

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)
        for widget in self.findChildren(pg.PlotWidget):
            style_plot(widget, pal)
        if self.ms_tab is not None:
            self.ms_tab.apply_theme(theme)
        self.theme_action.setText("Switch to &light mode" if theme == "dark" else "Switch to &dark mode")
        self.theme_action_tb.setText("Light mode" if theme == "dark" else "Dark mode")

    # ---- opening folders -------------------------------------------------

    def eventFilter(self, obj, event):
        if (event.type() == QEvent.MouseButtonDblClick
                and (self.session is None or self.session.is_empty)
                and not self._opening
                and QApplication.activeModalWidget() is None
                and isinstance(obj, QWidget) and self.isAncestorOf(obj)):
            self.choose_reorganized()
            return True
        return super().eventFilter(obj, event)

    def choose_reorganized(self):
        if self._opening:
            return
        self._opening = True
        try:
            start = self.settings.value("last_open_dir", "")
            if self.session and self.session.reorganized:
                start = str(self.session.reorganized)
            path = QFileDialog.getExistingDirectory(self, "Open project or reorganized folder", start)
        finally:
            self._opening = False
        if path:
            self.settings.setValue("last_open_dir", path)
            self.open_reorganized(path)

    def reload_current(self):
        if self.session is not None and self.session.reorganized is not None:
            self.load_session(self.session.reorganized)

    def open_reorganized(self, path):
        path = Path(path)
        # Accept either the reorganized dir or a project dir containing it.
        if (path / "files.tsv").exists():
            reorganized = path
        elif (path / "searches" / "reorganized" / "files.tsv").exists():
            reorganized = path / "searches" / "reorganized"
        else:
            QMessageBox.warning(self, "Not a project folder",
                                f"{path}\n\nExpected a reorganized folder (files.tsv) or a "
                                "project containing searches/reorganized.")
            return
        self.load_session(reorganized)

    def load_session(self, reorganized):
        try:
            self.session = ViewerSession(
                reorganized=reorganized, distribution_db=self.distribution_db,
                centroid_dir=self.centroid_dir, profile_dir=self.profile_dir)
        except Exception as exc:
            QMessageBox.critical(self, "Failed to open folder", str(exc))
            return

        db = None
        db_path = self.distribution_db or find_distributions_db(reorganized)
        if db_path and Path(db_path).exists():
            db = DistributionsDB(db_path)

        setup_path = find_experimental_setup(reorganized)
        self.experimental = ExperimentalSetup.load(setup_path) if setup_path else ExperimentalSetup([])

        self.ms_tab = MSViewerTab(self.session, distributions_db=db,
                                  xics_ppm=self.xics_ppm, xics_rt_window=self.rt_spin.value(),
                                  theme=self.theme)

        tabs = QTabWidget()
        tabs.addTab(self.ms_tab, "MS viewing")
        tabs.addTab(self._placeholder("Protein viewing",
                    "Whole-protein sequences with peptide coverage coloured by q-value "
                    "(shared q-value colour scale). Single-file or verticalized side-by-side "
                    "across files. — staged, see ARCHITECTURE.md"), "Proteins")
        tabs.addTab(self._placeholder("File-by-file comparison",
                    f"Quantitative comparison across files: time series + differential "
                    f"expression. Reads experimental-setup "
                    f"({'loaded: ' + str(len(self.experimental.rows)) + ' rows' if not self.experimental.is_empty() else 'not found'}). "
                    "— staged, see ARCHITECTURE.md"), "File comparison")
        tabs.addTab(self._placeholder("Motif quantification",
                    "Time series + DE at the motif level; proteins grouped by shared skeleton "
                    "motif, with include/exclude refinement saved back to a motif-sets folder. "
                    "— staged, see ARCHITECTURE.md"), "Motifs")
        self.setCentralWidget(tabs)
        # Restore the dock arrangement after the tab is laid out (sizes depend on
        # the final widget geometry), so opening a file doesn't reset it.
        from PySide6.QtCore import QTimer
        QTimer.singleShot(0, self.restore_layout)

        self.apply_theme(self.theme)
        self.reload_action.setEnabled(not self.session.is_empty)
        title = "no folder open (double-click to open)" if self.session.is_empty else str(reorganized)
        self.setWindowTitle(f"Motif Quantification Viewer — {title}")

    def _placeholder(self, title, text):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        heading = QLabel(title)
        heading.setStyleSheet("font-size: 16px; font-weight: bold;")
        body = QLabel(text)
        body.setWordWrap(True)
        layout.addWidget(heading)
        layout.addWidget(body)
        layout.addStretch(1)
        return widget

    # ---- dock layout persistence ----------------------------------------

    # Bump when the default dock arrangement changes so a stale saved layout
    # doesn't override the new default (the user can still rearrange + it saves).
    LAYOUT_VERSION = 3

    def save_layout(self):
        self.settings.setValue("window_geometry", self.saveGeometry())
        if self.ms_tab is not None:
            self.settings.setValue("ms_tab_state", self.ms_tab.saveState())
            self.settings.setValue("ms_tab_layout_version", self.LAYOUT_VERSION)

    def restore_geometry(self):
        geom = self.settings.value("window_geometry")
        if geom is not None:
            try:
                self.restoreGeometry(geom)
            except Exception:
                pass

    def restore_layout(self):
        if self.ms_tab is None:
            return
        version = self.settings.value("ms_tab_layout_version")
        state = self.settings.value("ms_tab_state")
        if state is not None and str(version) == str(self.LAYOUT_VERSION):
            try:
                self.ms_tab.restoreState(state)
            except Exception:
                pass

    def closeEvent(self, event):
        self.save_layout()
        super().closeEvent(event)
