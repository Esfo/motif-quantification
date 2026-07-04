from pathlib import Path

from PySide6.QtCore import QEvent, QSettings, Qt
from PySide6.QtWidgets import (
    QApplication,
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
    from .proteins_tab import ProteinsTab
    from .distributions_db import DistributionsDB
    from .experimental import ExperimentalSetup
except ImportError:
    from session import ViewerSession
    from theming import palette, style_plot
    from ms_viewer_tab import MSViewerTab
    from proteins_tab import ProteinsTab
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

    # The top controls toolbar (± m/z, ± RT spinboxes and navigation-history
    # buttons) was removed to reclaim vertical space: the ± m/z / ± RT controls
    # are gone, and navigation history moved into panel 1's button bar (left of
    # the charge-search arrows). Theme toggle lives in panel 1's bar; layout
    # reset stays in the View menu.

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
        if getattr(self, "proteins_tab", None) is not None:
            self.proteins_tab.apply_theme(theme)
        self.theme_action.setText("Switch to &light mode" if theme == "dark" else "Switch to &dark mode")

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
                                  xics_ppm=self.xics_ppm, xics_rt_window=self.xics_rt_window,
                                  theme=self.theme)
        self.ms_tab.on_theme_toggle = self.toggle_theme   # panel-1 Light/Dark button

        # Preserve whichever tab the user is on across open/reload (don't jump
        # back to MS Data).
        prev_index = 0
        old_central = self.centralWidget()
        if isinstance(old_central, QTabWidget):
            prev_index = old_central.currentIndex()

        tabs = QTabWidget()
        self.tabs = tabs
        # 'loading' now shows as a badge on top of panels 1/2 and as a row within
        # Table 1 (driven inside MSViewerTab), not in the tab bar.
        tabs.addTab(self.ms_tab, "MS Data")
        self.proteins_tab = ProteinsTab(self.session, theme=self.theme)
        self.proteins_tab.on_navigate_to_ms = self._navigate_to_ms
        self.proteins_tab.on_theme_toggle = self.toggle_theme
        tabs.addTab(self.proteins_tab, "Proteins")
        tabs.addTab(self._placeholder("File-by-file comparison",
                    f"Quantitative comparison across files: time series + differential "
                    f"expression. Reads experimental-setup "
                    f"({'loaded: ' + str(len(self.experimental.rows)) + ' rows' if not self.experimental.is_empty() else 'not found'}). "
                    "— staged, see ARCHITECTURE.md"), "Quantitative Comparisons")
        tabs.addTab(self._placeholder("Motif quantification",
                    "Time series + DE at the motif level; proteins grouped by shared skeleton "
                    "motif, with include/exclude refinement saved back to a motif-sets folder. "
                    "— staged, see ARCHITECTURE.md"), "Motifs")
        if 0 <= prev_index < tabs.count():
            tabs.setCurrentIndex(prev_index)
        self.setCentralWidget(tabs)
        # Restore the dock arrangement after the tab is laid out (sizes depend on
        # the final widget geometry), so opening a file doesn't reset it.
        from PySide6.QtCore import QTimer
        QTimer.singleShot(0, self.restore_layout)

        self.apply_theme(self.theme)
        self.reload_action.setEnabled(not self.session.is_empty)
        title = "no folder open (double-click to open)" if self.session.is_empty else str(reorganized)
        self.setWindowTitle(f"Motif Quantification Viewer — {title}")

    def _navigate_to_ms(self, filename, protein_id, peptide_plain):
        """Proteins-tab → MS Data tab jump: switch to MS Data and focus the given
        identification (double-clicking a peptide in the protein sequence)."""
        if self.ms_tab is None:
            return
        self.tabs.setCurrentIndex(self.tabs.indexOf(self.ms_tab))
        try:
            self.ms_tab.focus_identification(filename, protein_id, peptide_plain)
        except Exception:
            pass

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
        if getattr(self, "proteins_tab", None) is not None:
            h_state, v_state = self.proteins_tab.splitter_states()
            self.settings.setValue("proteins_h_state", h_state)
            self.settings.setValue("proteins_v_state", v_state)

    def restore_geometry(self):
        self._geometry_restored = False
        geom = self.settings.value("window_geometry")
        if geom is not None:
            try:
                if self.restoreGeometry(geom):
                    self._geometry_restored = True
            except Exception:
                pass

    def showEvent(self, event):
        # Re-apply the saved dock layout once the window is actually shown: a
        # restoreState done before the first show can be discarded as the nested
        # tab/dock geometry settles, which left panels back at their defaults.
        super().showEvent(event)
        if not getattr(self, "_shown_once", False):
            self._shown_once = True
            from PySide6.QtCore import QTimer
            QTimer.singleShot(0, self.restore_layout)

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
        if getattr(self, "proteins_tab", None) is not None \
                and str(version) == str(self.LAYOUT_VERSION):
            try:
                self.proteins_tab.restore_splitter_states(
                    self.settings.value("proteins_h_state"),
                    self.settings.value("proteins_v_state"))
            except Exception:
                pass

    def closeEvent(self, event):
        self.save_layout()
        super().closeEvent(event)
