from pathlib import Path

from PySide6.QtCore import QEvent, Qt
from PySide6.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QWidget,
)

import pyqtgraph as pg

try:
    from .session import ViewerSession
    from .theming import palette, style_plot
    from .workspace import Workspace
except ImportError:
    from session import ViewerSession
    from theming import palette, style_plot
    from workspace import Workspace


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
        self.workspace = None
        self._opening = False
        self.theme = "dark"

        self.build_menu()
        self.build_toolbar()

        QApplication.instance().installEventFilter(self)

        self.load_session(reorganized)

    # ---- chrome ----------------------------------------------------------

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

    def build_toolbar(self):
        bar = self.addToolBar("controls")
        bar.setMovable(False)

        bar.addWidget(QLabel(" ± m/z "))
        self.mz_spin = QDoubleSpinBox()
        self.mz_spin.setRange(0.1, 25.0)
        self.mz_spin.setDecimals(2)
        self.mz_spin.setSingleStep(0.5)
        self.mz_spin.setValue(2.5)
        self.mz_spin.valueChanged.connect(self.on_mz_changed)
        bar.addWidget(self.mz_spin)

        bar.addWidget(QLabel("  ± RT "))
        self.rt_spin = QDoubleSpinBox()
        self.rt_spin.setRange(0.02, 10.0)
        self.rt_spin.setDecimals(2)
        self.rt_spin.setSingleStep(0.1)
        self.rt_spin.setValue(self.xics_rt_window)
        self.rt_spin.valueChanged.connect(self.on_rt_changed)
        bar.addWidget(self.rt_spin)

        bar.addSeparator()
        reset_action = bar.addAction("Reset zoom")
        reset_action.setShortcut("Ctrl+0")
        reset_action.triggered.connect(self.on_reset_zoom)

        self.sync_action = bar.addAction("Sync: on")
        self.sync_action.setCheckable(True)
        self.sync_action.setChecked(True)
        self.sync_action.toggled.connect(self.on_sync_toggled)

        bar.addSeparator()
        self.theme_action = bar.addAction("Light mode")
        self.theme_action.setShortcut("Ctrl+T")
        self.theme_action.triggered.connect(self.toggle_theme)

    # ---- toolbar handlers (guard against no-workspace) -------------------

    def on_mz_changed(self, value):
        if self.workspace is not None:
            self.workspace.set_mz_half(value)

    def on_rt_changed(self, value):
        if self.workspace is not None:
            self.workspace.set_rt_half(value)

    def on_reset_zoom(self):
        if self.workspace is not None:
            self.workspace.reset_zoom()

    def on_sync_toggled(self, on):
        self.sync_action.setText("Sync: on" if on else "Sync: off")
        if self.workspace is not None:
            self.workspace.set_sync(on)

    def toggle_theme(self):
        self.apply_theme("light" if self.theme == "dark" else "dark")

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)

        for widget in self.findChildren(pg.PlotWidget):
            style_plot(widget, pal)

        if self.workspace is not None:
            self.workspace.apply_theme(theme)

        self.theme_action.setText("Light mode" if theme == "dark" else "Dark mode")

    # ---- opening folders -------------------------------------------------

    def eventFilter(self, obj, event):
        if (
            event.type() == QEvent.MouseButtonDblClick
            and (self.session is None or self.session.is_empty)
            and not self._opening
            and QApplication.activeModalWidget() is None
            and isinstance(obj, QWidget)
            and self.isAncestorOf(obj)
        ):
            self.choose_reorganized()
            return True

        return super().eventFilter(obj, event)

    def choose_reorganized(self):
        if self._opening:
            return

        self._opening = True
        try:
            start = ""
            if self.session is not None and self.session.reorganized is not None:
                start = str(self.session.reorganized)
            path = QFileDialog.getExistingDirectory(self, "Open reorganized search folder", start)
        finally:
            self._opening = False

        if path:
            self.open_reorganized(path)

    def reload_current(self):
        if self.session is not None and self.session.reorganized is not None:
            self.load_session(self.session.reorganized)

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

        self.workspace = Workspace(
            session=self.session,
            xics_ppm=self.xics_ppm,
            xics_rt_window=self.rt_spin.value(),
            theme=self.theme,
        )
        self.workspace.set_sync(self.sync_action.isChecked())
        self.setCentralWidget(self.workspace)

        self.apply_theme(self.theme)
        self.reload_action.setEnabled(not self.session.is_empty)

        if self.session.is_empty:
            self.setWindowTitle("Motif Quantification Viewer — no folder open (double-click to open)")
        else:
            self.setWindowTitle(f"Motif Quantification Viewer — {reorganized}")
