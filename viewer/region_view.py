"""3D / region viewer for MS1 profile (and centroid) data.

Driven by a selected Sage ID (or a specific MS1 distribution), this shows a
*bounded* RT x m/z window -- never the whole spectrum -- two ways:

  * a 2D heatmap (RT vs m/z, intensity as color) as the always-on overview;
  * a 3D view where the measured points are drawn individually and an
    interpolated surface is laid over them to show peak shape.

Extraction runs on a worker thread so reading multi-GB profile files never
freezes the UI. Nothing here fits or models peaks; it only bins and displays
measured intensities (the surface is plain grid interpolation).
"""

import traceback

import numpy as np
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

import pyqtgraph as pg

try:
    from .theming import palette, style_gl, style_plot
except ImportError:
    from theming import palette, style_gl, style_plot

try:
    import pyqtgraph.opengl as gl

    HAVE_GL = True
except Exception:
    gl = None
    HAVE_GL = False


# Hard caps so a careless window can't hang the reader.
MAX_HALF_MZ = 25.0
MAX_HALF_RT = 10.0
MAX_SCATTER_POINTS = 30000


def _colormap():
    try:
        return pg.colormap.get("viridis")
    except Exception:
        return pg.colormap.get("CET-L9")


class RegionWorker(QThread):
    done = Signal(object)

    def __init__(self, store, params):
        super().__init__()
        self.store = store
        self.params = params

    def run(self):
        try:
            self.store.load_metadata()
            region = self.store.extract_region(**self.params)
            self.done.emit(region)
        except Exception as exc:
            self.done.emit({"error": f"{exc}\n{traceback.format_exc()}"})


class RegionView(QWidget):
    def __init__(self):
        super().__init__()

        self.centroid_store = None
        self.profile_store = None
        self.label = ""
        self.pending = False
        self.worker = None
        self.theme = "dark"

        self.source_combo = QComboBox()
        self.source_combo.addItem("centroid", "centroid")
        self.source_combo.addItem("profile", "profile")

        self.mz_center = self._spin(0.0, 5000.0, 4, 500.0, step=0.5)
        self.mz_half = self._spin(0.05, MAX_HALF_MZ, 3, 3.0, step=0.5)
        self.rt_center = self._spin(0.0, 1000.0, 3, 0.0, step=0.1)
        self.rt_half = self._spin(0.02, MAX_HALF_RT, 3, 0.5, step=0.1)

        self.bins = QSpinBox()
        self.bins.setRange(50, 2000)
        self.bins.setValue(400)
        self.bins.setSingleStep(50)

        self.log_check = QCheckBox("log")
        self.log_check.setChecked(True)
        self.points_check = QCheckBox("points")
        self.points_check.setChecked(True)
        self.surface_check = QCheckBox("surface")
        self.surface_check.setChecked(True)

        self.render_button = QPushButton("render")
        self.render_button.clicked.connect(self.render_region)

        self.status = QLabel("select a Sage ID (Spectra tab) or double-click a distribution candidate")
        self.status.setWordWrap(True)

        self.controls = QWidget()
        controls = self.controls
        c = QHBoxLayout(controls)
        c.setContentsMargins(4, 4, 4, 4)
        c.addWidget(QLabel("source"))
        c.addWidget(self.source_combo)
        c.addWidget(QLabel("m/z"))
        c.addWidget(self.mz_center)
        c.addWidget(QLabel("±"))
        c.addWidget(self.mz_half)
        c.addWidget(QLabel("RT"))
        c.addWidget(self.rt_center)
        c.addWidget(QLabel("±"))
        c.addWidget(self.rt_half)
        c.addWidget(QLabel("bins"))
        c.addWidget(self.bins)
        c.addWidget(self.log_check)
        c.addWidget(self.points_check)
        c.addWidget(self.surface_check)
        c.addWidget(self.render_button)
        c.addStretch(1)

        self.heatmap = pg.PlotWidget()
        self.heatmap.setLabel("bottom", "RT", units="min")
        self.heatmap.setLabel("left", "m/z")
        self.image = pg.ImageItem()
        self.image.setColorMap(_colormap())
        self.heatmap.addItem(self.image)
        self.colorbar = pg.ColorBarItem(colorMap=_colormap())
        self.colorbar.setImageItem(self.image, insert_in=self.heatmap.getPlotItem())

        if HAVE_GL:
            self.gl_view = gl.GLViewWidget()
            self.gl_view.setCameraPosition(distance=3.2, elevation=30, azimuth=-60)
            # GLSurfacePlotItem crashes in paint() if it is drawn before any z
            # grid is set (vertexNormals on None). Seed it with a valid 2x2 grid
            # and keep both items hidden until the first real region renders.
            self.surface = gl.GLSurfacePlotItem(
                x=np.array([0.0, 1.0], dtype=np.float32),
                y=np.array([0.0, 1.0], dtype=np.float32),
                z=np.zeros((2, 2), dtype=np.float32),
                shader="shaded",
                smooth=True,
            )
            self.surface.setVisible(False)
            self.scatter = gl.GLScatterPlotItem(pos=np.zeros((1, 3), dtype=np.float32), size=3.0)
            self.scatter.setVisible(False)
            self.gl_view.addItem(self.surface)
            self.gl_view.addItem(self.scatter)
            view_3d = self.gl_view
        else:
            self.gl_view = None
            view_3d = QLabel("3D view needs pyqtgraph OpenGL (pip install PyOpenGL).\nThe heatmap above shows the same region.")
            view_3d.setAlignment(Qt.AlignCenter)
            view_3d.setWordWrap(True)

        plots = QSplitter(Qt.Vertical)
        plots.addWidget(self.heatmap)
        plots.addWidget(view_3d)
        plots.setStretchFactor(0, 2)
        plots.setStretchFactor(1, 5)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(controls)
        layout.addWidget(self.status)
        layout.addWidget(plots, stretch=1)

        self.apply_theme("dark")

    def set_controls_visible(self, visible):
        # Hide the internal control row when driven from the workspace toolbar,
        # but keep the status line for render feedback.
        self.controls.setVisible(visible)

    def configure_window(self, mz_center, mz_half, rt_center, rt_half):
        for spin, value in (
            (self.mz_center, mz_center),
            (self.mz_half, mz_half),
            (self.rt_center, rt_center),
            (self.rt_half, rt_half),
        ):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

    def _spin(self, low, high, decimals, value, step=1.0):
        spin = QDoubleSpinBox()
        spin.setRange(low, high)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(value)
        return spin

    def apply_theme(self, theme):
        self.theme = theme
        pal = palette(theme)
        style_plot(self.heatmap, pal)

        if self.gl_view is not None:
            style_gl(self.gl_view, pal)

    # ---- targeting -------------------------------------------------------

    def set_target(self, centroid_store, profile_store, mz_min, mz_max, rt_start, rt_end, label=""):
        self.centroid_store = centroid_store
        self.profile_store = profile_store
        self.label = label

        self.mz_center.setValue((mz_min + mz_max) / 2.0)
        self.mz_half.setValue(min(MAX_HALF_MZ, max(0.1, (mz_max - mz_min) / 2.0)))
        self.rt_center.setValue(max(0.0, (rt_start + rt_end) / 2.0))
        self.rt_half.setValue(min(MAX_HALF_RT, max(0.05, (rt_end - rt_start) / 2.0)))

        has_profile = profile_store is not None
        self.source_combo.model().item(1).setEnabled(has_profile)

        if not has_profile and self.source_combo.currentData() == "profile":
            self.source_combo.setCurrentIndex(0)

        self.pending = True
        self.status.setText(
            f"{label}: ready — m/z {self.mz_center.value():.3f} ± {self.mz_half.value():.2f}, "
            f"RT {self.rt_center.value():.3f} ± {self.rt_half.value():.2f}. Switch here or press render."
        )

    def render_if_pending(self):
        if self.pending:
            self.pending = False
            self.render_region()

    # ---- rendering -------------------------------------------------------

    def current_store(self):
        if self.source_combo.currentData() == "profile":
            return self.profile_store

        return self.centroid_store

    def render_region(self):
        if self.worker is not None and self.worker.isRunning():
            return

        store = self.current_store()
        source = self.source_combo.currentData()

        if store is None:
            self.status.setText(f"no {source} mzML available for this selection")
            return

        mz_min = self.mz_center.value() - self.mz_half.value()
        mz_max = self.mz_center.value() + self.mz_half.value()
        rt_start = max(0.0, self.rt_center.value() - self.rt_half.value())
        rt_end = self.rt_center.value() + self.rt_half.value()

        self.render_button.setEnabled(False)
        self.status.setText(
            f"rendering {source}: m/z {mz_min:.3f}–{mz_max:.3f}, RT {rt_start:.3f}–{rt_end:.3f}…"
        )

        self._win = (mz_min, mz_max, rt_start, rt_end)
        params = dict(
            mz_min=mz_min,
            mz_max=mz_max,
            rt_start=rt_start,
            rt_end=rt_end,
            mz_bins=self.bins.value(),
            mode=source,
        )

        self.worker = RegionWorker(store, params)
        self.worker.done.connect(self.on_done)
        self.worker.start()

    def on_done(self, region):
        self.render_button.setEnabled(True)

        if region is None:
            self.status.setText("region extraction returned nothing")
            return

        if "error" in region:
            self.status.setText(f"region extraction failed: {region['error']}")
            return

        self.update_views(region)

    def update_views(self, region):
        mz_min, mz_max, rt_start, rt_end = self._win
        rts = region["rts"]
        mz_grid = region["mz_grid"]
        z = region["z"]

        if z.size == 0 or rts.size == 0:
            self.status.setText(f"no MS1 scans in RT {rt_start:.3f}–{rt_end:.3f} for this source")
            self.image.clear()
            return

        display = np.log1p(z) if self.log_check.isChecked() else z

        # ImageItem is (x=rows=RT, y=cols=m/z).
        self.image.setImage(display, autoLevels=True)
        rt_span = max(rts[-1] - rts[0], 1e-6) if rts.size > 1 else max(rt_end - rt_start, 1e-6)
        self.image.setRect(pg.QtCore.QRectF(rts[0], mz_min, rt_span, mz_max - mz_min))
        self.colorbar.setLevels((float(display.min()), float(display.max())))

        kind = "log(1+I)" if self.log_check.isChecked() else "intensity"
        self.status.setText(
            f"{self.label}  |  {self.source_combo.currentData()}  |  {rts.size} scans × "
            f"{mz_grid.size} m/z bins  |  peak {float(z.max()):.4g}  |  color: {kind}"
        )

        self.update_3d(rts, mz_grid, display)

    def update_3d(self, rts, mz_grid, display):
        if not HAVE_GL:
            return

        pal = palette(self.theme)
        zmax = float(display.max()) or 1.0
        zr = (display / zmax).astype(np.float32)

        # Normalize axes to a stable unit cube regardless of absolute ranges.
        x = np.linspace(-1.0, 1.0, display.shape[0]).astype(np.float32)
        y = np.linspace(-1.0, 1.0, display.shape[1]).astype(np.float32)

        if self.surface_check.isChecked():
            try:
                self.surface.setData(x=x, y=y, z=zr)
                self.surface.setColor((*[c / 255.0 for c in pal["surface"]], 0.55))
                self.surface.show()
            except Exception:
                pass
        else:
            self.surface.hide()

        if self.points_check.isChecked():
            # One point per grid cell above a small threshold, downsampled.
            xi, yi = np.meshgrid(np.arange(display.shape[0]), np.arange(display.shape[1]), indexing="ij")
            mask = zr > 0.02
            xs = x[xi[mask]]
            ys = y[yi[mask]]
            zs = zr[mask]

            if xs.size > MAX_SCATTER_POINTS:
                step = int(np.ceil(xs.size / MAX_SCATTER_POINTS))
                xs, ys, zs = xs[::step], ys[::step], zs[::step]

            pts = np.column_stack([xs, ys, zs]).astype(np.float32)
            color = (*[c / 255.0 for c in pal["points"]], 0.9)
            try:
                self.scatter.setData(pos=pts, color=color, size=3.0)
                self.scatter.show()
            except Exception:
                pass
        else:
            self.scatter.hide()
