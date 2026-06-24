"""3D / region viewer for MS1 profile (and centroid) data.

Given a Sage-validated ID, this shows the surrounding RT x m/z x intensity
region two ways:

  * a 2D heatmap (RT vs m/z, intensity as color) -- the primary, always-on
    view, which stays readable no matter how dense the profile data is;
  * an optional 3D surface over the same gridded region, for inspecting peak
    shape, drawn with OpenGL when it is available.

Both are fed by MzmlStore.extract_region, which rasterizes the ragged
per-scan arrays onto a shared m/z grid. Nothing here fits or models peaks; it
only displays measured intensities.
"""

import traceback

import numpy as np
from PySide6.QtCore import Qt
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
    import pyqtgraph.opengl as gl

    HAVE_GL = True
except Exception:
    gl = None
    HAVE_GL = False


def _colormap():
    try:
        return pg.colormap.get("viridis")
    except Exception:
        return pg.colormap.get("CET-L9")


class RegionView(QWidget):
    def __init__(self):
        super().__init__()

        self.centroid_store = None
        self.profile_store = None
        self.label = ""
        self.last = None
        self.pending = False

        self.source_combo = QComboBox()
        self.source_combo.addItem("centroid", "centroid")
        self.source_combo.addItem("profile", "profile")

        self.mz_min = self._spin(100.0, 5000.0, 4, 400.0)
        self.mz_max = self._spin(100.0, 5000.0, 4, 410.0)
        self.rt_min = self._spin(0.0, 1000.0, 3, 0.0)
        self.rt_max = self._spin(0.0, 1000.0, 3, 5.0)

        self.bins = QSpinBox()
        self.bins.setRange(50, 4000)
        self.bins.setValue(600)
        self.bins.setSingleStep(50)

        self.log_check = QCheckBox("log intensity")
        self.log_check.setChecked(True)

        self.render_button = QPushButton("render region")
        self.render_button.clicked.connect(self.render_region)

        self.status = QLabel("select a Sage ID in the Spectra tab, or set a window and render")
        self.status.setWordWrap(True)

        controls = QWidget()
        controls_layout = QHBoxLayout(controls)
        controls_layout.setContentsMargins(4, 4, 4, 4)
        controls_layout.addWidget(QLabel("source"))
        controls_layout.addWidget(self.source_combo)
        controls_layout.addWidget(QLabel("m/z"))
        controls_layout.addWidget(self.mz_min)
        controls_layout.addWidget(QLabel("–"))
        controls_layout.addWidget(self.mz_max)
        controls_layout.addWidget(QLabel("RT"))
        controls_layout.addWidget(self.rt_min)
        controls_layout.addWidget(QLabel("–"))
        controls_layout.addWidget(self.rt_max)
        controls_layout.addWidget(QLabel("bins"))
        controls_layout.addWidget(self.bins)
        controls_layout.addWidget(self.log_check)
        controls_layout.addWidget(self.render_button)
        controls_layout.addStretch(1)

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
            self.gl_view.setCameraPosition(distance=3.0)
            self.surface = gl.GLSurfacePlotItem(
                shader="heightColor",
                computeNormals=False,
                smooth=False,
            )
            self.gl_view.addItem(self.surface)
            grid = gl.GLGridItem()
            grid.scale(0.1, 0.1, 0.1)
            self.gl_view.addItem(grid)
            view_3d = self.gl_view
        else:
            view_3d = QLabel(
                "3D surface needs pyqtgraph OpenGL (PyOpenGL). The 2D heatmap "
                "above shows the same region."
            )
            view_3d.setAlignment(Qt.AlignCenter)
            view_3d.setWordWrap(True)

        plots = QSplitter(Qt.Vertical)
        plots.addWidget(self.heatmap)
        plots.addWidget(view_3d)
        plots.setStretchFactor(0, 3)
        plots.setStretchFactor(1, 4)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(controls)
        layout.addWidget(self.status)
        layout.addWidget(plots, stretch=1)

    def _spin(self, low, high, decimals, value):
        spin = QDoubleSpinBox()
        spin.setRange(low, high)
        spin.setDecimals(decimals)
        spin.setValue(value)
        return spin

    def set_target(self, centroid_store, profile_store, mz_min, mz_max, rt_start, rt_end, label=""):
        self.centroid_store = centroid_store
        self.profile_store = profile_store
        self.label = label

        self.mz_min.setValue(float(mz_min))
        self.mz_max.setValue(float(mz_max))
        self.rt_min.setValue(max(0.0, float(rt_start)))
        self.rt_max.setValue(float(rt_end))

        has_profile = profile_store is not None
        self.source_combo.setItemData(1, "profile")
        self.source_combo.model().item(1).setEnabled(has_profile)

        if not has_profile and self.source_combo.currentData() == "profile":
            self.source_combo.setCurrentIndex(0)

        # Arm rather than render: extracting a region reads the whole file, so
        # we defer until this tab is actually shown (see render_if_pending).
        self.pending = True
        self.status.setText(
            f"{label}: ready — m/z {mz_min:.3f}–{mz_max:.3f}, RT {rt_start:.3f}–{rt_end:.3f}. "
            "Switch here or press 'render region'."
        )

    def render_if_pending(self):
        if self.pending:
            self.pending = False
            self.render_region()

    def current_store(self):
        if self.source_combo.currentData() == "profile":
            return self.profile_store

        return self.centroid_store

    def render_region(self):
        store = self.current_store()
        source = self.source_combo.currentData()

        if store is None:
            self.status.setText(f"no {source} mzML available for this selection")
            return

        mz_min = self.mz_min.value()
        mz_max = self.mz_max.value()
        rt_start = self.rt_min.value()
        rt_end = self.rt_max.value()

        if mz_max <= mz_min or rt_end <= rt_start:
            self.status.setText("window is empty: need mz_max > mz_min and rt_max > rt_min")
            return

        try:
            store.load_metadata()
            region = store.extract_region(
                mz_min=mz_min,
                mz_max=mz_max,
                rt_start=rt_start,
                rt_end=rt_end,
                mz_bins=self.bins.value(),
                mode=source,
            )
        except Exception as exc:
            self.status.setText(f"region extraction failed: {exc}\n{traceback.format_exc()}")
            return

        self.last = region
        self.update_views(region, mz_min, mz_max, rt_start, rt_end)

    def update_views(self, region, mz_min, mz_max, rt_start, rt_end):
        rts = region["rts"]
        mz_grid = region["mz_grid"]
        z = region["z"]

        if z.size == 0 or rts.size == 0:
            self.status.setText(
                f"no MS1 scans in RT {rt_start:.3f}–{rt_end:.3f} for this source"
            )
            self.image.clear()
            return

        display = np.log1p(z) if self.log_check.isChecked() else z

        # ImageItem is column-major in (x, y); rows are RT (x), cols are m/z (y).
        self.image.setImage(display, autoLevels=True)

        rt_span = max(rts[-1] - rts[0], 1e-6) if rts.size > 1 else 1.0
        scale_x = rt_span / max(display.shape[0], 1)
        scale_y = (mz_max - mz_min) / max(display.shape[1], 1)

        self.image.setRect(pg.QtCore.QRectF(rts[0], mz_min, rt_span, mz_max - mz_min))
        self.colorbar.setLevels((float(display.min()), float(display.max())))

        intensity_kind = "log(1+I)" if self.log_check.isChecked() else "intensity"
        peak = float(z.max())
        self.status.setText(
            f"{self.label}  |  {self.source_combo.currentData()}  |  "
            f"{rts.size} scans x {mz_grid.size} m/z bins  |  peak {peak:.4g}  |  color: {intensity_kind}"
        )

        self.update_surface(rts, mz_grid, display, rt_start, rt_end, mz_min, mz_max)

    def update_surface(self, rts, mz_grid, display, rt_start, rt_end, mz_min, mz_max):
        if not HAVE_GL:
            return

        # Normalize axes into a unit-ish cube so the camera framing is stable
        # regardless of the absolute RT / m/z / intensity magnitudes.
        x = np.linspace(-1.0, 1.0, display.shape[0])
        y = np.linspace(-1.0, 1.0, display.shape[1])

        zmax = float(display.max()) or 1.0
        zr = (display / zmax).astype(np.float32)

        try:
            self.surface.setData(x=x, y=y, z=zr)
        except Exception:
            pass
