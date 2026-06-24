import numpy as np
import pyqtgraph as pg


def clear_plot(plot):
    plot.clear()
    plot.showGrid(x=True, y=True, alpha=0.25)


def plot_spectrum(plot, mzs, intensities, title="", mz_label="m/z"):
    clear_plot(plot)
    plot.setTitle(title)
    plot.setLabel("bottom", mz_label)
    plot.setLabel("left", "intensity")

    mzs = np.asarray(mzs, dtype=np.float64)
    intensities = np.asarray(intensities, dtype=np.float64)

    if mzs.size == 0:
        return

    x = np.empty(mzs.size * 3, dtype=np.float64)
    y = np.empty(mzs.size * 3, dtype=np.float64)

    x[0::3] = mzs
    x[1::3] = mzs
    x[2::3] = np.nan

    y[0::3] = 0.0
    y[1::3] = intensities
    y[2::3] = np.nan

    plot.plot(x, y, pen=pg.mkPen(width=1))


def plot_points(plot, mzs, intensities, title=""):
    clear_plot(plot)
    plot.setTitle(title)
    plot.setLabel("bottom", "m/z")
    plot.setLabel("left", "intensity")

    mzs = np.asarray(mzs, dtype=np.float64)
    intensities = np.asarray(intensities, dtype=np.float64)

    if mzs.size == 0:
        return

    plot.plot(mzs, intensities, pen=None, symbol="o", symbolSize=5)


def add_profile_line(plot, mzs, intensities):
    mzs = np.asarray(mzs, dtype=np.float64)
    intensities = np.asarray(intensities, dtype=np.float64)

    if mzs.size == 0:
        return

    order = np.argsort(mzs)
    plot.plot(mzs[order], intensities[order], pen=pg.mkPen(width=1))


def plot_traces(plot, rts, traces, targets, title=""):
    clear_plot(plot)
    plot.setTitle(title)
    plot.setLabel("bottom", "RT", units="min")
    plot.setLabel("left", "extracted intensity")

    rts = np.asarray(rts, dtype=np.float64)

    if rts.size == 0:
        return

    for target, trace in zip(targets, traces):
        trace = np.asarray(trace, dtype=np.float64)

        if trace.size != rts.size:
            continue

        item = plot.plot(rts, trace, pen=pg.mkPen(width=2), name=f"{target:.4f}")
        item.setToolTip(f"{target:.4f}")
