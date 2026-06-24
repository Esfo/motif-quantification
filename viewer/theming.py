"""Light/dark theming for the viewer's pyqtgraph plots and GL views.

Toggled live from a toolbar button; styling is applied by walking the current
plot widgets rather than relying on global config, so it switches instantly.
"""

import pyqtgraph as pg

DARK = {
    "bg": "#101216",
    "fg": "#d8d8d8",
    "gl_bg": (16, 18, 22),
    "accent": "#4c72b0",
    "surface": (76, 114, 176),
    "points": (240, 240, 240),
}

LIGHT = {
    "bg": "#fafafa",
    "fg": "#202020",
    "gl_bg": (245, 245, 245),
    "accent": "#3060c0",
    "surface": (48, 96, 192),
    "points": (30, 30, 30),
}


def palette(name):
    return LIGHT if name == "light" else DARK


def style_plot(widget, pal):
    widget.setBackground(pal["bg"])

    plot_item = widget.getPlotItem() if hasattr(widget, "getPlotItem") else widget

    for name in ("left", "bottom", "right", "top"):
        axis = plot_item.getAxis(name)

        if axis is not None:
            axis.setPen(pg.mkPen(pal["fg"]))
            axis.setTextPen(pg.mkPen(pal["fg"]))


def style_gl(glview, pal):
    try:
        glview.setBackgroundColor(pal["gl_bg"])
    except Exception:
        pass
