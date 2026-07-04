"""Tab 2 — Protein viewing.

Shows identified protein sequences residue-by-residue (N→C), with the tryptic
peptides the search *attempted* drawn as outlined rectangles whose background is
coloured by q-value (best/green → worst/red on the shared q-value gradient).

Layout mirrors the MS Data tab: a file selector and a single protein list on the
left; the right side is a vertical split — panel 1 (top) over panel 2 (bottom),
a plain horizontal divider between them. Panel 1's header bar carries the **All**
button, the accepted-**FDR** spin box (governs the protein list), and a **colour
bar** legend mapping q-value/FDR → colour.

* **Panel 1** renders the selected protein horizontally, wrapping. Each base
  tryptic peptide (the protease's cut segments) is an outlined rectangle; its
  fill is the best q-value of any peptide identified over those residues in the
  current file. Segments too short/long to be searched carry no outline (plain
  letters). **All** combines identifications across every file. Double-clicking a
  peptide jumps to the MS Data tab focused on that identification (in All mode,
  the best file for that peptide).
* **Panel 2** renders the same protein *vertically* (N→C, top→bottom), one column
  per file, so the identification can be compared across files. Clicking a column
  loads that file into panel 1. It **zooms** (±/Ctrl-wheel): zooming out shrinks
  the residues to filled colour bands while the file names stay full size.

Sequences come from the project FASTA (``ViewerSession.protein_sequence``); the
reorganized tables carry none. Digestion parameters live in ``session.py`` and
mirror the Sage enzyme config in execution.xsh.
"""

import math

from PySide6.QtCore import Qt, QRectF, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QLinearGradient, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

try:
    from .session import digest_peptides, cleavage_sites, DIGEST_MIN_LEN, DIGEST_MAX_LEN
except ImportError:
    from session import digest_peptides, cleavage_sites, DIGEST_MIN_LEN, DIGEST_MAX_LEN


# q-value → colour on the shared gradient: best (low q) is green, worst (high q)
# is red, through yellow. Mapped on log10(q) between 1e-4 and 0.05 so the useful
# FDR band spans the full green→red sweep. q is None (identified, no q reported)
# counts as the best (green).
_Q_LOG_LO = -4.0
_Q_LOG_HI = math.log10(0.05)
_Q_LIN_MAX = 0.05


def q_to_t(q, mode="log"):
    """Map a q-value to the 0..1 colour position. ``log`` spreads the useful
    FDR band (1e-4…0.05) evenly; ``lin`` is linear in q over 0…0.05."""
    if q is None:
        return 0.0
    q = max(float(q), 1e-6)
    if mode == "lin":
        return min(1.0, max(0.0, q / _Q_LIN_MAX))
    return min(1.0, max(0.0, (math.log10(q) - _Q_LOG_LO) / (_Q_LOG_HI - _Q_LOG_LO)))


def color_from_t(t, alpha=190):
    # green (t=0, best) → yellow (t=0.5) → red (t=1, worst)
    r = int(255 * min(1.0, t * 2.0))
    g = int(255 * min(1.0, (1.0 - t) * 2.0))
    return QColor(r, g, 40, alpha)


def q_color(q, alpha=190, mode="log"):
    return color_from_t(q_to_t(q, mode), alpha)


THEMES = {
    "dark": {"fg": "#e6e6e6", "bg": "#101216", "muted": "#8a8f98", "panel": "#16181d"},
    "light": {"fg": "#202020", "bg": "#fafafa", "muted": "#606060", "panel": "#ffffff"},
}


def theme_is_dark(theme):
    return theme != "light"


# Protein-list sort metrics. Each is computed per (protein, file); the list sorts
# on the value, then the asc/desc toggle chooses the direction. "descending"
# means best/largest first for every metric (confidence uses -q so the most
# confident protein sorts first when descending).
SORT_OPTIONS = [
    "% Coverage",
    "Protein Length",
    "Total Identified Peptides",
    "FDR",
    "Spectral Count (PSMs)",
]


def base_segments(seq):
    """The fully-cleaved (0 missed cleavage) tryptic segments tiling ``seq``.

    Returns ``[(start, end, searchable), ...]`` covering the whole protein;
    ``searchable`` is False for segments outside the search's length bounds (the
    peptides the search never attempted → drawn plain, no rectangle)."""
    if not seq:
        return []
    sites = cleavage_sites(seq)
    out = []
    for i in range(len(sites) - 1):
        a, b = sites[i], sites[i + 1]
        searchable = DIGEST_MIN_LEN <= (b - a) <= DIGEST_MAX_LEN
        out.append((a, b, searchable))
    return out


def identified_spans(seq, qmap):
    """Occurrences of every identified peptide on ``seq``: ``[(start,end,plain,q)]``.

    Each in-silico peptide present in ``qmap`` is located on the protein (a
    peptide can repeat), yielding the residue spans and their q-values that drive
    both the colouring and the double-click → MS Data navigation."""
    spans = []
    if not qmap:
        return spans
    for _s, _e, pep in digest_peptides(seq):
        if pep not in qmap:
            continue
        q = qmap[pep]
        start = 0
        while True:
            idx = seq.find(pep, start)
            if idx < 0:
                break
            spans.append((idx, idx + len(pep), pep, q))
            start = idx + 1
    return spans


def residue_q(seq, qmap):
    """Per-residue best (min) q from identified peptides, and a hit mask.

    ``qmap`` maps plain peptide → q (or None). Residues covered by an
    identification carry its q (best kept on overlap); ``hit`` marks any covered
    residue so a None-q (unscored) identification still reads as identified."""
    n = len(seq)
    qs = [None] * n
    hit = [False] * n
    for start, end, _pep, q in identified_spans(seq, qmap):
        for i in range(start, end):
            hit[i] = True
            if q is not None and (qs[i] is None or q < qs[i]):
                qs[i] = q
    return qs, hit


def _segment_fill(qs, hit, a, b, mode="log"):
    """(fill_colour_or_None, covered_bool) for the residue range [a, b)."""
    seg_q = None
    covered = False
    for i in range(a, b):
        if hit[i]:
            covered = True
            if qs[i] is not None and (seg_q is None or qs[i] < seg_q):
                seg_q = qs[i]
    return (q_color(seg_q, mode=mode) if covered else None), covered


def _row_runs(a, b, cols):
    """Split residue range [a, b) into per-wrap-row contiguous runs:
    ``(row, col_start, count)`` — one rectangle per visual line."""
    i = a
    while i < b:
        row = i // cols
        col = i % cols
        count = min(b, (row + 1) * cols) - i
        yield row, col, count
        i += count


def _short_name(filename):
    name = filename or ""
    for suffix in (".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name if len(name) <= 40 else name[:38] + "…"


class VerticalColorBar(QWidget):
    """A vertical q-value/FDR → colour legend (matching ``q_color``) with a
    percentage y-axis. Best (green, low FDR) at the top, worst (red) at the
    bottom. Adaptive to light/dark. Occupies its own column beside panel 1."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._fg = QColor("#e6e6e6")
        self._bg = QColor("#101216")
        self._mode = "log"
        self.setFixedWidth(74)

    def set_theme(self, fg, bg):
        self._fg = QColor(fg)
        self._bg = QColor(bg)
        self.update()

    def set_color_mode(self, mode):
        self._mode = mode
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), self._bg)
        pad_top, pad_bot = 18, 12
        bar_w = 16
        x = 8
        h = self.height() - pad_top - pad_bot
        if h <= 0:
            painter.end()
            return
        # vertical gradient: top = best (green, 0.01%), bottom = worst (red, 5%)
        grad = QLinearGradient(0, pad_top, 0, pad_top + h)
        for k in range(21):
            t = k / 20.0
            grad.setColorAt(t, color_from_t(t, alpha=255))
        painter.fillRect(QRectF(x, pad_top, bar_w, h), grad)
        painter.setPen(QPen(self._fg, 1))
        painter.drawRect(QRectF(x, pad_top, bar_w, h))

        # title + percentage ticks (positions follow the lin/log mode)
        font = QFont()
        font.setPointSize(7)
        painter.setFont(font)
        painter.setPen(QPen(self._fg, 1))
        painter.drawText(QRectF(0, 2, self.width(), 14), Qt.AlignHCenter, "FDR")
        for q, label in ((1e-4, "0.01%"), (1e-3, "0.1%"), (1e-2, "1%"), (5e-2, "5%")):
            t = q_to_t(q, self._mode)
            y = pad_top + t * h
            painter.drawLine(int(x + bar_w), int(y), int(x + bar_w + 3), int(y))
            painter.drawText(QRectF(x + bar_w + 4, y - 7, self.width() - x - bar_w - 4, 14),
                             Qt.AlignLeft | Qt.AlignVCenter, label)
        painter.end()


class HorizontalSequenceView(QWidget):
    """Panel 1: the protein spelled out horizontally, wrapping, with peptide
    rectangles coloured by q-value. Double-click emits the residue index."""

    residue_double_clicked = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._seq = ""
        self._segments = []
        self._qs = []
        self._hit = []
        self._fg = QColor("#e6e6e6")
        self._panel_bg = QColor("#101216")
        self._muted = QColor("#808080")
        self._color_mode = "log"
        self.scroll_area = None
        self._pan = None
        self.setMinimumHeight(80)
        self.setCursor(Qt.OpenHandCursor)
        self._metrics()

    def set_color_mode(self, mode):
        self._color_mode = mode
        self.update()

    def _metrics(self):
        font = QFont("monospace")
        font.setStyleHint(QFont.Monospace)
        font.setPointSize(11)
        self._font = font
        fm = QFontMetrics(font)
        self.cw = max(12, fm.horizontalAdvance("W") + 6)
        self.ch = fm.height() + 8
        self.margin = 10

    def set_theme(self, fg, bg):
        self._fg = QColor(fg)
        self._panel_bg = QColor(bg)
        f = QColor(fg)
        f.setAlpha(110)
        self._muted = f
        self.update()

    def set_protein(self, seq, qmap):
        self._seq = seq or ""
        self._segments = base_segments(self._seq)
        self._qs, self._hit = residue_q(self._seq, qmap or {})
        self._recompute_height()
        self.update()

    def _cols(self):
        usable = max(1, self.width() - 2 * self.margin)
        return max(1, usable // self.cw)

    def _recompute_height(self):
        n = len(self._seq)
        if n == 0:
            self.setMinimumHeight(80)
            return
        rows = math.ceil(n / self._cols())
        self.setMinimumHeight(rows * self.ch + 2 * self.margin)

    def resizeEvent(self, event):
        self._recompute_height()
        super().resizeEvent(event)

    def _residue_at(self, x, y):
        cols = self._cols()
        col = int((x - self.margin) // self.cw)
        row = int((y - self.margin) // self.ch)
        if col < 0 or col >= cols or row < 0:
            return -1
        i = row * cols + col
        return i if 0 <= i < len(self._seq) else -1

    def mouseDoubleClickEvent(self, event):
        pos = event.position() if hasattr(event, "position") else event
        i = self._residue_at(pos.x(), pos.y())
        if i >= 0:
            self.residue_double_clicked.emit(i)

    # click-and-drag scrolls the (vertical) sequence, like grabbing the page
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and self.scroll_area is not None:
            gp = event.globalPosition() if hasattr(event, "globalPosition") else None
            y = gp.y() if gp is not None else event.globalY()
            self._pan = (y, self.scroll_area.verticalScrollBar().value())
            self.setCursor(Qt.ClosedHandCursor)

    def mouseMoveEvent(self, event):
        if self._pan is not None and self.scroll_area is not None:
            gp = event.globalPosition() if hasattr(event, "globalPosition") else None
            y = gp.y() if gp is not None else event.globalY()
            y0, v0 = self._pan
            self.scroll_area.verticalScrollBar().setValue(int(v0 - (y - y0)))

    def mouseReleaseEvent(self, event):
        self._pan = None
        self.setCursor(Qt.OpenHandCursor)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, False)
        painter.fillRect(self.rect(), self._panel_bg)
        if not self._seq:
            painter.setPen(self._fg)
            painter.drawText(self.rect(), Qt.AlignCenter,
                             "select a protein to view its sequence")
            painter.end()
            return

        cols = self._cols()
        painter.setFont(self._font)

        # 1) peptide rectangles (background fill + outline), behind the letters.
        # EVERY tryptic peptide gets the SAME outline; identified ones are filled
        # by q-value, the rest are just outlined (no fill).
        for a, b, _searchable in self._segments:
            fill, _covered = _segment_fill(self._qs, self._hit, a, b, self._color_mode)
            for row, col, count in _row_runs(a, b, cols):
                x = self.margin + col * self.cw
                y = self.margin + row * self.ch
                rect = QRectF(x, y + 2, count * self.cw, self.ch - 4)
                if fill is not None:
                    painter.fillRect(rect, fill)
                painter.setPen(QPen(self._fg, 1))
                painter.drawRect(rect)

        # 2) residue letters on top (clip to the exposed region for long proteins).
        exposed = event.rect()
        painter.setPen(self._fg)
        for i, chr_ in enumerate(self._seq):
            row = i // cols
            col = i % cols
            x = self.margin + col * self.cw
            y = self.margin + row * self.ch
            if y > exposed.bottom() or y + self.ch < exposed.top():
                continue
            painter.drawText(QRectF(x, y, self.cw, self.ch), Qt.AlignCenter, chr_)
        painter.end()


class VerticalMultiFileView(QWidget):
    """Panel 2: the protein spelled top→bottom (N→C), one column per file, each
    coloured by that file's identification. Zoomable — zooming out drops the
    letters to filled colour bands while the file names stay full size. Clicking
    a column emits its filename so panel 1 can switch to it."""

    file_clicked = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.HEADER_H = 116
        self._seq = ""
        self._segments = []
        self._files = []          # [(filename, qs, hit)]
        self._fg = QColor("#e6e6e6")
        self._panel_bg = QColor("#101216")
        self._muted = QColor("#808080")
        self._color_mode = "log"
        self._current = None
        self._scale = 1.0
        self._label_overhang = 40
        self.scroll_area = None
        self._pan = None          # (start_x, start_y, h0, v0)
        self._dragged = False
        self._metrics()
        self.setCursor(Qt.OpenHandCursor)

    def _metrics(self):
        font = QFont("monospace")
        font.setStyleHint(QFont.Monospace)
        font.setPointSize(10)
        self._font = font
        fm = QFontMetrics(font)
        self.base_ch = fm.height() + 2
        self.base_col_w = max(20, fm.horizontalAdvance("W") + 14)
        self.gap = 10
        self.margin = 10
        # file-name label font stays fixed regardless of zoom
        self._label_font = QFont(font)
        self._label_font.setPointSize(9)

    @property
    def ch(self):
        # zoom changes ONLY the vertical residue size; allow sub-pixel height so
        # the whole protein can shrink to fit.
        return max(0.1, self.base_ch * self._scale)

    @property
    def col_w(self):
        # columns FIT the panel width (zoom never changes width): divide the
        # usable width (minus room for the slanted labels) evenly across files.
        n = len(self._files)
        if n <= 0:
            return self.base_col_w
        # reserve the slanted-label overhang on BOTH sides so the centred block's
        # rightmost label never clips.
        avail = self.width() - 2 * self.margin - 2 * self._label_overhang - (n - 1) * self.gap
        return max(6.0, avail / n)

    def set_theme(self, fg, bg):
        self._fg = QColor(fg)
        self._panel_bg = QColor(bg)
        f = QColor(fg)
        f.setAlpha(110)
        self._muted = f
        self.update()

    def set_color_mode(self, mode):
        self._color_mode = mode
        self.update()

    def fit_to_height(self, viewport_h):
        """Scale so the whole protein column fits in ``viewport_h`` pixels."""
        n = len(self._seq)
        if n <= 0 or viewport_h <= 0:
            return
        avail = max(1, viewport_h - self.HEADER_H - self.margin)
        self.set_scale(avail / (n * self.base_ch))

    def set_protein(self, seq, files_qmaps, current=None):
        """``files_qmaps`` = ``[(filename, qmap), ...]`` over every file."""
        self._seq = seq or ""
        self._segments = base_segments(self._seq)
        self._current = current
        self._files = []
        for filename, qmap in files_qmaps:
            qs, hit = residue_q(self._seq, qmap or {})
            self._files.append((filename, qs, hit))
        # Header tall enough for the longest 45°-slanted file name (so names are
        # never cut off at the top).
        fm = QFontMetrics(self._label_font)
        maxw = max((fm.horizontalAdvance(_short_name(f)) for f, _q, _h in self._files),
                   default=0)
        # header tall enough, and right padding wide enough, for the longest
        # 45°-slanted file name (so names are never cut off top or right).
        self.HEADER_H = int(16 + maxw * 0.7071)
        self._label_overhang = int(maxw * 0.7071) + 8
        self._recompute_size()
        self.update()

    def set_current(self, current):
        self._current = current
        self.update()

    def set_scale(self, scale):
        # allow very small scales so a whole (long) protein can be zoomed to fit
        self._scale = min(4.0, max(0.01, scale))
        self._recompute_size()
        self.update()

    def zoom(self, factor):
        self.set_scale(self._scale * factor)

    def _col_pitch(self):
        return self.col_w + self.gap

    def _content_x0(self):
        # centre the block of columns within the panel (the right label-overhang
        # gap is balanced with an equal gap on the left).
        content = len(self._files) * self.col_w + max(0, len(self._files) - 1) * self.gap
        return max(self.margin, (self.width() - content) / 2.0)

    def _recompute_size(self):
        n = len(self._seq)
        # Width fits the panel (columns divide the viewport width), so keep the
        # minimum width tiny and let widgetResizable stretch us to the viewport;
        # only the height (from the zoomable residue size) drives scrolling.
        height = self.HEADER_H + n * self.ch + self.margin
        self.setMinimumSize(0, int(height))
        self.setMaximumWidth(16777215)

    def wheelEvent(self, event):
        # Ctrl-wheel zooms; a plain wheel scrolls the enclosing scroll area.
        if event.modifiers() & Qt.ControlModifier:
            self.zoom(1.15 if event.angleDelta().y() > 0 else 1 / 1.15)
            event.accept()
        else:
            event.ignore()

    def _global_xy(self, event):
        gp = event.globalPosition() if hasattr(event, "globalPosition") else None
        return (gp.x(), gp.y()) if gp is not None else (event.globalX(), event.globalY())

    def mousePressEvent(self, event):
        if event.button() != Qt.LeftButton:
            return
        self._dragged = False
        gx, gy = self._global_xy(event)
        if self.scroll_area is not None:
            self._pan = (gx, gy,
                         self.scroll_area.horizontalScrollBar().value(),
                         self.scroll_area.verticalScrollBar().value())
        self._press_x = event.position().x() if hasattr(event, "position") else event.x()
        self.setCursor(Qt.ClosedHandCursor)

    def mouseMoveEvent(self, event):
        if self._pan is None or self.scroll_area is None:
            return
        gx, gy = self._global_xy(event)
        x0, y0, h0, v0 = self._pan
        if abs(gx - x0) + abs(gy - y0) > 4:
            self._dragged = True
        self.scroll_area.horizontalScrollBar().setValue(int(h0 - (gx - x0)))
        self.scroll_area.verticalScrollBar().setValue(int(v0 - (gy - y0)))

    def mouseReleaseEvent(self, event):
        self.setCursor(Qt.OpenHandCursor)
        self._pan = None
        # a click (no drag) selects the column under the cursor → loads that file
        if not self._dragged and self._files:
            idx = int((self._press_x - self._content_x0()) // self._col_pitch())
            if 0 <= idx < len(self._files):
                self.file_clicked.emit(self._files[idx][0])

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), self._panel_bg)
        if not self._seq or not self._files:
            painter.setPen(self._fg)
            painter.drawText(self.rect(), Qt.AlignCenter,
                             "select a protein to compare across files")
            painter.end()
            return

        exposed = event.rect()
        pitch = self._col_pitch()
        ch = self.ch
        col_w = self.col_w
        base_x0 = self._content_x0()
        draw_letters = ch >= 11 and col_w >= 11
        for fidx, (filename, qs, hit) in enumerate(self._files):
            x0 = base_x0 + fidx * pitch

            if filename == self._current:
                painter.fillRect(QRectF(x0 - 3, 0, col_w + 6, self.height()),
                                 QColor(120, 150, 210, 45))

            # 45°-slanted file label, fixed size (independent of zoom)
            painter.save()
            painter.translate(x0 + col_w / 2.0, self.HEADER_H - 8)
            painter.rotate(-45)
            painter.setPen(self._fg)
            painter.setFont(self._label_font)
            painter.drawText(0, 0, _short_name(filename))
            painter.restore()

            painter.setFont(self._font)
            for a, b, searchable in self._segments:
                fill, covered = _segment_fill(qs, hit, a, b, self._color_mode)
                y_top = self.HEADER_H + a * ch
                y_bot = self.HEADER_H + b * ch
                if y_bot < exposed.top() or y_top > exposed.bottom():
                    continue
                # Keep a positive height even when zoomed way out, so every
                # peptide's outline stays visible (a tiny inset only when there's
                # room). Otherwise (b-a)*ch - 2 could go negative → nothing drawn.
                seg_h = (b - a) * ch
                inset = 1 if seg_h > 5 else 0
                rect = QRectF(x0, y_top + inset, col_w, max(1.0, seg_h - 2 * inset))
                if covered:
                    painter.fillRect(rect, fill)
                # every peptide gets the SAME outline; identified ones are filled
                # with their FDR colour, the rest are just outlined (no fill).
                if seg_h >= 3:
                    painter.setPen(QPen(self._fg, 1))
                    painter.drawRect(rect)

            if draw_letters:
                for i, chr_ in enumerate(self._seq):
                    y = self.HEADER_H + i * ch
                    if y > exposed.bottom() or y + ch < exposed.top():
                        continue
                    painter.setPen(self._fg)
                    painter.drawText(QRectF(x0, y, col_w, ch), Qt.AlignCenter, chr_)
        painter.end()


class ProteinsTab(QWidget):
    """Tab 2 assembled: file selector + protein list on the left; panel 1 (with
    the All button, FDR spin box and colour-bar legend) over panel 2 on the right.

    ``on_navigate_to_ms(filename, protein_id, peptide_plain)`` is set by the main
    window to jump to the MS Data tab on a panel-1 double-click.
    """

    def __init__(self, session, theme="dark"):
        super().__init__()
        self.session = session
        self.theme = theme
        self.current_file = None
        self.current_protein = None
        self._combined = False   # panel 1 "All" mode
        self._color_mode = "log"  # q-value colour scale: "log" or "lin"
        self._sort_desc = True
        self._len_cache = {}
        self._cov_cache = {}
        self._panel1_spans = []  # identified spans currently shown in panel 1
        self.on_navigate_to_ms = None
        self.on_theme_toggle = None

        self._build()
        self._populate_files()
        self._apply_theme(theme)
        self._refresh_protein_list()

    # ---- construction ----------------------------------------------------

    def _build(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)

        # left column: file selector + protein list
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(2, 2, 2, 2)

        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self._on_file_changed)
        left_layout.addWidget(self.file_combo)

        # Sort-by row: "Sort By" + metric dropdown + asc/desc toggle button.
        sort_row = QHBoxLayout()
        sort_row.setContentsMargins(0, 0, 0, 0)
        sort_row.addWidget(QLabel("Sort By"))
        self.sort_combo = QComboBox()
        for label in SORT_OPTIONS:
            self.sort_combo.addItem(label)
        self.sort_combo.currentIndexChanged.connect(self._on_sort_changed)
        sort_row.addWidget(self.sort_combo, stretch=1)
        self.sort_dir_btn = QPushButton("▼")
        self.sort_dir_btn.setFixedWidth(30)
        self.sort_dir_btn.setToolTip("Descending / ascending")
        self._sort_desc = True
        self.sort_dir_btn.clicked.connect(self._on_sort_dir_toggle)
        sort_row.addWidget(self.sort_dir_btn)
        left_layout.addLayout(sort_row)

        # Protein table: the sorted value on the left, protein name on the right.
        self.protein_table = QTableWidget(0, 2)
        self.protein_table.setHorizontalHeaderLabels([SORT_OPTIONS[0], "Protein"])
        self.protein_table.verticalHeader().setVisible(False)
        self.protein_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.protein_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.protein_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.protein_table.setShowGrid(False)
        self.protein_table.setWordWrap(False)
        hdr = self.protein_table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(1, QHeaderView.Stretch)
        self.protein_table.setStyleSheet(
            "QTableWidget::item:selected,"
            "QTableWidget::item:selected:!active"
            " { background-color: #2f6fb3; color: white; }"
        )
        self.protein_table.itemSelectionChanged.connect(self._on_protein_selected)
        left_layout.addWidget(self.protein_table, stretch=1)

        # right side: panel 1 over panel 2
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)

        # panel-1 header bar: All | FDR% | title | colour bar
        bar = QHBoxLayout()
        bar.setContentsMargins(4, 2, 4, 2)
        self.all_btn = QPushButton("All")
        self.all_btn.setToolTip("Combine identifications across every file into a "
                                "single consensus colouring")
        self.all_btn.setCheckable(True)
        self.all_btn.setFixedWidth(48)
        self.all_btn.clicked.connect(self._on_all_clicked)
        bar.addWidget(self.all_btn)

        # FDR acceptance criteria (percent) — to the right of All. Governs the
        # protein list; defaults to the search's q_max. Scrollable spin box.
        self.fdr_edit = QDoubleSpinBox()
        self.fdr_edit.setDecimals(2)
        self.fdr_edit.setRange(0.0, 100.0)
        self.fdr_edit.setSingleStep(0.1)
        self.fdr_edit.setValue(self._default_fdr_percent())
        self.fdr_edit.setFixedWidth(70)
        self.fdr_edit.setToolTip("Proteins are listed when identified at or below "
                                 "this FDR (percent)")
        self.fdr_edit.valueChanged.connect(self._on_fdr_changed)
        bar.addWidget(self.fdr_edit)
        bar.addWidget(QLabel("% FDR"))

        # q-value colour-scale toggle (like MS Data's Lin/Log colour button), but
        # for the FDR colours.
        self.colormode_btn = QPushButton("Log Color")
        self.colormode_btn.setFixedWidth(84)
        self.colormode_btn.setToolTip("Toggle the FDR colour scale between "
                                      "logarithmic and linear in q-value")
        self.colormode_btn.clicked.connect(self._on_colormode_toggle)
        bar.addWidget(self.colormode_btn)

        # Light/Dark toggle — to the right of the FDR bits. Calls back into the
        # main window so the theme stays in sync across tabs.
        self.theme_btn = QPushButton("Light Mode" if theme_is_dark(self.theme) else "Dark Mode")
        self.theme_btn.setFixedWidth(90)
        self.theme_btn.clicked.connect(lambda: self.on_theme_toggle and self.on_theme_toggle())
        bar.addWidget(self.theme_btn)

        self.p1_title = QLabel("")
        bar.addSpacing(10)
        bar.addWidget(self.p1_title)
        bar.addStretch(1)
        right_layout.addLayout(bar)

        self.panel1 = HorizontalSequenceView()
        self.panel1.residue_double_clicked.connect(self._on_panel1_double_click)
        self.p1_scroll = QScrollArea()
        self.p1_scroll.setWidgetResizable(True)
        self.p1_scroll.setWidget(self.panel1)
        self.p1_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.panel1.scroll_area = self.p1_scroll

        # Panel 1 gets the sequence on the left and the vertical q-value/FDR colour
        # bar in its own column on the right (no overlap with the sequence).
        p1_row = QWidget()
        p1_h = QHBoxLayout(p1_row)
        p1_h.setContentsMargins(0, 0, 0, 0)
        p1_h.setSpacing(0)
        p1_h.addWidget(self.p1_scroll, stretch=1)
        self.color_bar = VerticalColorBar()
        p1_h.addWidget(self.color_bar)

        self.panel2 = VerticalMultiFileView()
        self.panel2.file_clicked.connect(self._on_panel2_file)
        self.p2_scroll = QScrollArea()
        self.p2_scroll.setWidgetResizable(True)
        self.p2_scroll.setWidget(self.panel2)
        self.panel2.scroll_area = self.p2_scroll

        # Vertical splitter (panel 1 over panel 2) — resizable + persisted, so the
        # panel sizes are remembered like the MS Data tab's docks.
        self.v_splitter = QSplitter(Qt.Vertical)
        self.v_splitter.setHandleWidth(2)   # a plain horizontal divider, like MS Data
        self.v_splitter.addWidget(p1_row)
        self.v_splitter.addWidget(self.p2_scroll)
        self.v_splitter.setStretchFactor(0, 1)
        self.v_splitter.setStretchFactor(1, 1)
        right_layout.addWidget(self.v_splitter, stretch=1)

        # Horizontal splitter (left list column | right panels) — also resizable
        # and persisted.
        self.h_splitter = QSplitter(Qt.Horizontal)
        self.h_splitter.setHandleWidth(2)
        left.setMinimumWidth(160)
        self.h_splitter.addWidget(left)
        self.h_splitter.addWidget(right)
        self.h_splitter.setStretchFactor(0, 0)
        self.h_splitter.setStretchFactor(1, 1)
        self.h_splitter.setSizes([280, 900])
        root.addWidget(self.h_splitter)

    def _default_fdr_percent(self):
        q_max = self.session.summary().get("q_max") if self.session else None
        try:
            return max(0.0, float(q_max) * 100.0)
        except (TypeError, ValueError):
            return 1.0

    def _on_colormode_toggle(self):
        self._color_mode = "lin" if self._color_mode == "log" else "log"
        self.colormode_btn.setText("Lin Color" if self._color_mode == "lin" else "Log Color")
        self.panel1.set_color_mode(self._color_mode)
        self.panel2.set_color_mode(self._color_mode)
        self.color_bar.set_color_mode(self._color_mode)

    # ---- layout persistence ---------------------------------------------
    # The panel sizes are persisted by the main window via h_splitter/v_splitter
    # (mirrors how the MS Data dock layout is saved/restored).

    def splitter_states(self):
        return self.h_splitter.saveState(), self.v_splitter.saveState()

    def restore_splitter_states(self, h_state, v_state):
        if h_state is not None:
            self.h_splitter.restoreState(h_state)
        if v_state is not None:
            self.v_splitter.restoreState(v_state)

    # ---- data ------------------------------------------------------------

    def _all_files(self):
        return [row.get("filename", "") for row in self.session.files()
                if row.get("filename")]

    def _populate_files(self):
        self.file_combo.blockSignals(True)
        self.file_combo.clear()
        for name in self._all_files():
            self.file_combo.addItem(name, name)
        self.file_combo.blockSignals(False)
        self.current_file = self.file_combo.currentData()

    def _fdr_threshold(self):
        return max(0.0, self.fdr_edit.value()) / 100.0

    def _identified_proteins(self, filename):
        """Proteins identified at/below the FDR in ``filename``, ordered by the
        currently-selected sort metric and direction."""
        thr = self._fdr_threshold()
        out = []
        for row in self.session.file_proteins(filename or ""):
            pid = row.get("protein_id", "")
            if not pid:
                continue
            q = _safe_float(row.get("protein_q"))
            if q is not None and q > thr:
                continue
            out.append((pid, row))
        return self._sort_proteins(out, filename)

    # ---- sorting ---------------------------------------------------------

    def _protein_len(self, pid):
        if pid not in self._len_cache:
            seq = self.session.protein_sequence(pid)
            self._len_cache[pid] = len(seq) if seq else 0
        return self._len_cache[pid]

    def _protein_coverage(self, pid, filename):
        key = (pid, filename)
        if key not in self._cov_cache:
            seq = self.session.protein_sequence(pid)
            if not seq:
                self._cov_cache[key] = 0.0
            else:
                _qs, hit = residue_q(seq, self.session.peptide_q_for_file(filename or ""))
                self._cov_cache[key] = 100.0 * sum(hit) / len(seq)
        return self._cov_cache[key]

    def _metric_raw(self, metric, pid, row, filename):
        """The metric's actual value (for display), or None when unavailable."""
        if metric == "% Coverage":
            return self._protein_coverage(pid, filename)
        if metric == "Protein Length":
            return self._protein_len(pid)
        if metric == "Total Identified Peptides":
            return _safe_int(row.get("n_peptides"))
        if metric == "FDR":
            return _safe_float(row.get("protein_q"))
        if metric == "Spectral Count (PSMs)":
            return _safe_int(row.get("n_psms"))
        return None

    def _format_metric(self, metric, val):
        if val is None:
            return "—"
        if metric == "% Coverage":
            return f"{val:.1f}%"
        if metric == "FDR":
            return f"{val * 100:.2f}%"
        return str(int(val))

    def _metric_value(self, metric, pid, row, filename):
        """Sort key; larger = 'better' so descending puts the best first."""
        raw = self._metric_raw(metric, pid, row, filename)
        if metric == "FDR":                 # lower FDR is better
            return -(raw if raw is not None else 1.0)
        return raw if raw is not None else 0

    def _sort_proteins(self, items, filename):
        metric = self.sort_combo.currentText() if hasattr(self, "sort_combo") else SORT_OPTIONS[0]
        keyed = [(self._metric_value(metric, pid, row, filename), pid, row)
                 for pid, row in items]
        try:
            keyed.sort(key=lambda t: t[0], reverse=self._sort_desc)
        except TypeError:
            keyed.sort(key=lambda t: str(t[0]), reverse=self._sort_desc)
        return [(pid, row) for _v, pid, row in keyed]

    # ---- events ----------------------------------------------------------

    def _on_sort_changed(self):
        self._refresh_protein_list(scroll_top=True)

    def _on_sort_dir_toggle(self):
        self._sort_desc = not self._sort_desc
        self.sort_dir_btn.setText("▼" if self._sort_desc else "▲")
        self._refresh_protein_list(scroll_top=True)

    def _on_file_changed(self):
        self.current_file = self.file_combo.currentData()
        self._refresh_protein_list()

    def _on_fdr_changed(self):
        self._refresh_protein_list()

    def _refresh_protein_list(self, scroll_top=False):
        prev = self.current_protein
        metric = self.sort_combo.currentText()
        rows = self._identified_proteins(self.current_file)

        self.protein_table.blockSignals(True)
        self.protein_table.setHorizontalHeaderLabels([metric, "Protein"])
        self.protein_table.setRowCount(len(rows))
        restore_row = -1
        for i, (pid, row) in enumerate(rows):
            val = self._format_metric(metric, self._metric_raw(metric, pid, row, self.current_file))
            val_item = QTableWidgetItem(val)
            val_item.setData(Qt.UserRole, pid)
            name_item = QTableWidgetItem(pid)
            name_item.setData(Qt.UserRole, pid)
            self.protein_table.setItem(i, 0, val_item)
            self.protein_table.setItem(i, 1, name_item)
            if prev is not None and pid == prev:
                restore_row = i
        self.protein_table.blockSignals(False)

        if restore_row >= 0:
            self.protein_table.selectRow(restore_row)
        # On a sort change the user wants the view to jump back to the top,
        # regardless of which protein stays selected.
        if scroll_top:
            self.protein_table.scrollToTop()
        elif restore_row >= 0:
            self.protein_table.scrollToItem(self.protein_table.item(restore_row, 0))
        self._update_panels()

    def _on_protein_selected(self):
        items = self.protein_table.selectedItems()
        self.current_protein = items[0].data(Qt.UserRole) if items else None
        self._update_panels()

    def _on_all_clicked(self):
        self._combined = self.all_btn.isChecked()
        self._update_panel1()

    def _on_panel2_file(self, filename):
        # clicking a panel-2 column brings that file up in panel 1
        self._combined = False
        self.all_btn.setChecked(False)
        idx = self.file_combo.findData(filename)
        if idx >= 0:
            self.file_combo.setCurrentIndex(idx)   # triggers list refresh
        else:
            self.current_file = filename
            self._update_panels()

    def _on_panel1_double_click(self, residue):
        """Double-clicking a peptide in panel 1 → jump to the MS Data tab for
        that identification. In All mode, pick the best file for the peptide."""
        if not self.current_protein or self.on_navigate_to_ms is None:
            return
        # identified peptides covering this residue; pick the best q.
        covering = [s for s in self._panel1_spans if s[0] <= residue < s[1]]
        if not covering:
            return
        covering.sort(key=lambda s: (s[3] is None, s[3] if s[3] is not None else 0.0))
        _a, _b, plain, _q = covering[0]
        filename = self.current_file
        if self._combined:
            filename = self._best_file_for_peptide(plain) or filename
        self.on_navigate_to_ms(filename, self.current_protein, plain)

    def _best_file_for_peptide(self, plain):
        best_file, best_q = None, None
        for filename in self._all_files():
            qmap = self.session.peptide_q_for_file(filename)
            if plain not in qmap:
                continue
            q = qmap[plain]
            if best_file is None or (q is not None and (best_q is None or q < best_q)):
                best_file, best_q = filename, q
        return best_file

    # ---- rendering -------------------------------------------------------

    def _combined_qmap(self):
        """Best q per identified peptide across every file (the All view)."""
        combined = {}
        for filename in self._all_files():
            for pep, q in self.session.peptide_q_for_file(filename).items():
                if pep not in combined:
                    combined[pep] = q
                elif q is not None and (combined[pep] is None or q < combined[pep]):
                    combined[pep] = q
        return combined

    def _update_panels(self):
        self._update_panel1()
        self._update_panel2()

    def _update_panel1(self):
        seq = self.session.protein_sequence(self.current_protein) if self.current_protein else None
        if not seq:
            self.panel1.set_protein("", {})
            self._panel1_spans = []
            if self.current_protein and not self.session.has_sequences:
                self.p1_title.setText("no FASTA found — cannot show sequences")
            else:
                self.p1_title.setText("")
            return
        if self._combined:
            qmap = self._combined_qmap()
            self.p1_title.setText(f"{self.current_protein}   ·   All files ({len(seq)} aa)")
        else:
            qmap = self.session.peptide_q_for_file(self.current_file or "")
            self.p1_title.setText(f"{self.current_protein}   ·   "
                                  f"{_short_name(self.current_file or '')} ({len(seq)} aa)")
        self._panel1_spans = identified_spans(seq, qmap)
        self.panel1.set_protein(seq, qmap)

    def _update_panel2(self):
        seq = self.session.protein_sequence(self.current_protein) if self.current_protein else None
        if not seq:
            self.panel2.set_protein("", [])
            return
        files_qmaps = [(f, self.session.peptide_q_for_file(f)) for f in self._all_files()]
        self.panel2.set_protein(seq, files_qmaps, current=self.current_file)

    # ---- theme -----------------------------------------------------------

    def apply_theme(self, theme):
        self._apply_theme(theme)

    def _apply_theme(self, theme):
        self.theme = theme
        pal = THEMES.get(theme, THEMES["dark"])
        fg, bg = pal["fg"], pal["bg"]
        # Theme ONLY the sequence panels (+ their viewports) and the colour bar.
        # The left file-list column keeps the default app palette so it does not
        # turn dark — the theming is "just for the panels".
        self.panel1.set_theme(fg, bg)
        self.panel2.set_theme(fg, bg)
        self.color_bar.set_theme(fg, bg)
        # Panels get the same 1px frame the MS Data plots show (dark in light
        # mode, light in dark mode = the fg colour).
        for area in (self.p1_scroll, self.p2_scroll):
            area.setStyleSheet(
                f"QScrollArea {{ background: {bg}; border: 1px solid {fg}; }}")
            area.viewport().setStyleSheet(f"background: {bg};")
        self.color_bar.setStyleSheet(f"border: 1px solid {fg}; border-left: none;")
        if getattr(self, "theme_btn", None) is not None:
            self.theme_btn.setText("Light Mode" if theme_is_dark(theme) else "Dark Mode")


def _safe_float(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value):
    f = _safe_float(value)
    return int(f) if f is not None else 0
