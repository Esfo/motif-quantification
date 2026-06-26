#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication

try:
    from .main_window import MainWindow
except ImportError:
    from main_window import MainWindow


def parse_args():
    parser = argparse.ArgumentParser(
        description="Open the motif-quantification MS evidence viewer."
    )

    parser.add_argument(
        "--reorganized",
        type=Path,
        default=None,
        help="optional searches/reorganized directory; if omitted, open it from File ▸ Open in the app",
    )

    parser.add_argument(
        "--distribution-db",
        type=Path,
        default=None,
        help="optional distributions SQLite database",
    )

    parser.add_argument(
        "--centroid-dir",
        type=Path,
        default=None,
        help="optional directory containing centroid mzML files; defaults to manifest mzml_dir",
    )

    parser.add_argument(
        "--profile-dir",
        type=Path,
        default=None,
        help="optional directory containing matching profile mzML files",
    )

    parser.add_argument(
        "--xics-ppm",
        type=float,
        default=10.0,
        help="ppm tolerance for MS1 isotope XIC extraction",
    )

    parser.add_argument(
        "--xics-rt-window",
        type=float,
        default=0.8,
        help="RT window in minutes on each side of the selected ID",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Let Ctrl+C in the terminal actually quit: restore the default SIGINT
    # behavior and wake the Qt loop periodically so Python can run the handler.
    import signal

    from PySide6.QtCore import QTimer

    signal.signal(signal.SIGINT, signal.SIG_DFL)

    app = QApplication(sys.argv)
    app.setApplicationName("Motif Quantification Viewer")

    heartbeat = QTimer()
    heartbeat.start(200)
    heartbeat.timeout.connect(lambda: None)

    window = MainWindow(
        reorganized=args.reorganized,
        distribution_db=args.distribution_db,
        centroid_dir=args.centroid_dir,
        profile_dir=args.profile_dir,
        xics_ppm=args.xics_ppm,
        xics_rt_window=args.xics_rt_window,
    )

    # Only fall back to a default size if no saved geometry was restored --
    # otherwise this would clobber the remembered window/dock layout every launch.
    if not getattr(window, "_geometry_restored", False):
        window.resize(1600, 950)
    window.show()

    raise SystemExit(app.exec())


if __name__ == "__main__":
    main()
