"""Tab 4 — Motif Quantification.

Differential-expression view at the **motif** level: proteins are grouped by the
skeleton motif they share (``index-motifs.py``) and quantified as the SUM of
their member proteins' AUC quantities (computed off-GUI in ``quantify.py`` and
read from ``motif-sets/``).

It is deliberately the Quantitative-Comparisons tab with a different feature
vocabulary: ``MotifTab`` subclasses ``QuantTab`` and swaps in a motif model, a
"motif"/"observed" column pair, and a **minimum observed proteins** control in
place of the Peptides/Proteins switch. All of the faceting, fold-change scatter,
nested pivot, normalization, replicate handling, and theming are inherited
unchanged, so the two tabs stay in lockstep.
"""

from PySide6.QtWidgets import QLabel, QSpinBox, QVBoxLayout
from PySide6.QtCore import Qt

try:
    from .quant_tab import QuantTab
    from .motif_quant_model import MotifQuantModel
except ImportError:
    from quant_tab import QuantTab
    from motif_quant_model import MotifQuantModel


class MotifTab(QuantTab):
    SETTINGS_KEY = "motif_state"
    _default_level = "motif"

    def _make_model(self):
        return MotifQuantModel(self.session)

    # ---- availability: needs the motif-sets/ folder from quantify.py -----

    def _build_ui(self):
        self.min_observed = int(self._saved.get("min_observed", 2) or 2)
        if not self.model.is_available():
            # No motif quantities on disk yet — explain how to produce them
            # instead of rendering an empty DE view.
            self._active = False
            outer = QVBoxLayout(self)
            outer.setContentsMargins(12, 12, 12, 12)
            msg = QLabel(
                "No motif quantities found.\n\n"
                "Run the motif grouping stage to produce them:\n"
                "    python quantify.py --project <PROJECT> "
                "--motif-index <MOTIF INDEX>\n\n"
                "(the run_quantify stage in execution.xsh), which writes "
                "motif-sets/ next to searches/ and distributions/.")
            msg.setAlignment(Qt.AlignCenter)
            msg.setWordWrap(True)
            outer.addWidget(msg)
            return
        super()._build_ui()

    # ---- feature vocabulary ---------------------------------------------

    def _feature_column_label(self):
        return "motif"

    def _second_column_label(self):
        return "observed"

    def _second_column_value(self, feat):
        return str(self.model.observed_count(feat))

    def _visible_features(self):
        matrix = self._matrix()
        feats = [f for f in matrix.keys()
                 if self.model.observed_count(f) >= self.min_observed]
        return sorted(feats)

    # ---- min-observed control (replaces the Peptides/Proteins switch) ----

    def _build_level_controls(self, bar):
        bar.addWidget(QLabel("Min observed proteins:"))
        self.min_spin = QSpinBox()
        self.min_spin.setMinimum(2)
        self.min_spin.setMaximum(max(2, self.model.max_observed()))
        self.min_spin.setValue(max(2, self.min_observed))
        self.min_spin.setToolTip(
            "Only show motifs grouping at least this many quantified proteins. "
            "2 keeps every genuine multi-protein group.")
        self.min_spin.valueChanged.connect(self._on_min_observed_changed)
        bar.addWidget(self.min_spin)

    def _on_min_observed_changed(self, value):
        self.min_observed = int(value)
        self.selected_feature = None
        self._refresh_table()
        self._refresh_fold_change()
        self._auto_select_first()
        self._save_state()

    # ---- persist the extra control --------------------------------------

    def _save_state(self):
        super()._save_state()
        if getattr(self, "_restoring", False) or not self._active:
            return
        import json
        raw = self.settings.value(self.SETTINGS_KEY)
        try:
            state = json.loads(raw) if raw else {}
        except Exception:
            state = {}
        state["min_observed"] = self.min_observed
        self.settings.setValue(self.SETTINGS_KEY, json.dumps(state))
