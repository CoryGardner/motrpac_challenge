"""tfp — the tissue-fingerprint pipeline library: io, splits (frozen), models, conformal, transfer, batch, plots, report.
The evaluation rules every phase follows are in docs/EVALUATION_RULES.md."""
from . import config, conformal, io, models, plots, report, splits  # noqa: F401

__all__ = ["batch", "cli", "config", "conformal", "io", "models", "plots", "report", "splits", "transfer"]
