"""results_root(): where the read side (site export, tests, figures) finds the result tables."""
from pathlib import Path

from tfp import config as C


def _complete(root: Path) -> Path:
    (root / "06_conformal" / "TRNSCRPT").mkdir(parents=True)
    (root / "06_conformal" / "TRNSCRPT" / "coverage.csv").write_text("a\n1\n")
    return root


def test_explicit_root_wins(tmp_path, monkeypatch):
    monkeypatch.setenv("TFP_RESULTS", str(tmp_path / "env"))
    assert C.results_root(tmp_path / "given") == tmp_path / "given"


def test_env_root_beats_results_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("TFP_RESULTS", str(tmp_path / "env"))
    monkeypatch.setattr(C, "RESULTS_DIR", _complete(tmp_path / "results"))
    assert C.results_root() == tmp_path / "env"


def test_complete_results_dir_beats_frozen(tmp_path, monkeypatch):
    monkeypatch.delenv("TFP_RESULTS", raising=False)
    monkeypatch.setattr(C, "RESULTS_DIR", _complete(tmp_path / "results"))
    monkeypatch.setattr(C, "FROZEN_DIR", tmp_path / "results_frozen")
    assert C.results_root() == tmp_path / "results"


def test_incomplete_results_dir_falls_back_to_frozen(tmp_path, monkeypatch):
    monkeypatch.delenv("TFP_RESULTS", raising=False)
    (tmp_path / "results").mkdir()
    monkeypatch.setattr(C, "RESULTS_DIR", tmp_path / "results")
    monkeypatch.setattr(C, "FROZEN_DIR", tmp_path / "results_frozen")
    assert C.results_root() == tmp_path / "results_frozen"
