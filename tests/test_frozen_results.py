"""results_frozen/: the committed snapshot of every result file the site export reads. Its MANIFEST.json must
match the files byte for byte, and when a complete results/ run is present the snapshot must equal it."""
import hashlib
import json
from pathlib import Path

import pytest

from tfp import config as C

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "results_frozen"


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _manifest():
    assert (FROZEN / "MANIFEST.json").exists(), "results_frozen/MANIFEST.json missing"
    return json.loads((FROZEN / "MANIFEST.json").read_text())


def test_manifest_matches_the_files():
    m = _manifest()
    assert m["n_files"] == len(m["files"]) and m["n_files"] > 100
    on_disk = {str(p.relative_to(FROZEN)) for p in FROZEN.rglob("*") if p.is_file() and p.name != "MANIFEST.json"}
    listed = {f["path"] for f in m["files"]}
    assert on_disk == listed, f"extra: {sorted(on_disk - listed)[:5]}, missing: {sorted(listed - on_disk)[:5]}"
    for f in m["files"]:
        p = FROZEN / f["path"]
        assert p.stat().st_size == f["bytes"], f["path"]
        assert _sha(p) == f["sha256"], f["path"]


def test_manifest_records_its_origin():
    m = _manifest()
    for k in ("snapshot_date", "frozen_from_commit", "produced_by", "commands", "versions", "data"):
        assert m.get(k), k


@pytest.mark.skipif(not (C.RESULTS_DIR / "06_conformal" / "TRNSCRPT" / "coverage.csv").exists(), reason="no complete results/ run")
def test_snapshot_equals_the_live_results():
    m = _manifest()
    for f in m["files"]:
        live = C.RESULTS_DIR / f["path"]
        assert live.exists(), f"{f['path']} absent from results/"
        assert _sha(live) == f["sha256"], f"{f['path']} differs between results/ and results_frozen/ (re-run make freeze-results)"
