#!/usr/bin/env python
"""Phase 33 — freeze the result files the site export reads into results_frozen/ (committed), with a manifest.

The snapshot holds exactly: every file listed as a source in site/data/provenance.json, every file of the
--save-scores regeneration runs under results/31_site_regen/ (per-sample class probabilities, calibration scores,
class lists, recalibration thresholds), the published tables those runs are compared against by
tests/test_regen_scores.py, and the two identifiability side files the export reads without provenance.
No PNGs and no logs. With the snapshot, `make site` and every test run on a machine without data/ or results/.

Usage: python scripts/33_freeze_results.py [--source results] [--dest results_frozen] [--verify]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

from tfp import config as C

ROOT = C.ROOT
EXTRA_READS = [  # read by scripts/30_export_site_data.py without a provenance entry
    "16_identifiability/layers.json", "16_identifiability/NOTES.md",
]
REGEN_RUNS = ["06_conformal/TRNSCRPT", "08_shift_k20", "08_shift_k50", "12_bodymap", "13_gtex"]
PUBLISHED_COUNTERPART = {"06_conformal/TRNSCRPT": "06_conformal/TRNSCRPT", "12_bodymap": "12_bodymap", "13_gtex": "13_gtex"}
COMMANDS = {  # the regeneration and identifiability runs the snapshot came from (results/31_site_regen/logs/*.log)
    "06": "python scripts/06_conformal_certify.py --assay TRNSCRPT --alpha 0.1 --delta 0.1 --one-per-animal --cal-animals 22 --n-repeats 20 --conditional --grid 10,15,20,30,50,100 --out results/31_site_regen/06_conformal/TRNSCRPT --save-scores",
    "08_k20": "python scripts/08_shift_tests.py --assay TRNSCRPT --k 20 --recalibrate-target 3,5 --out results/31_site_regen/08_shift_k20 --save-scores",
    "08_k50": "python scripts/08_shift_tests.py --assay TRNSCRPT --k 50 --recalibrate-target 3,5 --out results/31_site_regen/08_shift_k50 --save-scores",
    "12": "python scripts/12_bodymap_validate.py --out results/31_site_regen/12_bodymap --save-scores",
    "13": "python scripts/13_gtex_transfer.py --out results/31_site_regen/13_gtex --save-scores",
    "16": "MOTRPAC_NO_REPORT=1 python scripts/16_identifiability.py --bridge",
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def wanted_files(source: Path) -> list[str]:
    prov = json.loads((ROOT / "site" / "data" / "provenance.json").read_text())
    files = {s.split("/", 1)[1] for s in prov["_meta"]["sources"] if s.startswith("results/")}
    for e in prov["entries"]:
        if e.get("file"):
            files.add(e["file"].split("/", 1)[1])
        for f in e.get("files") or []:
            files.add(f.split("/", 1)[1])
    for t in prov.get("tables", []):
        files.add(t["file"].split("/", 1)[1])
    files |= set(EXTRA_READS)
    for run in REGEN_RUNS:
        d = source / "31_site_regen" / run
        for f in sorted(d.rglob("*")):
            if f.is_file() and f.suffix in (".csv", ".json", ".md", ".txt"):
                relf = str(f.relative_to(source))
                files.add(relf)
                if run in PUBLISHED_COUNTERPART and f.suffix == ".csv":
                    pub = source / PUBLISHED_COUNTERPART[run] / f.name
                    if pub.exists():
                        files.add(str(pub.relative_to(source)))
    return sorted(files)


def freeze(source: Path, dest: Path) -> dict:
    files = wanted_files(source)
    missing = [f for f in files if not (source / f).exists()]
    if missing:
        sys.exit(f"missing under {source}: {missing[:5]} ...")
    if dest.exists():
        shutil.rmtree(dest)
    entries = []
    for relf in files:
        src, dst = source / relf, dest / relf
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        entries.append({"path": relf, "bytes": src.stat().st_size, "sha256": sha256(src)})
    try:
        commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        commit = "n/a"
    import numpy, pandas, scipy, sklearn, platform
    m = {"snapshot_date": date.today().isoformat(), "frozen_from_commit": commit,
         "what": "every result file scripts/30_export_site_data.py reads, the --save-scores regeneration runs and the published tables tests/test_regen_scores.py compares them against; derived statistics only, no raw data",
         "produced_by": {"phases_02_14": "make all bodymap gtex (2026-09-18); conformal-quantile correction reruns of phases 06, 08, 12-14 (2026-09-25)",
                         "identifiability": "make identifiability with --bridge (2026-09-26)", "regen_scores": "make regen-scores (2026-09-26/27)"},
         "commands": COMMANDS,
         "versions": {"python": platform.python_version(), "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__, "scikit-learn": sklearn.__version__},
         "data": {"motrpac": "MotrpacRatTraining6moData 2.0.0 (portal release c1.0, rn6; GitHub commit f831a4f), exported 2026-09-17",
                  "bodymap": "rat BodyMap GSE53960 via bodymapRat 1.28.0 (ExperimentHub), 2026-09-17",
                  "gtex": "GTEx v8 open-access files of 2017-06-05 (sha256 in the README), downloaded 2026-09-17/18"},
         "n_files": len(entries), "bytes": int(sum(e["bytes"] for e in entries)), "files": entries}
    (dest / "MANIFEST.json").write_text(json.dumps(m, indent=1))
    print(f"froze {len(entries)} files, {m['bytes'] / 1e6:.2f} MB → {dest} (commit {commit})")
    return m


def verify(dest: Path) -> int:
    m = json.loads((dest / "MANIFEST.json").read_text())
    listed = {e["path"] for e in m["files"]}
    on_disk = {str(p.relative_to(dest)) for p in dest.rglob("*") if p.is_file() and p.name != "MANIFEST.json"}
    bad = []
    for e in m["files"]:
        p = dest / e["path"]
        if not p.exists():
            bad.append(("missing", e["path"]))
        elif p.stat().st_size != e["bytes"] or sha256(p) != e["sha256"]:
            bad.append(("changed", e["path"]))
    bad += [("extra", f) for f in sorted(on_disk - listed)]
    for kind, f in bad:
        print(f"  {kind}: {f}")
    print(f"{'ok' if not bad else 'FAILED'}: {len(listed)} files listed, {len(bad)} problems")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description="freeze or verify the results snapshot")
    ap.add_argument("--source", default=str(C.RESULTS_DIR))
    ap.add_argument("--dest", default=str(C.FROZEN_DIR))
    ap.add_argument("--verify", action="store_true", help="check the snapshot against its manifest and exit")
    a = ap.parse_args()
    if a.verify:
        sys.exit(verify(Path(a.dest)))
    freeze(Path(a.source), Path(a.dest))


if __name__ == "__main__":
    main()
