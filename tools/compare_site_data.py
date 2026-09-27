#!/usr/bin/env python
"""Compare two site/data directories (e.g. the committed one and a fresh export from results_frozen/).
Every JSON must be equal apart from its _meta block; manifest.json is compared on the result files it lists
(paths and sha256) and on the regeneration runs, not on versions, hashes or dates. Exit 1 on any difference.

Usage: python tools/compare_site_data.py site/data /tmp/site_data
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def strip(obj):
    if isinstance(obj, dict):
        return {k: strip(v) for k, v in obj.items() if k != "_meta"}
    if isinstance(obj, list):
        return [strip(v) for v in obj]
    return obj


def diff(a, b, path="", out=None, tol=1e-9):
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: only in {'B' if k not in a else 'A'}")
            else:
                diff(a[k], b[k], f"{path}/{k}", out, tol)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            diff(x, y, f"{path}[{i}]", out, tol)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if not math.isclose(float(a), float(b), abs_tol=tol, rel_tol=0):
            out.append(f"{path}: {a} vs {b}")
    elif a != b:
        out.append(f"{path}: {a!r} vs {b!r}")
    return out


def main(a_dir: Path, b_dir: Path) -> int:
    names = sorted({p.name for p in a_dir.glob("*.json")} | {p.name for p in b_dir.glob("*.json")})
    problems = []
    for name in names:
        pa, pb = a_dir / name, b_dir / name
        if not pa.exists() or not pb.exists():
            problems.append(f"{name}: only in {'B' if not pa.exists() else 'A'}")
            continue
        a, b = json.loads(pa.read_text()), json.loads(pb.read_text())
        if name == "manifest.json":
            a = {"phases": a.get("phases"), "regeneration_runs": a.get("regeneration_runs"), "phases_used": a.get("phases_used")}
            b = {"phases": b.get("phases"), "regeneration_runs": b.get("regeneration_runs"), "phases_used": b.get("phases_used")}
        d = diff(strip(a), strip(b), name)
        problems += d[:20] + ([f"{name}: ... {len(d) - 20} more"] if len(d) > 20 else [])
    for p in problems:
        print("  ", p)
    print("identical apart from _meta" if not problems else f"{len(problems)} differences")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
