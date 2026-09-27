"""The README's numbered text is generated (tools/readme_blocks.py) from site/data, whose values carry provenance; the
Key results table equals scripts/30_export_site_data.py --readme-table; the sections moved to docs/REPRODUCIBILITY.md
are there."""
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_readme_blocks_match_the_generator():
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "readme_blocks.py"), "--check"], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stdout + r.stderr


def test_key_results_table_equals_the_site_export_generator():
    import json
    spec = importlib.util.spec_from_file_location("export30", ROOT / "scripts" / "30_export_site_data.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    table = mod.readme_table(json.loads((ROOT / "site" / "data" / "provenance.json").read_text())["entries"])
    readme = (ROOT / "README.md").read_text()
    assert table in readme
    assert "All-gene simple baselines on the same folds" in table


def test_moved_sections_live_in_reproducibility():
    doc = (ROOT / "docs" / "REPRODUCIBILITY.md").read_text()
    for h in ("## Setup", "## Data: sources, access and provenance", "## Inputs and outputs", "## Methods and provenance",
              "## Validation", "## Example runs", "## Repository map"):
        assert h in doc, h
    assert len(re.findall(r"\b[0-9a-f]{64}\b", doc)) >= 20, "the sha256 blocks moved with the data section"
    readme = (ROOT / "README.md").read_text()
    assert len(readme.splitlines()) <= 190
    for s in ("NOTICE.md", "docs/REPRODUCIBILITY.md", "evaluation must prevent animal-level data leakage"):
        assert s in readme, s
