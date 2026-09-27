#!/usr/bin/env python
"""Prepare the rat BodyMap (GEO GSE53960) for external validation.

--source geo: the 320 per-sample files in data/external/bodymap/raw/ (GEO supplementary archive):
  AceView gene symbols with the authors' expression values (not raw counts; one file per biological
  sample, technical runs already merged). Writes bodymap_expr_geo.csv (genes × samples) and
  bodymap_meta_geo.csv. The animal identifier is the replicate index within organ × stage × sex
  (GEO carries no animal id); animal_id = <sex>_<stage>_<replicate>.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

from tfp import config as C

ORGAN = {"Adr": "Adrenal", "Brn": "Brain", "Hrt": "Heart", "Kdn": "Kidney", "Lng": "Lung", "Lvr": "Liver",
         "Msc": "Muscle", "Spl": "Spleen", "Thm": "Thymus", "Tst": "Testes", "Utr": "Uterus"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=str(C.EXTERNAL_DIR / "bodymap" / "raw"))
    ap.add_argument("--out", default=str(C.EXTERNAL_DIR))
    args = ap.parse_args()
    raw = Path(args.raw)
    files = sorted(raw.glob("GSM*_SEQC_*.txt.gz"))
    cols, meta = [], []
    for f in files:
        name = f.name.split("_SEQC_", 1)[1].replace(".txt.gz", "")
        organ, sex, stage, rep = name.split("_")
        with gzip.open(f, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", index_col=0)
        s = df.iloc[:, 0].astype(float)
        s.name = name
        cols.append(s)
        meta.append({"sample": name, "gsm": f.name.split("_")[0], "organ": ORGAN.get(organ, organ), "organ_code": organ,
                     "sex": {"F": "female", "M": "male"}[sex], "stage_weeks": int(stage), "replicate": int(rep),
                     "animal_id": f"{sex}_{int(stage)}_{rep}", "source": "GEO GSE53960 AceView expression"})
    X = pd.concat(cols, axis=1)
    X.index.name = "feature_ID"
    out = Path(args.out)
    X.to_csv(out / "bodymap_expr_geo.csv")
    pd.DataFrame(meta).to_csv(out / "bodymap_meta_geo.csv", index=False)
    print(f"wrote {X.shape[0]} genes × {X.shape[1]} samples; organs {pd.DataFrame(meta)['organ'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
