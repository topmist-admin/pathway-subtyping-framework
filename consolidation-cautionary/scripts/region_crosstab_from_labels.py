#!/usr/bin/env python3
"""
Layer-A exact recompute of the Figure-1 headline (manuscript §3.2).

The schizophrenia primary partition (39/92/10) is a DEPOSITED artifact (framework
v0.3.0, seed 42). Its region and diagnosis cross-tabulations are pure functions of
the deposited labels x sample metadata -- no clustering, no seed, no framework
version. This script recomputes them and must reproduce the manuscript EXACTLY:

  subtype x region (3-level AnCg/DLPFC/nAcc): chi2 = 125.12, p = 4.30e-26,
      Cramer's V = 0.666 (uncorrected) / 0.660 (Bergsma-corrected)
  subtype x diagnosis: chi2 = 1.79, p = 0.408, V = 0.113
  composition: S0 = 2/0/37, S1 = 45/47/0, S2 = 1/1/8
  (sanity) the BUGGY 2-level region coding gives only chi2 = 34.6, V = 0.50 -- if
      you get this, you used the wrong (nAcc->ACC merged) column.

Input: the deposited partition table (sample_metadata_with_subtypes.csv) carrying
per-sample `subtype`, the 3-level `brain region`, and `diagnosis`/`clinical diagnosis`.
No framework dependency; requires only pandas + scipy.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


def cv(ct, bc=False):
    chi2, p, _, _ = chi2_contingency(ct.values, correction=False)
    n = ct.values.sum(); r, k = ct.shape
    if bc:
        phi2 = chi2 / n
        phi2c = max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
        rc = r - (r - 1) ** 2 / (n - 1); kc = k - (k - 1) ** 2 / (n - 1)
        return chi2, p, np.sqrt(phi2c / max(min(kc - 1, rc - 1), 1e-12))
    return chi2, p, np.sqrt((chi2 / n) / (min(r, k) - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", required=True,
                    help="deposited partition CSV (sample_metadata_with_subtypes.csv)")
    args = ap.parse_args()
    df = pd.read_csv(args.labels, index_col=0)
    reg = "brain region" if "brain region" in df.columns else "brain_region"
    dcol = "diagnosis" if "diagnosis" in df.columns else "clinical diagnosis"
    print(f"n = {len(df)} | subtypes {df['subtype'].value_counts().sort_index().to_dict()}")
    print(f"region levels: {df[reg].value_counts().to_dict()}")

    # Warn BEFORE printing any numbers. This check used to sit at the end of main(),
    # which meant a reader saw the (wrong) 2-level chi-square and V first and the
    # caveat last -- exactly the wrong order for a trap this easy to fall into.
    if df[reg].nunique() < 3:
        print(
            "\n*** WARNING: region column has only "
            f"{df[reg].nunique()} levels -- this is the buggy 2-level coding. ***\n"
            "*** The numbers below will NOT match the manuscript. Use the 3-level  ***\n"
            "*** `brain region` re-derived from raw GEO.                           ***\n"
        )

    print("\n=== SUBTYPE x REGION (use the 3-level coding) ===")
    ctr = pd.crosstab(df["subtype"], df[reg]); print(ctr)
    chi2, p, v = cv(ctr); _, _, vc = cv(ctr, True)
    print(f"chi2={chi2:.2f}  p={p:.3e}  V(uncorr)={v:.4f}  V(Bergsma)={vc:.4f}")

    # The 2-level coding, computed rather than asserted. The manuscript quotes a
    # Bergsma-corrected V for this merged coding so that the drop from the 3-level
    # headline is a like-for-like comparison; before 2026-08-06 only the UNCORRECTED
    # 0.50 was recorded anywhere in the deposit, so the corrected figure the paper
    # actually quotes had no source a referee could check. Derived here from the same
    # deposited labels, by the same estimator, in the same run.
    if df[reg].nunique() >= 3:
        merged = df[reg].replace({"nAcc": "AnCg", "nACC": "AnCg", "nucleus accumbens": "AnCg"})
        if merged.nunique() == 2:
            print("\n=== SUBTYPE x REGION, BUGGY 2-LEVEL CODING (nAcc merged into AnCg) ===")
            print("    Reported for contrast only -- this is NOT the manuscript headline.")
            ct2 = pd.crosstab(df["subtype"], merged); print(ct2)
            chi22, p2, v2 = cv(ct2); _, _, vc2 = cv(ct2, True)
            print(f"chi2={chi22:.2f}  p={p2:.3e}  V(uncorr)={v2:.4f}  V(Bergsma)={vc2:.4f}")
            print(f"    manuscript quotes chi2=34.6, V(Bergsma)=0.48, V(uncorr)=0.50 -> "
                  f"match: {abs(chi22 - 34.6) < 0.5 and abs(vc2 - 0.48) < 0.01 and abs(v2 - 0.50) < 0.01}")
        else:
            print(f"\n(2-level contrast skipped: merge produced {merged.nunique()} levels, expected 2)")

    print("\n=== SUBTYPE x DIAGNOSIS ===")
    ctd = pd.crosstab(df["subtype"], df[dcol]); print(ctd)
    chi2d, pdv, vd = cv(ctd)
    print(f"chi2={chi2d:.2f}  p={pdv:.3e}  V={vd:.4f}")

    ok = abs(chi2 - 125.12) < 0.5 and abs(v - 0.666) < 0.01
    print(f"\nMatches manuscript §3.2 headline (125.12 / V=0.666): {ok}")
    if df[reg].nunique() < 3:
        print("WARNING: region column has <3 levels -- this is the buggy 2-level coding; "
              "use the 3-level `brain region` re-derived from raw GEO.")


if __name__ == "__main__":
    main()
