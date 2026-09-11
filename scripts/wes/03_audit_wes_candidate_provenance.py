#!/usr/bin/env python3
"""Audit provenance of the four WES candidates recovered outside the integrated peptide-pair landscape.

This supporting audit reads the final WES candidate table, extracts the four
target peptides, and searches text-like project files for exact peptide
occurrences. Paths are supplied explicitly; no user-specific directories are
hard-coded.
"""
import argparse
import os
import re
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--wes-table", required=True, help="Final WES candidate TSV (Supplementary Table S1-compatible).")
ap.add_argument("--project-root", required=True, help="Project directory to search for provenance text hits.")
ap.add_argument("--out-dir", required=True, help="Directory for audit outputs.")
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
report_path = os.path.join(args.out_dir, "WES4_PROVENANCE_FULL_RECORDS.tsv")
hits_path = os.path.join(args.out_dir, "WES4_PROJECT_TEXT_HITS.tsv")

TARGETS = {
    "RYFDEPVEL",
    "AAAHIHRYL",
    "RAAPHHISL",
    "REPDLVLRL"
}

# ============================================================
# 1. READ FINAL WES TABLE
# ============================================================

wes = pd.read_csv(
    args.wes_table,
    sep="\t",
    low_memory=False
)

def norm(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())

lookup = {norm(c): c for c in wes.columns}

pep_col = None

for candidate in [
    "Peptide sequence",
    "Peptide"
]:
    k = norm(candidate)
    if k in lookup:
        pep_col = lookup[k]
        break

if pep_col is None:
    raise RuntimeError(
        "Peptide column not found in WES table"
    )

wes["_pep_clean"] = (
    wes[pep_col]
    .astype(str)
    .str.upper()
    .str.replace(r"[^A-Z]", "", regex=True)
)

q = wes[
    wes["_pep_clean"].isin(TARGETS)
].copy()

if len(q) != 4:
    raise RuntimeError(
        f"Expected 4 WES records; found {len(q)}"
    )

q.drop(
    columns=["_pep_clean"],
    inplace=True
)

q.to_csv(
    str(report_path),
    sep="\t",
    index=False
)

print("")
print("=" * 78)
print("WES 4 FULL FINAL-TABLE RECORDS")
print("=" * 78)

for _, row in q.iterrows():

    print("")
    print("-" * 78)

    for col in q.columns:

        val = row[col]

        if pd.isna(val):
            continue

        s = str(val).strip()

        if not s:
            continue

        print(f"{col}: {s}")


# ============================================================
# 2. SEARCH TEXT-LIKE PROJECT FILES FOR EXACT PEPTIDES
#
# Purpose:
# determine provenance/location of these four sequences.
#
# Skip generated Figure 3 rebuild folders to avoid finding only
# the audit products we just created.
# ============================================================

allowed_ext = {
    ".tsv",
    ".csv",
    ".txt"
}

skip_fragments = [
    os.path.normcase(
        os.path.join(
            "MANUSCRIPT_ASSEMBLY",
            "01_MAIN_FIGURES",
            "FIG03",
            "PANEL_A_REBUILD"
        )
    )
]

hits = []

files_scanned = 0

for dirpath, dirnames, filenames in os.walk(args.project_root):

    norm_dir = os.path.normcase(dirpath)

    if any(
        frag in norm_dir
        for frag in skip_fragments
    ):
        dirnames[:] = []
        continue

    for fn in filenames:

        ext = os.path.splitext(fn)[1].lower()

        if ext not in allowed_ext:
            continue

        path = os.path.join(
            dirpath,
            fn
        )

        try:
            size = os.path.getsize(path)
        except OSError:
            continue

        # Avoid giant accidental text files.
        if size > 500_000_000:
            continue

        files_scanned += 1

        try:

            with open(
                path,
                "r",
                encoding="utf-8",
                errors="ignore"
            ) as fh:

                for lineno, line in enumerate(
                    fh,
                    start=1
                ):

                    upper = line.upper()

                    found = [
                        pep
                        for pep in TARGETS
                        if pep in upper
                    ]

                    for pep in found:

                        hits.append({
                            "Peptide": pep,
                            "File": path,
                            "Line": lineno,
                            "Text": line.strip()[:1000]
                        })

        except Exception:
            continue


hit_df = pd.DataFrame(
    hits,
    columns=[
        "Peptide",
        "File",
        "Line",
        "Text"
    ]
)

hit_df.to_csv(
    str(hits_path),
    sep="\t",
    index=False
)


print("")
print("=" * 78)
print("PROJECT PROVENANCE SEARCH")
print("=" * 78)

print("")
print(
    "Text files scanned:",
    files_scanned
)

for pep in sorted(TARGETS):

    h = hit_df[
        hit_df["Peptide"] == pep
    ]

    print("")
    print("-" * 78)
    print(pep)
    print("Hits:", len(h))

    if len(h):

        for _, r in h.iterrows():

            print(
                r["File"]
            )

            print(
                "  line",
                int(r["Line"])
            )

            print(
                " ",
                r["Text"][:300]
            )


print("")
print("=" * 78)
print("AUDIT COMPLETE")
print("=" * 78)

print("")
print(
    "Full WES records:",
    str(report_path)
)

print(
    "Project hits:",
    str(hits_path)
)