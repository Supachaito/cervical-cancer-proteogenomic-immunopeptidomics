#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAPER 3 — SCREEN ALL 54 WT-ABSENT STRG HLA-I LIGANDS FOR
STRINGTIE TRANSCRIPT-STRUCTURE NOVELTY

Why this step?
--------------
The current flagship AKSDVTNQQY comes from STRG.7884.3, and the original
HeLa StringTie GTF contains:
    reference_id "ENST00000497707"
    ref_gene_name "COMMD6"

StringTie documents that transcripts lacking a reference_id can be treated
as novel transcript structures relative to the guide annotation supplied by -G.

Therefore, instead of chasing one candidate, screen all 54 WT-absent
STRG-derived HLA-I ligands systematically.

READS
-----
PAPER3_ROOT/PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv and cell-line-specific
assembled GTF files. PAPER3_WES_ROOT can point to an external analysis tree
containing recovered GTF inputs.

WRITES
------
PAPER3_OUTPUT_ROOT/STRG54_TRANSCRIPT_NOVELTY_SCREEN (or under PAPER3_ROOT).
"""

from __future__ import annotations

import os

import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("PAPER3_ROOT", Path.cwd())).expanduser().resolve()
EXTERNAL_ROOT = Path(os.environ.get("PAPER3_WES_ROOT", ROOT / "external" / "WES_Project_Analysis")).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get("PAPER3_OUTPUT_ROOT", ROOT / "results")).expanduser().resolve()
OUTDIR = OUTPUT_ROOT / "STRG54_TRANSCRIPT_NOVELTY_SCREEN"
OUTDIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_SEARCH = [
    ROOT / "PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv",
    ROOT / "INPUT" / "PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv",
]

GTF_CANDIDATES = {
    "C33A": [
        ROOT / "TRANSCRIPT_STRUCTURE" / "C33A" / "C33A_assembled.gtf",
        EXTERNAL_ROOT / "05_Paper3_Rerun" / "Figure2_event_flow" / "01_input" / "gtf" / "C33A_assembled.gtf",
    ],
    "HELA": [
        ROOT / "TRANSCRIPT_STRUCTURE" / "HeLa" / "HELA_assembled.gtf",
        EXTERNAL_ROOT / "05_Paper3_Rerun" / "Figure2_event_flow" / "01_input" / "gtf" / "HELA_assembled.gtf",
    ],
    "SIHA": [
        ROOT / "TRANSCRIPT_STRUCTURE" / "SiHa" / "SIHA_assembled.gtf",
        EXTERNAL_ROOT / "05_Paper3_Rerun" / "Figure2_event_flow" / "01_input" / "gtf" / "SIHA_assembled.gtf",
    ],
}


def norm_cell(x: str) -> str:
    x = str(x).strip().upper()
    if x == "C33A":
        return "C33A"
    if x == "HELA":
        return "HELA"
    if x == "SIHA":
        return "SIHA"
    return x


def parse_attrs(text: str) -> dict:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', text):
        out[m.group(1)] = m.group(2)
    for item in text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out.setdefault(k.strip(), v.strip())
    return out


def first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def read_gtf_header_command(gtf: Path) -> str:
    try:
        with gtf.open("rt", encoding="utf-8", errors="ignore") as fh:
            for _ in range(20):
                line = fh.readline()
                if not line:
                    break
                if line.startswith("#"):
                    return line.strip()
    except Exception:
        pass
    return ""


def parse_selected_transcripts(gtf: Path, target_ids: set[str]):
    """
    Exact transcript_id parsing: no substring matching.
    Returns dict transcript_id -> transcript metadata.
    """
    meta = {}

    with gtf.open("rt", encoding="utf-8", errors="ignore") as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                continue

            a = parse_attrs(cols[8])
            tid = a.get("transcript_id", "")
            if tid not in target_ids:
                continue

            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue

            rec = meta.setdefault(
                tid,
                {
                    "transcript_id": tid,
                    "seqname": cols[0],
                    "strand": cols[6],
                    "gene_id": a.get("gene_id", ""),
                    "reference_ids": set(),
                    "ref_gene_ids": set(),
                    "ref_gene_names": set(),
                    "transcript_start": None,
                    "transcript_end": None,
                    "exons": [],
                    "line_numbers": [],
                },
            )

            rec["line_numbers"].append(line_no)

            if a.get("reference_id"):
                rec["reference_ids"].add(a["reference_id"])
            if a.get("ref_gene_id"):
                rec["ref_gene_ids"].add(a["ref_gene_id"])
            if a.get("ref_gene_name"):
                rec["ref_gene_names"].add(a["ref_gene_name"])

            if cols[2].lower() == "transcript":
                rec["transcript_start"] = start
                rec["transcript_end"] = end

            if cols[2].lower() == "exon":
                rec["exons"].append((start, end))

    # Finalize.
    for tid, rec in meta.items():
        rec["exons"] = sorted(set(rec["exons"]))
        if rec["transcript_start"] is None and rec["exons"]:
            rec["transcript_start"] = min(x[0] for x in rec["exons"])
            rec["transcript_end"] = max(x[1] for x in rec["exons"])

    return meta


def junctions_from_exons(exons):
    exons = sorted(exons)
    out = []
    for left, right in zip(exons[:-1], exons[1:]):
        # Junction breakpoint representation used in the previous deep mapping:
        # left exon end -- (right exon start - 1)
        out.append((left[1], right[0] - 1))
    return out


def boolish(v):
    if pd.isna(v):
        return False
    return str(v).strip().lower() in {"true", "1", "yes", "y"}


# ------------------------------------------------------------------
# 1) Load the 54-candidate table
# ------------------------------------------------------------------

candidate_file = first_existing(CANDIDATE_SEARCH)
if candidate_file is None:
    raise FileNotFoundError(
        "Could not find PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv in ROOT or ROOT/INPUT."
    )

cand = pd.read_csv(candidate_file, sep="\t", dtype=str)

required = {
    "CellLine",
    "Peptide",
    "Protein_ID",
    "STRG_transcript_id",
    "allele",
    "Binding_Class",
    "rep_count",
}
missing = required - set(cand.columns)
if missing:
    raise RuntimeError(f"Candidate table missing required columns: {sorted(missing)}")

cand["CellLine_norm"] = cand["CellLine"].map(norm_cell)

print("=" * 88)
print("PAPER 3 — 54 STRG HLA TRANSCRIPT-NOVELTY SCREEN")
print("=" * 88)
print("Candidate table:", candidate_file)
print("Candidates:", len(cand))
print()


# ------------------------------------------------------------------
# 2) Load each cell-line-specific assembled GTF exactly once
# ------------------------------------------------------------------

gtf_inventory = []
meta_by_cell_tx = {}

for cell in ("C33A", "HELA", "SIHA"):
    tids = set(
        cand.loc[
            cand["CellLine_norm"].eq(cell),
            "STRG_transcript_id"
        ].dropna().astype(str)
    )

    gtf = first_existing(GTF_CANDIDATES[cell])

    if gtf is None:
        gtf_inventory.append({
            "CellLine": cell,
            "GTF": "",
            "Status": "MISSING",
            "TargetTranscriptCount": len(tids),
            "ResolvedTranscriptCount": 0,
            "StringTieCommandHeader": "",
        })
        continue

    print(f"{cell}: parsing {gtf}")
    meta = parse_selected_transcripts(gtf, tids)

    for tid, rec in meta.items():
        meta_by_cell_tx[(cell, tid)] = rec

    gtf_inventory.append({
        "CellLine": cell,
        "GTF": str(gtf),
        "Status": "OK",
        "TargetTranscriptCount": len(tids),
        "ResolvedTranscriptCount": len(meta),
        "StringTieCommandHeader": read_gtf_header_command(gtf),
    })


pd.DataFrame(gtf_inventory).to_csv(
    OUTDIR / "00_GTF_INVENTORY.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------------
# 3) Annotate all 54 candidates
# ------------------------------------------------------------------

rows = []

for _, r in cand.iterrows():
    cell = r["CellLine_norm"]
    tid = str(r["STRG_transcript_id"])
    m = meta_by_cell_tx.get((cell, tid))

    out = r.to_dict()

    if m is None:
        out.update({
            "GTF_transcript_found": False,
            "StringTie_reference_id": "",
            "StringTie_ref_gene_id": "",
            "StringTie_ref_gene_name": "",
            "Transcript_exon_count": "",
            "Transcript_exons": "",
            "Transcript_junctions": "",
            "Transcript_structure_class": "UNRESOLVED_GTF_TRANSCRIPT",
            "Novel_structure_vs_StringTie_guide": "",
            "Novelty_interpretation": (
                "Transcript not resolved in the cell-line-specific assembled GTF."
            ),
        })
    else:
        refs = sorted(m["reference_ids"])
        ref_gene_ids = sorted(m["ref_gene_ids"])
        ref_gene_names = sorted(m["ref_gene_names"])
        exons = m["exons"]
        juncs = junctions_from_exons(exons)

        has_ref = len(refs) > 0

        if has_ref:
            structure_class = "REFERENCE_LINKED_BY_STRINGTIE"
            novel = False
            interp = (
                "reference_id present in StringTie output; do not classify this "
                "transcript structure as novel relative to the -G guide annotation."
            )
        else:
            structure_class = "NOVEL_STRUCTURE_VS_STRINGTIE_GUIDE"
            novel = True
            interp = (
                "No reference_id in the resolved StringTie transcript; candidate "
                "novel transcript structure relative to the -G guide annotation."
            )

        out.update({
            "GTF_transcript_found": True,
            "StringTie_reference_id": ";".join(refs),
            "StringTie_ref_gene_id": ";".join(ref_gene_ids),
            "StringTie_ref_gene_name": ";".join(ref_gene_names),
            "Transcript_exon_count": len(exons),
            "Transcript_exons": ";".join(f"{a}-{b}" for a, b in exons),
            "Transcript_junctions": ";".join(f"{a}-{b}" for a, b in juncs),
            "Transcript_structure_class": structure_class,
            "Novel_structure_vs_StringTie_guide": novel,
            "Novelty_interpretation": interp,
        })

    rows.append(out)

ann = pd.DataFrame(rows)


# ------------------------------------------------------------------
# 4) Candidate priority among true transcript-structure-novel candidates
# ------------------------------------------------------------------

def priority(row):
    novel = str(row.get("Novel_structure_vs_StringTie_guide", "")).lower() == "true"
    if not novel:
        return "NOT_NOVEL_STRUCTURE"

    linked_novel_junction = boolish(
        row.get("Locus_overlap_with_rMATS_novelJunction", False)
    )

    binding = str(row.get("Binding_Class", ""))
    rep = pd.to_numeric(row.get("rep_count", None), errors="coerce")
    rep = 0 if pd.isna(rep) else int(rep)

    if linked_novel_junction and binding == "Strong" and rep >= 4:
        return "A_NOVEL_STRUCTURE_PLUS_RMATS_NOVELJUNCTION_STRONG"
    if linked_novel_junction:
        return "B_NOVEL_STRUCTURE_PLUS_RMATS_NOVELJUNCTION"
    if binding == "Strong" and rep >= 4:
        return "C_NOVEL_STRUCTURE_STRONG_HLA"
    return "D_NOVEL_STRUCTURE_OTHER"


ann["NovelStructurePriority"] = ann.apply(priority, axis=1)

ann.to_csv(
    OUTDIR / "01_ALL54_TRANSCRIPT_NOVELTY_ANNOTATED.tsv",
    sep="\t",
    index=False,
)

novel = ann[
    ann["Novel_structure_vs_StringTie_guide"].astype(str).str.lower().eq("true")
].copy()

priority_order = {
    "A_NOVEL_STRUCTURE_PLUS_RMATS_NOVELJUNCTION_STRONG": 1,
    "B_NOVEL_STRUCTURE_PLUS_RMATS_NOVELJUNCTION": 2,
    "C_NOVEL_STRUCTURE_STRONG_HLA": 3,
    "D_NOVEL_STRUCTURE_OTHER": 4,
}

novel["_priority"] = novel["NovelStructurePriority"].map(priority_order).fillna(99)
novel["_affinity"] = pd.to_numeric(
    novel.get("mhcflurry_affinity", pd.Series(index=novel.index, dtype=float)),
    errors="coerce",
)

novel = novel.sort_values(
    ["_priority", "_affinity", "CellLine", "Peptide"],
    ascending=[True, True, True, True],
).drop(columns=["_priority", "_affinity"], errors="ignore")

novel.to_csv(
    OUTDIR / "02_NOVEL_STRUCTURE_HLA_CANDIDATES.tsv",
    sep="\t",
    index=False,
)

reference_linked = ann[
    ann["Transcript_structure_class"].eq("REFERENCE_LINKED_BY_STRINGTIE")
].copy()

reference_linked.to_csv(
    OUTDIR / "03_REFERENCE_LINKED_HLA_CANDIDATES.tsv",
    sep="\t",
    index=False,
)

unresolved = ann[
    ann["Transcript_structure_class"].eq("UNRESOLVED_GTF_TRANSCRIPT")
].copy()

unresolved.to_csv(
    OUTDIR / "04_UNRESOLVED_GTF_TRANSCRIPTS.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------------
# 5) Explicit flagship audit
# ------------------------------------------------------------------

flagship = ann[ann["Peptide"].eq("AKSDVTNQQY")].copy()
flagship.to_csv(
    OUTDIR / "05_AKSDVTNQQY_STRINGTIE_REFERENCE_AUDIT.tsv",
    sep="\t",
    index=False,
)


# ------------------------------------------------------------------
# 6) Summary
# ------------------------------------------------------------------

n_total = len(ann)
n_resolved = int(ann["GTF_transcript_found"].astype(str).str.lower().eq("true").sum())
n_novel = len(novel)
n_ref = len(reference_linked)
n_unresolved = len(unresolved)

n_novel_rnj = int(
    novel.get(
        "Locus_overlap_with_rMATS_novelJunction",
        pd.Series(index=novel.index, dtype=str),
    ).map(boolish).sum()
) if len(novel) else 0

n_novel_strong = int(
    novel["Binding_Class"].eq("Strong").sum()
) if len(novel) else 0

n_priority_A = int(
    novel["NovelStructurePriority"].eq(
        "A_NOVEL_STRUCTURE_PLUS_RMATS_NOVELJUNCTION_STRONG"
    ).sum()
) if len(novel) else 0

summary = [
    "PAPER 3 — 54 STRG HLA TRANSCRIPT-NOVELTY SCREEN",
    "",
    f"Total WT-absent STRG HLA-I ligands: {n_total}",
    f"GTF transcripts resolved: {n_resolved}",
    f"Novel transcript structures vs StringTie -G guide: {n_novel}",
    f"Reference-linked transcript structures: {n_ref}",
    f"Unresolved transcripts: {n_unresolved}",
    "",
    f"Novel-structure candidates with rMATS novel-junction locus overlap: {n_novel_rnj}",
    f"Novel-structure candidates with Strong predicted HLA binding: {n_novel_strong}",
    f"Priority-A candidates (novel structure + rMATS novelJunction + Strong): {n_priority_A}",
    "",
    "FLAGSHIP REASSESSMENT:",
]

if flagship.empty:
    summary.append("AKSDVTNQQY not found in the 54-candidate table.")
else:
    x = flagship.iloc[0]
    summary += [
        f"AKSDVTNQQY source transcript: {x['STRG_transcript_id']}",
        f"StringTie reference_id: {x.get('StringTie_reference_id', '')}",
        f"Reference gene: {x.get('StringTie_ref_gene_name', '')}",
        f"Structure class: {x.get('Transcript_structure_class', '')}",
        "",
    ]

summary += [
    "DECISION:",
    "Use 02_NOVEL_STRUCTURE_HLA_CANDIDATES.tsv to choose the next mechanistic candidates.",
    "Do not use a reference-linked transcript as proof of 'novel splicing'.",
    "",
    "The already locked result remains valid:",
    "AKSDVTNQQY is a WT-absent HLA-I ligand directly spanning a significant alternative-splice event.",
]

(OUTDIR / "06_FINAL_NOVELTY_SCREEN_SUMMARY.txt").write_text(
    "\n".join(summary),
    encoding="utf-8",
)

print()
print("=" * 88)
print("DONE")
print("=" * 88)
print("\n".join(summary))
print()
print("Results:")
print(OUTDIR)
