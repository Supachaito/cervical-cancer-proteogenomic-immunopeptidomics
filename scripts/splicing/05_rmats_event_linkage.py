#!/usr/bin/env python3
"""
Paper 3 — HPV-associated rMATS splicing analysis
================================================

Purpose
-------
1) Read rMATS differential-splicing outputs for:
   - HeLa vs C33A
   - SiHa vs C33A
2) Filter significant events:
   FDR <= 0.05 and |IncLevelDifference| >= 0.10
3) Identify events shared between both HPV-positive cell lines:
   - exact same event coordinates + same direction
   - same gene + event type + same direction (supporting, less strict)
4) Intersect significant/shared events with the 54 WT-absent STRG-derived
   HLA-I ligand transcript loci.
5) Produce tables for manuscript-oriented interpretation.

IMPORTANT
---------
This is an association analysis, not a causal HPV experiment.

Direction assumes that in each rMATS run:
  group 1 = HeLa or SiHa
  group 2 = C33A

If this is not true, use --group1-is-hpv-positive no.  Event significance and
locus overlap remain valid, but direction will be reported as unverified.

This script DOES NOT modify H:. It only reads rMATS files and writes outputs
to the chosen C: output directory.
"""

from __future__ import annotations

import os

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")

COORD_COLUMNS = {
    "SE": [
        "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ],
    "A3SS": [
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
    ],
    "A5SS": [
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
    ],
    "MXE": [
        "1stExonStart_0base", "1stExonEnd",
        "2ndExonStart_0base", "2ndExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ],
    "RI": [
        "riExonStart_0base", "riExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ],
}


def yes_no(v: str) -> bool:
    s = str(v).strip().lower()
    if s in {"yes", "y", "true", "1"}:
        return True
    if s in {"no", "n", "false", "0"}:
        return False
    raise argparse.ArgumentTypeError("Use yes or no.")


def normalize_chr(v) -> str:
    if pd.isna(v):
        return ""
    s = str(v).strip()
    s = re.sub(r"^chr", "", s, flags=re.I)
    return s.upper()


def normalize_gene(v) -> str:
    if pd.isna(v):
        return ""
    s = str(v).strip().strip('"')
    return s


def resolve_candidate_file(candidate_arg: str) -> Path:
    if candidate_arg.lower() != "auto":
        p = Path(candidate_arg)
        if not p.exists():
            raise FileNotFoundError(f"Candidate table not found: {p}")
        return p

    name = "PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv"
    root = Path(os.environ.get("PAPER3_ROOT", Path.cwd())).expanduser().resolve()
    hits = []
    for d in (root, root / "input", root / "data"):
        p = d / name
        if p.exists():
            hits.append(p)

    if len(hits) == 1:
        return hits[0]
    if len(hits) == 0:
        raise FileNotFoundError(
            f"Could not find {name} under PAPER3_ROOT, input, or data.\n"
            "Place the candidate table there or pass "
            "--candidates with the exact path."
        )
    raise RuntimeError(
        "More than one candidate table was found. Do not guess which is current:\n"
        + "\n".join(str(x) for x in hits)
        + "\nPass --candidates with the exact path."
    )


def find_rmats_file(root: Path, event_type: str, preferred_mode: str) -> tuple[Path, str]:
    """
    Prefer the requested rMATS mode. If absent, fall back to the other mode.
    Never silently choose among duplicates.
    """
    modes = [preferred_mode, "JCEC" if preferred_mode == "JC" else "JC"]
    for mode in modes:
        basename = f"{event_type}.MATS.{mode}.txt"
        hits = list(root.rglob(basename))
        if len(hits) == 1:
            return hits[0], mode
        if len(hits) > 1:
            raise RuntimeError(
                f"Multiple {basename} files found under {root}; refusing to guess:\n"
                + "\n".join(str(x) for x in hits)
            )
    raise FileNotFoundError(
        f"No {event_type}.MATS.{preferred_mode}.txt or fallback file found under {root}"
    )


def event_key(df: pd.DataFrame, event_type: str) -> pd.Series:
    cols = COORD_COLUMNS[event_type]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{event_type}: missing coordinate columns required for cross-run event matching: "
            + ", ".join(missing)
        )

    parts = [
        pd.Series(event_type, index=df.index),
        df["_chr_norm"],
        df["strand"].astype(str),
    ]
    for c in cols:
        parts.append(pd.to_numeric(df[c], errors="coerce").astype("Int64").astype(str))

    out = parts[0].astype(str)
    for p in parts[1:]:
        out = out + "|" + p.astype(str)
    return out


def load_rmats_comparison(
    root: Path,
    comparison: str,
    preferred_mode: str,
    fdr_cutoff: float,
    dpsi_cutoff: float,
    group1_is_hpv_positive: bool,
) -> pd.DataFrame:
    tables = []

    for event_type in EVENT_TYPES:
        path, used_mode = find_rmats_file(root, event_type, preferred_mode)
        df = pd.read_csv(path, sep="\t", low_memory=False)

        required = ["FDR", "IncLevelDifference", "chr", "strand"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(missing)}"
            )

        df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
        df["IncLevelDifference"] = pd.to_numeric(
            df["IncLevelDifference"], errors="coerce"
        )

        df["_chr_norm"] = df["chr"].map(normalize_chr)
        df["_gene_id"] = (
            df["GeneID"].map(normalize_gene)
            if "GeneID" in df.columns
            else ""
        )
        if "geneSymbol" in df.columns:
            df["_gene_symbol"] = df["geneSymbol"].map(normalize_gene)
        else:
            df["_gene_symbol"] = ""

        # Prefer geneSymbol when available; otherwise use GeneID.
        df["_gene_key"] = np.where(
            df["_gene_symbol"].astype(str).str.len() > 0,
            df["_gene_symbol"],
            df["_gene_id"],
        )

        df["EventType"] = event_type
        df["Comparison"] = comparison
        df["rMATS_Mode"] = used_mode
        df["SourceFile"] = str(path)
        df["EventKey"] = event_key(df, event_type)

        # Genomic span for locus-level overlap.
        coords = COORD_COLUMNS[event_type]
        numeric_coords = df[coords].apply(pd.to_numeric, errors="coerce")
        df["EventStart"] = numeric_coords.min(axis=1)
        df["EventEnd"] = numeric_coords.max(axis=1)

        # rMATS IncLevelDifference = IncLevel1 - IncLevel2.
        # Direction is biologically labelled only if group order is verified.
        if group1_is_hpv_positive:
            df["Direction"] = np.where(
                df["IncLevelDifference"] > 0,
                "Higher inclusion in HPV+",
                np.where(
                    df["IncLevelDifference"] < 0,
                    "Lower inclusion in HPV+",
                    "No direction",
                ),
            )
            df["DirectionSign"] = np.sign(df["IncLevelDifference"]).astype(int)
        else:
            df["Direction"] = "Unverified group order"
            df["DirectionSign"] = 0

        df["Significant"] = (
            df["FDR"].le(fdr_cutoff)
            & df["IncLevelDifference"].abs().ge(dpsi_cutoff)
        )

        tables.append(df)

    return pd.concat(tables, ignore_index=True, sort=False)


def build_shared_exact(h: pd.DataFrame, s: pd.DataFrame, direction_verified: bool) -> pd.DataFrame:
    hsig = h[h["Significant"]].copy()
    ssig = s[s["Significant"]].copy()

    keep = [
        "EventKey", "EventType", "_chr_norm", "strand",
        "_gene_key", "_gene_symbol", "_gene_id",
        "EventStart", "EventEnd",
        "FDR", "IncLevelDifference", "DirectionSign", "Direction",
    ]

    hm = hsig[keep].rename(columns={
        "FDR": "FDR_HeLa_vs_C33A",
        "IncLevelDifference": "dPSI_HeLa_vs_C33A",
        "DirectionSign": "DirectionSign_HeLa",
        "Direction": "Direction_HeLa",
        "_gene_key": "Gene",
        "_gene_symbol": "GeneSymbol",
        "_gene_id": "GeneID",
    })

    sm = ssig[keep].rename(columns={
        "FDR": "FDR_SiHa_vs_C33A",
        "IncLevelDifference": "dPSI_SiHa_vs_C33A",
        "DirectionSign": "DirectionSign_SiHa",
        "Direction": "Direction_SiHa",
        "_gene_key": "Gene",
        "_gene_symbol": "GeneSymbol",
        "_gene_id": "GeneID",
    })

    # EventKey is coordinate-defined. Retain one row per event per comparison,
    # choosing the lowest FDR if a duplicate somehow occurs.
    hm = hm.sort_values("FDR_HeLa_vs_C33A").drop_duplicates("EventKey")
    sm = sm.sort_values("FDR_SiHa_vs_C33A").drop_duplicates("EventKey")

    shared = hm.merge(
        sm[[
            "EventKey", "FDR_SiHa_vs_C33A", "dPSI_SiHa_vs_C33A",
            "DirectionSign_SiHa", "Direction_SiHa"
        ]],
        on="EventKey",
        how="inner",
    )

    if direction_verified:
        shared["SameDirection"] = (
            shared["DirectionSign_HeLa"].eq(shared["DirectionSign_SiHa"])
            & shared["DirectionSign_HeLa"].ne(0)
        )
    else:
        shared["SameDirection"] = np.nan

    return shared


def build_shared_gene_eventtype(
    h: pd.DataFrame, s: pd.DataFrame, direction_verified: bool
) -> pd.DataFrame:
    hsig = h[h["Significant"]].copy()
    ssig = s[s["Significant"]].copy()

    if direction_verified:
        keys = ["_gene_key", "EventType", "DirectionSign"]
    else:
        keys = ["_gene_key", "EventType"]

    hsum = (
        hsig.groupby(keys, dropna=False)
        .agg(
            n_events_HeLa=("EventKey", "nunique"),
            minFDR_HeLa=("FDR", "min"),
            maxAbsDPSI_HeLa=("IncLevelDifference", lambda x: x.abs().max()),
        )
        .reset_index()
    )
    ssum = (
        ssig.groupby(keys, dropna=False)
        .agg(
            n_events_SiHa=("EventKey", "nunique"),
            minFDR_SiHa=("FDR", "min"),
            maxAbsDPSI_SiHa=("IncLevelDifference", lambda x: x.abs().max()),
        )
        .reset_index()
    )

    out = hsum.merge(ssum, on=keys, how="inner")
    out = out.rename(columns={"_gene_key": "Gene"})
    return out


def interval_overlap(a_start, a_end, b_start, b_end) -> bool:
    if any(pd.isna(x) for x in (a_start, a_end, b_start, b_end)):
        return False
    return float(a_start) <= float(b_end) and float(b_start) <= float(a_end)


def summarize_overlaps(
    candidates: pd.DataFrame,
    hela_events: pd.DataFrame,
    siha_events: pd.DataFrame,
    shared_exact: pd.DataFrame,
    shared_gene_evt: pd.DataFrame,
    direction_verified: bool,
) -> pd.DataFrame:
    out_rows = []

    shared_exact_same = shared_exact.copy()
    if direction_verified:
        shared_exact_same = shared_exact_same[shared_exact_same["SameDirection"].eq(True)]

    # Gene/event-type shared key set.
    if direction_verified and "DirectionSign" in shared_gene_evt.columns:
        shared_gene_keys = set(
            zip(
                shared_gene_evt["Gene"].astype(str),
                shared_gene_evt["EventType"].astype(str),
                shared_gene_evt["DirectionSign"].astype(int),
            )
        )
    else:
        shared_gene_keys = set(
            zip(
                shared_gene_evt["Gene"].astype(str),
                shared_gene_evt["EventType"].astype(str),
            )
        )

    for _, c in candidates.iterrows():
        row = c.to_dict()
        cell = str(c.get("CellLine", ""))
        if cell.lower() == "c33a":
            cell = "C33A"
        elif cell.lower() == "hela":
            cell = "HeLa"
        elif cell.lower() == "siha":
            cell = "SiHa"

        chrom = normalize_chr(c.get("STRG_chr", ""))
        st = pd.to_numeric(c.get("STRG_start", np.nan), errors="coerce")
        en = pd.to_numeric(c.get("STRG_end", np.nan), errors="coerce")
        strand = str(c.get("STRG_strand", ""))

        row["CellLine"] = cell

        if cell == "C33A":
            row.update({
                "n_significant_splicing_events_in_locus": np.nan,
                "Significant_EventTypes": "",
                "Significant_Genes": "",
                "min_FDR_significant_event": np.nan,
                "max_abs_dPSI_significant_event": np.nan,
                "n_shared_exact_same_direction_events_in_locus": np.nan,
                "SharedExact_EventTypes": "",
                "SharedExact_Genes": "",
                "n_shared_gene_eventtype_same_direction_support": np.nan,
                "HPV_associated_splicing_support_class":
                    "Reference comparator (C33A); not classified as HPV+",
            })
            out_rows.append(row)
            continue

        events = hela_events if cell == "HeLa" else siha_events
        sig = events[events["Significant"]].copy()

        # Same chromosome and strand, then interval overlap.
        loc = sig[
            sig["_chr_norm"].eq(chrom)
            & sig["strand"].astype(str).eq(strand)
        ].copy()
        if not loc.empty:
            loc = loc[
                loc.apply(
                    lambda r: interval_overlap(
                        st, en, r["EventStart"], r["EventEnd"]
                    ),
                    axis=1,
                )
            ]

        # Exact shared events in the locus.
        shared_loc = shared_exact_same[
            shared_exact_same["_chr_norm"].eq(chrom)
            & shared_exact_same["strand"].astype(str).eq(strand)
        ].copy()
        if not shared_loc.empty:
            shared_loc = shared_loc[
                shared_loc.apply(
                    lambda r: interval_overlap(
                        st, en, r["EventStart"], r["EventEnd"]
                    ),
                    axis=1,
                )
            ]

        # Less strict: event in this candidate locus belongs to a gene/event-type
        # pattern also significant in the other HPV+ comparison.
        shared_gene_support = []
        for _, e in loc.iterrows():
            if direction_verified:
                k = (str(e["_gene_key"]), str(e["EventType"]), int(e["DirectionSign"]))
            else:
                k = (str(e["_gene_key"]), str(e["EventType"]))
            if k in shared_gene_keys:
                shared_gene_support.append(e)

        if shared_gene_support:
            shared_gene_df = pd.DataFrame(shared_gene_support)
        else:
            shared_gene_df = pd.DataFrame()

        row["n_significant_splicing_events_in_locus"] = len(loc)
        row["Significant_EventTypes"] = ";".join(
            sorted(set(loc["EventType"].astype(str)))
        ) if len(loc) else ""
        row["Significant_Genes"] = ";".join(
            sorted(set(x for x in loc["_gene_key"].astype(str) if x))
        ) if len(loc) else ""
        row["min_FDR_significant_event"] = loc["FDR"].min() if len(loc) else np.nan
        row["max_abs_dPSI_significant_event"] = (
            loc["IncLevelDifference"].abs().max() if len(loc) else np.nan
        )

        row["n_shared_exact_same_direction_events_in_locus"] = len(shared_loc)
        row["SharedExact_EventTypes"] = ";".join(
            sorted(set(shared_loc["EventType"].astype(str)))
        ) if len(shared_loc) else ""
        row["SharedExact_Genes"] = ";".join(
            sorted(set(x for x in shared_loc["Gene"].astype(str) if x))
        ) if len(shared_loc) else ""

        row["n_shared_gene_eventtype_same_direction_support"] = len(shared_gene_df)

        if direction_verified and len(shared_loc) > 0:
            cls = "Shared exact rMATS event, same direction in HeLa and SiHa"
        elif len(shared_gene_df) > 0:
            if direction_verified:
                cls = "Shared gene+event-type, same direction; exact event differs"
            else:
                cls = "Shared gene+event-type; direction unverified"
        elif len(loc) > 0:
            cls = "Cell-line-specific significant rMATS event in STRG locus"
        else:
            cls = "No significant rMATS event in STRG locus"

        row["HPV_associated_splicing_support_class"] = cls
        out_rows.append(row)

    return pd.DataFrame(out_rows)


def main():
    ap = argparse.ArgumentParser(
        description="Paper 3: shared HPV-associated rMATS splicing analysis."
    )
    ap.add_argument(
        "--candidates",
        default="auto",
        help=(
            "Path to PAPER3_54_STRG_HLA_rMATS_locus_audit.tsv, "
            "or 'auto' to look under PAPER3_ROOT/input/data."
        ),
    )
    ap.add_argument(
        "--hela-root",
        default=str(Path(os.environ.get("PAPER3_WES_ROOT", "data")) / "04_Splicing_Analysis" / "rMATS_HELA_vs_C33A"),
    )
    ap.add_argument(
        "--siha-root",
        default=str(Path(os.environ.get("PAPER3_WES_ROOT", "data")) / "04_Splicing_Analysis" / "rMATS_SIHA_vs_C33A"),
    )
    ap.add_argument(
        "--outdir",
        default=str(Path(os.environ.get("PAPER3_OUTPUT_ROOT", "results")) / "PAPER3_RMATS_HPV_SHARED_ANALYSIS"),
    )
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--dpsi", type=float, default=0.10)
    ap.add_argument(
        "--mode",
        choices=["JC", "JCEC"],
        default="JC",
        help="Prefer JC (default, conservative) or JCEC. Falls back if preferred mode is absent.",
    )
    ap.add_argument(
        "--group1-is-hpv-positive",
        type=yes_no,
        default=True,
        help=(
            "yes if rMATS group1 is HeLa/SiHa and group2 is C33A. "
            "Default yes. Use no if group order is not verified."
        ),
    )
    args = ap.parse_args()

    candidates_path = resolve_candidate_file(args.candidates)
    hela_root = Path(args.hela_root)
    siha_root = Path(args.siha_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not hela_root.exists():
        raise FileNotFoundError(f"HeLa rMATS root not found: {hela_root}")
    if not siha_root.exists():
        raise FileNotFoundError(f"SiHa rMATS root not found: {siha_root}")

    candidates = pd.read_csv(candidates_path, sep="\t", low_memory=False)
    required_candidates = [
        "CellLine", "Peptide", "STRG_transcript_id",
        "STRG_chr", "STRG_start", "STRG_end", "STRG_strand"
    ]
    miss = [c for c in required_candidates if c not in candidates.columns]
    if miss:
        raise ValueError(
            "Candidate table is missing required columns: " + ", ".join(miss)
        )

    print("=" * 72)
    print("PAPER 3 — HPV-ASSOCIATED rMATS SPLICING ANALYSIS")
    print("=" * 72)
    print("Candidates:", candidates_path)
    print("HeLa rMATS:", hela_root)
    print("SiHa rMATS:", siha_root)
    print("Output:", outdir)
    print(f"Thresholds: FDR <= {args.fdr}; |dPSI| >= {args.dpsi}")
    print("Preferred rMATS mode:", args.mode)
    print("Direction verified:", args.group1_is_hpv_positive)
    print()

    hela = load_rmats_comparison(
        hela_root, "HeLa_vs_C33A", args.mode,
        args.fdr, args.dpsi, args.group1_is_hpv_positive
    )
    siha = load_rmats_comparison(
        siha_root, "SiHa_vs_C33A", args.mode,
        args.fdr, args.dpsi, args.group1_is_hpv_positive
    )

    hela_sig = hela[hela["Significant"]].copy()
    siha_sig = siha[siha["Significant"]].copy()

    shared_exact = build_shared_exact(
        hela, siha, args.group1_is_hpv_positive
    )
    if args.group1_is_hpv_positive:
        shared_exact_same = shared_exact[shared_exact["SameDirection"].eq(True)].copy()
    else:
        shared_exact_same = shared_exact.copy()

    shared_gene_evt = build_shared_gene_eventtype(
        hela, siha, args.group1_is_hpv_positive
    )

    hla = summarize_overlaps(
        candidates,
        hela,
        siha,
        shared_exact,
        shared_gene_evt,
        args.group1_is_hpv_positive,
    )

    # Save core outputs.
    hela_sig.to_csv(
        outdir / "01_significant_events_HeLa_vs_C33A.tsv",
        sep="\t", index=False
    )
    siha_sig.to_csv(
        outdir / "02_significant_events_SiHa_vs_C33A.tsv",
        sep="\t", index=False
    )
    shared_exact.to_csv(
        outdir / "03_shared_exact_significant_events.tsv",
        sep="\t", index=False
    )
    shared_exact_same.to_csv(
        outdir / "04_shared_exact_same_direction_events.tsv",
        sep="\t", index=False
    )
    shared_gene_evt.to_csv(
        outdir / "05_shared_gene_eventtype_patterns.tsv",
        sep="\t", index=False
    )
    hla.to_csv(
        outdir / "06_HLA54_rMATS_significant_splicing_classification.tsv",
        sep="\t", index=False
    )

    # Candidate-focused subset.
    hpv_hla = hla[hla["CellLine"].isin(["HeLa", "SiHa"])].copy()
    supported = hpv_hla[
        hpv_hla["HPV_associated_splicing_support_class"].ne(
            "No significant rMATS event in STRG locus"
        )
    ].copy()
    supported.to_csv(
        outdir / "07_HPVpos_HLA_candidates_with_significant_splicing_support.tsv",
        sep="\t", index=False
    )

    # Summary tables.
    event_summary = pd.concat([
        hela_sig.groupby("EventType").size().rename("HeLa_vs_C33A"),
        siha_sig.groupby("EventType").size().rename("SiHa_vs_C33A"),
    ], axis=1).fillna(0).astype(int).reset_index()

    class_summary = (
        hla.groupby(["CellLine", "HPV_associated_splicing_support_class"])
        .size()
        .rename("n_HLA_ligands")
        .reset_index()
    )

    summary_rows = [
        ["HLA candidates total", len(hla)],
        ["HPV+ HLA candidates total", len(hpv_hla)],
        ["Significant rMATS events HeLa vs C33A", len(hela_sig)],
        ["Significant rMATS events SiHa vs C33A", len(siha_sig)],
        ["Shared exact significant events", len(shared_exact)],
        [
            "Shared exact significant events, same direction",
            len(shared_exact_same) if args.group1_is_hpv_positive else np.nan,
        ],
        [
            "HPV+ HLA candidates with any significant splicing locus support",
            int((hpv_hla["n_significant_splicing_events_in_locus"] > 0).sum()),
        ],
        [
            "HPV+ HLA candidates with shared exact same-direction event support",
            int((hpv_hla["n_shared_exact_same_direction_events_in_locus"] > 0).sum())
            if args.group1_is_hpv_positive else np.nan,
        ],
        [
            "HPV+ HLA candidates with shared gene+event-type support",
            int((hpv_hla["n_shared_gene_eventtype_same_direction_support"] > 0).sum()),
        ],
    ]
    pd.DataFrame(summary_rows, columns=["Metric", "Value"]).to_csv(
        outdir / "08_summary.tsv", sep="\t", index=False
    )
    event_summary.to_csv(
        outdir / "09_significant_event_type_counts.tsv", sep="\t", index=False
    )
    class_summary.to_csv(
        outdir / "10_HLA_support_class_counts.tsv", sep="\t", index=False
    )

    interpretation = f"""PAPER 3 — rMATS HPV-associated splicing analysis

Thresholds
----------
FDR <= {args.fdr}
|IncLevelDifference| >= {args.dpsi}
Preferred rMATS quantification = {args.mode}

Biological interpretation
-------------------------
This analysis asks whether WT-absent STRG-derived HLA-I ligands in HPV-positive
cervical cancer models occur in transcript loci that also show significant
alternative-splicing changes relative to C33A.

It does NOT prove that HPV caused the splice event, because HeLa, SiHa and C33A
are non-isogenic cell lines.

It also does NOT prove that an individual HLA peptide spans the splice junction.
The HLA-to-rMATS connection here is locus-level unless exact peptide-to-junction
sequence reconstruction is performed separately.

Shared evidence hierarchy
-------------------------
1. Shared exact rMATS event, same direction in HeLa and SiHa
   = strongest HPV-associated splicing pattern in this analysis.
2. Shared gene+event-type, same direction; exact event differs
   = supportive but less specific.
3. Cell-line-specific significant rMATS event in STRG locus
   = cell-line-associated splicing support.
4. No significant rMATS event in STRG locus
   = no differential-splicing support under the chosen thresholds.

Direction note
--------------
group1_is_hpv_positive = {args.group1_is_hpv_positive}

When True, IncLevelDifference is interpreted as:
  group1 (HeLa or SiHa) minus group2 (C33A).

Before manuscript use, verify the original rMATS group ordering from the run
command/sample-list files. If it was reversed, rerun with:
  --group1-is-hpv-positive no
and do not interpret direction until group order is reconstructed.

Do NOT use raw event counts as proof that HPV-positive cells have globally more
splicing. A proper global enrichment claim requires an appropriate background
universe and matched transcript detectability/search space.

Recommended manuscript wording if shared support is observed
------------------------------------------------------------
"WT-absent STRG-derived HLA-I ligands were preferentially prioritized when their
source transcript loci also exhibited significant alternative-splicing changes
in HPV-positive cervical cancer models relative to C33A, including a subset
associated with concordant rMATS events across HeLa and SiHa."

Use "HPV-associated", not "HPV-induced".
"""
    (outdir / "README_INTERPRETATION.txt").write_text(
        interpretation, encoding="utf-8"
    )

    print("SIGNIFICANT EVENT COUNTS")
    print(event_summary.to_string(index=False))
    print()
    print("HLA SUPPORT CLASS COUNTS")
    print(class_summary.to_string(index=False))
    print()
    print("SUMMARY")
    for metric, value in summary_rows:
        print(f"- {metric}: {value}")
    print()
    print("DONE")
    print("Output:", outdir)


if __name__ == "__main__":
    main()
