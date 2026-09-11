#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAPER 3 — FINAL TARGETED NOVELTY CHECK FOR STRG.7884.3
======================================================

Goal
----
1) Reconstruct the exon/intron chain of HeLa STRG.7884.3 from HELA_assembled.gtf.
2) Find plausible human reference GTF files used/available in the Paper 3 project.
3) Compare STRG.7884.3 against reference transcripts at the same locus:
   - exact intron-chain match?
   - shares some splice junctions but has a different chain?
   - no matching reference transcript?
4) Search old rMATS commands/logs/scripts for group1/group2 order clues.

READ ONLY:
  PAPER3_ROOT
  PAPER3_WES_ROOT

WRITE:
  PAPER3_OUTPUT_ROOT/TARGET_STRG7884_3_REFERENCE_COMPARISON

This script does NOT modify H:.
"""

from __future__ import annotations

import os

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

TARGET = "STRG.7884.3"

WORKSPACE = Path(os.environ.get("PAPER3_ROOT", Path.cwd())).expanduser().resolve()
HROOT = Path(os.environ.get("PAPER3_WES_ROOT", WORKSPACE / "external" / "WES_Project_Analysis")).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get("PAPER3_OUTPUT_ROOT", WORKSPACE / "results")).expanduser().resolve()

OUTDIR = OUTPUT_ROOT / "TARGET_STRG7884_3_REFERENCE_COMPARISON"
OUTDIR.mkdir(parents=True, exist_ok=True)

HELA_GTF = (
    WORKSPACE / "TRANSCRIPT_STRUCTURE" / "HeLa" / "HELA_assembled.gtf"
)

TARGETED_ROOTS = [
    HROOT / "04_Splicing_Analysis",
    HROOT / "05_Paper3_Rerun",
    HROOT / "Project_Paper_3_Databases",
]

REFERENCE_NAME_HINTS = (
    "gencode",
    "grch38",
    "grch37",
    "hg38",
    "hg19",
    "annotation",
    "reference",
    "ref_",
    "genes.gtf",
    "homo_sapiens",
)

RUN_TEXT_EXTS = {".py", ".ps1", ".sh", ".bat", ".cmd", ".log", ".txt"}


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("rt", encoding="utf-8", errors="ignore")


def attrs(s: str) -> dict:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', s):
        out[m.group(1)] = m.group(2)
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out.setdefault(k.strip(), v.strip())
    return out


def norm_chr(s: str) -> str:
    s = str(s).strip()
    return re.sub(r"^chr", "", s, flags=re.I).upper()


def write_tsv(path: Path, rows: list[dict], columns: list[str]):
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def intron_chain(exons: list[tuple[int, int]], strand: str):
    if len(exons) < 2:
        return tuple()

    exons = sorted(exons)
    chain = []
    for a, b in zip(exons[:-1], exons[1:]):
        left_end = a[1]
        right_start = b[0]
        chain.append((left_end, right_start - 1))
    return tuple(chain)


def read_target_transcript(gtf: Path):
    if not gtf.exists():
        raise FileNotFoundError(f"HeLa assembled GTF not found: {gtf}")

    rows = []
    exons = []

    with gtf.open("rt", encoding="utf-8", errors="ignore") as fh:
        for n, line in enumerate(fh, start=1):
            if TARGET not in line:
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                continue
            a = attrs(cols[8])
            tid = a.get("transcript_id", "")
            if tid != TARGET:
                continue

            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue

            rec = {
                "LineNumber": n,
                "SeqID": cols[0],
                "Feature": cols[2],
                "Start": start,
                "End": end,
                "Strand": cols[6],
                "GeneID": a.get("gene_id", ""),
                "GeneName": a.get("gene_name", ""),
                "TranscriptID": tid,
                "Attributes": cols[8],
            }
            rows.append(rec)

            if cols[2].lower() == "exon":
                exons.append((start, end))

    if not exons:
        raise RuntimeError(f"No exon rows found for {TARGET} in {gtf}")

    chroms = {norm_chr(r["SeqID"]) for r in rows if r["Feature"].lower() == "exon"}
    strands = {r["Strand"] for r in rows if r["Feature"].lower() == "exon"}

    if len(chroms) != 1 or len(strands) != 1:
        raise RuntimeError("Target transcript has ambiguous chromosome/strand.")

    chrom = next(iter(chroms))
    strand = next(iter(strands))
    exons = sorted(set(exons))

    return {
        "rows": rows,
        "exons": exons,
        "chrom": chrom,
        "strand": strand,
        "start": min(x[0] for x in exons),
        "end": max(x[1] for x in exons),
        "introns": intron_chain(exons, strand),
    }


def looks_like_reference_gtf(path: Path) -> tuple[bool, dict]:
    """
    Inspect a small prefix rather than parsing the whole file.
    Reference GTF usually contains ENSG/ENST or known gene/transcript IDs
    and is not dominated by STRG transcript IDs.
    """
    stats = {
        "lines_checked": 0,
        "feature_lines": 0,
        "strg_lines": 0,
        "ens_lines": 0,
        "refseq_lines": 0,
    }

    try:
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                stats["lines_checked"] += 1
                if "\t" in line:
                    stats["feature_lines"] += 1
                if "STRG." in line:
                    stats["strg_lines"] += 1
                if "ENSG" in line or "ENST" in line:
                    stats["ens_lines"] += 1
                if "NM_" in line or "NR_" in line:
                    stats["refseq_lines"] += 1
                if stats["lines_checked"] >= 20000:
                    break
    except Exception:
        return False, stats

    ref_signal = stats["ens_lines"] + stats["refseq_lines"]
    looks_ref = (
        stats["feature_lines"] > 100
        and ref_signal > 10
        and stats["strg_lines"] < max(10, stats["feature_lines"] * 0.05)
    )
    return looks_ref, stats


def collect_gtf_candidates():
    hits = []

    for root in [WORKSPACE] + TARGETED_ROOTS:
        if not root.exists():
            continue

        for p in root.rglob("*"):
            if not p.is_file():
                continue

            low = p.name.lower()
            if not (
                low.endswith(".gtf")
                or low.endswith(".gtf.gz")
                or low.endswith(".gff")
                or low.endswith(".gff.gz")
                or low.endswith(".gff3")
                or low.endswith(".gff3.gz")
            ):
                continue

            # Skip cell-line assembled/transdecoder files as reference candidates.
            full = str(p).lower()
            if "transdecoder" in full:
                continue
            if "assembled.gtf" in low:
                continue

            hint = any(k in full for k in REFERENCE_NAME_HINTS)
            looks_ref, stats = looks_like_reference_gtf(p)

            if hint or looks_ref:
                hits.append({
                    "Path": str(p),
                    "NameHint": hint,
                    "LooksLikeReference": looks_ref,
                    **stats,
                })

    # Deduplicate exact path strings.
    uniq = {}
    for r in hits:
        uniq.setdefault(r["Path"].lower(), r)

    return list(uniq.values())


def parse_reference_locus(path: Path, target):
    """
    Parse only rows overlapping the target locus on same chromosome/strand.
    Group exons by transcript_id.
    """
    tx = defaultdict(lambda: {
        "exons": [],
        "gene_id": "",
        "gene_name": "",
        "transcript_id": "",
        "transcript_name": "",
    })

    tchr = target["chrom"]
    tstrand = target["strand"]
    tstart = target["start"]
    tend = target["end"]

    try:
        with open_text(path) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) < 9:
                    continue

                if norm_chr(cols[0]) != tchr:
                    continue
                if cols[6] != tstrand:
                    continue

                try:
                    start, end = int(cols[3]), int(cols[4])
                except ValueError:
                    continue

                if end < tstart or start > tend:
                    continue

                a = attrs(cols[8])
                tid = (
                    a.get("transcript_id", "")
                    or a.get("transcript", "")
                    or a.get("Parent", "")
                )

                if not tid:
                    continue

                rec = tx[tid]
                rec["transcript_id"] = tid
                rec["gene_id"] = a.get("gene_id", rec["gene_id"])
                rec["gene_name"] = a.get("gene_name", rec["gene_name"])
                rec["transcript_name"] = a.get(
                    "transcript_name", rec["transcript_name"]
                )

                if cols[2].lower() == "exon":
                    rec["exons"].append((start, end))
    except Exception:
        return []

    rows = []
    target_introns = set(target["introns"])
    target_exons = set(target["exons"])

    for tid, rec in tx.items():
        exons = sorted(set(rec["exons"]))
        if not exons:
            continue

        introns = set(intron_chain(exons, tstrand))

        exact_chain = tuple(sorted(introns)) == tuple(sorted(target_introns))
        shared_introns = len(introns & target_introns)
        union_introns = len(introns | target_introns)
        intron_jaccard = shared_introns / union_introns if union_introns else 1.0

        shared_exons = len(set(exons) & target_exons)
        union_exons = len(set(exons) | target_exons)
        exon_jaccard = shared_exons / union_exons if union_exons else 1.0

        rows.append({
            "ReferenceFile": str(path),
            "GeneID": rec["gene_id"],
            "GeneName": rec["gene_name"],
            "TranscriptID": tid,
            "TranscriptName": rec["transcript_name"],
            "ReferenceExonCount": len(exons),
            "TargetExonCount": len(target["exons"]),
            "ReferenceIntronCount": len(introns),
            "TargetIntronCount": len(target_introns),
            "ExactIntronChainMatch": exact_chain,
            "SharedIntrons": shared_introns,
            "IntronJaccard": round(intron_jaccard, 6),
            "SharedExactExons": shared_exons,
            "ExonJaccard": round(exon_jaccard, 6),
            "ReferenceExons": ";".join(f"{a}-{b}" for a,b in exons),
            "ReferenceIntrons": ";".join(f"{a}-{b}" for a,b in sorted(introns)),
        })

    return rows


def search_rmats_order_clues():
    clues = []

    patterns = (
        "rmats",
        "--b1",
        "--b2",
        "b1.txt",
        "b2.txt",
        "rMATS_HELA_vs_C33A",
        "rMATS_SIHA_vs_C33A",
    )

    for root in TARGETED_ROOTS:
        if not root.exists():
            continue

        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if p.suffix.lower() not in RUN_TEXT_EXTS:
                continue

            # Avoid giant result tables.
            try:
                if p.stat().st_size > 10 * 1024 * 1024:
                    continue
            except OSError:
                continue

            lowname = p.name.lower()
            if any(x in lowname for x in ("mats.", "fromgtf", "read_out", "rawinput")):
                continue

            try:
                text = p.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue

            low = text.lower()
            if not any(x.lower() in low for x in patterns):
                continue
            if not any(x in low for x in ("hela", "siha", "c33a", "c33")):
                continue

            # Keep only relevant lines plus nearby context.
            lines = text.splitlines()
            matched = []
            for i, line in enumerate(lines):
                l = line.lower()
                if (
                    "rmats" in l
                    or "--b1" in l
                    or "--b2" in l
                    or "b1.txt" in l
                    or "b2.txt" in l
                    or "hela" in l
                    or "siha" in l
                    or "c33a" in l
                ):
                    for j in range(max(0, i-2), min(len(lines), i+3)):
                        matched.append(lines[j])

            # De-duplicate while preserving order.
            seen = set()
            context = []
            for x in matched:
                if x not in seen:
                    seen.add(x)
                    context.append(x)

            clues.append({
                "File": str(p),
                "Context": "\n".join(context[:200]),
            })

    return clues


# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------

print("=" * 80)
print("PAPER 3 — FINAL NOVELTY CHECK: STRG.7884.3")
print("=" * 80)

target = read_target_transcript(HELA_GTF)

target_rows = [{
    "TranscriptID": TARGET,
    "Chr": target["chrom"],
    "Strand": target["strand"],
    "Start": target["start"],
    "End": target["end"],
    "ExonCount": len(target["exons"]),
    "Exons": ";".join(f"{a}-{b}" for a,b in target["exons"]),
    "IntronCount": len(target["introns"]),
    "Introns": ";".join(f"{a}-{b}" for a,b in target["introns"]),
}]

write_tsv(
    OUTDIR / "01_TARGET_STRG7884_3_STRUCTURE.tsv",
    target_rows,
    [
        "TranscriptID","Chr","Strand","Start","End",
        "ExonCount","Exons","IntronCount","Introns"
    ],
)

print("Target locus:")
print(target_rows[0])
print()

candidates = collect_gtf_candidates()

write_tsv(
    OUTDIR / "02_REFERENCE_GTF_CANDIDATES.tsv",
    candidates,
    [
        "Path","NameHint","LooksLikeReference","lines_checked",
        "feature_lines","strg_lines","ens_lines","refseq_lines"
    ],
)

print("Reference-like GTF candidates:", len(candidates))

comparisons = []

for i, r in enumerate(candidates, start=1):
    p = Path(r["Path"])
    print(f"[{i}/{len(candidates)}] comparing {p}")
    comparisons.extend(parse_reference_locus(p, target))

# Sort best matches first.
comparisons.sort(
    key=lambda r: (
        not bool(r["ExactIntronChainMatch"]),
        -float(r["IntronJaccard"]),
        -float(r["ExonJaccard"]),
    )
)

write_tsv(
    OUTDIR / "03_REFERENCE_TRANSCRIPT_COMPARISON.tsv",
    comparisons,
    [
        "ReferenceFile","GeneID","GeneName","TranscriptID","TranscriptName",
        "ReferenceExonCount","TargetExonCount",
        "ReferenceIntronCount","TargetIntronCount",
        "ExactIntronChainMatch","SharedIntrons","IntronJaccard",
        "SharedExactExons","ExonJaccard",
        "ReferenceExons","ReferenceIntrons"
    ],
)

best = comparisons[:50]
write_tsv(
    OUTDIR / "04_TOP_REFERENCE_MATCHES.tsv",
    best,
    [
        "ReferenceFile","GeneID","GeneName","TranscriptID","TranscriptName",
        "ReferenceExonCount","TargetExonCount",
        "ReferenceIntronCount","TargetIntronCount",
        "ExactIntronChainMatch","SharedIntrons","IntronJaccard",
        "SharedExactExons","ExonJaccard",
        "ReferenceExons","ReferenceIntrons"
    ],
)

order_clues = search_rmats_order_clues()
write_tsv(
    OUTDIR / "05_RMATS_GROUP_ORDER_CLUES.tsv",
    order_clues,
    ["File","Context"],
)

exact = [r for r in comparisons if r["ExactIntronChainMatch"]]
best_match = comparisons[0] if comparisons else None

summary = []
summary.append("PAPER 3 — FINAL STRG.7884.3 REFERENCE COMPARISON")
summary.append("")
summary.append(f"Target transcript: {TARGET}")
summary.append(f"Target locus: chr{target['chrom']}:{target['start']}-{target['end']} ({target['strand']})")
summary.append(f"Target exon count: {len(target['exons'])}")
summary.append(f"Target intron count: {len(target['introns'])}")
summary.append("")
summary.append(f"Reference-like GTF candidates inspected: {len(candidates)}")
summary.append(f"Reference transcripts overlapping target locus: {len(comparisons)}")
summary.append(f"Exact intron-chain matches: {len(exact)}")
summary.append(f"rMATS group-order clue files: {len(order_clues)}")
summary.append("")

if exact:
    summary.append("NOVELTY DECISION:")
    summary.append("An exact reference intron-chain match was found.")
    summary.append("Do NOT call STRG.7884.3 a novel splice isoform on intron-chain grounds.")
    summary.append("")
    for r in exact[:10]:
        summary.append(
            f"Exact match: {r['GeneName']} | {r['TranscriptID']} | {r['ReferenceFile']}"
        )
elif best_match:
    summary.append("NOVELTY DECISION:")
    summary.append("No exact reference intron-chain match was found among the recovered reference-like annotations.")
    summary.append("This supports a non-reference transcript/splice configuration, subject to confirming that the correct original reference GTF is among the inspected files.")
    summary.append("")
    summary.append(
        f"Best reference overlap: {best_match['GeneName']} | "
        f"{best_match['TranscriptID']} | intron Jaccard={best_match['IntronJaccard']} | "
        f"shared introns={best_match['SharedIntrons']}"
    )
else:
    summary.append("NOVELTY DECISION:")
    summary.append("No usable overlapping reference transcripts were recovered.")
    summary.append("Novelty cannot yet be classified from reference structure.")
    summary.append("")

if order_clues:
    summary.append("rMATS GROUP ORDER:")
    summary.append("Potential command/run-order clues were recovered.")
    summary.append("Inspect 05_RMATS_GROUP_ORDER_CLUES.tsv before final ΔPSI direction wording.")
else:
    summary.append("rMATS GROUP ORDER:")
    summary.append("No explicit group-order clue was recovered.")
    summary.append("Keep 'shared significant event' wording and avoid HPV+ higher/lower direction wording for now.")

summary.append("")
summary.append("CURRENT LOCKED CLAIM REGARDLESS OF NOVELTY RESULT:")
summary.append(
    "AKSDVTNQQY is a WT-absent HLA-I ligand that directly spans a significant alternative-splice event."
)

(OUTDIR / "06_FINAL_REFERENCE_COMPARISON_SUMMARY.txt").write_text(
    "\n".join(summary),
    encoding="utf-8"
)

print()
print("=" * 80)
print("DONE")
print("=" * 80)
print("\n".join(summary))
print()
print("Results:")
print(OUTDIR)
