#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAPER 3 — DEEP PEPTIDE-TO-SPLICE MAPPING
========================================

Workspace:
PAPER3_ROOT (environment variable; defaults to current working directory)

Goal
----
Take the highest-priority HLA-I candidates from the previous rMATS analysis:
  A) shared exact significant rMATS event, same direction in HeLa and SiHa
  B) shared gene + event type, same direction, exact event differs

Then determine how far the existing data can support:
  peptide -> STRG ORF -> transcript nucleotide interval -> splice event

The script also explicitly tests whether the linked significant rMATS event is
listed in rMATS fromGTF.novelJunction / fromGTF.novelSpliceSite outputs.

If a genomic STRG exon GTF/GFF is available in the workspace, the script will
attempt direct peptide-to-junction coordinate mapping. If not, it will STOP the
claim at locus/ORF level rather than guessing.

No web access. No H: access. It reads only the workspace.
"""

from __future__ import annotations

import os

import gzip
import io
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

ROOT = Path(os.environ.get("PAPER3_ROOT", Path.cwd())).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get("PAPER3_OUTPUT_ROOT", ROOT / "results")).expanduser().resolve()
OUTDIR = OUTPUT_ROOT / "DEEP_SPLICE_MAPPING_V3_EVENTTYPE_FIX"

FDR_CUTOFF = 0.05
DPSI_CUTOFF = 0.10
PREFERRED_MODE = "JC"

# Must match the previous analysis.
# HISTORICAL SETTING from the exploratory mapping script.
# The final manuscript does not use the sign of IncLevelDifference to infer
# biological direction because the original rMATS b1/b2 order was not recovered.
# Direct peptide-event reconstruction should therefore be interpreted independently
# of direction-sensitive labels produced by this script.
GROUP1_IS_HPV_POSITIVE = True

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


# ---------------------------------------------------------------------
# GENERIC HELPERS
# ---------------------------------------------------------------------

def normalize_chr(v) -> str:
    if pd.isna(v):
        return ""
    s = str(v).strip()
    s = re.sub(r"^chr", "", s, flags=re.I)
    return s.upper()


def clean_gene(v) -> str:
    if pd.isna(v):
        return ""
    return str(v).strip().strip('"')


def as_int(v):
    x = pd.to_numeric(v, errors="coerce")
    if pd.isna(x):
        return None
    return int(x)


def interval_overlap(a1, a2, b1, b2) -> bool:
    vals = (a1, a2, b1, b2)
    if any(pd.isna(v) for v in vals):
        return False
    return float(a1) <= float(b2) and float(b1) <= float(a2)


def all_occurrences(sequence: str, peptide: str) -> list[int]:
    """Return 1-based AA start positions."""
    sequence = sequence.replace("*", "").upper()
    peptide = peptide.upper()
    out = []
    start = 0
    while True:
        i = sequence.find(peptide, start)
        if i < 0:
            break
        out.append(i + 1)
        start = i + 1
    return out


def canonical_pair(a, b):
    a, b = int(a), int(b)
    return (a, b) if a <= b else (b, a)


# ---------------------------------------------------------------------
# RESOURCE RESOLVER
# ---------------------------------------------------------------------

class ResourceResolver:
    """
    Finds resources either as normal files or as members inside ZIP files.
    The workspace may contain recovered ZIP archives rather than extracted
    copies, so this avoids another manual extraction step.
    """

    def __init__(self, root: Path):
        self.root = root
        self.zips = [p for p in root.rglob("*.zip") if p.is_file()]

    def filesystem_hits(self, basename: str) -> list[Path]:
        return [p for p in self.root.rglob(basename) if p.is_file()]

    def zip_hits(self, basename: str) -> list[tuple[Path, str]]:
        hits = []
        b = basename.lower()
        for zp in self.zips:
            try:
                with zipfile.ZipFile(zp) as z:
                    for member in z.namelist():
                        member_base = member.replace("\\", "/").rstrip("/").split("/")[-1]
                        if member_base.lower() == b:
                            hits.append((zp, member))
            except zipfile.BadZipFile:
                continue
        return hits

    def read_all_text(self, basename: str) -> list[tuple[str, str]]:
        """
        Return list of (source_label, text) for every exact-basename hit.
        Duplicate contents are collapsed by exact text identity.
        """
        found = []

        for p in self.filesystem_hits(basename):
            try:
                if p.suffix.lower() == ".gz":
                    with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as fh:
                        txt = fh.read()
                else:
                    txt = p.read_text(encoding="utf-8", errors="ignore")
                found.append((str(p), txt))
            except Exception:
                pass

        for zp, member in self.zip_hits(basename):
            try:
                with zipfile.ZipFile(zp) as z:
                    raw = z.read(member)
                if member.lower().endswith(".gz"):
                    raw = gzip.decompress(raw)
                txt = raw.decode("utf-8", errors="ignore")
                found.append((f"{zp}::{member}", txt))
            except Exception:
                pass

        # collapse exact duplicate content
        uniq = {}
        for label, txt in found:
            uniq.setdefault(txt, label)
        return [(label, txt) for txt, label in uniq.items()]

    def read_unique_text(self, basename: str) -> tuple[str, str]:
        hits = self.read_all_text(basename)
        if not hits:
            raise FileNotFoundError(f"Could not find resource: {basename}")
        if len(hits) > 1:
            labels = "\n".join("  " + x[0] for x in hits)
            raise RuntimeError(
                f"Multiple non-identical copies found for {basename}; refusing to guess:\n{labels}"
            )
        return hits[0]

    def structure_sources(self) -> list[tuple[str, str]]:
        """
        Return candidate GTF/GFF/GFF3 text sources from filesystem and ZIPs.
        These are scanned later only for six STRG transcript IDs.
        """
        exts = (".gtf", ".gff", ".gff3", ".gtf.gz", ".gff.gz", ".gff3.gz")
        out = []

        for p in self.root.rglob("*"):
            if not p.is_file():
                continue
            low = p.name.lower()
            if not low.endswith(exts):
                continue
            try:
                if low.endswith(".gz"):
                    with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as fh:
                        txt = fh.read()
                else:
                    txt = p.read_text(encoding="utf-8", errors="ignore")
                out.append((str(p), txt))
            except Exception:
                pass

        for zp in self.zips:
            try:
                with zipfile.ZipFile(zp) as z:
                    for member in z.namelist():
                        low = member.lower()
                        if not low.endswith(exts):
                            continue
                        try:
                            raw = z.read(member)
                            if low.endswith(".gz"):
                                raw = gzip.decompress(raw)
                            txt = raw.decode("utf-8", errors="ignore")
                            out.append((f"{zp}::{member}", txt))
                        except Exception:
                            pass
            except zipfile.BadZipFile:
                pass

        # collapse identical text
        uniq = {}
        for label, txt in out:
            uniq.setdefault(txt, label)
        return [(label, txt) for txt, label in uniq.items()]


# ---------------------------------------------------------------------
# FASTA + STRG MAP
# ---------------------------------------------------------------------

def parse_fasta(text: str) -> dict[str, str]:
    seqs = {}
    header = None
    chunks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seqs[header] = "".join(chunks)
            header = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        seqs[header] = "".join(chunks)
    return seqs


def parse_strg_map(text: str) -> dict[str, dict]:
    """
    Old map TSV is effectively two columns but often has no formal header.
    Example:
    P3I_HELA_STRG_000001 <tab>
    STRG.1.1.p6 ... STRG.1.1:1-123(+)
    """
    out = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\r\n").split("\t", 1)
        if len(parts) < 2:
            continue

        p3i = parts[0].strip()
        desc = parts[1].strip()

        if not p3i.startswith("P3I_"):
            continue

        m_trans = re.search(r"\b(STRG\.\d+\.\d+)(?:\.p\d+)?\b", desc)
        m_orf = re.search(
            r"\b(STRG\.\d+\.\d+):(\d+)-(\d+)\(([+-])\)",
            desc
        )

        rec = {
            "P3I_ID": p3i,
            "Descriptor": desc,
            "STRG_transcript_from_map": m_trans.group(1) if m_trans else "",
            "ORF_nt_start": int(m_orf.group(2)) if m_orf else np.nan,
            "ORF_nt_end": int(m_orf.group(3)) if m_orf else np.nan,
            "ORF_orientation": m_orf.group(4) if m_orf else "",
        }
        out[p3i] = rec
    return out


def peptide_transcript_interval(
    aa_start: int,
    peptide_len: int,
    orf_start: int,
    orf_end: int,
    orf_orientation: str,
):
    aa_end = aa_start + peptide_len - 1

    if orf_orientation == "+":
        nt_start = orf_start + (aa_start - 1) * 3
        nt_end = nt_start + peptide_len * 3 - 1
    elif orf_orientation == "-":
        nt_start = orf_end - 3 * aa_end + 1
        nt_end = orf_end - 3 * (aa_start - 1)
    else:
        return None, None

    return int(min(nt_start, nt_end)), int(max(nt_start, nt_end))


# ---------------------------------------------------------------------
# rMATS
# ---------------------------------------------------------------------

def event_key(df: pd.DataFrame, event_type: str) -> pd.Series:
    cols = COORD_COLUMNS[event_type]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{event_type}: missing rMATS coordinate columns: {', '.join(missing)}"
        )

    out = (
        pd.Series(event_type, index=df.index).astype(str)
        + "|" + df["_chr_norm"].astype(str)
        + "|" + df["strand"].astype(str)
    )

    for c in cols:
        x = pd.to_numeric(df[c], errors="coerce").astype("Int64").astype(str)
        out = out + "|" + x

    return out


def load_rmats_table(
    resolver: ResourceResolver,
    comparison: str,
    event_type: str,
) -> pd.DataFrame:
    """
    Read MATS.JC first, fallback JCEC.
    Since both comparisons have identical basenames, source label is used to
    select the copy whose path/member contains the comparison name.
    """
    chosen = None
    used_mode = None

    for mode in (PREFERRED_MODE, "JCEC"):
        basename = f"{event_type}.MATS.{mode}.txt"
        hits = resolver.read_all_text(basename)

        comp_hits = [
            (label, txt) for label, txt in hits
            if comparison.lower() in label.lower()
        ]

        if len(comp_hits) == 1:
            chosen = comp_hits[0]
            used_mode = mode
            break
        elif len(comp_hits) > 1:
            labels = "\n".join("  " + h[0] for h in comp_hits)
            raise RuntimeError(
                f"Multiple non-identical {basename} resources for {comparison}:\n{labels}"
            )

    if chosen is None:
        raise FileNotFoundError(
            f"Could not find {event_type}.MATS.JC/JCEC for {comparison}"
        )

    label, txt = chosen
    df = pd.read_csv(io.StringIO(txt), sep="\t", low_memory=False)

    required = ["ID", "FDR", "IncLevelDifference", "chr", "strand"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")

    df["ID"] = df["ID"].astype(str)
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["IncLevelDifference"] = pd.to_numeric(
        df["IncLevelDifference"], errors="coerce"
    )
    df["_chr_norm"] = df["chr"].map(normalize_chr)
    df["_gene_symbol"] = (
        df["geneSymbol"].map(clean_gene)
        if "geneSymbol" in df.columns else ""
    )
    df["_gene_id"] = (
        df["GeneID"].map(clean_gene)
        if "GeneID" in df.columns else ""
    )
    df["_gene_key"] = np.where(
        pd.Series(df["_gene_symbol"]).astype(str).str.len() > 0,
        df["_gene_symbol"],
        df["_gene_id"],
    )

    df["EventType"] = event_type
    df["Comparison"] = comparison
    df["rMATS_Mode"] = used_mode
    df["SourceFile"] = label
    df["EventKey"] = event_key(df, event_type)

    numeric = df[COORD_COLUMNS[event_type]].apply(
        pd.to_numeric, errors="coerce"
    )
    df["EventStart"] = numeric.min(axis=1)
    df["EventEnd"] = numeric.max(axis=1)

    df["Significant"] = (
        df["FDR"].le(FDR_CUTOFF)
        & df["IncLevelDifference"].abs().ge(DPSI_CUTOFF)
    )

    if GROUP1_IS_HPV_POSITIVE:
        df["DirectionSign"] = np.sign(
            df["IncLevelDifference"]
        ).astype("Int64")
    else:
        df["DirectionSign"] = pd.Series(
            [pd.NA] * len(df), dtype="Int64"
        )

    # Novel-junction / novel-splice-site membership by rMATS event ID.
    for novelty_type, prefix in [
        ("IsNovelJunction", "fromGTF.novelJunction"),
        ("IsNovelSpliceSite", "fromGTF.novelSpliceSite"),
    ]:
        basename = f"{prefix}.{event_type}.txt"
        hits = resolver.read_all_text(basename)
        comp_hits = [
            (lab, tx) for lab, tx in hits
            if comparison.lower() in lab.lower()
        ]

        ids = set()
        if len(comp_hits) == 1:
            ndf = pd.read_csv(io.StringIO(comp_hits[0][1]), sep="\t")
            if "ID" in ndf.columns:
                ids = set(ndf["ID"].astype(str))
        elif len(comp_hits) > 1:
            # If duplicate non-identical copies exist, union is safer for
            # annotation; report source ambiguity separately later.
            for _, tx in comp_hits:
                ndf = pd.read_csv(io.StringIO(tx), sep="\t")
                if "ID" in ndf.columns:
                    ids.update(ndf["ID"].astype(str))

        df[novelty_type] = df["ID"].isin(ids)

    return df


def load_all_rmats(resolver: ResourceResolver, comparison: str) -> pd.DataFrame:
    return pd.concat(
        [load_rmats_table(resolver, comparison, et) for et in EVENT_TYPES],
        ignore_index=True,
        sort=False,
    )


def build_shared_exact(hela: pd.DataFrame, siha: pd.DataFrame) -> pd.DataFrame:
    h = hela[hela["Significant"]].copy()
    s = siha[siha["Significant"]].copy()

    h = h.sort_values("FDR").drop_duplicates("EventKey")
    s = s.sort_values("FDR").drop_duplicates("EventKey")

    hcols = [
        "EventKey", "EventType", "_chr_norm", "strand",
        "_gene_key", "_gene_symbol", "_gene_id",
        "EventStart", "EventEnd",
        "ID", "FDR", "IncLevelDifference", "DirectionSign",
        "IsNovelJunction", "IsNovelSpliceSite",
    ]

    hout = h[hcols].rename(columns={
        "_gene_key": "Gene",
        "_gene_symbol": "GeneSymbol",
        "_gene_id": "GeneID",
        "ID": "ID_HeLa",
        "FDR": "FDR_HeLa",
        "IncLevelDifference": "dPSI_HeLa",
        "DirectionSign": "DirectionSign_HeLa",
        "IsNovelJunction": "NovelJunction_HeLa",
        "IsNovelSpliceSite": "NovelSpliceSite_HeLa",
    })

    sout = s[[
        "EventKey", "ID", "FDR", "IncLevelDifference", "DirectionSign",
        "IsNovelJunction", "IsNovelSpliceSite"
    ]].rename(columns={
        "ID": "ID_SiHa",
        "FDR": "FDR_SiHa",
        "IncLevelDifference": "dPSI_SiHa",
        "DirectionSign": "DirectionSign_SiHa",
        "IsNovelJunction": "NovelJunction_SiHa",
        "IsNovelSpliceSite": "NovelSpliceSite_SiHa",
    })

    out = hout.merge(sout, on="EventKey", how="inner")

    if GROUP1_IS_HPV_POSITIVE:
        out["SameDirection"] = (
            out["DirectionSign_HeLa"].eq(out["DirectionSign_SiHa"])
            & out["DirectionSign_HeLa"].ne(0)
        )
    else:
        out["SameDirection"] = pd.NA

    out["NovelJunction_in_both_runs"] = (
        out["NovelJunction_HeLa"].fillna(False)
        & out["NovelJunction_SiHa"].fillna(False)
    )
    out["NovelJunction_in_either_run"] = (
        out["NovelJunction_HeLa"].fillna(False)
        | out["NovelJunction_SiHa"].fillna(False)
    )

    return out


# ---------------------------------------------------------------------
# EVENT JUNCTION COORDINATES
# ---------------------------------------------------------------------

def event_junction_pairs(row: pd.Series) -> set[tuple[int, int]]:
    """
    Junctions are represented as genomic intron boundary breakpoints:
      (left exon end, right exon start - 1)
    rMATS *_ES columns are 0-based exon starts and *_EE/end are exon ends,
    so they already align naturally to breakpoint representation.

    IMPORTANT:
    Deep-link rows store the linked event type as EventType_linked.
    Fall back to EventType only for compatibility with raw rMATS rows.
    """
    et = row.get("EventType_linked", row.get("EventType", ""))
    pairs = set()

    def add(a, b):
        if pd.isna(a) or pd.isna(b):
            return
        pairs.add(canonical_pair(int(a), int(b)))

    if et == "SE":
        add(row["upstreamEE"], row["exonStart_0base"])
        add(row["exonEnd"], row["downstreamES"])
        add(row["upstreamEE"], row["downstreamES"])

    elif et in ("A3SS", "A5SS"):
        long_s = int(row["longExonStart_0base"])
        long_e = int(row["longExonEnd"])
        short_s = int(row["shortES"])
        short_e = int(row["shortEE"])
        flank_s = int(row["flankingES"])
        flank_e = int(row["flankingEE"])

        # flanking exon is genomically upstream
        if flank_e <= min(long_s, short_s):
            add(flank_e, long_s)
            add(flank_e, short_s)
        # flanking exon is genomically downstream
        elif max(long_e, short_e) <= flank_s:
            add(long_e, flank_s)
            add(short_e, flank_s)
        else:
            # fallback: retain every plausible boundary pairing
            add(flank_e, long_s)
            add(flank_e, short_s)
            add(long_e, flank_s)
            add(short_e, flank_s)

    elif et == "MXE":
        add(row["upstreamEE"], row["1stExonStart_0base"])
        add(row["1stExonEnd"], row["downstreamES"])
        add(row["upstreamEE"], row["2ndExonStart_0base"])
        add(row["2ndExonEnd"], row["downstreamES"])

    elif et == "RI":
        add(row["upstreamEE"], row["downstreamES"])

    return pairs


# ---------------------------------------------------------------------
# GENOMIC STRG EXON STRUCTURE
# ---------------------------------------------------------------------

def parse_attrs(attr: str) -> dict[str, str]:
    out = {}

    # GTF form: key "value";
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        out[m.group(1)] = m.group(2)

    # GFF3 form: key=value;
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out.setdefault(k.strip(), v.strip())

    return out


def scan_exon_structures(
    sources: list[tuple[str, str]],
    transcript_ids: set[str],
):
    """
    Search only for candidate STRG transcript IDs.

    A valid genomic exon structure must have:
    - feature = exon
    - chromosome-like seqid rather than STRG transcript seqid
    - a transcript ID / Parent containing one of the target STRG IDs

    TransDecoder GFF3 files have seqid=STRG.* and are therefore intentionally
    NOT accepted as genomic exon structure.
    """
    by_transcript = defaultdict(list)
    source_inventory = []

    target_regex = re.compile(
        r"\b(" + "|".join(re.escape(x) for x in sorted(transcript_ids, key=len, reverse=True)) + r")\b"
    ) if transcript_ids else None

    for label, text in sources:
        matched_lines = 0
        accepted_exons = 0
        looks_transdecoder = False

        if target_regex is None:
            continue

        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            if not target_regex.search(line):
                continue

            matched_lines += 1
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = cols[:9]

            if seqid.startswith("STRG."):
                looks_transdecoder = True
                continue

            if feature.lower() != "exon":
                continue

            am = parse_attrs(attrs)
            candidate_text = " ".join([
                am.get("transcript_id", ""),
                am.get("Parent", ""),
                am.get("ID", ""),
                attrs,
            ])

            tm = target_regex.search(candidate_text)
            if not tm:
                continue

            tid = tm.group(1)

            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue

            by_transcript[tid].append({
                "source": label,
                "chr": normalize_chr(seqid),
                "start": s,
                "end": e,
                "strand": strand,
            })
            accepted_exons += 1

        source_inventory.append({
            "StructureSource": label,
            "TargetMatchingLines": matched_lines,
            "AcceptedGenomicExons": accepted_exons,
            "LooksLikeTransDecoderTranscriptCoordinateGFF": looks_transdecoder,
        })

    # Collapse duplicate exon records.
    for tid, exons in by_transcript.items():
        uniq = {}
        for x in exons:
            key = (x["chr"], x["start"], x["end"], x["strand"])
            uniq.setdefault(key, x)
        by_transcript[tid] = list(uniq.values())

    return by_transcript, pd.DataFrame(source_inventory)


def transcript_nt_to_genomic_segments(
    exons: list[dict],
    tx_start: int,
    tx_end: int,
):
    """
    Map 1-based transcript interval to genomic exon segments.
    Assumes assembled transcript follows its annotated transcript strand.
    """
    if not exons:
        return [], "NO_EXONS"

    chromosomes = {e["chr"] for e in exons}
    strands = {e["strand"] for e in exons}

    if len(chromosomes) != 1 or len(strands) != 1:
        return [], "AMBIGUOUS_EXON_STRUCTURE"

    strand = next(iter(strands))
    if strand not in ("+", "-"):
        return [], "UNKNOWN_TRANSCRIPT_STRAND"

    exons_ordered = sorted(
        exons,
        key=lambda e: e["start"],
        reverse=(strand == "-")
    )

    cursor = 1
    tx_exons = []
    for i, e in enumerate(exons_ordered, start=1):
        length = e["end"] - e["start"] + 1
        tx_s = cursor
        tx_e = cursor + length - 1
        tx_exons.append((i, e, tx_s, tx_e))
        cursor = tx_e + 1

    if tx_start < 1 or tx_end > cursor - 1:
        return [], "PEPTIDE_TX_INTERVAL_OUTSIDE_ANNOTATED_TRANSCRIPT"

    segs = []

    for exon_index, e, ex_tx_s, ex_tx_e in tx_exons:
        ov_s = max(tx_start, ex_tx_s)
        ov_e = min(tx_end, ex_tx_e)
        if ov_s > ov_e:
            continue

        offset_s = ov_s - ex_tx_s
        offset_e = ov_e - ex_tx_s

        if strand == "+":
            g_s = e["start"] + offset_s
            g_e = e["start"] + offset_e
        else:
            g_high = e["end"] - offset_s
            g_low = e["end"] - offset_e
            g_s, g_e = min(g_low, g_high), max(g_low, g_high)

        segs.append({
            "exon_index": exon_index,
            "tx_start": ov_s,
            "tx_end": ov_e,
            "genomic_start": g_s,
            "genomic_end": g_e,
            "chr": e["chr"],
            "strand": strand,
        })

    return segs, "OK"


def crossed_exon_junction_pairs(exons: list[dict], segments: list[dict]):
    """
    For a peptide spanning >=2 exons, return genomic intron breakpoint pairs
    crossed by the peptide.
    """
    if len(segments) < 2:
        return set()

    chromosomes = {e["chr"] for e in exons}
    strands = {e["strand"] for e in exons}
    if len(chromosomes) != 1 or len(strands) != 1:
        return set()

    strand = next(iter(strands))
    exons_ordered = sorted(
        exons,
        key=lambda e: e["start"],
        reverse=(strand == "-")
    )

    by_idx = {i + 1: e for i, e in enumerate(exons_ordered)}
    idxs = sorted({s["exon_index"] for s in segments})

    pairs = set()
    for a, b in zip(idxs[:-1], idxs[1:]):
        if b != a + 1:
            continue
        e1 = by_idx[a]
        e2 = by_idx[b]

        left, right = sorted(
            [e1, e2],
            key=lambda e: e["start"]
        )
        pairs.add(canonical_pair(left["end"], right["start"] - 1))

    return pairs


# ---------------------------------------------------------------------
# CORE CANDIDATE + EVENT LINKING
# ---------------------------------------------------------------------

EXACT_CLASS = (
    "Shared exact significant rMATS event, same direction in HeLa and SiHa"
)
GENE_CLASS_PREFIX = "Shared gene+event-type"


def load_classification(resolver: ResourceResolver) -> pd.DataFrame:
    _, txt = resolver.read_unique_text("06_HLA54_splicing_classification.tsv")
    return pd.read_csv(io.StringIO(txt), sep="\t", low_memory=False)


def select_core_candidates(class_df: pd.DataFrame) -> pd.DataFrame:
    exact = class_df[class_df["SplicingSupportClass"].eq(EXACT_CLASS)].copy()
    gene = class_df[
        class_df["SplicingSupportClass"].astype(str).str.startswith(GENE_CLASS_PREFIX)
    ].copy()

    core = pd.concat([exact, gene], ignore_index=True)
    core["CorePriorityClass"] = np.where(
        core["SplicingSupportClass"].eq(EXACT_CLASS),
        "A_EXACT_SHARED_SAME_DIRECTION",
        "B_SHARED_GENE_EVENTTYPE"
    )
    return core


def candidate_locus_events(
    candidate: pd.Series,
    events: pd.DataFrame,
) -> pd.DataFrame:
    chrom = normalize_chr(candidate.get("STRG_chr", ""))
    strand = str(candidate.get("STRG_strand", ""))
    st = pd.to_numeric(candidate.get("STRG_start"), errors="coerce")
    en = pd.to_numeric(candidate.get("STRG_end"), errors="coerce")

    if not chrom or pd.isna(st) or pd.isna(en):
        return events.iloc[0:0].copy()

    x = events[
        events["Significant"]
        & events["_chr_norm"].eq(chrom)
        & events["strand"].astype(str).eq(strand)
    ].copy()

    if x.empty:
        return x

    return x[
        x.apply(
            lambda r: interval_overlap(
                st, en, r["EventStart"], r["EventEnd"]
            ),
            axis=1,
        )
    ].copy()


def link_core_to_events(
    core: pd.DataFrame,
    hela: pd.DataFrame,
    siha: pd.DataFrame,
    shared: pd.DataFrame,
):
    rows = []

    shared_same = shared.copy()
    if GROUP1_IS_HPV_POSITIVE:
        shared_same = shared_same[shared_same["SameDirection"].eq(True)].copy()

    for _, c in core.iterrows():
        cell = str(c["CellLine"])
        if cell.lower() == "hela":
            own = hela
            other = siha
        elif cell.lower() == "siha":
            own = siha
            other = hela
        else:
            continue

        own_loc = candidate_locus_events(c, own)

        if c["CorePriorityClass"] == "A_EXACT_SHARED_SAME_DIRECTION":
            chrom = normalize_chr(c.get("STRG_chr", ""))
            strand = str(c.get("STRG_strand", ""))
            st = pd.to_numeric(c.get("STRG_start"), errors="coerce")
            en = pd.to_numeric(c.get("STRG_end"), errors="coerce")

            matches = shared_same[
                shared_same["_chr_norm"].eq(chrom)
                & shared_same["strand"].astype(str).eq(strand)
            ].copy()

            if not matches.empty:
                matches = matches[
                    matches.apply(
                        lambda r: interval_overlap(
                            st, en, r["EventStart"], r["EventEnd"]
                        ),
                        axis=1,
                    )
                ]

            for _, e in matches.iterrows():
                row = c.to_dict()
                row.update({
                    "LinkType": "EXACT_SHARED_EVENT",
                    "EventKey": e["EventKey"],
                    "EventType_linked": e["EventType"],
                    "Gene_linked": e["Gene"],
                    "FDR_HeLa": e["FDR_HeLa"],
                    "dPSI_HeLa": e["dPSI_HeLa"],
                    "FDR_SiHa": e["FDR_SiHa"],
                    "dPSI_SiHa": e["dPSI_SiHa"],
                    "NovelJunction_HeLa": e["NovelJunction_HeLa"],
                    "NovelJunction_SiHa": e["NovelJunction_SiHa"],
                    "NovelJunction_in_both_runs": e["NovelJunction_in_both_runs"],
                    "NovelJunction_in_either_run": e["NovelJunction_in_either_run"],
                    "NovelSpliceSite_HeLa": e["NovelSpliceSite_HeLa"],
                    "NovelSpliceSite_SiHa": e["NovelSpliceSite_SiHa"],
                })

                # Attach coordinate columns from own comparison event.
                own_match = own[own["EventKey"].eq(e["EventKey"])]
                if not own_match.empty:
                    oe = own_match.iloc[0]
                    for col in COORD_COLUMNS[e["EventType"]]:
                        row[col] = oe.get(col)
                    row["OwnComparisonEventID"] = oe.get("ID")
                    row["OwnComparisonNovelJunction"] = oe.get("IsNovelJunction")
                    row["OwnComparisonNovelSpliceSite"] = oe.get("IsNovelSpliceSite")

                rows.append(row)

        else:
            # Shared gene + event type + same direction, exact coordinates differ.
            for _, oe in own_loc.iterrows():
                other_sig = other[
                    other["Significant"]
                    & other["_gene_key"].astype(str).eq(str(oe["_gene_key"]))
                    & other["EventType"].astype(str).eq(str(oe["EventType"]))
                ].copy()

                if GROUP1_IS_HPV_POSITIVE:
                    other_sig = other_sig[
                        other_sig["DirectionSign"].eq(oe["DirectionSign"])
                    ]

                if other_sig.empty:
                    continue

                for _, xe in other_sig.iterrows():
                    row = c.to_dict()
                    row.update({
                        "LinkType": "SHARED_GENE_EVENTTYPE_EXACT_EVENT_DIFFERS",
                        "EventKey": oe["EventKey"],
                        "EventType_linked": oe["EventType"],
                        "Gene_linked": oe["_gene_key"],
                        "OwnComparisonEventID": oe["ID"],
                        "OwnComparisonFDR": oe["FDR"],
                        "OwnComparison_dPSI": oe["IncLevelDifference"],
                        "OwnComparisonNovelJunction": oe["IsNovelJunction"],
                        "OwnComparisonNovelSpliceSite": oe["IsNovelSpliceSite"],
                        "OtherComparisonEventID": xe["ID"],
                        "OtherComparisonFDR": xe["FDR"],
                        "OtherComparison_dPSI": xe["IncLevelDifference"],
                        "OtherComparisonNovelJunction": xe["IsNovelJunction"],
                        "OtherComparisonNovelSpliceSite": xe["IsNovelSpliceSite"],
                        "OtherEventKey": xe["EventKey"],
                    })
                    for col in COORD_COLUMNS[oe["EventType"]]:
                        row[col] = oe.get(col)
                    rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():
    if not ROOT.exists():
        raise FileNotFoundError(f"Workspace not found: {ROOT}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    resolver = ResourceResolver(ROOT)

    # 1) Core candidates from prior analysis
    class_df = load_classification(resolver)
    core = select_core_candidates(class_df)
    core.to_csv(
        OUTDIR / "01_CORE_candidates.tsv",
        sep="\t", index=False
    )

    # 2) Reload rMATS + novelty annotations
    hela = load_all_rmats(resolver, "rMATS_HELA_vs_C33A")
    siha = load_all_rmats(resolver, "rMATS_SIHA_vs_C33A")
    shared = build_shared_exact(hela, siha)

    links = link_core_to_events(core, hela, siha, shared)
    links.to_csv(
        OUTDIR / "02_CORE_candidate_to_rMATS_event_links.tsv",
        sep="\t", index=False
    )

    # 3) Load STRG protein DBs + maps from filesystem or recovered ZIP
    protein_sources = {}
    map_sources = {}

    for cell, prefix in [
        ("HeLa", "HELA"),
        ("SiHa", "SIHA"),
        ("C33A", "C33A"),
    ]:
        fasta_name = f"{prefix}_STRG_Personalized.fasta"
        map_name = f"{prefix}_STRG_Map.tsv"

        fasta_hits = resolver.read_all_text(fasta_name)
        map_hits = resolver.read_all_text(map_name)

        protein_sources[cell] = []
        for label, txt in fasta_hits:
            protein_sources[cell].append((label, parse_fasta(txt)))

        map_sources[cell] = []
        for label, txt in map_hits:
            map_sources[cell].append((label, parse_strg_map(txt)))

    # 4) Peptide -> ORF -> transcript nucleotide interval
    orf_rows = []

    for _, c in core.iterrows():
        cell = str(c["CellLine"])
        peptide = str(c["Peptide"]).upper()
        p3i = str(c["Protein_ID"])
        strg_id = str(c["STRG_transcript_id"])

        seq_records = []
        for label, db in protein_sources.get(cell, []):
            if p3i in db:
                seq_records.append((label, db[p3i].replace("*", "")))

        # collapse identical sequence
        seq_uniq = {}
        for label, seq in seq_records:
            seq_uniq.setdefault(seq, label)

        map_records = []
        for label, mp in map_sources.get(cell, []):
            if p3i in mp:
                map_records.append((label, mp[p3i]))

        # collapse identical descriptor
        map_uniq = {}
        for label, rec in map_records:
            map_uniq.setdefault(rec["Descriptor"], (label, rec))

        base = c.to_dict()

        if len(seq_uniq) == 0:
            base.update({
                "ORF_mapping_status": "NO_STRG_PROTEIN_SEQUENCE",
            })
            orf_rows.append(base)
            continue

        if len(seq_uniq) > 1:
            base.update({
                "ORF_mapping_status": "AMBIGUOUS_NONIDENTICAL_STRG_PROTEIN_SEQUENCES",
                "n_sequence_versions": len(seq_uniq),
            })
            orf_rows.append(base)
            continue

        seq, seq_label = next(iter(seq_uniq.items()))
        aa_positions = all_occurrences(seq, peptide)

        if len(map_uniq) == 0:
            base.update({
                "ORF_mapping_status": "NO_STRG_MAP_RECORD",
                "ProteinSequenceSource": seq_label,
                "Protein_length_aa": len(seq),
                "Peptide_AA_starts": ";".join(map(str, aa_positions)),
            })
            orf_rows.append(base)
            continue

        if len(map_uniq) > 1:
            base.update({
                "ORF_mapping_status": "AMBIGUOUS_NONIDENTICAL_STRG_MAP_RECORDS",
                "ProteinSequenceSource": seq_label,
                "Protein_length_aa": len(seq),
                "Peptide_AA_starts": ";".join(map(str, aa_positions)),
                "n_map_versions": len(map_uniq),
            })
            orf_rows.append(base)
            continue

        _, (map_label, rec) = next(iter(map_uniq.items()))

        if len(aa_positions) == 0:
            status = "PEPTIDE_NOT_FOUND_IN_ASSIGNED_STRG_PROTEIN"
        elif len(aa_positions) > 1:
            status = "PEPTIDE_OCCURS_MULTIPLE_TIMES_IN_ASSIGNED_STRG_PROTEIN"
        else:
            status = "OK"

        base.update({
            "ORF_mapping_status": status,
            "ProteinSequenceSource": seq_label,
            "STRGMapSource": map_label,
            "Protein_length_aa": len(seq),
            "Peptide_AA_starts": ";".join(map(str, aa_positions)),
            "Map_STRG_transcript": rec["STRG_transcript_from_map"],
            "Candidate_vs_Map_STRG_match":
                rec["STRG_transcript_from_map"] == strg_id,
            "ORF_nt_start": rec["ORF_nt_start"],
            "ORF_nt_end": rec["ORF_nt_end"],
            "ORF_orientation": rec["ORF_orientation"],
            "STRG_Map_descriptor": rec["Descriptor"],
        })

        if len(aa_positions) == 1 and not pd.isna(rec["ORF_nt_start"]):
            aa_start = aa_positions[0]
            aa_end = aa_start + len(peptide) - 1
            tx_s, tx_e = peptide_transcript_interval(
                aa_start,
                len(peptide),
                int(rec["ORF_nt_start"]),
                int(rec["ORF_nt_end"]),
                rec["ORF_orientation"],
            )
            base.update({
                "Peptide_AA_start": aa_start,
                "Peptide_AA_end": aa_end,
                "Peptide_transcript_nt_start": tx_s,
                "Peptide_transcript_nt_end": tx_e,
            })

        orf_rows.append(base)

    orf_df = pd.DataFrame(orf_rows)
    orf_df.to_csv(
        OUTDIR / "03_CORE_peptide_to_STRG_ORF_mapping.tsv",
        sep="\t", index=False
    )

    # 5) Search genomic exon structure — CELL-LINE SPECIFIC
    #
    # IMPORTANT:
    # STRG transcript IDs are local to each StringTie assembly and are NOT
    # globally unique across cell lines. Therefore STRG.14863.1 in HeLa must
    # never be merged with an identically named transcript from SiHa/C33A.
    #
    # Use only the explicitly copied cell-line-specific assembled GTFs.
    transcript_ids_by_cell = {
        cell: set(
            core.loc[
                core["CellLine"].astype(str).str.lower().eq(cell.lower()),
                "STRG_transcript_id"
            ].astype(str)
        )
        for cell in ("HeLa", "SiHa", "C33A")
    }

    structure_inventory_parts = []
    exon_by_cell_tx = {}

    gtf_by_cell = {
        "HeLa": ROOT / "TRANSCRIPT_STRUCTURE" / "HeLa" / "HELA_assembled.gtf",
        "SiHa": ROOT / "TRANSCRIPT_STRUCTURE" / "SiHa" / "SIHA_assembled.gtf",
        "C33A": ROOT / "TRANSCRIPT_STRUCTURE" / "C33A" / "C33A_assembled.gtf",
    }

    for cell, gtf_path in gtf_by_cell.items():
        tids = transcript_ids_by_cell.get(cell, set())
        if not tids:
            continue

        if not gtf_path.exists():
            structure_inventory_parts.append(pd.DataFrame([{
                "CellLine": cell,
                "StructureSource": str(gtf_path),
                "TargetMatchingLines": 0,
                "AcceptedGenomicExons": 0,
                "LooksLikeTransDecoderTranscriptCoordinateGFF": False,
                "Status": "MISSING_CELL_SPECIFIC_GTF",
            }]))
            continue

        txt = gtf_path.read_text(encoding="utf-8", errors="ignore")
        exon_map, inv = scan_exon_structures(
            [(str(gtf_path), txt)],
            tids
        )

        if inv.empty:
            inv = pd.DataFrame([{
                "StructureSource": str(gtf_path),
                "TargetMatchingLines": 0,
                "AcceptedGenomicExons": 0,
                "LooksLikeTransDecoderTranscriptCoordinateGFF": False,
            }])

        inv.insert(0, "CellLine", cell)
        inv["Status"] = "OK"
        structure_inventory_parts.append(inv)

        for tid, exons in exon_map.items():
            exon_by_cell_tx[(cell, tid)] = exons

    structure_inventory = (
        pd.concat(structure_inventory_parts, ignore_index=True, sort=False)
        if structure_inventory_parts else pd.DataFrame()
    )

    structure_inventory.to_csv(
        OUTDIR / "04_transcript_structure_source_inventory.tsv",
        sep="\t", index=False
    )

    # 6) Deep resolution: candidate + event link + peptide genomic mapping
    deep_rows = []

    # If one candidate has multiple linked events, retain every event row.
    if links.empty:
        links_iter = [None]
    else:
        links_iter = list(links.iterrows())

    for _, c in core.iterrows():
        p3i = str(c["Protein_ID"])
        peptide = str(c["Peptide"])
        txid = str(c["STRG_transcript_id"])

        om = orf_df[
            orf_df["Protein_ID"].astype(str).eq(p3i)
            & orf_df["Peptide"].astype(str).eq(peptide)
        ]
        omrow = om.iloc[0].to_dict() if not om.empty else {}

        candidate_links = links[
            links["Protein_ID"].astype(str).eq(p3i)
            & links["Peptide"].astype(str).eq(peptide)
        ] if not links.empty else pd.DataFrame()

        if candidate_links.empty:
            candidate_links = pd.DataFrame([{
                "Protein_ID": p3i,
                "Peptide": peptide,
                "LinkType": "NO_LINKED_EVENT_ROW",
            }])

        for _, lk in candidate_links.iterrows():
            row = c.to_dict()
            row.update(omrow)
            row.update(lk.to_dict())

            exons = exon_by_cell_tx.get((str(c["CellLine"]), txid), [])

            tx_s = as_int(row.get("Peptide_transcript_nt_start"))
            tx_e = as_int(row.get("Peptide_transcript_nt_end"))

            if not exons:
                row.update({
                    "Genomic_exon_mapping_status":
                        "NO_GENOMIC_STRG_EXON_STRUCTURE_IN_WORKSPACE",
                    "Peptide_genomic_segments": "",
                    "Peptide_crosses_exon_junction": np.nan,
                    "Crossed_junction_pairs": "",
                    "Event_junction_pairs": "",
                    "Direct_event_junction_match": False,
                    "FinalStructuralResolution":
                        "UNRESOLVED_LOCUS_ORF_LEVEL_ONLY",
                })
                deep_rows.append(row)
                continue

            if tx_s is None or tx_e is None:
                row.update({
                    "Genomic_exon_mapping_status":
                        "NO_UNIQUE_PEPTIDE_TRANSCRIPT_INTERVAL",
                    "Peptide_genomic_segments": "",
                    "Peptide_crosses_exon_junction": np.nan,
                    "Crossed_junction_pairs": "",
                    "Event_junction_pairs": "",
                    "Direct_event_junction_match": False,
                    "FinalStructuralResolution":
                        "UNRESOLVED_NO_UNIQUE_PEPTIDE_POSITION",
                })
                deep_rows.append(row)
                continue

            segments, map_status = transcript_nt_to_genomic_segments(
                exons, tx_s, tx_e
            )

            crossed = crossed_exon_junction_pairs(exons, segments)

            # Reconstruct event row from linked coordinate columns.
            et = row.get("EventType_linked", "")
            event_pairs = set()
            if et in EVENT_TYPES:
                try:
                    event_pairs = event_junction_pairs(pd.Series(row))
                except Exception:
                    event_pairs = set()

            direct = bool(crossed & event_pairs)

            seg_text = ";".join(
                f'{s["chr"]}:{s["genomic_start"]}-{s["genomic_end"]}'
                for s in segments
            )
            crossed_text = ";".join(
                f"{a}-{b}" for a, b in sorted(crossed)
            )
            event_text = ";".join(
                f"{a}-{b}" for a, b in sorted(event_pairs)
            )

            own_novel = bool(row.get("OwnComparisonNovelJunction", False))
            novel_both = bool(row.get("NovelJunction_in_both_runs", False))

            if direct and novel_both:
                final_class = "DIRECT_SHARED_NOVEL_JUNCTION_PEPTIDE"
            elif direct and own_novel:
                final_class = "DIRECT_NOVEL_JUNCTION_PEPTIDE"
            elif direct:
                final_class = "DIRECT_SPLICE_EVENT_JUNCTION_PEPTIDE_NOT_NOVEL_ANNOTATED"
            elif len(segments) > 0:
                final_class = "PEPTIDE_GENOMICALLY_MAPPED_BUT_NOT_ACROSS_LINKED_EVENT_JUNCTION"
            else:
                final_class = "UNRESOLVED_GENOMIC_MAPPING"

            row.update({
                "Genomic_exon_mapping_status": map_status,
                "n_genomic_exons_for_STRG": len(exons),
                "Peptide_genomic_segments": seg_text,
                "Peptide_crosses_exon_junction": len(crossed) > 0,
                "Crossed_junction_pairs": crossed_text,
                "Event_junction_pairs": event_text,
                "Direct_event_junction_match": direct,
                "FinalStructuralResolution": final_class,
            })

            deep_rows.append(row)

    deep = pd.DataFrame(deep_rows)
    deep.to_csv(
        OUTDIR / "05_CORE_peptide_to_splice_resolution.tsv",
        sep="\t", index=False
    )

    # 7) Strongest candidates
    if not deep.empty:
        strongest = deep[
            deep["FinalStructuralResolution"].astype(str).isin([
                "DIRECT_SHARED_NOVEL_JUNCTION_PEPTIDE",
                "DIRECT_NOVEL_JUNCTION_PEPTIDE",
                "DIRECT_SPLICE_EVENT_JUNCTION_PEPTIDE_NOT_NOVEL_ANNOTATED",
            ])
        ].copy()
    else:
        strongest = deep.copy()

    strongest.to_csv(
        OUTDIR / "06_STRONGEST_direct_splice_candidates.tsv",
        sep="\t", index=False
    )

    # 8) Candidate-level compact summary
    compact_rows = []

    for _, c in core.iterrows():
        p3i = str(c["Protein_ID"])
        pep = str(c["Peptide"])

        d = deep[
            deep["Protein_ID"].astype(str).eq(p3i)
            & deep["Peptide"].astype(str).eq(pep)
        ].copy()

        final_classes = sorted(
            set(d["FinalStructuralResolution"].dropna().astype(str))
        ) if not d.empty else []

        novel_own = False
        novel_both = False
        if not d.empty:
            if "OwnComparisonNovelJunction" in d.columns:
                novel_own = d["OwnComparisonNovelJunction"].fillna(False).astype(bool).any()
            if "NovelJunction_in_both_runs" in d.columns:
                novel_both = d["NovelJunction_in_both_runs"].fillna(False).astype(bool).any()

        direct = (
            d["Direct_event_junction_match"].fillna(False).astype(bool).any()
            if (not d.empty and "Direct_event_junction_match" in d.columns)
            else False
        )

        compact_rows.append({
            "CellLine": c["CellLine"],
            "Peptide": pep,
            "Protein_ID": p3i,
            "STRG_transcript_id": c["STRG_transcript_id"],
            "allele": c.get("allele", ""),
            "Binding_Class": c.get("Binding_Class", ""),
            "rep_count": c.get("rep_count", ""),
            "CorePriorityClass": c["CorePriorityClass"],
            "n_linked_event_rows": len(d),
            "Any_linked_event_is_novelJunction_in_candidate_comparison": novel_own,
            "Any_shared_exact_event_is_novelJunction_in_both_runs": novel_both,
            "Direct_peptide_to_event_junction_coordinate_match": direct,
            "FinalStructuralResolution":
                ";".join(final_classes),
        })

    compact = pd.DataFrame(compact_rows)
    compact.to_csv(
        OUTDIR / "07_FINAL_candidate_summary.tsv",
        sep="\t", index=False
    )

    # 9) Overall interpretation summary
    n_core = len(core)
    n_exact = int((core["CorePriorityClass"] == "A_EXACT_SHARED_SAME_DIRECTION").sum())
    n_gene = int((core["CorePriorityClass"] == "B_SHARED_GENE_EVENTTYPE").sum())

    n_novel_own = int(
        compact[
            "Any_linked_event_is_novelJunction_in_candidate_comparison"
        ].fillna(False).astype(bool).sum()
    )

    n_novel_both = int(
        compact[
            "Any_shared_exact_event_is_novelJunction_in_both_runs"
        ].fillna(False).astype(bool).sum()
    )

    n_direct = int(
        compact[
            "Direct_peptide_to_event_junction_coordinate_match"
        ].fillna(False).astype(bool).sum()
    )

    n_with_genomic_exons = sum(
        1
        for _, c in core.iterrows()
        if exon_by_cell_tx.get(
            (str(c["CellLine"]), str(c["STRG_transcript_id"]))
        )
    )

    summary = pd.DataFrame([
        ["Core HLA candidates", n_core],
        ["Exact shared same-direction candidates", n_exact],
        ["Shared gene+event-type candidates", n_gene],
        ["Candidates with linked rMATS novelJunction in own comparison", n_novel_own],
        ["Exact-shared candidates with novelJunction in both rMATS runs", n_novel_both],
        ["Core STRG transcripts with genomic exon structure found", n_with_genomic_exons],
        ["Candidates with direct peptide-to-event-junction coordinate match", n_direct],
    ], columns=["Metric", "Value"])

    summary.to_csv(
        OUTDIR / "08_FINAL_INTERPRETATION_SUMMARY.tsv",
        sep="\t", index=False
    )

    # 10) README
    no_exon = n_with_genomic_exons == 0

    readme = f"""
PAPER 3 — DEEP PEPTIDE-TO-SPLICE MAPPING

Core candidates selected:
  total = {n_core}
  exact shared same-direction = {n_exact}
  shared gene+event-type = {n_gene}

Novel-junction annotation:
  candidates with linked novelJunction in own comparison = {n_novel_own}
  exact-shared candidates novelJunction in both runs = {n_novel_both}

Genomic STRG exon structure:
  core transcripts with usable genomic exon structure = {n_with_genomic_exons}

Direct peptide-to-event junction:
  candidates with direct coordinate match = {n_direct}

Interpretation rules
--------------------
DIRECT_SHARED_NOVEL_JUNCTION_PEPTIDE
  Strongest structural category: the peptide spans an exon junction whose
  coordinates match the linked shared rMATS event, and that event is annotated
  as novelJunction in both comparisons.

DIRECT_NOVEL_JUNCTION_PEPTIDE
  Peptide spans the linked event junction and the event is novelJunction in
  the candidate cell-line comparison.

DIRECT_SPLICE_EVENT_JUNCTION_PEPTIDE_NOT_NOVEL_ANNOTATED
  Peptide spans the linked rMATS event junction, but the event was not flagged
  by fromGTF.novelJunction.

UNRESOLVED_LOCUS_ORF_LEVEL_ONLY
  We have peptide -> STRG protein/ORF -> transcript locus -> significant rMATS
  event, but no genomic exon-by-exon STRG structure was available to prove the
  peptide crosses the event junction.

Important:
  This analysis supports HPV-associated splicing, not HPV causation.
  HeLa, SiHa and C33A are non-isogenic.

  A significant rMATS event at the same locus is NOT, by itself, proof that the
  HLA peptide was generated by that splice event.

Genomic exon structure missing entirely: {no_exon}
""".strip()

    (OUTDIR / "README_DEEP_MAPPING.txt").write_text(
        readme, encoding="utf-8"
    )

    print("=" * 78)
    print("PAPER 3 — DEEP PEPTIDE-TO-SPLICE MAPPING (V2 CELL-SPECIFIC GTF)")
    print("=" * 78)
    print()
    print("CORE CANDIDATES")
    print(compact.to_string(index=False))
    print()
    print("FINAL SUMMARY")
    print(summary.to_string(index=False))
    print()
    print("Results:")
    print(OUTDIR)
    print("=" * 78)


if __name__ == "__main__":
    main()
