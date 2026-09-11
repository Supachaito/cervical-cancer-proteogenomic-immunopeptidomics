#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PAPER 3 — FINAL RESOLVER FOR STRG.7884.3 / rMATS RUN PROVENANCE

Purpose
-------
A) Recover any existing reference/classification mapping for STRG.7884.3
   from high-value Paper 3 rescue / annotation files.
B) Locate GffCompare-like outputs (.tmap/.tracking/.refmap/gffcompare/gffcmp).
C) Recover the original rMATS run provenance:
   --gtf, --b1, --b2, and sample order, including PowerShell/WSL history when possible.

READ ONLY:
  PAPER3_ROOT
  PAPER3_WES_ROOT

WRITE:
  PAPER3_OUTPUT_ROOT/TARGET_STRG7884_3_FINAL_RESOLVER
"""

from __future__ import annotations

import csv
import os
import re
import subprocess
from pathlib import Path

TARGET = "STRG.7884.3"

WORKSPACE = Path(os.environ.get("PAPER3_ROOT", Path.cwd())).expanduser().resolve()
HROOT = Path(os.environ.get("PAPER3_WES_ROOT", WORKSPACE / "external" / "WES_Project_Analysis")).expanduser().resolve()
OUTPUT_ROOT = Path(os.environ.get("PAPER3_OUTPUT_ROOT", WORKSPACE / "results")).expanduser().resolve()
OUTDIR = OUTPUT_ROOT / "TARGET_STRG7884_3_FINAL_RESOLVER"
ENABLE_SHELL_HISTORY = os.environ.get("PAPER3_ENABLE_SHELL_HISTORY", "0").strip().lower() in {"1", "true", "yes"}
OUTDIR.mkdir(parents=True, exist_ok=True)

HIGH_VALUE_DIRS = [
    HROOT / "05_Paper3_Rerun" / "_REFMAP_STRG_Rescue_WSL",
    HROOT / "05_Paper3_Rerun" / "_TRUE_STRG_AnnotationSearch_WSL",
    HROOT / "05_Paper3_Rerun" / "_TRUE_STRG_MAPPING_ONLY_WSL",
    HROOT / "05_Paper3_Rerun" / "_REAL_STRG_MAPPING_SEARCH_V2_WSL",
    HROOT / "05_Paper3_Rerun" / "_STRG_MappingSourceInspection_WSL",
    HROOT / "05_Paper3_Rerun" / "Personalized_DB_System" / "01_Novel_Splicing_DB",
    HROOT / "Project_Paper_3_Databases" / "Personalized_DB_Rebuild" / "01_STRG_Splicing",
    HROOT / "05_Paper3_Rerun" / "Figure2_event_flow",
]

RMATS_DIR = HROOT / "04_Splicing_Analysis"

TARGET_EXTS = {
    ".tsv", ".csv", ".txt", ".gtf", ".gff", ".gff3",
    ".tmap", ".tracking", ".refmap", ".loci"
}

RUN_EXTS = {
    ".py", ".ps1", ".sh", ".bash", ".bat", ".cmd",
    ".log", ".txt", ".md", ".yaml", ".yml", ".json"
}

GFFCOMPARE_NAME_RE = re.compile(
    r"(gffcompare|gffcmp|\.tmap$|\.tracking$|\.refmap$|annotated.*\.gtf$|merged.*\.gtf$)",
    re.I,
)

RMATS_TOKEN_RE = re.compile(
    r"(rmats|--gtf|--b1|--b2|b1\.txt|b2\.txt|hela|siha|c33a)",
    re.I,
)


def write_tsv(path: Path, rows: list[dict], columns: list[str]):
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def safe_size(path: Path):
    try:
        return path.stat().st_size
    except OSError:
        return None


def read_text(path: Path, max_bytes: int | None = None):
    try:
        if max_bytes is not None:
            size = safe_size(path)
            if size is not None and size > max_bytes:
                return None
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None


def exact_target_hits(root: Path):
    hits = []
    if not root.exists():
        return hits

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in TARGET_EXTS:
            continue

        size = safe_size(p)
        if size is not None and size > 250 * 1024 * 1024:
            continue

        text = read_text(p)
        if text is None or TARGET not in text:
            continue

        lines = text.splitlines()
        header = lines[0] if lines else ""

        for n, line in enumerate(lines, start=1):
            if TARGET not in line:
                continue
            hits.append({
                "File": str(p),
                "LineNumber": n,
                "Header": header[:20000],
                "Line": line[:50000],
            })

    return hits


def find_gffcompare_like_files(root: Path):
    rows = []
    if not root.exists():
        return rows

    for p in root.rglob("*"):
        if not p.is_file():
            continue

        low = p.name.lower()
        if (
            GFFCOMPARE_NAME_RE.search(low)
            or p.suffix.lower() in {".tmap", ".tracking", ".refmap"}
        ):
            rows.append({
                "File": str(p),
                "SizeBytes": safe_size(p) or "",
            })

    return rows


def target_hits_in_files(file_rows):
    hits = []
    for r in file_rows:
        p = Path(r["File"])
        text = read_text(p, max_bytes=500 * 1024 * 1024)
        if text is None or TARGET not in text:
            continue

        lines = text.splitlines()
        header = lines[0] if lines else ""

        for n, line in enumerate(lines, start=1):
            if TARGET in line:
                hits.append({
                    "File": str(p),
                    "LineNumber": n,
                    "Header": header[:20000],
                    "Line": line[:50000],
                })
    return hits


def find_named_b1_b2(root: Path):
    rows = []
    if not root.exists():
        return rows

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if p.name.lower() not in {"b1.txt", "b2.txt"}:
            continue

        text = read_text(p, max_bytes=5 * 1024 * 1024)
        rows.append({
            "File": str(p),
            "Name": p.name,
            "Content": (text or "")[:50000],
        })
    return rows


def scan_rmats_run_files(root: Path):
    rows = []
    if not root.exists():
        return rows

    for p in root.rglob("*"):
        if not p.is_file():
            continue

        low = p.name.lower()

        if low.startswith("fromgtf.") or ".mats." in low:
            continue

        if p.suffix.lower() not in RUN_EXTS:
            continue

        text = read_text(p, max_bytes=20 * 1024 * 1024)
        if not text or not RMATS_TOKEN_RE.search(text):
            continue

        matched = []
        lines = text.splitlines()

        for i, line in enumerate(lines):
            l = line.lower()
            if (
                "--gtf" in l
                or "--b1" in l
                or "--b2" in l
                or "b1.txt" in l
                or "b2.txt" in l
                or "rmats.py" in l
                or "rmats " in l
            ):
                for j in range(max(0, i - 3), min(len(lines), i + 4)):
                    matched.append(lines[j])

        if not matched:
            continue

        seen = set()
        context = []
        for x in matched:
            if x not in seen:
                seen.add(x)
                context.append(x)

        rows.append({
            "File": str(p),
            "Context": "\n".join(context[:400]),
        })

    return rows


def powershell_history_clues():
    appdata = os.environ.get("APPDATA", "")
    candidates = []
    if appdata:
        candidates.append(
            Path(appdata)
            / "Microsoft"
            / "Windows"
            / "PowerShell"
            / "PSReadLine"
            / "ConsoleHost_history.txt"
        )

    rows = []

    for p in candidates:
        if not p.exists():
            continue
        text = read_text(p, max_bytes=100 * 1024 * 1024)
        if not text:
            continue

        lines = text.splitlines()
        matched = []

        for i, line in enumerate(lines):
            l = line.lower()
            if (
                "rmats" in l
                or "--gtf" in l
                or "--b1" in l
                or "--b2" in l
            ):
                for j in range(max(0, i - 3), min(len(lines), i + 4)):
                    matched.append(lines[j])

        if matched:
            rows.append({
                "Source": str(p),
                "Context": "\n".join(matched[-1000:]),
            })

    return rows


def wsl_history_clues():
    rows = []
    shell_cmd = (
        'for f in ~/.bash_history /root/.bash_history /home/*/.bash_history; do '
        '[ -f "$f" ] || continue; '
        'echo "### HISTORY_FILE=$f"; '
        'grep -Ein -C 3 "rmats|--gtf|--b1|--b2|b1\\.txt|b2\\.txt" "$f" 2>/dev/null | tail -n 1500; '
        'done'
    )

    try:
        proc = subprocess.run(
            ["wsl.exe", "bash", "-lc", shell_cmd],
            capture_output=True,
            text=True,
            timeout=60,
            errors="ignore",
        )
        out = (proc.stdout or "") + "\n" + (proc.stderr or "")
        if out.strip():
            rows.append({
                "Source": "WSL bash histories",
                "Context": out[:200000],
            })
    except Exception as e:
        rows.append({
            "Source": "WSL bash histories",
            "Context": f"WSL history scan unavailable: {type(e).__name__}: {e}",
        })

    return rows


print("=" * 90)
print("PAPER 3 — FINAL RESOLVER")
print("=" * 90)
print("Target:", TARGET)
print()

target_hits = []
for root in HIGH_VALUE_DIRS:
    print("Scanning target in:", root)
    target_hits.extend(exact_target_hits(root))

uniq = {}
for r in target_hits:
    uniq[(r["File"].lower(), r["LineNumber"])] = r
target_hits = list(uniq.values())

write_tsv(
    OUTDIR / "01_STRG7884_3_HIGH_VALUE_HITS.tsv",
    target_hits,
    ["File", "LineNumber", "Header", "Line"],
)

gff_files = []
for root in [
    HROOT / "05_Paper3_Rerun",
    HROOT / "Project_Paper_3_Databases",
    HROOT / "04_Splicing_Analysis",
]:
    print("Inventory GffCompare-like files:", root)
    gff_files.extend(find_gffcompare_like_files(root))

uniq = {}
for r in gff_files:
    uniq[r["File"].lower()] = r
gff_files = list(uniq.values())

write_tsv(
    OUTDIR / "02_GFFCOMPARE_LIKE_FILE_INVENTORY.tsv",
    gff_files,
    ["File", "SizeBytes"],
)

gff_target_hits = target_hits_in_files(gff_files)

write_tsv(
    OUTDIR / "03_STRG7884_3_GFFCOMPARE_HITS.tsv",
    gff_target_hits,
    ["File", "LineNumber", "Header", "Line"],
)

print("Searching b1/b2 files...")
b_files = find_named_b1_b2(RMATS_DIR)

write_tsv(
    OUTDIR / "04_RMATS_B1_B2_FILES.tsv",
    b_files,
    ["File", "Name", "Content"],
)

print("Searching rMATS run scripts/logs...")
rmats_run = scan_rmats_run_files(RMATS_DIR)

write_tsv(
    OUTDIR / "05_RMATS_RUN_PROVENANCE.tsv",
    rmats_run,
    ["File", "Context"],
)

print("Searching PowerShell history..." if ENABLE_SHELL_HISTORY else "Skipping PowerShell history (disabled by default).")
ps_hist = powershell_history_clues() if ENABLE_SHELL_HISTORY else []

write_tsv(
    OUTDIR / "06_POWERSHELL_RMATS_HISTORY.tsv",
    ps_hist,
    ["Source", "Context"],
)

print("Searching WSL bash history..." if ENABLE_SHELL_HISTORY else "Skipping WSL history (disabled by default).")
wsl_hist = wsl_history_clues() if ENABLE_SHELL_HISTORY else []

write_tsv(
    OUTDIR / "07_WSL_RMATS_HISTORY.tsv",
    wsl_hist,
    ["Source", "Context"],
)

summary = [
    "PAPER 3 — FINAL RESOLVER SUMMARY",
    "",
    f"Target transcript: {TARGET}",
    f"High-value exact target hits: {len(target_hits)}",
    f"GffCompare-like files inventoried: {len(gff_files)}",
    f"Target hits inside GffCompare-like files: {len(gff_target_hits)}",
    f"b1/b2 files recovered: {len(b_files)}",
    f"rMATS run-provenance files with command clues: {len(rmats_run)}",
    f"PowerShell history sources with rMATS clues: {len(ps_hist)}",
    f"WSL history result blocks: {len(wsl_hist)}",
    "",
    "INTERPRETATION:",
]

if gff_target_hits:
    summary += [
        "Direct STRG.7884.3 hits were found in GffCompare-like output(s).",
        "Inspect 03_STRG7884_3_GFFCOMPARE_HITS.tsv first.",
    ]
elif target_hits:
    summary += [
        "STRG.7884.3 was found in high-value rescue/annotation files.",
        "Inspect 01_STRG7884_3_HIGH_VALUE_HITS.tsv; header and exact row are preserved together.",
    ]
else:
    summary += [
        "No high-value rescue/classification hit was recovered for STRG.7884.3.",
    ]

if b_files or rmats_run or ps_hist:
    summary += [
        "Potential rMATS group-order provenance was recovered.",
        "Inspect outputs 04–06 for --b1/--b2/--gtf and sample-list contents.",
    ]
else:
    summary += [
        "No explicit Windows/project rMATS group-order provenance was recovered.",
        "Inspect 07_WSL_RMATS_HISTORY.tsv for WSL command-history evidence.",
    ]

summary += [
    "",
    "LOCKED BIOLOGICAL RESULT (unchanged):",
    "AKSDVTNQQY is a WT-absent HLA-I ligand that directly spans a significant alternative-splice event.",
    "",
    "Do not upgrade to 'novel-splicing-derived' until reference novelty is resolved.",
]

(OUTDIR / "08_FINAL_RESOLVER_SUMMARY.txt").write_text(
    "\n".join(summary),
    encoding="utf-8",
)

print()
print("=" * 90)
print("DONE")
print("=" * 90)
print("\n".join(summary))
print()
print("Results:")
print(OUTDIR)
