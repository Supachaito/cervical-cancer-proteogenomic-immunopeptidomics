#!/usr/bin/env python3
"""Build cell-line-specific WES-derived variant protein FASTA databases.

The script expects three input FASTA files in --base-dir:
  C33A_MUTANTS_RAW_INDEX.fasta
  HELA_MUTANTS_RAW_INDEX.fasta
  SIHA_MUTANTS_RAW_INDEX.fasta

Outputs are written to --out-dir with stable Paper 3 variant identifiers and
an accompanying ID-to-original-header map.
"""
from __future__ import annotations
import argparse
from pathlib import Path

TARGET_FILES = {
    "C33A": "C33A_MUTANTS_RAW_INDEX.fasta",
    "HELA": "HELA_MUTANTS_RAW_INDEX.fasta",
    "SIHA": "SIHA_MUTANTS_RAW_INDEX.fasta",
}

def build_database(in_path: Path, out_fasta: Path, out_map: Path, cell: str) -> int:
    count = 0
    seq: list[str] = []
    header_original = ""
    with in_path.open('r', encoding='utf-8', errors='ignore') as f_in, \
         out_fasta.open('w', encoding='utf-8') as f_out, \
         out_map.open('w', encoding='utf-8') as m_out:
        for line in f_in:
            if line.startswith('>'):
                if seq:
                    count += 1
                    new_id = f"P3I_{cell}_MUT_{count:06d}"
                    f_out.write(f">{new_id} {header_original}\n{''.join(seq)}\n")
                    m_out.write(f"{new_id}\t{header_original}\n")
                header_original = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            count += 1
            new_id = f"P3I_{cell}_MUT_{count:06d}"
            f_out.write(f">{new_id} {header_original}\n{''.join(seq)}\n")
            m_out.write(f"{new_id}\t{header_original}\n")
    return count

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-dir', required=True, help='Directory containing the three *_MUTANTS_RAW_INDEX.fasta inputs.')
    ap.add_argument('--out-dir', default='', help='Output directory. Default: <base-dir>/Personalized_DB_Rebuild/02_Mutant_Variants')
    args = ap.parse_args()
    base_dir = Path(args.base_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else base_dir / 'Personalized_DB_Rebuild' / '02_Mutant_Variants'
    out_dir.mkdir(parents=True, exist_ok=True)
    for cell, file_name in TARGET_FILES.items():
        in_path = base_dir / file_name
        if not in_path.exists():
            raise FileNotFoundError(f"Required input not found: {in_path}")
        out_fasta = out_dir / f"{cell}_MUT_Personalized.fasta"
        out_map = out_dir / f"{cell}_MUT_Map.tsv"
        count = build_database(in_path, out_fasta, out_map, cell)
        print(f"{cell}: {count:,} sequences -> {out_fasta}")

if __name__ == '__main__':
    main()
