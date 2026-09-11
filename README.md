# Cervical cancer proteogenomic immunopeptidomics

Code accompanying the manuscript:

**Proteogenomic immunopeptidomics uncovers altered-host HLA-I ligands from genomic variation and alternative splicing in cervical cancer**

## Scope

This repository contains a curated, path-sanitized set of custom scripts used during the Paper 3 analyses. The release focuses on WES-derived sequence-space construction, transcript/STRG candidate assessment, rMATS event linkage, peptide-to-splicing follow-up, and supporting provenance checks. Large primary data files are not distributed here.

The code was cleaned for public release by replacing workstation-specific absolute paths with command-line arguments or environment variables. Scientific thresholds and candidate-level logic were otherwise retained unless explicitly noted below.

## Data availability

- Public HLA-I immunopeptidomics: PRIDE/ProteomeXchange **PXD028738**
- In-house global proteomics: ProteomeXchange **PXD083029**

Additional RNA-seq/WES-derived inputs used by individual scripts are described by each script and are not bundled in this repository.

## Environment variables

Scripts that originated from the local analysis workspace now use these optional variables:

- `PAPER3_ROOT` — root of the Paper 3 working directory; defaults to the current working directory.
- `PAPER3_WES_ROOT` — external WES/rMATS analysis root.
- `PAPER3_OUTPUT_ROOT` — output root; defaults to `<PAPER3_ROOT>/results` or `results`, depending on the script.
- `PAPER3_SEARCH_ROOTS` — optional `os.pathsep`-separated directories for the VEP/immunopeptidome convenience search.
- `PAPER3_ENABLE_SHELL_HISTORY` — default `0`. Set to `1` only if reproducing the historical local provenance search; shell-history scanning is disabled in the public release by default.

See `config.example.env` for an example.

## Script map

### WES / variant branch

1. `scripts/wes/01_build_wes_variant_database.py`  
   Rebuilds cell-line-specific variant FASTA files and ID maps from the WES-derived mutant sequence inputs.

2. `scripts/wes/02_compare_vep_immunopeptidome.py`  
   Summarizes protein-altering VEP calls and compares mutated-gene sets with immunopeptidome-associated gene sets.

3. `scripts/wes/03_audit_wes_candidate_provenance.py`  
   Supporting audit for the four WES candidates recovered outside the integrated peptide-pair landscape.

### Transcript / splicing branch

4. `scripts/splicing/04_screen_strg_transcript_novelty.py`  
   Screens the 54 WT-absent STRG HLA-I candidates against StringTie transcript annotations.

5. `scripts/splicing/05_rmats_event_linkage.py`  
   Applies the manuscript rMATS event filters (`FDR <= 0.05`, `|DeltaPSI| >= 0.10`) and links significant events to the STRG candidate loci.

6. `scripts/splicing/06_map_peptides_to_splice_events_legacy.py`  
   Historical deep mapping used to resolve peptide/STRG/event relationships. **Important:** this script contains a historical direction-sensitive setting. The final manuscript does not infer biological direction from the sign of `IncLevelDifference`, because the original rMATS b1/b2 order could not be independently recovered. Direct structural reconstruction should be interpreted independently of those historical direction labels.

### Supporting provenance analyses

7. `scripts/supporting/07_resolve_commd6_rmats_provenance.py`  
   Searches recovered project artifacts for STRG.7884.3 and rMATS run-provenance evidence. Shell-history scanning is disabled by default in this public version.

8. `scripts/supporting/08_compare_commd6_strg_to_reference.py`  
   Compares the STRG.7884.3 transcript structure with reference transcript candidates and recovers rMATS group-order clues from project artifacts.

## Installation

Python 3.10+ is recommended.

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
# source .venv/bin/activate

pip install -r requirements.txt
```

## Reproducibility notes

This repository does not bundle raw MS files, FASTQ/BAM files, VEP outputs, assembled GTF files, rMATS output directories, or large intermediate search databases. Scripts expect those inputs to be supplied through documented arguments/environment variables.

The public immunopeptidomics spectra were generated independently of the in-house sequencing/global-proteome experiments. The repository therefore supports the computational reanalysis and candidate-provenance steps; it should not be interpreted as a prospectively matched multi-omic workflow.

## Release note

Before making the repository public, verify software authorship, choose an institutional/IP-compatible license, and create a tagged release (for example `v1.0.0`). That release can then be archived in Zenodo to obtain a DOI for the manuscript Code availability statement.
