import os, re, csv, argparse, datetime

TAB = chr(9)
SAMPLES = ["C33A","SiHa","HeLa"]
PROTEIN_ALTERING = {
    "missense_variant","frameshift_variant","stop_gained","stop_lost","start_lost",
    "inframe_insertion","inframe_deletion","protein_altering_variant"
}

def norm_sample(x):
    if x is None: return None
    t = str(x).strip().lower()
    if "siha" in t: return "SiHa"
    if "hela" in t: return "HeLa"
    if "c33" in t:  return "C33A"
    return None

def norm_gene(g):
    if g is None: return None
    s = str(g).strip()
    if s == "" or s.upper() in {"NA","N/A","NONE","-","."}: return None
    # keep first token if multi-gene
    s = re.split(r"[;,|\s]+", s)[0]
    return s.upper()

def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=TAB)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

def parse_symbol_from_extra(extra):
    if not extra: return None
    m = re.search(r"(?:^|;)SYMBOL=([^;]+)", extra)
    return m.group(1) if m else None

def detect_vep_header(path):
    # return (header_cols, sep)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip():
                continue
            s = line.rstrip("\n")
            if ("Uploaded_variation" in s) and ("Consequence" in s):
                s2 = s.lstrip("#")
                if TAB in s2:
                    return [c.strip() for c in s2.split(TAB)], TAB
                else:
                    return [c.strip() for c in re.split(r"\s+", s2.strip())], None
    return None, None

def parse_vep(path, sample):
    if not os.path.exists(path):
        return {"sample":sample, "status":"MISSING_FILE", "total":0, "protein":0, "syn":0, "genes":set()}

    header, sep = detect_vep_header(path)
    if header is None:
        return {"sample":sample, "status":"NO_HEADER", "total":0, "protein":0, "syn":0, "genes":set()}

    idx = {c:i for i,c in enumerate(header)}
    idx_cons  = idx.get("Consequence", -1)
    idx_sym   = idx.get("SYMBOL", -1)
    idx_gene  = idx.get("Gene", -1)
    idx_extra = idx.get("Extra", -1)

    total = protein = syn = 0
    genes = set()

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            s = line.rstrip("\n")
            parts = s.split(TAB) if sep == TAB else re.split(r"\s+", s.strip())
            if idx_cons < 0 or len(parts) <= idx_cons:
                continue

            cons = parts[idx_cons]
            total += 1
            terms = cons.split("&") if cons else []

            if "synonymous_variant" in terms:
                syn += 1

            if any(t in PROTEIN_ALTERING for t in terms):
                protein += 1
                g = None
                if idx_sym >= 0 and len(parts) > idx_sym and parts[idx_sym].strip():
                    g = parts[idx_sym]
                elif idx_gene >= 0 and len(parts) > idx_gene and parts[idx_gene].strip():
                    g = parts[idx_gene]
                elif idx_extra >= 0 and len(parts) > idx_extra:
                    g = parse_symbol_from_extra(parts[idx_extra])
                g = norm_gene(g)
                if g:
                    genes.add(g)

    return {"sample":sample, "status":"OK", "total":total, "protein":protein, "syn":syn, "genes":genes}

def detect_sep(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline()
    return TAB if first.count(TAB) >= first.count(",") else ","

def choose_col(headers, exact, contains):
    hmap = {h.lower(): h for h in headers}
    for c in exact:
        if c in headers: return c
        if c.lower() in hmap: return hmap[c.lower()]
    for h in headers:
        lh = h.lower()
        if any(k in lh for k in contains):
            return h
    return None

def load_immuno_gene_sets(path):
    sep = detect_sep(path)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f, delimiter=sep)
        headers = r.fieldnames if r.fieldnames else []
        sample_col = choose_col(headers,
                                ["Sample","CellLine","Condition","Group","Source","SampleID","sample","cell_line"],
                                ["sample","cell","condition","group"])
        gene_col = choose_col(headers,
                              ["Gene","GeneSymbol","SYMBOL","GeneName","gene","gene_symbol","Genes"],
                              ["gene","symbol"])

        if gene_col is None:
            raise RuntimeError("NO_GENE_COLUMN in immunopeptidomics file")

        sets = {s:set() for s in SAMPLES}
        n = 0
        inferred = norm_sample(os.path.basename(path)) if sample_col is None else None

        for row in r:
            n += 1
            if sample_col is None:
                s = inferred
            else:
                s = norm_sample(row.get(sample_col))
            if s is None:
                continue
            g = norm_gene(row.get(gene_col))
            if g:
                sets[s].add(g)

    return sets, sample_col, gene_col, n

def find_immuno_candidates():
    # Optional convenience search. To avoid machine-specific paths, roots are
    # supplied through PAPER3_SEARCH_ROOTS (os.pathsep-separated). If unset,
    # only the current working directory is searched.
    env_roots = os.environ.get("PAPER3_SEARCH_ROOTS", "").strip()
    roots = [p for p in env_roots.split(os.pathsep) if p] if env_roots else [os.getcwd()]
    kws = ("immuno","hla","ligand","peptide","neoantigen","classi","classii")
    exts = (".tsv",".csv",".txt")
    out = []
    for root in roots:
        if not os.path.exists(root):
            continue
        root_depth = len(os.path.normpath(root).split(os.sep))
        for d, _, files in os.walk(root):
            depth = len(os.path.normpath(d).split(os.sep)) - root_depth
            if depth > 6:
                continue
            for fn in files:
                lf = fn.lower()
                if lf.endswith(exts) and any(k in lf for k in kws):
                    p = os.path.join(d, fn)
                    try:
                        sz = os.path.getsize(p)
                        mt = datetime.datetime.fromtimestamp(os.path.getmtime(p)).isoformat(sep=" ", timespec="seconds")
                    except:
                        sz, mt = 0, ""
                    out.append((p, sz, mt))
    # unique by path
    seen = set()
    uniq = []
    for p,sz,mt in out:
        if p in seen: 
            continue
        seen.add(p)
        uniq.append((p,sz,mt))
    return uniq[:60]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--immuno", default="")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 1) Parse VEP
    stats_rows = []
    gene_rows = []
    vep_sets = {}

    for s in SAMPLES:
        path = os.path.join(args.vep_dir, f"{s}.vep.txt")
        st = parse_vep(path, s)
        vep_sets[s] = st["genes"]
        stats_rows.append([s, st["total"], st["protein"], st["syn"], st["status"]])
        for g in sorted(st["genes"]):
            gene_rows.append([s, g])

    write_tsv(os.path.join(args.out_dir, "vep_protein_altering_stats.tsv"),
              ["Sample","TotalRows","ProteinAltering","Synonymous","Status"], stats_rows)

    write_tsv(os.path.join(args.out_dir, "vep_mutated_genes_by_sample.tsv"),
              ["Sample","Gene"], gene_rows)

    # 2) If immuno missing -> write candidates and stop
    immuno = args.immuno.strip()
    if immuno == "" or (not os.path.exists(immuno)):
        cands = find_immuno_candidates()
        cand_path = os.path.join(args.out_dir, "immuno_candidates.tsv")
        cand_rows = [[p,sz,mt] for p,sz,mt in cands]
        write_tsv(cand_path, ["Path","SizeBytes","ModifiedTime"], cand_rows)
        print("IMMUNO_NOT_SET")
        print("CANDIDATES=" + cand_path)
        print("VEP_STATS=" + os.path.join(args.out_dir, "vep_protein_altering_stats.tsv"))
        print("VEP_GENES=" + os.path.join(args.out_dir, "vep_mutated_genes_by_sample.tsv"))
        return

    # 3) Load immunopep + overlap
    imm_sets, sample_col, gene_col, nrows = load_immuno_gene_sets(immuno)

    im_rows = []
    for s in SAMPLES:
        for g in sorted(imm_sets[s]):
            im_rows.append([s,g])
    write_tsv(os.path.join(args.out_dir, "immunopep_detected_genes_by_sample.tsv"),
              ["Sample","Gene"], im_rows)

    sum_rows = []
    for s in SAMPLES:
        mut = vep_sets[s]
        imm = imm_sets[s]
        ov = mut & imm
        pct = (len(ov) / len(mut) * 100.0) if len(mut) > 0 else 0.0
        sum_rows.append([s, len(mut), len(imm), len(ov), round(pct,3)])

    mut_all = set().union(*[vep_sets[s] for s in SAMPLES])
    imm_all = set().union(*[imm_sets[s] for s in SAMPLES])
    ov_all = mut_all & imm_all
    pct_all = (len(ov_all) / len(mut_all) * 100.0) if len(mut_all) > 0 else 0.0
    sum_rows.append(["ALL_UNION", len(mut_all), len(imm_all), len(ov_all), round(pct_all,3)])

    write_tsv(os.path.join(args.out_dir, "gene_overlap_summary.tsv"),
              ["Sample","MutatedGenes(VEP)","DetectedGenes(Immunopep)","OverlapGenes","OverlapPct_of_Mutated"], sum_rows)

    hpv_mut = (vep_sets["SiHa"] | vep_sets["HeLa"]) - vep_sets["C33A"]
    hpv_imm = (imm_sets["SiHa"] | imm_sets["HeLa"])
    prio = sorted(hpv_mut & hpv_imm)

    pr_rows = []
    for g in prio:
        pr_rows.append([
            g,
            int(g in vep_sets["SiHa"]), int(g in vep_sets["HeLa"]), int(g in vep_sets["C33A"]),
            int(g in imm_sets["SiHa"]), int(g in imm_sets["HeLa"]), int(g in imm_sets["C33A"])
        ])

    write_tsv(os.path.join(args.out_dir, "HPVplus_priority_genes.tsv"),
              ["Gene","In_VEP_SiHa","In_VEP_HeLa","In_VEP_C33A","In_Immuno_SiHa","In_Immuno_HeLa","In_Immuno_C33A"], pr_rows)

    runlog = os.path.join(args.out_dir, "RUNLOG_compare_vep_immunopep_v2.txt")
    with open(runlog, "w", encoding="utf-8") as f:
        f.write(f"IMMUNO_FILE={immuno}\n")
        f.write(f"IMMUNO_ROWS={nrows}\n")
        f.write(f"IMMUNO_SAMPLE_COL={sample_col}\n")
        f.write(f"IMMUNO_GENE_COL={gene_col}\n")
        f.write("DONE\n")

    print("DONE")
    print("OUT_DIR=" + args.out_dir)
    print("SUMMARY=" + os.path.join(args.out_dir, "gene_overlap_summary.tsv"))
    print("PRIORITY=" + os.path.join(args.out_dir, "HPVplus_priority_genes.tsv"))
    print("RUNLOG=" + runlog)

if __name__ == "__main__":
    main()
