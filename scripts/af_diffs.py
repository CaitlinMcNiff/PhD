#!/usr/bin/env python3
"""
Pairs EUR/EAS variants with same phenotype & same gene, computes AF difference, and
outputs pairs with |AF_EUR - AF_EAS| > 0.5.

Produces two files:
1) CGC-filtered:
   cancer_type, gene, chr, start, AF_EUR, AF_EAS, AF_diff
2) Unfiltered (no CGC constraint):
   cancer_type, gene, chr, start, AF_EUR, AF_EAS, AF_diff
"""

import re
import sys
from pathlib import Path
import pandas as pd

# ---File Paths & Config---------
BASE = Path("/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/variant_selection")
GLOB = "*_cancer_closest_genes.bed"
OUTFILE = BASE / "both_pops_same_gene_pheno_AFdiff.tsv"             # CGC-filtered
OUTFILE_ALL = BASE / "both_pops_same_gene_pheno_AFdiff_noCGC.tsv"   # no CGC
AF_COL_1BASED = 16                                                  # AF is column 16 (1-based)
CGC_PATH = BASE / "CGC_EUR_EAS_overlap_genes.txt"
GTf_ATTR_REGEX = re.compile(r'gene_name\s+"([^"]+)"')
# -------------------------------

# fixed indices for the bed-like file
IDX = {"chr": 0, "start": 1, "end": 2, "rsid": 3, "phenotype": 4}
AF_IDX0 = AF_COL_1BASED - 1  # 0-based

def map_cancer_type(phenotype_raw: str) -> str:
    if not isinstance(phenotype_raw, str):
        return "Other/General"
    p = phenotype_raw.lower()
    if "prostate" in p: return "Prostate"
    if "breast" in p: return "Breast"
    if any(k in p for k in ["colorectal","colon","rectal","rectum","crc"]): return "Colorectal"
    if any(k in p for k in ["ovarian","ovary"]): return "Ovarian"
    if any(k in p for k in ["nsclc","non-small cell lung","small-cell lung","sclc","lung","bronchus"]): return "Lung"
    if any(k in p for k in ["skin","keratinocyte","non-melanoma","melanoma"]): return "Skin"
    if any(k in p for k in ["gastric","stomach","cardia"]): return "Gastric"
    if "pancrea" in p: return "Pancreatic"
    if any(k in p for k in ["cervical","cervix"]): return "Cervical"
    if "endometrial" in p: return "Endometrial"
    if any(k in p for k in ["uterine","leiomyoma"]): return "Uterine"
    if "bladder" in p: return "Bladder"
    if any(k in p for k in ["kidney","renal","renal pelvis"]): return "Kidney"
    if "thyroid" in p: return "Thyroid"
    if any(k in p for k in ["esophageal","oesophageal","esophagus","oesophagus"]): return "Esophageal"
    if any(k in p for k in ["testicular","germ cell"]): return "Testicular"
    if any(k in p for k in ["liver","hepatic","intrahepatic bile duct","bile duct","gallbladder","cholangiocarcinoma"]): return "Hepatobiliary"
    if any(k in p for k in ["brain","nervous system","glioma","glioblastoma"]): return "Brain/CNS"
    if any(k in p for k in ["oropharynx","hypopharynx","larynx","pharynx","oral cavity","tongue","head and neck","salivary"]): return "Head & Neck"
    if p.strip() in {"cancer"} or "cancer (" in p or "multiple cancers" in p or "hormone-sensitive cancer" in p:
        return "Other/General"
    return "Other/General"

def infer_population_from_name(path: Path) -> str:
    name = path.name.lower()
    if "eur" in name: return "EUR"
    if "eas" in name: return "EAS"
    return "UNK"

def try_find_gene_column(df: pd.DataFrame):
    lower_names = {str(c).lower(): i for i, c in enumerate(df.columns)}
    for key in ["gene", "gene_symbol", "genesymbol", "symbol"]:
        if key in lower_names:
            return lower_names[key]
    return None

def read_cgc_symbols(path: Path) -> set:
    df = pd.read_csv(path, sep=None, engine="python")
    col = "GeneSymbol" if "GeneSymbol" in df.columns else df.columns[0]
    syms = df[col].astype(str).str.strip().str.strip('"').str.upper()
    return set(s for s in syms if s)

def load_one(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, engine="python")
    ncol = df.shape[1]
    core = pd.DataFrame({
        "chr":   df.iloc[:, IDX["chr"]].astype(str),
        "start": df.iloc[:, IDX["start"]].astype(str),
        "end":   df.iloc[:, IDX["end"]].astype(str),
        "rsid":  df.iloc[:, IDX["rsid"]].fillna(""),
        "phenotype": df.iloc[:, IDX["phenotype"]].fillna("").astype(str).str.strip(),
    })

    if AF_IDX0 >= ncol:
        raise ValueError(f"{path.name}: AF column index {AF_IDX0} out of bounds for {ncol} columns")
    core["AF"] = pd.to_numeric(df.iloc[:, AF_IDX0], errors="coerce")

    gene_col_idx = try_find_gene_column(df)
    if gene_col_idx is not None:
        gene = df.iloc[:, gene_col_idx].astype(str)
    else:
        gene = pd.Series([""] * len(df), dtype=str)
        picked_idx = None
        for j in range(min(ncol, 40)):
            s = df.iloc[:, j].astype(str)
            if s.str.contains(r'gene_name\s+"[^"]+"', regex=True).any():
                picked_idx = j
                g = s.str.extract(GTf_ATTR_REGEX, expand=False)
                gene = g.fillna("").astype(str)
                break
        if picked_idx is None:
            sys.stderr.write(f"[WARN] {path.name}: no gene column found; gene will be empty.\n")

    core["gene"] = gene.str.strip().str.strip('"').str.upper().fillna("")
    core["population"] = infer_population_from_name(path)
    core["phenotype_norm"] = core["phenotype"].str.lower()
    core["cancer_type"] = core["phenotype_norm"].apply(map_cancer_type)
    core["var_key"] = core["rsid"].where(core["rsid"].ne("") & core["rsid"].notna(),
                                         core["chr"] + ":" + core["start"] + "-" + core["end"])
    return core

def format_output(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Turn a merged (EUR/EAS) DF into:
    cancer_type, gene, chr, start, AF_EUR, AF_EAS, AF_diff
    """
    out = pd.DataFrame({
        "cancer_type": merged["cancer_type_EUR"],
        "gene":        merged["gene"],
        "chr":         merged["chr_EUR"],
        "start":       merged["start_EUR"],
        "AF_EUR":      merged["AF_EUR"],
        "AF_EAS":      merged["AF_EAS"],
        "AF_diff":     merged["AF_diff"],
    })
    out = out.sort_values(by=["cancer_type","gene","AF_diff"], ascending=[True, True, False])
    return out

def main():
    # Load CGC
    if not Path(CGC_PATH).exists():
        print(f"[ERROR] CGC file not found: {CGC_PATH}", file=sys.stderr)
        sys.exit(1)
    cgc = read_cgc_symbols(CGC_PATH)
    if not cgc:
        print(f"[ERROR] No CGC symbols parsed from: {CGC_PATH}", file=sys.stderr)
        sys.exit(1)

    files = sorted(BASE.glob(GLOB))
    if not files:
        print(f"[ERROR] No files matched {BASE / GLOB}", file=sys.stderr)
        sys.exit(2)

    tables = []
    for f in files:
        try:
            t = load_one(f)
            tables.append(t)
        except Exception as e:
            print(f"[ERROR] Failed on {f.name}: {e}", file=sys.stderr)
            sys.exit(3)

    all_df = pd.concat(tables, ignore_index=True)
    all_df = all_df[all_df["population"].isin(["EUR","EAS"])].copy()

    # ==============================================================
    # 1) CGC-FILTERED PATH
    # ==============================================================
    cgc_df = all_df[all_df["gene"].isin(cgc)].copy()
    print(f"[INFO] CGC filter kept {len(cgc_df)}/{len(all_df)} rows.")

    cgc_df = cgc_df.drop_duplicates(subset=["population","var_key","gene","phenotype_norm"]).reset_index(drop=True)
    eur_cgc = cgc_df[cgc_df["population"]=="EUR"].copy()
    eas_cgc = cgc_df[cgc_df["population"]=="EAS"].copy()

    merged_cgc = eur_cgc.merge(
        eas_cgc,
        on=["var_key","gene","phenotype_norm"],
        how="inner",
        suffixes=("_EUR","_EAS")
    )
    merged_cgc["AF_diff"] = (merged_cgc["AF_EUR"] - merged_cgc["AF_EAS"]).abs()
    before = len(merged_cgc)
    merged_cgc = merged_cgc[merged_cgc["AF_diff"] > 0.5].copy()
    print(f"[INFO] (CGC) Kept {len(merged_cgc)}/{before} pairs with |AF_EUR - AF_EAS| > 0.5")

    out_cgc = format_output(merged_cgc)
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)
    out_cgc.to_csv(OUTFILE, sep="\t", index=False)
    print(f"[OK] Wrote CGC-filtered: {OUTFILE}")

    # ==============================================================
    # 2) UNFILTERED (NO CGC) PATH
    # ==============================================================
    all_df_nocgc = all_df.drop_duplicates(subset=["population","var_key","gene","phenotype_norm"]).reset_index(drop=True)
    eur_all = all_df_nocgc[all_df_nocgc["population"]=="EUR"].copy()
    eas_all = all_df_nocgc[all_df_nocgc["population"]=="EAS"].copy()

    merged_all = eur_all.merge(
        eas_all,
        on=["var_key","gene","phenotype_norm"],
        how="inner",
        suffixes=("_EUR","_EAS")
    )
    merged_all["AF_diff"] = (merged_all["AF_EUR"] - merged_all["AF_EAS"]).abs()
    before_all = len(merged_all)
    merged_all = merged_all[merged_all["AF_diff"] > 0.5].copy()
    print(f"[INFO] (noCGC) Kept {len(merged_all)}/{before_all} pairs with |AF_EUR - AF_EAS| > 0.5")

    out_all = format_output(merged_all)
    out_all.to_csv(OUTFILE_ALL, sep="\t", index=False)
    print(f"[OK] Wrote unfiltered: {OUTFILE_ALL}")

if __name__ == "__main__":
    main()
