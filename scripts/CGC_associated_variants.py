#!/usr/bin/env python3
# Finds EUR-unique and EAS-unique variants per shared phenotype, keeping only variants whose closest gene is in the CGC overlap list, and also groups phenotypes into broader cancer types for cross-phenotype merging.

import re
import pandas as pd
from pathlib import Path

# ---------- EDIT THESE ----------
BASE = Path("/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/variant_selection")
EUR_PATH = BASE / "EUR_cancer_closest_genes.bed"
EAS_PATH = BASE / "EAS_cancer_closest_genes.bed"
CGC_PATH = BASE / "CGC_EUR_EAS_overlap_genes.txt"
OUTDIR   = BASE
MAX_DIST = 10000   # set to None to skip distance filtering
# -------------------------------

def read_cgc_symbols(path: Path) -> set:
    df = pd.read_csv(path, sep=None, engine="python")  # header "GeneSymbol"
    col = "GeneSymbol" if "GeneSymbol" in df.columns else df.columns[0]
    syms = df[col].astype(str).str.strip().str.strip('"').str.upper()
    return set(s for s in syms if s != "")

# fixed column indices for your bedtools-closest outputs
IDX = {
    "chr": 0, "start": 1, "end": 2, "rsid": 3, "phenotype": 4,
    "gtf_attrs": 22, "dist": 23
}

gene_pat = re.compile(r'gene_name\s+"([^"]+)"')

def map_cancer_type(phenotype_raw: str) -> str:
    """Map free-text phenotype to a broad cancer type."""
    if not isinstance(phenotype_raw, str):
        return "Other/General"
    p = phenotype_raw.lower()

    # Order matters: check more specific/site terms before broad ones
    if "prostate" in p:
        return "Prostate"
    if "breast" in p:
        return "Breast"
    if any(k in p for k in ["colorectal", "colon", "rectal", "rectum", "crc"]):
        return "Colorectal"
    if any(k in p for k in ["ovarian", "ovary"]):
        return "Ovarian"
    if any(k in p for k in ["nsclc", "non-small cell lung", "small-cell lung", "sclc", "lung", "bronchus"]):
        return "Lung"
    if any(k in p for k in ["skin", "keratinocyte", "non-melanoma", "melanoma"]):
        return "Skin"
    if any(k in p for k in ["gastric", "stomach", "cardia"]):
        return "Gastric"
    if "pancrea" in p:
        return "Pancreatic"
    if any(k in p for k in ["cervical", "cervix"]):
        return "Cervical"
    if "endometrial" in p:
        return "Endometrial"
    if any(k in p for k in ["uterine", "leiomyoma"]):
        return "Uterine"
    if "bladder" in p:
        return "Bladder"
    if any(k in p for k in ["kidney", "renal", "renal pelvis"]):
        return "Kidney"
    if "thyroid" in p:
        return "Thyroid"
    if any(k in p for k in ["esophageal", "oesophageal", "esophagus", "oesophagus"]):
        return "Esophageal"
    if any(k in p for k in ["testicular", "germ cell"]):
        return "Testicular"
    if any(k in p for k in ["liver", "hepatic", "intrahepatic bile duct", "bile duct", "gallbladder", "cholangiocarcinoma"]):
        return "Hepatobiliary"
    if any(k in p for k in ["brain", "nervous system", "glioma", "glioblastoma"]):
        return "Brain/CNS"
    if any(k in p for k in ["oropharynx", "hypopharynx", "larynx", "pharynx", "oral cavity", "tongue", "head and neck", "salivary"]):
        return "Head & Neck"

    # Broad/unspecific buckets
    if p.strip() in {"cancer"} or "cancer (" in p or "multiple cancers" in p or "hormone-sensitive cancer" in p:
        return "Other/General"

    return "Other/General"

def load_pop_table(path: Path, pop: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, engine="python")
    # extract gene_name from attributes
    attrs = df.iloc[:, IDX["gtf_attrs"]].fillna("")
    gene = attrs.str.extract(gene_pat, expand=False).fillna("").str.upper()

    out = pd.DataFrame({
        "chr":   df.iloc[:, IDX["chr"]],
        "start": df.iloc[:, IDX["start"]],
        "end":   df.iloc[:, IDX["end"]],
        "rsid":  df.iloc[:, IDX["rsid"]].fillna(""),
        "phenotype": df.iloc[:, IDX["phenotype"]].fillna("").str.strip(),
        "gene":  gene,
    })

    # distance (numeric if present)
    if IDX["dist"] < df.shape[1]:
        out["dist"] = pd.to_numeric(df.iloc[:, IDX["dist"]], errors="coerce")
    else:
        out["dist"] = pd.NA

    # stable variant key (prefer rsID)
    out["var_key"] = out["rsid"].mask(out["rsid"].eq("") | out["rsid"].isna(),
                                      out["chr"] + ":" + out["start"] + "-" + out["end"])
    out["population"] = pop

    # normalized helpers
    out["phenotype_norm"] = out["phenotype"].str.lower()
    out["cancer_type"] = out["phenotype_norm"].apply(map_cancer_type)

    return out

def unique_by_pheno(dfA, dfB, pop):
    bmap = dfB.groupby("phenotype_norm")["var_key"].apply(set).to_dict()
    parts = []
    for ph, sub in dfA.groupby("phenotype_norm"):
        others = bmap.get(ph, set())
        parts.append(sub[~sub["var_key"].isin(others)])
    out = pd.concat(parts, ignore_index=True) if parts else dfA.iloc[0:0].copy()
    out["population"] = pop  # overwrite or create safely
    return out

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cgc = read_cgc_symbols(CGC_PATH)

    # 1) Load population tables
    eur = load_pop_table(EUR_PATH, "EUR")
    eas = load_pop_table(EAS_PATH, "EAS")

    # 2) Keep only CGC-nearest genes you provided
    eur = eur[eur["gene"].isin(cgc)].copy()
    eas = eas[eas["gene"].isin(cgc)].copy()

    # 3) Optional distance cutoff
    if MAX_DIST is not None:
        if eur["dist"].notna().any(): eur = eur[eur["dist"].le(MAX_DIST).fillna(False)]
        if eas["dist"].notna().any(): eas = eas[eas["dist"].le(MAX_DIST).fillna(False)]

    # 4) Shared phenotypes + dedup (by phenotype)
    shared_pheno = sorted(set(eur["phenotype_norm"]).intersection(set(eas["phenotype_norm"])))
    eur = eur[eur["phenotype_norm"].isin(shared_pheno)].drop_duplicates(subset=["phenotype_norm","var_key"])
    eas = eas[eas["phenotype_norm"].isin(shared_pheno)].drop_duplicates(subset=["phenotype_norm","var_key"])

    # 5) Population-unique by phenotype
    eur_u = unique_by_pheno(eur, eas, "EUR")
    eas_u = unique_by_pheno(eas, eur, "EAS")

    # 6) Build EUR–EAS pair table on SAME phenotype & SAME CGC gene
    both_same_cgc = (
        eur[["phenotype","phenotype_norm","cancer_type","gene","rsid","chr","start","end","dist"]]
        .merge(
            eas[["phenotype","phenotype_norm","cancer_type","gene","rsid","chr","start","end","dist"]],
            on=["phenotype_norm","gene"],  # strict: same phenotype text (normalized) & same CGC gene
            suffixes=("_EUR","_EAS")
        )
    )

    # 7) NEW FILTER: drop rows where EUR and EAS coords are identical
    same_coords = (
        (both_same_cgc["chr_EUR"]   == both_same_cgc["chr_EAS"]) &
        (both_same_cgc["start_EUR"] == both_same_cgc["start_EAS"]) &
        (both_same_cgc["end_EUR"]   == both_same_cgc["end_EAS"])
    )
    dropped_n = int(same_coords.sum())
    both_same_cgc = both_same_cgc[~same_coords].copy()
    print(f"[INFO] Dropped {dropped_n} pairs with identical EUR/EAS coordinates.")

    # 8) Recompute compact summary AFTER filtering (phenotype level)
    triplet = (
        both_same_cgc.groupby(["phenotype_norm","gene"])
        .size()
        .reset_index(name="n_variant_pairs")
        .rename(columns={"phenotype_norm":"phenotype"})
    )

    # 9) ALSO: Merge by CANCER TYPE (cross-phenotype grouping)
    #    This allows different phenotype phrasings that map to the same site to align.
    eur_ct = eur.drop_duplicates(subset=["cancer_type","var_key"])
    eas_ct = eas.drop_duplicates(subset=["cancer_type","var_key"])

    both_same_cgc_ct = (
        eur_ct[["cancer_type","gene","rsid","chr","start","end","dist"]]
        .merge(
            eas_ct[["cancer_type","gene","rsid","chr","start","end","dist"]],
            on=["cancer_type","gene"],
            suffixes=("_EUR","_EAS")
        )
    )

    # Drop identical-coordinate pairs here too
    same_coords_ct = (
        (both_same_cgc_ct["chr_EUR"]   == both_same_cgc_ct["chr_EAS"]) &
        (both_same_cgc_ct["start_EUR"] == both_same_cgc_ct["start_EAS"]) &
        (both_same_cgc_ct["end_EUR"]   == both_same_cgc_ct["end_EAS"])
    )
    dropped_ct = int(same_coords_ct.sum())
    both_same_cgc_ct = both_same_cgc_ct[~same_coords_ct].copy()
    print(f"[INFO] (CancerType) Dropped {dropped_ct} pairs with identical EUR/EAS coordinates.")

    triplet_ct = (
        both_same_cgc_ct.groupby(["cancer_type","gene"])
        .size()
        .reset_index(name="n_variant_pairs")
    )

    # 10) Save
    cols = ["population","phenotype","cancer_type","gene","rsid","chr","start","end","dist"]
    eur_u[cols].to_csv(OUTDIR/"unique_EUR_by_phenotype.tsv", sep="\t", index=False)
    eas_u[cols].to_csv(OUTDIR/"unique_EAS_by_phenotype.tsv", sep="\t", index=False)

    # Original phenotype-level paired file (now includes cancer_type columns for convenience)
    both_same_cgc.to_csv(OUTDIR/"both_pops_same_cgc_by_phenotype.tsv", sep="\t", index=False)
    triplet.to_csv(OUTDIR/"triplet_overlap.tsv", sep="\t", index=False)

    # NEW cancer-type-level paired file
    both_same_cgc_ct.to_csv(OUTDIR/"both_pops_same_cgc_by_cancertype.tsv", sep="\t", index=False)
    triplet_ct.to_csv(OUTDIR/"triplet_overlap_by_cancertype.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
