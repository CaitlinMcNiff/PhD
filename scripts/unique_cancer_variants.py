# Does the same thing as the CGC_associated_variants..py script, but does not limit the variants to those found in genes on the CGC gene list

import re
import pandas as pd
from pathlib import Path

# ---------- EDIT THESE ----------
BASE = Path("/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/variant_selection")
EUR_PATH = BASE / "EUR_cancer_closest_genes.bed"
EAS_PATH = BASE / "EAS_cancer_closest_genes.bed"
OUTDIR   = BASE
MAX_DIST = 10000   # set to None to skip distance filtering
# -------------------------------

IDX = {
    "chr": 0, "start": 1, "end": 2, "rsid": 3, "phenotype": 4,
    "gtf_attrs": 22, "dist": 23
}

gene_pat = re.compile(r'gene_name\s+"([^"]+)"')

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

def load_pop_table(path: Path, pop: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, engine="python")
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
    out["population"] = pop
    return out

def _variant_label(row: pd.Series) -> str:
    rsid = (row.get("rsid") or "").strip()
    if rsid:
        return rsid
    return f"{row.get('chr')}:{row.get('start')}-{row.get('end')}"

def _agg_unique_join(values) -> str:
    return ",".join(sorted({str(v) for v in values if pd.notna(v) and str(v) != ""}))

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # 1) Load population tables
    eur = load_pop_table(EUR_PATH, "EUR")
    eas = load_pop_table(EAS_PATH, "EAS")

    # Copies for population-specific derivation (no CGC filtering here)
    eur_all = eur.copy()
    eas_all = eas.copy()

    # (Optional) distance cutoff
    if MAX_DIST is not None:
        if eur_all["dist"].notna().any():
            eur_all = eur_all[eur_all["dist"].le(MAX_DIST).fillna(False)]
        if eas_all["dist"].notna().any():
            eas_all = eas_all[eas_all["dist"].le(MAX_DIST).fillna(False)]

    # Only compare phenotypes present in BOTH pops
    shared_pheno_all = sorted(set(eur_all["phenotype_norm"]).intersection(set(eas_all["phenotype_norm"])))
    eur_all = eur_all[eur_all["phenotype_norm"].isin(shared_pheno_all)].drop_duplicates(subset=["phenotype_norm", "var_key"])
    eas_all = eas_all[eas_all["phenotype_norm"].isin(shared_pheno_all)].drop_duplicates(subset=["phenotype_norm", "var_key"])

    # Population-specific variants by phenotype (NO CGC filter)
    eur_u_noCGC = unique_by_pheno(eur_all, eas_all, "EUR")
    eas_u_noCGC = unique_by_pheno(eas_all, eur_all, "EAS")

    # Save original unique-by-phenotype lists (unchanged)
    cols_noCGC = ["population","phenotype","cancer_type","gene","rsid","chr","start","end","dist"]
    eur_u_noCGC[cols_noCGC].to_csv(OUTDIR / "unique_EUR_by_phenotype_noCGC.tsv", sep="\t", index=False)
    eas_u_noCGC[cols_noCGC].to_csv(OUTDIR / "unique_EAS_by_phenotype_noCGC.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------------
    # Population-specific gene–cancer_type pairs
    # ---------------------------------------------------------------------
    eur_pairs_ps = (eur_u_noCGC.loc[eur_u_noCGC["gene"].ne(""), ["gene", "cancer_type"]]
                                .drop_duplicates()
                                .sort_values(["gene", "cancer_type"])
                                .reset_index(drop=True))
    eas_pairs_ps = (eas_u_noCGC.loc[eas_u_noCGC["gene"].ne(""), ["gene", "cancer_type"]]
                                .drop_duplicates()
                                .sort_values(["gene", "cancer_type"])
                                .reset_index(drop=True))

    # Intersection keys (pairs present in both pops' pop-specific sets)
    both_pairs_keys = (eur_pairs_ps.merge(eas_pairs_ps, on=["gene", "cancer_type"], how="inner")
                                   .sort_values(["gene", "cancer_type"])
                                   .reset_index(drop=True))

    # Save per-pop pair lists (unchanged)
    eur_pairs_ps.to_csv(OUTDIR / "gene_cancer_pairs_popSpecific_EUR.tsv", sep="\t", index=False)
    eas_pairs_ps.to_csv(OUTDIR / "gene_cancer_pairs_popSpecific_EAS.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------------
    # NEW: add VARIANT DATA to the "both" file
    # ---------------------------------------------------------------------
    # Prepare variant rows (population-specific only) for EUR/EAS, limited to pairs in both
    eur_ps = eur_u_noCGC.copy()
    eas_ps = eas_u_noCGC.copy()

    # Add a human-friendly variant label (rsid or chr:start-end)
    eur_ps["variant_label"] = eur_ps.apply(_variant_label, axis=1)
    eas_ps["variant_label"] = eas_ps.apply(_variant_label, axis=1)

    # Keep only entries whose (gene, cancer_type) are in the intersection
    eur_both = eur_ps.merge(both_pairs_keys, on=["gene", "cancer_type"], how="inner")
    eas_both = eas_ps.merge(both_pairs_keys, on=["gene", "cancer_type"], how="inner")

    # Long / tidy table (one row per variant)
    both_long = pd.concat([
        eur_both.assign(population="EUR")[["gene","cancer_type","population","rsid","chr","start","end","dist","var_key","variant_label"]],
        eas_both.assign(population="EAS")[["gene","cancer_type","population","rsid","chr","start","end","dist","var_key","variant_label"]],
    ], ignore_index=True).sort_values(["gene","cancer_type","population","variant_label"])
    both_long.to_csv(OUTDIR / "gene_cancer_pairs_popSpecific_both_long.tsv", sep="\t", index=False)

    # Aggregated / wide table: variant counts + lists per pop, plus phenotype summaries
    eur_agg = (eur_both.groupby(["gene","cancer_type"])
                      .agg(
                          EUR_variant_count=("var_key","nunique"),
                          EUR_variants=("variant_label", _agg_unique_join)
                      ).reset_index())

    eas_agg = (eas_both.groupby(["gene","cancer_type"])
                      .agg(
                          EAS_variant_count=("var_key","nunique"),
                          EAS_variants=("variant_label", _agg_unique_join)
                      ).reset_index())

    both_agg = (both_pairs_keys
                .merge(eur_agg, on=["gene","cancer_type"], how="left")
                .merge(eas_agg, on=["gene","cancer_type"], how="left")
                .sort_values(["gene","cancer_type"])
                .reset_index(drop=True))

    # Write the aggregated "both" file (now includes variant data)
    both_agg.to_csv(OUTDIR / "gene_cancer_pairs_popSpecific_both.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
