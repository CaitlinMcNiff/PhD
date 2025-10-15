# Convert COSMIC CGC file to a 1-based BED-style file: Chr, Start, Stop, gene_symbol

import pandas as pd
import re

# --------- EDIT THESE ---------
INPUT_PATH  = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/variant_selection/Cosmic_CGC.tsv" 
OUTPUT_PATH = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/variant_selection/CGC_genes.bed"
# ------------------------------

df = pd.read_csv(INPUT_PATH, sep="\t", dtype=str)

rows = []
for _, row in df.iterrows():
    symbol = (row["GENE_SYMBOL"] or "").strip()
    locstr = (row["Genome Location"] or "").strip()
    tier   = (row["Tier"] or "").strip()
    if not symbol or not locstr:
        continue

    # Support multiple locations separated by comma/semicolon
    for part in re.split(r"[;,]", locstr):
        part = part.strip()
        m = re.match(r'(?:chr)?([0-9]{1,2}|X|Y)\s*:\s*(\d+)\s*-\s*(\d+)', part, flags=re.IGNORECASE)
        if not m:
            continue
        chrom, s, e = m.group(1).upper(), int(m.group(2)), int(m.group(3))

        rows.append((chrom, s, e, symbol, tier))


# Write BED-style (no header)
out = pd.DataFrame(rows, columns=["Chr", "Start", "Stop", "gene_symbol", "Tier"])
out.to_csv(OUTPUT_PATH, sep="\t", index=False, header=False)
print(f"[OK] Wrote {len(out)} rows to {OUTPUT_PATH}")
