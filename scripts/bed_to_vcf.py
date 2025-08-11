import os
import glob

input_dir = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/gwas_1000_genomes"   
output_dir = "/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/preliminary_exploration/vep"
os.makedirs(output_dir, exist_ok=True)


vcf_header = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]

def process_bed_to_vcf(in_path, out_path):
    """
    Expects lines like:
    0:CHROM 1:START 2:END 3:ID 4:TRAIT 5:... 6:PVAL ... 13:REF 14:ALT ...
    Example:
    1  832873  832874  rs2977608  Major depressive disorder  0.259029  8E-6  1.0729614  31969693  1  832873  832874  .  A  C  0.71  1
    """
    seen = set()
    with open(out_path, "w") as out:
        for line in vcf_header:
            out.write(line + "\n")

        with open(in_path, "rt") as f:
            for raw in f:
                s = raw.strip()
                if not s or s.startswith("#"):
                    continue
                cols = s.split("\t")
                # require minimum columns for CHR, START, END, ID, REF, ALT
                if len(cols) < 15:
                    continue

                try:
                    chrom = cols[0]
                    start_0based = int(cols[1])    # BED start (0-based)
                    pos = start_0based + 1         # VCF POS (1-based)
                    rsid = cols[3] if cols[3] not in ("", ".") else "."
                    ref = cols[13]
                    alt = cols[14]

                    # Form the VCF line
                    vcf_line = f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t.\t."

                    # De-duplicate within this file
                    if vcf_line not in seen:
                        seen.add(vcf_line)
                        out.write(vcf_line + "\n")
                except Exception:
                    # Skip malformed rows gracefully
                    continue

def main():
    bed_paths = glob.glob(os.path.join(input_dir, "*_all_gwas.bed"))
    if not bed_paths:
        print(f"No files matching *_all.bed found in {input_dir}")
        return

    for bed in bed_paths:
        base = os.path.splitext(os.path.basename(bed))[0]
        out_vcf = os.path.join(output_dir, f"{base}.vcf")
        process_bed_to_vcf(bed, out_vcf)
        print(f"Wrote: {out_vcf}")

if __name__ == "__main__":
    main()
