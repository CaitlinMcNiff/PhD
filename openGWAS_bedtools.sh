#$ -cwd  #current dir
#$ -V
#$ -l h_rt=2:00:00   #runtime
#$ -l h_vmem=12G
#$ -l rl9=true

# load bedtools module
module load roslin/bedtools/2.29.2

sort -k1,1 -k2,2n EUR_variants.vcf > EUR_variants_sorted.bed
bedtools intersect -a openGWAS_genome_wide_sig_sorted.bed -b EUR_variants_sorted.bed -sorted > EUR_openGWAS.bed

#for pop_vcf in *_variants.vcf; do
    
    #vcf_base=$(basename "$pop_vcf" _variants.vcf)
    # bedtools intersect for each of the superpopulations for all GWAS hits
    #oG_file=${vcf_base}_openGWAS.bed
    #bedtools intersect -a openGWAS_genome_wide_sig.bed -b "$pop_vcf" > "$oG_file"

#done