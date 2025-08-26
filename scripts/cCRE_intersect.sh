#$ -wd  /exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/scripts
#$ -V
#$ -l h_rt=1:00:00   #runtime
#$ -l h_vmem=8G
#$ -l rl9=true

module load roslin/bedtools/2.31.1

bedtools intersect -a ../gwas_1000_genomes/AFR_all_gwas.bed -b ../preliminary_exploration/cCRE/encode_ccres_hg38.bed -wo > ../preliminary_exploration/cCRE/AFR_all_ccres.bed

bedtools intersect -a ../gwas_1000_genomes/AMR_all_gwas.bed -b ../preliminary_exploration/cCRE/encode_ccres_hg38.bed -wo > ../preliminary_exploration/cCRE/AMR_all_ccres.bed

bedtools intersect -a ../gwas_1000_genomes/EAS_all_gwas.bed -b ../preliminary_exploration/cCRE/encode_ccres_hg38.bed -wo > ../preliminary_exploration/cCRE/EAS_all_ccres.bed

bedtools intersect -a ../gwas_1000_genomes/EUR_all_gwas.bed -b ../preliminary_exploration/cCRE/encode_ccres_hg38.bed -wo > ../preliminary_exploration/cCRE/EUR_all_ccres.bed

bedtools intersect -a ../gwas_1000_genomes/SAS_all_gwas.bed -b ../preliminary_exploration/cCRE/encode_ccres_hg38.bed -wo > ../preliminary_exploration/cCRE/SAS_all_ccres.bed