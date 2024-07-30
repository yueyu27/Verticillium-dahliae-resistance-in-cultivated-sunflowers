
#########################
#
# Part 2: LD clumping (SNP + PAV) and Gene annotation
#
#########################



#==================
#
#      PAV
#
#==================

#----------- LD clumping for 135 PAVs

#ALL PAV BED FILE: /SundanceScratch/yueyu_SD/1_PhD/GWAS_jun/PAV/SAM.5maf.selected231.raw.tped

#Convert tped to bed file

#cd /SundanceScratch/yueyu_SD/1_PhD/GWAS_jun/PAV
#plink --tfile SAM.5maf.selected231.raw --make-bed --out SAM.5maf.selected231.raw.2024Mar27

#Perform LD clumping using the --clump option in PLINK. 
#This command will identify LD blocks around each 
#significant SNP and collapse overlapping blocks

# Update your PLINK command with the new .fam file
cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/July_PAV

BED="/SundanceScratch/yueyu_SD/1_PhD/GWAS_jun/PAV/SAM.5maf.selected231.raw.2024Mar27"
GENE="/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/HA412_Gene_list/USE_IN_PLINK.txt"
ASSOC="PAV_boxcox_above4.404_135PAVs.txt"

plink --bfile $BED \
      --clump $ASSOC \
      --clump-snp-field PAV \
      --clump-field P \
      --clump-kb 250 \
      --clump-r2 0.4 \
      --allow-extra-chr \
      --clump-range $GENE \
      --clump-range-border 5 \
      --out PAV_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene


#plink.clumped.ranges

     CHR      Chromosome code
     SNP      Index SNP per clump
     P        p-value 
     N        Number of clumped SNPs
     POS      Genomic co-ordinates
     KB       kb span of clumped SNPs 
     RANGES   List of ranges/genes that intersect the clumped region


wc -l PAV_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene.clumped.ranges
# Resulted in 107 HAP for PAV = 108 - header

#EXTRACT ALL GENES INTO NEW TXT
awk -F' ' '{print $7}' PAV_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene.clumped.ranges | sed 's/[][]//g' | tr ',' '\n' | sed '/^\s*$/d' > GENES_in_Top135_PAVs.txt

#34 genes 

# ----------- SAVE TO PLOT HAP OUT --------
awk -F' ' '{print $5}' PAV_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene.clumped.ranges | tail -n 107 > forexcel_manplot_PAV_Hap.txt




#==================
#
#      SNP
#
#==================

#cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data
#plink --vcf sam.pop.231.vcf.gz --set-missing-var-ids @:# --maf 0.03 --allow-extra-chr --make-bed --output-missing-genotype 0 --out SNP/sam.pop.231.VCF.plink.2024Mar

BED="/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/SNP/sam.pop.231.VCF.plink.2024Mar"
GENE="/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/HA412_Gene_list/HA412.gene.list.Just.Gene.2024July18.fromYue.txt"
ASSOC="SNP_boxcox_above4.232_60PAVs.txt"

plink --bfile $BED \
      --clump $ASSOC \
      --clump-snp-field SNP \
      --clump-field P \
      --clump-kb 250 \
      --clump-r2 0.4 \
      --allow-extra-chr \
      --clump-range $GENE \
      --clump-range-border 5 \
      --out SNP_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene



#=========== Haploblock conclusion ==========
# Resulted in 41 HAP for SNP = 42 - header
# Resulted in 107 HAP for PAV = 108 - header

# PAV + SNP = 41 + 107 = 148 Haploblocks



#EXTRACT ALL GENES INTO NEW TXT
awk -F' ' '{print $7}' SNP_LD_CLUMP_250kb_r0.4_flank5kbp_with_gene.clumped.ranges | sed 's/[][]//g' | tr ',' '\n' | sed '/^\s*$/d' > GENES_in_Top60_SNPs.txt

#19 genes 





#======================================
#
#     MERGE > Candidate Genes (51)
#
#======================================

cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/July_SNP

cat GENES_19_in_Top60_SNPs.txt /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/July_PAV/GENES_34_in_Top135_PAVs.txt | sort | uniq > MERGED_51_Unique_Genes_from_top60SNP_top135PAV_2024July18.txt

wc -l MERGED_51_Unique_Genes_from_top60SNP_top135PAV_2024July18.txt

# Merged Candidate Genes: 51 = 19 + 34 - 2 (2 overlaps both found in PAV: Ha412HOChr04g0149911 & Ha412HOChr04g0149921)

#Ha412HOChr04g0149911	4	12506690	12515912		PAV_Hap_11	chr4:12515100..12515100	0.001		PAV_Hap_12	chr4:12517700..12517700 0.001
#Ha412HOChr04g0149921	4	12519457	12521978		PAV_Hap_11	chr4:12515100..12515100	0.001		PAV_Hap_12	chr4:12517700..12517700 0.001




#==============================================
#
#     MERGE > Candidate Genes >> extract from gff3 file
#
#==============================================

GENE="/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/HA412_Gene_list/HA412.gene.list.Gene.PFAM.IPR.GO.2024July18.fromYue.nocontig.txt"
head $GENE    # Metadata with all gene information

head MERGED_51_Unique_Genes_from_top60SNP_top135PAV_2024July18.txt  #CANDIDATE 51 GENE NAMES

#------ Extract 51 candidate gene's info from the metadata file 

awk 'NR==FNR {names[$1]; next} $2 in names' MERGED_51_Unique_Genes_from_top60SNP_top135PAV_2024July18.txt $GENE  > MERGED_51_Genes_from_top60SNP_top135PAV_GENE_PFAM_IPR_GO_2024July18.txt

wc -l MERGED_51_Genes_from_top60SNP_top135PAV_GENE_PFAM_IPR_GO_2024July18.txt

# 51 unique genes with functions!! 

le MERGED_51_Genes_from_top60SNP_top135PAV_GENE_PFAM_IPR_GO_2024July18.txt


#END







