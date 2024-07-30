#########################
#
# Part 1: Run EMMAX for SNP and PAV
#
#########################


#-----------------------
#
#Step 1: Transform Phenoytpe data
#Step 2: Prepare VCF (filter 231 samples)
#        both: SNP & PAV
#Step 3: Prepare for EMMAX (P/Q matrix)
#        both: SNP & PAV
#Step 4: Run EMMAX
#        both: SNP & PAV
#Step 5: Visualize (Manhattan plot)
#
#-----------------------


#-----------------------
#Step 1: Transform Phenoytpe data
#-----------------------


#---------------------- in R -------------------
##GWAS check for Shapori-Wilk test
#packages
library("dplyr")
library("ggpubr")

getwd()
setwd("/Users/yueyu/Desktop/Jun_GWAS/Raw_Data_Jun_2024March")
mydata <- read.delim("RAW_Phenotype_DICorr_2024March_withNA.txt", sep = '\t', header = T, check.names = FALSE)

#Noted as ND
mydata_nona <- mydata[mydata$DIcorr != "ND",]
na_samples <- mydata[mydata$DIcorr == "ND",]
#SAVE NA samples names
write.table(na_samples, file = "NA_26sample_names.txt", sep = "\t", row.names = F, col.names=F, quote = FALSE)


mydata_nona$DIcorr <- as.numeric(mydata_nona$DIcorr)
hist(mydata_nona$DIcorr)
#Smaller breaks
hist(mydata_nona$DIcorr,breaks = 30)

ggdensity(mydata_nona$DIcorr, 
          main = "Density plot of DIcorr",
          xlab = "sample_name")

ggqqplot(mydata_nona$DIcorr)
shapiro.test(mydata_nona$DIcorr)
#p-value = 6.622e-07, smaller than 0.05, reject null hypo (same as normal)
#this means that the current data is not normally distributed 


#---boxcox transform
#install.packages("EnvStats")
library("EnvStats")

mydata_nona[mydata_nona ==0] <- 0.0001
#SIGNIFICANT 0 CHANGED TO 0.0001 a very small number

DIcorr_bc <- mydata_nona$DIcorr;str(DIcorr_bc)

bc <- boxcox(DIcorr_bc,optimize = TRUE)
bc$lambda #lambda is 0.6806267

#added col with labmda transform
mydata_nona$boxcox <- mydata_nona$DIcorr^(bc$lambda)
hist(mydata_nona$boxcox)
ggqqplot(mydata_nona$boxcox)
shapiro.test(mydata_nona$boxcox)
#W = 0.98113, p-value = 0.003464

#ADD sample name to match VCF file names
mydata_nona$Sample <- gsub("1-","",mydata_nona$Number)

#Save all
write.table(mydata_nona, file = "FULL.txt", sep = "\t", row.names = F, col.names=T, quote = FALSE)
#Save df
write.table(mydata_nona[,c(4,4,3)], file = "Phenotype_DICorr_2024March_232lines_noNA_boxcox.txt", sep = "\t", row.names = F, col.names=F, quote = FALSE)
#Save sample names
write.table(mydata_nona[,4], file = "SampleNames_232lines_2024Mar.txt", sep = "\t", row.names = F, col.names=F, quote = FALSE)

#---------------------- in R (END) -------------------








#-----------------------
#Step 2: Prepare VCF file
#-----------------------


#####################
# Step 2.1 : for SNPs
#####################

screen -r vcf

cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data


#-----------------
#
# PCA + KINSHIP WITHOUT (SAM_289) >> Thus only 231 lines
#
#-----------------

#Phenotype data 

#Use SNP PCA and Kinship matrix
#subset 231 samples that have germinated
bcftools view /moonriseNFS/VCF_HA412/Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99.vcf.gz -S SampleNames_231lines_NoSAM289_2024Mar.txt | less > sam.pop.231.vcf

bgzip sam.pop.231.vcf
tabix sam.pop.231.vcf.gz


#find a list of sample names in the 
bcftools query -l sam.pop.231.vcf.gz | wc -l

#vcf >> plink text 
plink --vcf sam.pop.231.vcf.gz --set-missing-var-ids @:# --maf 0.03 --allow-extra-chr --recode12 --output-missing-genotype 0 --transpose --out sam.pop.231.plink

rm sam.pop.231.plink.log sam.pop.231.plink.nosex 

#Edit the ID colimn in tped and map files
sed -i 's/^Ha412HOChr//g' sam.pop.231.plink.tped
sed -i 's/^Ha412HOChr//g' sam.pop.231.plink.map

#LD Pruning with PLINK  >> for PCA 
#(this $NAME here is for the name before all the .map/.tped stuff): sam.pop.232.plink
plink --tfile sam.pop.231.plink --allow-extra-chr --indep-pairphase 50kb 50 0.2 --out LD_pruned_231_lines

#PCA calculation with plink
plink --tfile sam.pop.231.plink --allow-extra-chr --extract LD_pruned_231_lines.prune.in --pca 3 --out PCA

#organize into .cov file, adding 1 into data? Why?
awk '{print $1,"\t",$2,"\t",1,"\t",$3,"\t",$4,"\t",$5}' PCA.eigenvec > PCA.cov


#-------INFO: EMMAX software saved at ~/software/EMMAX
#kinship calculation with EMMAX
#EMMAX kinship input is the prefix of the .tped files
~/software/EMMAX/emmax-kin-intel64 -v -s -d 10 sam.pop.231.plink


#$PHENOTYPE_FILE #see EMMAX wiki page
#FINAL EMMAX RUN
~/software/EMMAX/emmax-intel64 -v -d 10 -t sam.pop.231.plink -p Phenotype_DICorr_2024March_231lines_noNA_noSAM289_boxcox.txt -k sam.pop.231.plink.aIBS.kinf -c PCA.cov -o /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/EMMAX_output_SNP_20240325




#####################
# Step 2.2 : for PAVs
#####################


#Run EMMAX 
cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data

~/software/EMMAX/emmax-intel64 -v -d 10 -t /SundanceScratch/yueyu_SD/1_PhD/GWAS_jun/PAV/SAM.5maf.selected231.raw -p Phenotype_DICorr_2024March_231lines_noNA_noSAM289_boxcox.txt -k sam.pop.231.plink.aIBS.kinf -c PCA.cov -o /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/EMMAX_output_PAV_20240325

rm EMMAX_output_PAV_20240325.log EMMAX_output_PAV_20240325.reml

wc -l EMMAX_output_PAV_20240325.ps
#7,541,946 PAVs



#END
