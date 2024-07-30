

#########################
#
# Part 3: rrBLUP for predictive ability of top PAVs and SNPs
#
#########################



###################################

# rrBLUP + top PAVs

##################################

cd /SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/PAV

#--------Extract top 1000 PAV for input BLUP

R

library(rrBLUP)


#Load tped file of top 1000 PAVs

all <- read.delim("PAV/PAV_boxcox_top1000_name_subset_GT_2024Mar26.tped", check.names = FALSE, header = F, sep = "\t")
dim(all)
all <- all[,-seq(6,466, by = 2)]
all <- all[,-c(1,3,4)]

#Load names 
name <- read.delim("SampleNames_231lines_NoSAM289_2024Mar.txt", check.names = FALSE, header = F, sep = "\t")

colnames(all)[1] <- c("PAV")
colnames(all)[2:232] <- name[,1]

# Create a function to convert genotypes to {-1, 0, 1}
convert_genotype_PAV <- function(gt) {
  if (gt == "0/0") {
    return(-1)
  } else if (gt == "0/1" || gt == "1/0") {
    return(0)
  } else if (gt == "1/1") {
    return(1)
  } else {
    return(NA)
  }
}

# Apply the conversion function to each element of the genotype matrix
all[,2:232] <- apply(all[,2:232], c(1, 2), convert_genotype_PAV)

dim(all)
all[1:10,1:10]


#==================
#
# Loop PAV: 5-1000, step 5
#
#==================

#-------- Step 0: Load data
#Load Phenotype data:
#Load VD BOXCOX data
vd <- read.delim("Phenotype_DICorr_2024March_231lines_noNA_noSAM289_boxcox.txt", check.names = FALSE, header = F, sep = "\t", stringsAsFactors = T)
rownames(vd) <- vd[,1]
vd <- vd[,2:3]
head(vd)

#Load GENOTYPE data
dim(all)
all[1:10,1:10]


#-------- Step 1: Top 1000 SNPs, rank by P value
#Load RANKING data
rank <- read.delim("PAV/PAV_boxcox_top1000_name.txt", check.names = FALSE, header = F, sep = "\t", stringsAsFactors = F)
colnames(rank)[1] <- "PAV"
head(rank);dim(rank)

#--Set 70% for training
frac.train <- 0.7

#---LOOP START
#-------- Step 2: subset G set, P stay same
#-------- Step 3: test and validate set, save prediction accuracy
#-------- Step 3.1: Make new df to save
perm <- data.frame(matrix(nrow = length(seq(5, 1000, by = 5)), ncol = 500))
#-------- Step 3.2: seperate test and validate
#-------- Step 3.3: loop for 500 times, take mean save in new df

set.seed(777)

n <- 5
r <- 1

for (n in seq(5, 1000, by = 5)){
  nrow <- n/5
  print(nrow)
  need <- rank$PAV[1:n]
  subset_all <- all[all$PAV %in% need, ]
  all2 <- t(subset_all[,2:232])
  colnames(all2) <- subset_all$PAV

  A <- A.mat(all2, min.MAF=0.05, impute.method = "mean")

  n.lines <- round(frac.train*nrow(all2))

  for (r in 1:500){

    train <- sort(sample(1:nrow(all2), n.lines))
    y.pred <- as.matrix(vd[-train,2])

    vd2 <- vd
    vd2[-train,2] <- NA 
    minimal <- data.frame(y=vd2[,2],gid=rownames(all2))

    model <- kin.blup(minimal,K=A,geno="gid",pheno="y")

    perm[nrow,r] <- round(cor(model$g[-train],y.pred[,1]),2)

  }

rm(all2, A)
}

perm[1:20,1:10]
dim(perm)

#each column is permutation
#each row is steps from 5, 10, 15... to 1000

#Calculate mean
perm$average <- rowMeans(perm)
perm[1:20,498:501]

#Add col for SNP steps
perm$step <- seq(5,1000, by = 5)
perm[1:20,498:502]

write.table(perm, file = "PAV/PAV_boxcox_kinBLUP_231SAM_PAV5to1000step5_PERM500_2024Mar26.txt", sep = "\t", row.names = F, col.names= T, quote = F)



#------------- Plot the thershold graph: scatter plot

toplot <- perm[,c(502,501)]
dim(toplot)
head(toplot)


#------------- output to publish -----------

png(filename = "July_PAV/PAV_BOXCOX_thershold_2024July18_top130SNPs_Publish.png", width = 1200, height = 900, res = 200) 

ggplot(toplot, aes(x = step, y =  average)) +
  geom_point(size = 1.5) +
  #geom_vline(xintercept = 362, linetype = "dashed", color = "black", linewidth = 0.5) +  # do not plot, make it thin if needed
  geom_vline(xintercept = 135, linetype = "dashed", color = "black", linewidth = 1.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  labs( x = "Number of top ranked PAVs",
        y = "Predictive ability") +
  theme_classic() +
  theme(
        axis.text = element_text(size = 13,colour = "black"),
        axis.title = element_text(size = 16,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 5), colour = "black") 
        )


dev.off()














###################################

# rrBLUP + top SNPs

##################################

R

library(rrBLUP)

#Load tped file of top 1000 SNPs
all <- read.delim("SNP_boxcox_top1000_name_for_kinBLUP_231SAM_2024Mar.txt", check.names = FALSE, header = F, sep = "\t")
dim(all)
all <- all[,-c(2,3,235)]

#Load names 
name <- read.delim("/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/SampleNames_231lines_NoSAM289_2024Mar.txt", check.names = FALSE, header = F, sep = "\t")
colnames(all)[1] <- c("SNP")
colnames(all)[2:232] <- name[,1]

# Create a function to convert genotypes to {-1, 0, 1}
convert_genotype <- function(gt) {
  if (gt == "0|0") {
    return(-1)
  } else if (gt == "0|1" || gt == "1|0") {
    return(0)
  } else if (gt == "1|1") {
    return(1)
  } else {
    return(NA)
  }
}

# Apply the conversion function to each element of the genotype matrix
all[,2:232] <- apply(all[,2:232], c(1, 2), convert_genotype)


#==================
#
# Loop PAV: 5-1000, step 5
#
#==================

#-------- Step 0: Load data
#Load Phenotype data:
#Load VD BOXCOX data
vd <- read.delim("/SundanceScratch/yueyu_SD/3_PhD/GWAS_2024Mar/raw_data/Phenotype_DICorr_2024March_231lines_noNA_noSAM289_boxcox.txt", check.names = FALSE, header = F, sep = "\t", stringsAsFactors = T)
rownames(vd) <- vd[,1]
vd <- vd[,2:3]
head(vd)

#Load GENOTYPE data
dim(all)
all[1:10,1:10]



#-------- Step 1: Top 1000 SNPs, rank by P value
#Load RANKING data
rank <- read.delim("SNP_boxcox_top1000_name.txt", check.names = FALSE, header = F, sep = "\t", stringsAsFactors = F)
colnames(rank)[1] <- "SNP"
head(rank);dim(rank)

#--Set 70% for training
frac.train <- 0.7

#---LOOP START
#-------- Step 2: subset G set, P stay same
#-------- Step 3: test and validate set, save prediction accuracy
#-------- Step 3.1: Make new df for saving
perm <- data.frame(matrix(nrow = length(seq(5, 1000, by = 5)), ncol = 500))
#-------- Step 3.2: seperate test and validate
#-------- Step 3.3: loop for 500 times, take mean save in new df

set.seed(778)

n <- 5
r <- 1

for (n in seq(5, 1000, by = 5)){
  nrow <- n/5
  print(nrow)
  need <- rank$SNP[1:n]
  subset_all <- all[all$SNP %in% need, ]
  all2 <- t(subset_all[,2:232])
  colnames(all2) <- subset_all$SNP

  A <- A.mat(all2, min.MAF=0.05, impute.method = "mean")

  n.lines <- round(frac.train*nrow(all2))

  for (r in 1:500){

    train <- sort(sample(1:nrow(all2), n.lines))
    y.pred <- as.matrix(vd[-train,2])

    vd2 <- vd
    vd2[-train,2] <- NA 
    minimal <- data.frame(y=vd2[,2],gid=rownames(all2))

    model <- kin.blup(minimal,K=A,geno="gid",pheno="y")

    perm[nrow,r] <- round(cor(model$g[-train],y.pred[,1]),2)

  }

rm(all2, A)
}

perm[1:20,1:10]
dim(perm)


#each column is permutation
#each row is steps from 5, 10, 15... to 1000

#Calculate mean
perm$average <- rowMeans(perm)
perm[1:20,498:501]

#Add col for SNP steps
perm$step <- seq(5,1000, by = 5)
perm[1:20,498:502]

write.table(perm, file = "SNP_boxcox_kinBLUP_231SAM_PAV5to1000step5_PERM500_2024Mar27.txt", sep = "\t", row.names = F, col.names= T, quote = F)


toplot <- perm[,c(502,501)]
dim(toplot)
head(toplot)


#-------------- output to publish -----------

png(filename = "July_SNP/SNP_BOXCOX_thershold_2024July18_top60SNPs_Publish.png", width = 1200, height = 900, res = 200) 

ggplot(toplot, aes(x = step, y =  average)) +
  geom_point(size = 1.5) +
  #geom_vline(xintercept = 362, linetype = "dashed", color = "black", linewidth = 0.5) +  # do not plot, make it thin if needed
  geom_vline(xintercept = 60, linetype = "dashed", color = "black", linewidth = 1.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  labs( x = "Number of top ranked SNPs",
        y = "Predictive ability") +
  theme_classic() +
  theme(
        axis.text = element_text(size = 13,colour = "black"),
        axis.title = element_text(size = 16,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 5), colour = "black") 
        )


dev.off()



#END