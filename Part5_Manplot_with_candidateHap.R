

#########################
#
# Part 6: Plot Manhattan plot with candidate haploblocks above thershold
#
#########################


################
#  PAV
################

#Open R
R

#Set dir
#Load packages
library(data.table)
library(tidyverse)
library(ggplot2)


#--------------Checking this right now

#I used 7,541,946 PAVs
ben_threshold <- -log10(0.1/7541946); ben_threshold
ben_threshold.5 <- -log10(0.05/7541946); ben_threshold.5


a <- fread("EMMAX_output_PAV_20240325.ps", header =F)
colnames(a)[1] <- "PAV"
colnames(a)[4] <- "P"
dim(a)

a2 <- a[ ,c(1,4)]
rm(a)

#PAV use this
a3 <- a2 %>% separate(PAV,into = c("CHR","BP"), sep = "_" , remove = F)

# Change 01 02 03 ... to 1 2 3 
a3$CHR <- sprintf("%d", as.numeric(a3$CHR))

a3[,3] <- as.numeric(unlist(a3[,3])) 
a3[,2] <- as.numeric(unlist(a3[,2]))
#use UNLIST!!



#============================

#  ADD Selected region PAV

#============================

region <- fread("July_PAV/manplot_PAV_Hap.txt", header = T) 
#colnames(region) <- c("CHR","START","END")

region$CHR <- gsub("Ha412HOChr0","",region$CHR)
region$CHR <- gsub("Ha412HOChr","",region$CHR)
region$CHR <- as.numeric(region$CHR)

gwasResults <- a3

#Position of the Chr
don <- gwasResults %>%
  
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

head(don)
dim(don)

#define X axis
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#Transform P 
don$logP <- (-log10(don$P))



# Create a new column within_region and initialize it with "No"
don$within_region <- "NO"

# Loop through each row in the don data frame

r <- 1

for (r in 1:nrow(region)){

		CHR <- region$CHR[r]; CHR 
		START <- region$START[r]; START
		END <- region$END[r]; END

		T_row_idc <- which(don$CHR == CHR & don$BP >= START & don$BP <= END)

		#Col 8th is "within_region" for YES or NO
		don[T_row_idc,8] <- "YES"

		print(r)
 
  }




#Check number of YES and NOs
table(don$within_region)

     NO     YES
7540337    1609


# Assuming CHR is a character variable
don$point_color <- ifelse(don$within_region == "YES", "black", 
                    ifelse(as.numeric(don$CHR) %% 2 == 1, "cornflowerblue", "darkorange1"))

#Need to be a factor to be plotted!!
don$point_color <- factor(don$point_color,levels = c("cornflowerblue", "darkorange1", "black"))



#PLOT -logP > 2.5 points

ggplot(don, aes(x = BPcum, y = logP, color = point_color)) +
  geom_point(alpha = 0.8, size = 2.7) +
  scale_color_manual(values = c("cornflowerblue", "darkorange1", "black")) +
  geom_hline(yintercept = 4.404, linetype = "dashed", color = "black", linewidth = 1.6) +
  scale_x_continuous(labels = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(limits = c(2.5, NA)) +
  theme_bw() +
  theme(
    legend.position = "none",
   panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),  #axis text and title colour black
    axis.title = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(margin = margin(t = 15), colour = "black"), # Adjust t (top margin) as needed
    axis.title.y = element_text(margin = margin(r = 15), colour = "black")  # Adjust r (right margin) as needed
  ) +
  labs(
    y = "-log10 P-Value",
    x = "Chromosome")
  
ggsave("July_PAV/pav_manplot_above4.404_135PAV_p25_blue_orange_Publish_July23.png",width=22,height=8)






################
#  SNP
################

#Open R
R

#Set dir
#Load packages
library(data.table)
library(tidyverse)
library(ggplot2)


#--------------Checking this right now

#I used 3,699,248 SNPs, threshold is even bigger!
ben_threshold <- -log10(0.1/3699248); ben_threshold
ben_threshold.5 <- -log10(0.05/3699248); ben_threshold.5


a <- fread("EMMAX_output_SNP_20240325.ps", header =F)
colnames(a)[1] <- "SNP"
colnames(a)[4] <- "P"
dim(a)

a2 <- a[ ,c(1,4)]
rm(a)

#SNP use this
a3 <- a2 %>% separate(SNP,into = c("CHR","BP"), sep = ":" , remove = F)
a3[,3] <- as.numeric(unlist(a3[,3])) #use UNLIST!!

rm(a2)


#RUNNING HERE

#============================

#  ADD Selected region SNP

#============================

region <- fread("July_SNP/manplot_SNP_Hap.txt", header = T) 
#colnames(region) <- c("CHR","START","END")


region$CHR <- gsub("Ha412HOChr0","",region$CHR)
region$CHR <- gsub("Ha412HOChr","",region$CHR)
region$CHR <- as.numeric(region$CHR)




gwasResults <- a3

#Position of the Chr
don <- gwasResults %>%
  
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

head(don)
dim(don)

#Transform P 
don$logP <- (-log10(don$P))


# Create a new variable indicating whether each SNP is in a significant region


# Create a new column within_region and initialize it with "No"
don$within_region <- "NO"

# Loop through each row in the don data frame

r <- 1

for (r in 1:nrow(region)){

		CHR <- region$CHR[r]; CHR 
		START <- region$START[r]; START
		END <- region$END[r]; END

		T_row_idc <- which(don$CHR == region$CHR[r] & don$BP >= region$START[r] & don$BP <= region$END[r])

		#Col 8th is "within_region" for YES or NO
		don[T_row_idc,8] <- "YES"

		print(r)
 
  }


#Check number of YES and NOs
table(don$within_region)
#   NO     YES
#3699031     217


don$CHR <- gsub("Ha412HOChr0","",don$CHR)
don$CHR <- gsub("Ha412HOChr","",don$CHR)
don$CHR <- as.numeric(don$CHR)
class(don$CHR)
unique(don$CHR)

#define X axis (after CHR is numeric)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


# Assuming CHR is a character variable

don$point_color <- ifelse(don$within_region == "YES", "black", 
                    ifelse(as.numeric(don$CHR) %% 2 == 1, "cornflowerblue", "darkorange1"))

#Need to be a factor to be plotted!!
don$point_color <- factor(don$point_color,levels = c("cornflowerblue", "darkorange1", "black"))



#PLOT -logP > 2.5 points
ggplot(don, aes(x = BPcum, y = logP, color = point_color)) +
  geom_point(alpha = 0.8, size = 2.7) +                                               # Points bigger 
  scale_color_manual(values = c("cornflowerblue", "darkorange1", "black")) +
  geom_hline(yintercept = 4.232, linetype = "dashed", color = "black", linewidth = 1.6) + # line width not too thick
  scale_x_continuous(labels = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(limits = c(2.5, NA)) +
  theme_bw() +
  theme(
    legend.position = "none",
   panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),  #axis text and title colour black
    axis.title = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(margin = margin(t = 15), colour = "black"), # Adjust t (top margin) as needed
    axis.title.y = element_text(margin = margin(r = 15), colour = "black")  # Adjust r (right margin) as needed
  ) +
  labs(
    y = "-log10 P-Value",
    x = "Chromosome")
  
ggsave("July_SNP/snp_manplot_above4.232_60SNPs_p25_blue_orange_Publish_July23.png",width=22,height=8)

#END