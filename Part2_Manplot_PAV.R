#########################
#
# Part 2: Plot Manhattan plot (PAV)
#
#########################

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

#Load PAV output
a <- fread("EMMAX_output_PAV_20240325.ps", header =F)
colnames(a)[1] <- "PAV"
colnames(a)[4] <- "P"
dim(a)

a2 <- a[ ,c(1,4)]

#PAV use this
a3 <- a2 %>% separate(PAV,into = c("CHR","BP"), sep = "_" , remove = F)
a3[,3] <- as.numeric(unlist(a3[,3])) #use UNLIST!!

rm(a2)

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


#define X axis
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don$logP <- (-log10(don$P))
dim(don)
head(don)

#Plot all points above 2.5
ggplot(don,aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size= 1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  geom_hline(yintercept=as.numeric(ben_threshold.5), linetype="dashed",color = "red", size=0.8) +
  geom_hline(yintercept=as.numeric(ben_threshold), linetype="dashed",color = "blue", size=0.8) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(limits = c(2.5, NA)) +
  
  # Custom the theme:
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank()
  )+
  labs(
    y = "-log10 P-Value",
    x = "PAV_Position")

#SAVE BOXCOX
ggsave("PAV/PAV.boxcox.ManPlot.2024March25.png",width=22,height=8)



#============================
#  Extract top 1000 and above thershold 
#============================
don2 <- don[order(don$logP,decreasing = T),]
head(don2)

don_top1000 <- don2[1:1000,]
write.table(don_top1000, file = "PAV/PAV_boxcox_top1000.txt", sep = "\t", row.names = F, col.names=T, quote = F)




#============================
#  Plot QQ-PLOT
#============================

# Extract the p-values
p_values <- a$P

# Calculate the expected p-values under the null hypothesis
expected_p_values <- -log10(ppoints(length(p_values)))
observed_p_values <- -log10(sort(p_values))

# Create a dataframe for plotting
qq_data <- data.frame(expected = expected_p_values, observed = observed_p_values)

# Generate the QQ plot
ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "QQ Plot",
       x = "Expected -log10(p-values)",
       y = "Observed -log10(p-values)") +
  theme_minimal()

ggsave("EMMAX_qqplot/PAV.boxcox.afterEMMAX.qqplot.png",width=22,height=8)


