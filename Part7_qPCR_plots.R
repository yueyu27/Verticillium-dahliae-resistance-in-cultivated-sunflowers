
#########################
#
# Part 7: Plot q-PCR results
#
#########################



########################
#####
#####
#####   Summarized: Line graphs
#####
#####
#########################



# ----- for all five lines used in the qPCR experiemnt 

cd /Users/yueyu/Desktop/0_ManuscriptPre_July18/qPCR_plot

R

library(ggplot2)
library(dplyr)
library(reshape2)

# Table with all q_PCR raw results, per tested GENE
a <- read.csv("ALL_input_Chr08g0327051.txt", sep = '\t', header = T, check.names = F)  

a <- a[,1:5]
dim(a)
head(a)
lapply(a,class)


#--------------------
#  IMMUNE
#--------------------

I <- a %>% filter(Type == "I")

head(I)
dim(I)

# Transform the data to long format
data_long <- melt(I[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

data_long$Days <- as.numeric(gsub("d", "", data_long$Days))

# Calculate summary statistics
summary_data <- data_long %>%
  group_by(Days, Group) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Lower = quantile(Value, 0.25, na.rm = TRUE),
    Upper = quantile(Value, 0.75, na.rm = TRUE)
  )

plot <- ggplot(data = summary_data, aes(x = Days, y = Median, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(x = "Treatment Time (days)", y = "Gene expression level") +
  scale_color_manual(values = c("Control" = "grey60", "Treatment" = "lightsteelblue3")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits =c(0,10))  + # Adjust the 'by' parameter as needed
  scale_x_continuous(breaks = c(0,1,3,5,7), limits =c(-0.2,7.3))  +
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black"),
        legend.text = element_text(size = 13), # Increase legend text size
        legend.title = element_text(size = 15) # Increase legend title size
        )

ggsave("Summerizeall5samples_IMMUNE_lineplot_2024July29.png", plot = plot, width = 8, height = 5, dpi = 200)



#--------------------
#  High Resistant
#--------------------

R <- a %>% filter(Type == "R")

head(R)
dim(R)

# Transform the data to long format
data_long <- melt(R[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")


data_long$Days <- as.numeric(gsub("d", "", data_long$Days))

# Calculate summary statistics
summary_data <- data_long %>%
  group_by(Days, Group) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Lower = quantile(Value, 0.25, na.rm = TRUE),
    Upper = quantile(Value, 0.75, na.rm = TRUE)
  )

plot <- ggplot(data = summary_data, aes(x = Days, y = Median, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(x = "Treatment Time (days)", y = "Gene expression level") +
  scale_color_manual(values = c("Control" = "grey60", "Treatment" = "lightsalmon2")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits =c(0,10))  + # Adjust the 'by' parameter as needed
  scale_x_continuous(breaks = c(0,1,3,5,7), limits =c(-0.2,7.3))  +
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black"),
        legend.text = element_text(size = 13), # Increase legend text size
        legend.title = element_text(size = 15) # Increase legend title size
        )


ggsave("Summerizeall5samples_HighR_lineplot_2024July29.png", plot = plot, width = 8, height = 5, dpi = 200)





#--------------------
#  Medium Susceptible
#--------------------

S <- a %>% filter(Type == "S")

head(S)
dim(S)

# Transform the data to long format
data_long <- melt(S[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

data_long$Days <- as.numeric(gsub("d", "", data_long$Days))

# Calculate summary statistics
summary_data <- data_long %>%
  group_by(Days, Group) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Lower = quantile(Value, 0.25, na.rm = TRUE),
    Upper = quantile(Value, 0.75, na.rm = TRUE)
  )

plot <- ggplot(data = summary_data, aes(x = Days, y = Median, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(x = "Treatment Time (days)", y = "Gene expression level") +
  scale_color_manual(values = c("Control" = "grey60", "Treatment" = "lightgoldenrod")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits =c(0,10))  + # Adjust the 'by' parameter as needed
  scale_x_continuous(breaks = c(0,1,3,5,7), limits =c(-0.2,7.3))  +
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black"),
        legend.text = element_text(size = 13), # Increase legend text size
        legend.title = element_text(size = 15) # Increase legend title size
        )


ggsave("Summerizeall5samples_MedSuscep_lineplot_2024July29.png", plot = plot, width = 8, height = 5, dpi = 200)






########################
#####
#####
#####   Individual: bar plots
#####
#####
#########################


R

library(ggplot2)
library(dplyr)
library(reshape2)

a <- read.csv("ALL_input_Chr08g0327051.txt", sep = '\t', header = T, check.names = F)  
a <- a[,1:5]

a$Days <- as.numeric(gsub("d", "", a$Days))

a$Days <- as.factor(a$Days)

dim(a)
head(a)
lapply(a,class)



#--------------------
#  IMMUNE SAMPLE
#--------------------

test <- a %>% filter(Type == "I") %>% filter(SAM == 177)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightsteelblue3")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM177.png", plot = plot, width = 8, height = 5, dpi = 200)



#--------------------
#  IMMUNE SAMPLE
#--------------------

test <- a %>% filter(Type == "I") %>% filter(SAM == 136)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightsteelblue3")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM136.png", plot = plot, width = 8, height = 5, dpi = 200)





#--------------------
#  HR SAMPLE
#--------------------

test <- a %>% filter(Type == "R") %>% filter(SAM == 14)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightsalmon2")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM014.png", plot = plot, width = 8, height = 5, dpi = 200)



#--------------------
#  HR SAMPLE
#--------------------

test <- a %>% filter(Type == "R") %>% filter(SAM == 10)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightsalmon2")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM010.png", plot = plot, width = 8, height = 5, dpi = 200)



#--------------------
#  S SAMPLE
#--------------------

test <- a %>% filter(Type == "S") %>% filter(SAM == 85)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightgoldenrod")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM086.png", plot = plot, width = 8, height = 5, dpi = 200)






#--------------------
#  S SAMPLE
#--------------------

test <- a %>% filter(Type == "S") %>% filter(SAM == 175)

head(test)
dim(test)

# Transform the data to long format
data_long <- melt(test[,3:5], id.vars = "Days", variable.name = "Group", value.name = "Value")

plot <- ggplot(data_long, aes(x = Days, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Control" = "grey", "Treatment" = "lightgoldenrod")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14), limits =c(0,14))  + 
  theme_classic() +
  theme(
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(margin = margin(t = 13), colour = "black"), 
        axis.title.y = element_text(margin = margin(r = 13), colour = "black")
        )

ggsave("I_SAM175.png", plot = plot, width = 8, height = 5, dpi = 200)





#END