#split GangSTR output
library(dplyr)
library(plyr)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(grid)

GangSTR_output_final <- mutate(GangSTR_output_final, A2 = gsub(":.*", "", GangSTR_output_final$V3))
GangSTR_output_final <- mutate(GangSTR_output_final, REPCI.1 = gsub(".*:", "", GangSTR_output_final$V3))
GangSTR_output_final <- mutate(GangSTR_output_final, GT = gsub(":.*", "", GangSTR_output_final$V2))
GangSTR_output_final <- mutate(GangSTR_output_final, DP = gsub("(:[^:]+):.*", "\\1", GangSTR_output_final$V2))
GangSTR_output_final <- mutate(GangSTR_output_final, DP = gsub(".*:", "", GangSTR_output_final$DP))
GangSTR_output_final <- mutate(GangSTR_output_final, A1 = gsub(".*:", "", GangSTR_output_final$V2))
GangSTR_output_final <- mutate(GangSTR_output_final, Q = gsub("^.*:(.*):.*$", "", GangSTR_output_final$V2))
GangSTR_output_final <- mutate(GangSTR_output_final, Q = gsub("(.+?)(\\.*:)", "\\1", GangSTR_output_final$V2))
GangSTR_output_final <- mutate(GangSTR_output_final, test = gsub("*.:(+[:^]:)", "\\1", GangSTR_output_final$V2))
GangSTR_output_final$test <- gsub("^.{0,5}", "", GangSTR_output_final$test)
GangSTR_output_final <- mutate(GangSTR_output_final, Q = gsub("(:[^:]+):.*", "\\1", GangSTR_output_final$test))
GangSTR_output_final <- mutate(GangSTR_output_final, Q = gsub(".*:", "", GangSTR_output_final$Q))

#merge csv as dataframes by LP numbers
GangSTR_output_final
PCR_vs_WGS <- HTT_PCR_vs_WGS

names(GangSTR_output_final)[1] <- "LP_number"
names(GangSTR_output_final)[2] <- "GangSTR_ouput_1"
names(GangSTR_output_final)[3] <- "GangSTR_ouput_2"
names(GangSTR_output_final)[4] <- "REPCI.2"
GangSTR_output_final
#df as merged dataframe
df_join <- join_all(list(GangSTR_output_final, PCR_vs_WGS), by="LP_number", type = "left")
df_join
#edit dataframe to keep only needed columns
df_complete = select(df_join, LP_number, min.PCR.a1, maxPCR.a2, A1, A2, min.EHv255.a1, max.EHv255.a2, min.EHv312.a1, max.EHv312.a2)
df_complete

#each allele its own row by method
df_analysis <- data.frame(df_complete[1], stack(df_complete[2:3]))
df_analysis1 <- data.frame(df_complete[1], stack(df_complete[4:5]))
df_analysis2 <- data.frame(df_complete[1], stack(df_complete[6:7]))
df_analysis3 <- data.frame(df_complete[1], stack(df_complete[8:9]))
names(df_analysis)[2] <- "PCR"
names(df_analysis1)[2] <- "GangSTR"
names(df_analysis2)[2] <- "EHv255"
names(df_analysis3)[2] <- "EHv312"
#select only LP once and then each method from the new data frames
df_analysis = select(df_analysis, LP_number, PCR)
df_analysis1 = select(df_analysis1, GangSTR)
df_analysis2 = select(df_analysis2, EHv255)
df_analysis3 = select(df_analysis3, EHv312)
#create index in a new column to left_join by index NOT LP
df_analysis$index <- seq.int(nrow(df_analysis))
df_analysis1$index <- seq.int(nrow(df_analysis1))
df_analysis2$index <- seq.int(nrow(df_analysis2))
df_analysis3$index <- seq.int(nrow(df_analysis3))
#left_join so each allele is stacked on top of each other in each method column
df_analysis_join <- left_join(df_analysis, df_analysis1, by = "index")
df_analysis_join1 <- left_join(df_analysis_join, df_analysis2, by = "index")
df_analysis_join2 <- left_join(df_analysis_join1, df_analysis3, by = "index")
df_analysis_join3 <- left_join(df_analysis_join2, df_analysis3, by = "index")
df_analysis_complete = select(df_analysis_join3, LP_number, PCR, GangSTR, EHv255, EHv312.x)

DF <- data.frame(df_analysis_complete)
write.csv(df_analysis_complete, "/mnt/ovd/slaveserver/o1611312021gBQgK/sharedFolder_5e06ef7d85ae976dbb6b59b83c35d909/shared_allGeCIPs/VG_DH_AT/Valentina_GangSTR/GangSTR/GangSTR_output//df_analysis_complete.csv")

#graphs to compare EH and GangSTR against PCR sizes
Gang <- DF %>%
  select(PCR, GangSTR) %>%
  gather(key = "variable", value = "value", -PCR)
head(Gang)
Gang$value <- as.numeric(as.character(Gang$value))
ggplot(data = Gang, aes(x=PCR, y=value))+
  geom_point(color="blue")+
  labs(subtitle = paste("Number of alleles:", nrow(Gang)))+
  xlab("PCR expansion size")+
  ylab("GangSTR prediction")+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  xlim(5, 125)+
  ylim(5, 125)+
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

LP_GangSTR <- Gang %>% select(PCR, value)
counts_GangSTR <- ddply(LP_GangSTR, .(LP_GangSTR$PCR, LP_GangSTR$value), nrow)
names(counts_GangSTR) <- c("PCR", "GangSTR", "Frequency")
counts_GangSTR$Less_frequency <- counts_GangSTR$Frequency/4
counts_complete_frequency_GangSTR <- counts_GangSTR$Less_frequency


ggplot(data = counts_GangSTR, aes(x=PCR, y=GangSTR))+
  geom_point(color="blue", size = counts_complete_frequency_GangSTR, shape=20)+
  labs(subtitle = paste("Number of alleles:", nrow(Gang)))+
  xlab("PCR expansion size")+
  ylab("GangSTR prediction")+
  xlim(5, 125)+
  ylim(5, 125)+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = df_analysis_complete, aes(x=PCR, y=EHv255))+
  geom_point(color="dark green")+
  labs(subtitle = paste("Number of alleles:", nrow(df_analysis_complete)))+
  xlab("PCR expansion size")+
  ylab("EHv255 prediction")+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  xlim(5, 125)+
  ylim(5, 125)+
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

counts_EH2 <- ddply(df_analysis_complete, .(df_analysis_complete$PCR, df_analysis_complete$EHv255), nrow)
names(counts_EH2) <- c("PCR", "EHv2", "Frequency")
counts_EH2$Less_frequency <- counts_EH2$Frequency/4
counts_complete_frequency_EH2 <- counts_EH2$Less_frequency

ggplot(data = counts_EH2, aes(x=PCR, y=EHv2))+
  geom_point(color="dark green", size = counts_complete_frequency_EH2, shape=20)+
  labs(subtitle = paste("Number of alleles:", nrow(df_analysis_complete)))+
  xlab("PCR expansion size")+
  ylab("EHv255 prediction")+
  xlim(5, 125)+
  ylim(5, 125)+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = df_analysis_complete, aes(x=PCR, y=EHv312.x))+
  geom_point(color="purple")+
  labs(subtitle = paste("Number of alleles:", nrow(df_analysis_complete)))+
  xlab("PCR expansion size")+
  ylab("EHv312 prediction")+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  xlim(5, 125)+
  ylim(5, 125)+
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

counts_EH3 <- ddply(df_analysis_complete, .(df_analysis_complete$PCR, df_analysis_complete$EHv312.x), nrow)
names(counts_EH3) <- c("PCR", "EHv3", "Frequency")
counts_EH3$Less_frequency <- counts_EH3$Frequency/4
counts_complete_frequency_EH3 <- counts_EH3$Less_frequency

ggplot(data = counts_EH3, aes(x=PCR, y=EHv3))+
  geom_point(color="purple", size = counts_complete_frequency_EH3, shape=20)+
  labs(subtitle = paste("Number of alleles:", nrow(df_analysis_complete)))+
  xlab("PCR expansion size")+
  ylab("EHv312 prediction")+
  xlim(5, 125)+
  ylim(5, 125)+
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  geom_vline(xintercept = 36, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 36, colour = "red", linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


to_plot = ggplot() +
  geom_point(data = AR, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles)) +
  xlim(5,55) +
  ylim(5,55) +
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +
  coord_equal() +
  geom_vline(xintercept = 34, colour = "red", linetype="dotted") +
  geom_hline(yintercept = 34, colour = "red", linetype="dotted") +
  labs(title = "",
       y = "EH sizes (#repeats)",
       x = "PCR sizes (#repeats)") +
  theme_light() +

#graph to compare EH and GangSTR against PCR sizes - changed discrete y axis to continuous data
test <- DF %>%
  select(PCR, GangSTR, EHv255, EHv312.x) %>%
  gather(key = "variable", value = "value", -PCR)
head(test)
test$value <- as.numeric(as.character(test$value))
ggplot(test, aes(x = PCR, y = value))+
  geom_point(aes(color = variable))+
  labs(subtitle = paste("Number of observations:", nrow(test)))+
  scale_color_manual(values = c("blue", "red", "green"))+
  xlab("PCR expansion size")+
  ylab("Predicted expansion size")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

HTT_PCR_vs_WGS1 <- HTT_PCR_vs_WGS %>% select(LP_number, min.EHv255.a1, max.EHv255.a2)
HTT_PCR_vs_WGS1 <- cbind(HTT_PCR_vs_WGS1[1], stack(HTT_PCR_vs_WGS1[2:3]))
names(HTT_PCR_vs_WGS1)[names(HTT_PCR_vs_WGS1) == 'values'] <- 'EHv2_size'
HTT_PCR_vs_WGS1 <- HTT_PCR_vs_WGS1 %>% mutate(id = row_number())
names(HTT_PCR_vs_WGS1)[names(HTT_PCR_vs_WGS1) == 'id'] <- 'X'
names(Allele_numeric_table)[names(Allele_numeric_table) == 'LP_number.x'] <- 'LP_number'
updateallele <- left_join(Allele_numeric_table, df_analysis_complete, by='X')
updateallele1 <- updateallele %>% select(X, LP_number.x, PCR_size, GangSTR, EHv255, EHv3_size)
updateallele1 <- left_join(updateallele1, HTT_PCR_vs_WGS1, by='X')
updateallele1 <- updateallele1 %>% select(LP_number, PCR_size, GangSTR, EHv2_size, EHv3_size)
updateallele_below60 <- updateallele %>% filter(PCR_size < 60, EHv3_size < 60)

#Prediction vs PCR visualisation
ggplot(data = updateallele1, aes(x=PCR_size, y=GangSTR)) +
  geom_count(colour = 'green', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,100) +
  ylim(0,100) +
  xlab('PCR expansion size') +
  ylab('GangSTR prediction') +
  ggtitle('Number of alleles = 154') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))
ggplot(data = updateallele1, aes(x=PCR_size, y=EHv2_size)) +
  geom_count(colour = 'purple', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +  
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,100) +
  ylim(0,100) +
  xlab('PCR expansion size') +
  ylab('EHv2 prediction') +
  ggtitle('Number of alleles = 154') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))
ggplot(data = updateallele1, aes(x=PCR_size, y=EHv3_size)) +
  geom_count(colour = 'blue', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,100) +
  ylim(0,100) +
  xlab('PCR expansion size') +
  ylab('EHv3 prediction') +
  ggtitle('Number of alleles = 154') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))

ggplot(data = updateallele1, aes(x=PCR_size, y=GangSTR)) +
  geom_count(colour = 'green', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,60) +
  ylim(0,60) +
  xlab('PCR expansion size') +
  ylab('GangSTR prediction') +
  ggtitle('Number of alleles = 150') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))
ggplot(data = updateallele1, aes(x=PCR_size, y=EHv2_size)) +
  geom_count(colour = 'purple', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +  
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,60) +
  ylim(0,60) +
  xlab('PCR expansion size') +
  ylab('EHv2 prediction') +
  ggtitle('Number of alleles = 150') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))
ggplot(data = updateallele1, aes(x=PCR_size, y=EHv3_size)) +
  geom_count(colour = 'blue', alpha = 0.4) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = 35, color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 35, color = 'red', linetype = 'dashed') +
  xlim(0,60) +
  ylim(0,60) +
  xlab('PCR expansion size') +
  ylab('EHv3 prediction') +
  ggtitle('Number of alleles = 150') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  guides(color = guide_legend(title = 'number_of_alleles'))+
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))

updateallele$PCR_EH_difference <- updateallele$PCR_size - updateallele$EHv3_size
updateallele$PCR_EH_difference.abs <- abs(updateallele$PCR_EH_difference)
updateallele$Spanning[updateallele$Spanning > 8] <- '>8'
updateallele[complete.cases(updateallele), ]
updateallele_below36 <- updateallele %>% filter(PCR_size < 36, EHv3_size < 36)

ggplot(data = updateallele, aes(x=factor(Spanning, levels = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '>8')), y=PCR_EH_difference.abs)) +
  geom_count(colour = 'red') +
  theme_bw() +
  xlim(0, 1, 2, 3, 4, 5, 6, 7, 8, '>8') +
  xlab('No. of Spanning reads') +
  ylab('Difference between EHv3 and PCR') +
  ggtitle('Number of alleles = 154') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))

ggplot(data = updateallele_below36, aes(x=factor(Spanning, levels = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '>8')), y=PCR_EH_difference.abs)) +
  geom_count(colour = 'blue') +
  theme_bw() +
  xlim(0, 1, 2, 3, 4, 5, 6, 7, 8, '>8') +
  xlab('No. of Spanning reads') +
  ylab('Difference between EHv3 and PCR (<36)') +
  ggtitle('Number of alleles = 140') +
  guides(color = guide_legend(override.aes = list(size = 15))) +
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 10))

#rmse for each prediction method
G <- sqrt(mean((Gang$PCR - Gang$value)^2, na.rm = TRUE))
#7.119
E2 <- sqrt(mean((df_analysis_complete$PCR - df_analysis_complete$EHv255)^2, na.rm = TRUE))
#4.249
E3 <- sqrt(mean((df_analysis_complete$PCR - df_analysis_complete$EHv312.x)^2, na.rm = TRUE))
#4.644

#rmse without two large expansions
DF_1 <- DF[-c(51, 77, 128, 154),]
E2_2 <- sqrt(mean((DF_1$PCR - DF_1$EHv255)^2, na.rm = TRUE))
#3.073
E3_2 <- sqrt(mean((DF_1$PCR - DF_1$EHv312.x)^2, na.rm = TRUE))
#2.949

#rmse without expansions above 40
DF_2 <- DF[-c(51, 77, 128, 154, 16, 93, 21, 98, 42, 119, 49, 126, 54, 131, 74, 151, 75, 152),]
E2_3 <- sqrt(mean((DF_2$PCR - DF_2$EHv255)^2, na.rm = TRUE))
#0.248
E3_3 <- sqrt(mean((DF_2$PCR - DF_2$EHv312.x)^2, na.rm = TRUE))
#0.421

#rmse without expansions above 40 GangSTR
DF_G <- Gang[-c(154, 152, 131, 128, 126, 93),]
G_1 <- sqrt(mean((DF_G$PCR - DF_G$value)^2, na.rm = TRUE))
#6.731

