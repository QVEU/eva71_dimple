#Loading libraries required for generation of dataframes for indel paper analysis.
library(ggplot2)
library(tidyverse)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggridges)
library(ineq)
library(RColorBrewer)
library(stringr)
library(gglorenz)
library(readr)
library(scales)

#Reading in dataframes need for figure generation and analysis
allStickles_final_input <- read_csv("allStickles_final_input_handle.csv") 
alldel_final_input <- read_csv("alldel_final_input.csv") 
Capsid_P2_Scores_shared <- read_csv("Capsid_P2_Scores_shared.csv") 
Replication_P2_Scores_shared <- read_csv("Replication_P2_Scores_shared.csv") 
Scores_Insertional_Handle_Fullproteome <- read_csv("Scores_Insertional_Handle_Fullproteome.csv") 
Scores_Deletions_Fullproteome <- read_csv("Scores_Deletions_Fullproteome.csv") 
Scores_Deletions1AA_Fullproteome_1AA <- read_csv("Scores_Deletions1AA_Fullproteome_1AA.csv") 
Capsid_P2_Deletions_shared <- read_csv("Capsid_P2_Deletions_shared.csv") 
Replication_P2_Deletion_shared <- read_csv("Replication_P2_Deletion_shared.csv") 
Fullproteome_input_DMS_Enrich2_long <- read_csv("Fullproteome_input_DMS_Enrich2_long.csv") 
Fullproteome_P2_DMS_Enrich2_long <- read_csv("Fullproteome_P2_DMS_Enrich2_long.csv") 
Capsid_P2_DMS_Enrich2_shared <- read_csv("Capsid_P2_DMS_Enrich2_shared.csv") 
Replication_P2_DMS_Enrich2_shared <- read_csv("Replication_P2_DMS_Enrich2_shared.csv") 
allStickles_AA_final_input <- read_csv("allStickles_AA_final_input.csv") 
Fullproteome_1AA_Enrich2_long <- read_csv("Fullproteome_1AA_Enrich2_long.csv") 
Competition_Pool_Enrich2_Scores <- read_tsv("main_identifiers_scores_Competition_Pooled_Experiment.tsv", col_names = TRUE)
Scores_Competition_input <- read_tsv("Competition_Pool_Input.tsv", col_names = TRUE)
EV_Features <- read_csv("EV71_4643_Features.csv")

#Figure 1. Deep Insertion, Deletion, and Mutational Scanning of the EV-A71 proteome

#Bar plot of insertional handle input library 
G_input_handle<-ggplot(allStickles_final_input)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=5),alpha=rep(times = 1,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(indel,n,fill=Seq),binwidth = 5,fun = function(X) {log10(sum(X))} )+scale_fill_manual(values = "#8BC3E4")+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(77,157,235,314,392,470,548,626,706,784,862,940,1018,1096,1174,1252,1330,1409,1488,1566,1645,1724,1802,1881,1958,2037,2115,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10) 
G_input_handle
#Save figure
ggsave(plot = G_input_handle,filename = "meanPlot_inserts_input_handle.pdf",width =5,height=2)

#Ridgeline plot of insertional handle input library
Sublib <- c(0,77,157,235,314,392,470,548,626,706,784,862,940,1018,1096,1174,1252,1330,1409,1488,1566,1645,1724,1802,1881,1958,2037,2115,2193)
Sublib_label  <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28")
allStickles_final_input$cat <- cut(allStickles_final_input$indel, Sublib, Sublib_label)
allStickles_final_input_filter<- allStickles_final_input %>%group_by(cat)%>% mutate(Median=median(n+1))
Ridgeplot_input_lib <- ggplot(allStickles_final_input_filter, aes(x = n+1, y = cat)) + geom_density_ridges2(aes(fill=Median))+scale_y_discrete(expand = c(0, 0), breaks=c(1, 10, 20, 30)) + scale_x_log10(lim = c(50, 1000))+ylab("Sublibraries")+ xlab("Count")
Ridgeplot_input_lib
#Save figure
ggsave(plot = Ridgeplot_input_lib,filename = "Ridge_plot.pdf",device = pdf,width =4,height=4)

#Lorenz curve and Gini Coefficient of insertional handle input library

# Calculate Lorenz curve points for all of the dataset
lorenz_data_master <- ineq::Lc(allStickles_final_input$n)
# Calculate Gini coefficient
gini_coef_master <- ineq::Gini(allStickles_final_input$n)
cat("Gini Coefficient_master:", gini_coef_master, "\n")
# Create a dataframe for Lorenz curve plotting
lorenz_df_master <- data.frame(x = lorenz_data_master$p, y = lorenz_data_master$L)
# Calculate Gini coefficient for each category
gini_results <- tapply(allStickles_final_input_filter$n, allStickles_final_input_filter$cat, Gini)
# Convert results to a data frame
gini_df <- data.frame(cat = names(gini_results), gini = gini_results)
gini_df
min_value <- min(gini_df$gini)
max_value <- max(gini_df$gini)
cat("Minimum value:", min_value, "\n")
cat("Maximum value:", max_value, "\n")

# Create a function to calculate Lorenz curve
calculate_lorenz <- function(n) {
  lorenz_data <- Lc(n)
  lorenz_df <- data.frame(x = lorenz_data$p, y = lorenz_data$L)
  return(lorenz_df)
}

# Calculate Lorenz curves for each category
lorenz_results <- by(allStickles_final_input_filter$n, allStickles_final_input_filter$cat, calculate_lorenz)

# Convert the results to a list of data frames
lorenz_list <- lapply(lorenz_results, as.data.frame)

# Combine the list of data frames into a single data frame
lorenz_df <- do.call(rbind, lorenz_list)

# Add a column for category
lorenz_df$category <- rep(names(lorenz_results), sapply(lorenz_results, nrow))

#Bind master_df with sub-libraries
loren<-bind_rows(lorenz_df_master, lorenz_df)

#Create a graph to plot Lorenz curves 
Lorenz_Curve_sublib <- ggplot(loren, aes(x = x, y = y, col=factor(category, levels=c(1:28) ))) +
  geom_line() +geom_abline(intercept = 0, slope = 1, linetype = "dashed") +labs(title = "Gini Coefficient=0.22(0.08-0.34)", x = "Cumulative Proportion of Variants",y = "Cumulative Proportion of Reads")+theme_classic()
#Create a new graph to plot Lorenz curve sub-library and the master df with a dotted line
Lorenz_Curve_Final  <- Lorenz_Curve_sublib+geom_line(data = loren[is.na(loren$category), ], color = "black", size = 1,linetype="dotted")+theme(legend.position = "none")
Lorenz_Curve_Final

#Save figure
ggsave(plot = Lorenz_Curve_Final,filename = "Lorenz_Curve_Final.pdf",device = pdf,width =4,height=4)

#Statistics included in the Figure 2 description
# Count positions with count more than 1
positions_above_1_insertional_handle <- sum(allStickles_final_input$n > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_insertional_handle/2193, "\n")


#Bar plot of deletion 3 bp input library 
alldel_final_input_3bp  <- filter(alldel_final_input,Seq=="3bp")
G_input_del<-ggplot(alldel_final_input_3bp)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=5),alpha=rep(times = 1,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(indel,n,fill=Seq),binwidth = 5,fun = function(X) {log10(sum(X))} )+scale_fill_manual(values = "#EB8F93")+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10)
G_input_del
ggsave(plot = G_input_del,filename = "meanPlot_input_3bp.pdf",width = 5,height=2)

#Ridgeline plot of deletion 3 bp input library 
Sublib <- c(0,47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193)
Sublib_label  <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45")
alldel_final_input_3bp$cat <- cut(alldel_final_input_3bp$indel, Sublib, Sublib_label)
alldel_final_input_filter<-alldel_final_input_3bp %>%group_by(cat)%>%mutate(Median=median(n+1))
Ridgeplot_input_lib_del <- ggplot(alldel_final_input_filter, aes(x = n+1, y = cat)) + geom_density_ridges2(aes(fill=Median))+scale_y_discrete(expand = c(0, 0), breaks=c(1, 10, 20, 30, 40, 50)) + scale_x_log10(lim = c(50, 3000))+ylab("Sublibraries")+ xlab("Count")
Ridgeplot_input_lib_del
ggsave(plot = Ridgeplot_input_lib_del,filename = "Ridge_plot_del.pdf",device = pdf,width = 4,height=4)

#Lorenz curve and Gini Coefficient of deletion 3 bp input library 

# Calculate Lorenz curve points for all of the dataset
lorenz_data_master <- ineq::Lc(alldel_final_input_3bp$n)

# Calculate Gini coefficient
gini_coef_master <- ineq::Gini(alldel_final_input_3bp$n)
cat("Gini Coefficient_master:", gini_coef_master, "\n")

# Create a dataframe for Lorenz curve plotting
lorenz_df_master <- data.frame(x = lorenz_data_master$p, y = lorenz_data_master$L)

# Calculate Gini coefficient for each category
gini_results <- tapply(alldel_final_input_3bp$n, alldel_final_input_3bp$cat, Gini)

# Convert results to a data frame
gini_df <- data.frame(cat = names(gini_results), gini = gini_results)
gini_df
min_value <- min(gini_df$gini)
max_value <- max(gini_df$gini)
cat("Minimum value:", min_value, "\n")
cat("Maximum value:", max_value, "\n")

# Create a function to calculate Lorenz curve
calculate_lorenz <- function(n) {
  lorenz_data <- Lc(n)
  lorenz_df <- data.frame(x = lorenz_data$p, y = lorenz_data$L)
  return(lorenz_df)
}

# Calculate Lorenz curves for each category
lorenz_results <- by(alldel_final_input_3bp$n, alldel_final_input_3bp$cat, calculate_lorenz)

# Convert the results to a list of data frames
lorenz_list <- lapply(lorenz_results, as.data.frame)

# Combine the list of data frames into a single data frame
lorenz_df <- do.call(rbind, lorenz_list)

# Add a column for category
lorenz_df$category <- rep(names(lorenz_results), sapply(lorenz_results, nrow))

#Bind master_df with sub-libraries
loren<-bind_rows(lorenz_df_master, lorenz_df)

#Create a graph to plot Lorenz curves NA=grey=Masterdf)
Lorenz_Curve_sublib <- ggplot(loren, aes(x = x, y = y, col=factor(category, levels=c(1:45) ))) +
  geom_line() +geom_abline(intercept = 0, slope = 1, linetype = "dashed") +labs(title = "Gini Coefficient=0.36(0.06-0.46)", x = "Cumulative Proportion of Variants",y = "Cumulative Proportion of Reads")+theme_classic()

#Create a new graph to plot Lorenz curve sub-library and the master df with a dotted line
Lorenz_Curve_Final_del  <- Lorenz_Curve_sublib+geom_line(data = loren[is.na(loren$category), ], color = "black", size = 1,linetype="dotted") +theme(legend.position = "none")
Lorenz_Curve_Final_del

#Save Figure
ggsave(plot = Lorenz_Curve_Final_del,filename = "Lorenz_Curve_Final_del.pdf",device = pdf,width =4,height=4)

#Statistics included in the Figure 2 description
#Count positions with count more than 1
positions_above_1_deletion_3bp <- sum(alldel_final_input_3bp$n > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_deletion_3bp/2193, "\n")


#Bar plot of DMS input library 
DMS_plot_input <- ggplot(Fullproteome_input_DMS_Enrich2_long)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=7),alpha=rep(times = 1,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(position,count,fill="#3c7430"),binwidth = 5,fun = function(X) {log10(sum(X))} )+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(51,102,153,204,255,306,357,408,459,510,561,612,662,712,762,812,862,916,970,1024,1078,1132,1186,1239,1292,1345,1398,1451,1504,1557,1610,1663,1716,1769,1822,1875,1928,1981,2034,2087,2140,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10)+scale_fill_manual(values = "#3c7430")
ggsave(plot = DMS_plot_input,filename = "DMS_plot_input.pdf",width = 5,height=2)
DMS_plot_input

#DMS boxplot
Boxplot_DMS_Aminoacidcounts <- ggplot(Fullproteome_input_DMS_Enrich2_long, aes(y = count, x = factor(AminoAcid, levels = unique(AminoAcid)))) +ylim(0,400)+geom_boxplot(outlier.shape = NA)+theme_classic()+labs(y="Median Variant Count",x="")
Boxplot_DMS_Aminoacidcounts

#Save figure
ggsave(plot = Boxplot_DMS_Aminoacidcounts,filename = "Boxplot_DMS_Aminoacidcounts.pdf",width = 5,height=2)

#Statistics included in the Figure 2 description
#Count positions with count more than 1
Fullproteome_input_DMS_Enrich2_long_NAreplaced0 <- Fullproteome_input_DMS_Enrich2_long %>% replace(is.na(.), 0)
positions_above_1_DMS <- sum(Fullproteome_input_DMS_Enrich2_long_NAreplaced0$count > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_DMS/41667, "\n")



#Supplementary Figure 2. Deep Deletion Scanning of the EV-A71 proteome

#Bar plot of deletion 6 bp input library 
alldel_final_input_6bp  <- filter(alldel_final_input,Seq=="6bp")
G_input_del_6bp<-ggplot(alldel_final_input_6bp)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=5),alpha=rep(times = 1,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(indel,n,fill=Seq),binwidth = 5,fun = function(X) {log10(sum(X))} )+scale_fill_manual(values = "#EB8F93")+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10)
G_input_del_6bp
#Save figure
ggsave(plot = G_input_del_6bp,filename = "meanPlot_input_6bp.pdf",width = 5,height=2)

#Ridgeline plot of deletion 6 bp input library 
Sublib <- c(0,47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193)
Sublib_label  <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45")
alldel_final_input_6bp$cat <- cut(alldel_final_input_6bp$indel, Sublib, Sublib_label)
alldel_final_input_6bp_filter<-alldel_final_input_6bp %>%group_by(cat)%>%mutate(Median=median(n+1))
Ridgeplot_input_6bp_lib_del <- ggplot(alldel_final_input_6bp_filter, aes(x = n+1, y = cat)) + geom_density_ridges2(aes(fill=Median))+scale_y_discrete(expand = c(0, 0), breaks=c(1, 10, 20, 30, 40, 50)) + scale_x_log10(lim = c(50, 3000))+ylab("Sublibraries")+ xlab("Count")
Ridgeplot_input_6bp_lib_del
#Save figure
ggsave(plot = Ridgeplot_input_6bp_lib_del,filename = "Ridge_plot_del_6bp.pdf",device = pdf,width = 4,height=4)

#Lorenz curve and Gini Coefficient of deletion 6 bp input library 

# Calculate Lorenz curve points for all of the dataset
lorenz_data_master <- ineq::Lc(alldel_final_input_6bp_filter$n)

# Calculate Gini coefficient
gini_coef_master <- ineq::Gini(alldel_final_input_6bp_filter$n)
cat("Gini Coefficient_master:", gini_coef_master, "\n")

# Create a dataframe for Lorenz curve plotting
lorenz_df_master <- data.frame(x = lorenz_data_master$p, y = lorenz_data_master$L)

# Calculate Gini coefficient for each category
gini_results <- tapply(alldel_final_input_6bp_filter$n, alldel_final_input_6bp_filter$cat, Gini)

# Convert results to a data frame
gini_df <- data.frame(cat = names(gini_results), gini = gini_results)
gini_df
min_value <- min(gini_df$gini)
max_value <- max(gini_df$gini)
cat("Minimum value:", min_value, "\n")
cat("Maximum value:", max_value, "\n")

# Create a function to calculate Lorenz curve
calculate_lorenz <- function(n) {
  lorenz_data <- Lc(n)
  lorenz_df <- data.frame(x = lorenz_data$p, y = lorenz_data$L)
  return(lorenz_df)
}

# Calculate Lorenz curves for each category
lorenz_results <- by(alldel_final_input_6bp_filter$n, alldel_final_input_6bp_filter$cat, calculate_lorenz)

# Convert the results to a list of data frames
lorenz_list <- lapply(lorenz_results, as.data.frame)

# Combine the list of data frames into a single data frame
lorenz_df <- do.call(rbind, lorenz_list)

# Add a column for category
lorenz_df$category <- rep(names(lorenz_results), sapply(lorenz_results, nrow))

#Bind master_df with sub-libraries
loren<-bind_rows(lorenz_df_master, lorenz_df)

#Create a graph to plot Lorenz curves 
Lorenz_Curve_sublib_6bp <- ggplot(loren, aes(x = x, y = y, col=factor(category, levels=c(1:45) ))) +
  geom_line() +geom_abline(intercept = 0, slope = 1, linetype = "dashed") +labs(title = "Gini Coefficient=0.34(0.06-0.47)", x = "Cumulative Proportion of Variants",y = "Cumulative Proportion of Reads")+theme_classic()

#Create a new graph to plot Lorenz curve sub-library and the master df with a dotted line
Lorenz_Curve_Final_del_6bp  <- Lorenz_Curve_sublib_6bp+geom_line(data = loren[is.na(loren$category), ], color = "black", size = 1,linetype="dotted")+theme(legend.position = "none")
Lorenz_Curve_Final_del_6bp

#Save figure
ggsave(plot = Lorenz_Curve_Final_del_6bp,filename = "Lorenz_Curve_Final_del_6bp.pdf",device = pdf,width = 4,height=4)

#Statistics included in the Figure 2 description
#Count positions with count more than 1
positions_above_1_deletion_6bp <- sum(alldel_final_input_6bp$n > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_deletion_6bp/2193, "\n")

#Bar plot of deletion 9 bp input library 
alldel_final_input_9bp  <- filter(alldel_final_input,Seq=="9bp")
G_input_del_9bp<-ggplot(alldel_final_input_9bp)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=5),alpha=rep(times = 1,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(indel,n,fill=Seq),binwidth = 5,fun = function(X) {log10(sum(X))} )+scale_fill_manual(values = "#EB8F93")+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10)
G_input_del_9bp

#Save figure
ggsave(plot = G_input_del_9bp,filename = "meanPlot_input_9bp.pdf",width = 5,height=2)

#Ridgeline plot of deletion 9 bp input library 
Sublib <- c(0,47,96,144,192,240,288,336,384,432,480,528,576,624,672,720,768,815,862,911,961,1011,1061,1111,1161,1211,1261,1310,1359,1408,1457,1506,1555,1604,1653,1702,1751,1800,1848,1898,1948,1997,2046,2095,2144,2193)
Sublib_label  <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45")
alldel_final_input_9bp$cat <- cut(alldel_final_input_9bp$indel, Sublib, Sublib_label)
alldel_final_input_9bp_filter<-alldel_final_input_9bp %>%group_by(cat)%>%mutate(Median=median(n+1))
Ridgeplot_input_9bp_lib_del <- ggplot(alldel_final_input_9bp_filter, aes(x = n+1, y = cat)) + geom_density_ridges2(aes(fill=Median))+scale_y_discrete(expand = c(0, 0), breaks=c(1, 10, 20, 30, 40, 50)) + scale_x_log10(lim = c(50, 3000))+ylab("Sublibraries")+ xlab("Count")
Ridgeplot_input_9bp_lib_del
ggsave(plot = Ridgeplot_input_9bp_lib_del,filename = "Ridge_plot_del_9bp.pdf",device = pdf,width = 4,height=4)

#Lorenz curve and Gini Coefficient of deletion 9 bp input library 

# Calculate Lorenz curve points for all of the dataset
lorenz_data_master <- ineq::Lc(alldel_final_input_9bp_filter$n)

# Calculate Gini coefficient
gini_coef_master <- ineq::Gini(alldel_final_input_9bp_filter$n)
cat("Gini Coefficient_master:", gini_coef_master, "\n")

# Create a dataframe for Lorenz curve plotting
lorenz_df_master <- data.frame(x = lorenz_data_master$p, y = lorenz_data_master$L)

# Calculate Gini coefficient for each category
gini_results <- tapply(alldel_final_input_9bp_filter$n, alldel_final_input_9bp_filter$cat, Gini)

# Convert results to a data frame
gini_df <- data.frame(cat = names(gini_results), gini = gini_results)
gini_df
min_value <- min(gini_df$gini)
max_value <- max(gini_df$gini)
cat("Minimum value:", min_value, "\n")
cat("Maximum value:", max_value, "\n")

# Create a function to calculate Lorenz curve
calculate_lorenz <- function(n) {
  lorenz_data <- Lc(n)
  lorenz_df <- data.frame(x = lorenz_data$p, y = lorenz_data$L)
  return(lorenz_df)
}

# Calculate Lorenz curves for each category
lorenz_results <- by(alldel_final_input_9bp_filter$n, alldel_final_input_9bp_filter$cat, calculate_lorenz)

# Convert the results to a list of data frames
lorenz_list <- lapply(lorenz_results, as.data.frame)

# Combine the list of data frames into a single data frame
lorenz_df <- do.call(rbind, lorenz_list)

# Add a column for category
lorenz_df$category <- rep(names(lorenz_results), sapply(lorenz_results, nrow))

#Bind master_df with sub-libraries
loren<-bind_rows(lorenz_df_master, lorenz_df)

#Create a graph to plot Lorenz curves NA=grey=Masterdf)
Lorenz_Curve_sublib_9bp <- ggplot(loren, aes(x = x, y = y, col=factor(category, levels=c(1:45) ))) +
  geom_line() +geom_abline(intercept = 0, slope = 1, linetype = "dashed") +labs(title = "Gini Coefficient=0.34(0.06-0.48)", x = "Cumulative Proportion of Variants",y = "Cumulative Proportion of Reads")+theme_classic()

#Create a new graph to plot Lorenz curve sub-library and the master df with a dotted line
Lorenz_Curve_Final_del_9bp  <- Lorenz_Curve_sublib_9bp+geom_line(data = loren[is.na(loren$category), ], color = "black", size = 1,linetype="dotted")+theme(legend.position = "none")
Lorenz_Curve_Final_del_9bp

#Save figure
ggsave(plot = Lorenz_Curve_Final_del_9bp,filename = "Lorenz_Curve_Final_del_9bp.pdf",device = pdf,width = 4,height=4)

#Statistics included in the Figure 2 description
#Count positions with count more than 1
positions_above_1_deletion_9bp <- sum(alldel_final_input_9bp$n > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_deletion_9bp/2193, "\n")

#Boxplot of different deletion sizes
Boxplot_Delsize_input <- ggplot(alldel_final_input, aes(y=n,x=Seq))+geom_boxplot(outlier.shape = NA)+ylim(0,500)+ylab("Median Variant Count")+xlab("")+theme_classic()
Boxplot_Delsize_input

#Save figure
ggsave(plot = Boxplot_Delsize_input,filename = "Boxplot_Delsize_input.pdf",device = pdf,width = 1,height=1)

#Scatter plots for different sizes of deletions

#Tidying dataframe for 3bp, creating a column 3bp with input counts
alldel_final_input_filtered_3bp <- filter(alldel_final_input,Seq=="3bp") 
alldel_final_input_filtered_3bp <-  pivot_wider(alldel_final_input_filtered_3bp,names_from = Seq, values_from =n)
alldel_final_input_filtered_3bp <- select(alldel_final_input_filtered_3bp,indel,"3bp")

#Tidying dataframe for 6bp, creating a column 6bp with input counts
alldel_final_input_filtered_6bp <- filter(alldel_final_input,Seq=="6bp") 
alldel_final_input_filtered_6bp <-  pivot_wider(alldel_final_input_filtered_6bp,names_from = Seq, values_from =n)
alldel_final_input_filtered_6bp <- select(alldel_final_input_filtered_6bp,indel,"6bp")

#Tidying dataframe for 9bp, creating a column 9bp with input counts
alldel_final_input_filtered_9bp <- filter(alldel_final_input,Seq=="9bp") 
alldel_final_input_filtered_9bp <-  pivot_wider(alldel_final_input_filtered_9bp,names_from = Seq, values_from =n)
alldel_final_input_filtered_9bp <- select(alldel_final_input_filtered_9bp,indel,"9bp")

#Merging all three dataframes
alldel_final_input_filtered_3_6_bp <- merge(alldel_final_input_filtered_3bp,alldel_final_input_filtered_6bp,by="indel",all.x=TRUE,all.y=TRUE)
alldel_final_input_filtered_3_9_bp <- merge(alldel_final_input_filtered_3bp,alldel_final_input_filtered_9bp,by="indel",all.x=TRUE,all.y=TRUE)
alldel_final_input_filtered_6_9_bp <- merge(alldel_final_input_filtered_6bp,alldel_final_input_filtered_9bp,by="indel",all.x=TRUE,all.y=TRUE)


#Plotting 1AA against 2AA deletions and then using a linear model to obtain correlation, R-squared
Scatter_input_1AA_2AA <- ggplot(alldel_final_input_filtered_3_6_bp, aes(x = `3bp`, y = `6bp`)) + geom_point(size = 0.05,color="#EB8F93") +scale_x_log10()+scale_y_log10()+
labs(x = "1 AA", y = "2 AA") +theme_classic()
Scatter_input_1AA_2AA
ggsave(plot = Scatter_input_1AA_2AA,filename = "Scatter_input_1AA_2AA.pdf",device = pdf,width = 2,height=2)
model_1AA_2AA <- lm(`3bp` ~ `6bp`, data = alldel_final_input_filtered_3_6_bp)
summary_model_1AA_2AA <- summary(model_1AA_2AA)
summary_model_1AA_2AA

#Plotting 1AA against 3AA deletions and then using a linear model to obtain correlation, R-squared
Scatter_input_1AA_3AA <- ggplot(alldel_final_input_filtered_3_9_bp, aes(x = `3bp`, y = `9bp`)) + geom_point(size = 0.05,color="#EB8F93") +scale_x_log10()+scale_y_log10()+
labs(x = "1 AA", y = "3 AA") +theme_classic()
Scatter_input_1AA_3AA
ggsave(plot = Scatter_input_1AA_3AA,filename = "Scatter_input_1AA_3AA.pdf",device = pdf,width = 2,height=2)
model_1AA_3AA <- lm(`3bp` ~ `9bp`, data = alldel_final_input_filtered_3_9_bp)
summary_model_1AA_3AA <- summary(model_1AA_3AA)
summary_model_1AA_3AA


#Plotting 1AA against 3AA deletions and then using a linear model to obtain correlation, R-squared
Scatter_input_2AA_3AA <- ggplot(alldel_final_input_filtered_6_9_bp, aes(x = `6bp`, y = `9bp`)) + geom_point(size = 0.05,color="#EB8F93") +scale_x_log10()+scale_y_log10()+
labs(x = "2 AA", y = "3 AA") +theme_classic()
Scatter_input_2AA_3AA
ggsave(plot = Scatter_input_2AA_3AA,filename = "Scatter_input_2AA_3AA.pdf",device = pdf,width = 2,height=2)
model_2AA_3AA <- lm(`6bp` ~ `9bp`, data = alldel_final_input_filtered_6_9_bp)
summary_model_2AA_3AA <- summary(model_2AA_3AA)
summary_model_2AA_3AA

#Supplementary Figure 3. Deep Insertion Scanning of the EV-A71 proteome

#Bar plot of 5AA input library 
G_input_AA<-ggplot(allStickles_AA_final_input)+
geom_rect(data=EV_Features,aes(xmin=(Start-746)/3,xmax=(End-746)/3,ymin=0,ymax=5),alpha=rep(times = 20,c(0.15,0,0.15,0,0.15,0,0.15,0,0.15,0.15,0)))+
stat_summary_bin(show.legend = F,geom = "bar",aes(indel,n),binwidth = 5,fun = function(X) {log10(sum(X))} )+facet_grid(factor(Seq,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"))~.)+ facet_wrap(~Seq, ncol = 5)+
theme_classic()+ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+geom_vline(xintercept=c(77,157,235,314,392,470,548,626,706,784,862,940,1018,1096,1174,1252,1330,1409,1488,1566,1645,1724,1802,1881,1958,2037,2115,2193),linetype="dashed",colour="black",size=0.5,alpha=0.10)

#Save figure
ggsave(plot = G_input_AA,filename = "G_input_AA.pdf",width = 15,height=10)

#Print the result
positions_above_1_AA <- sum(allStickles_AA_final_input$n > 1)
cat("Number of positions with count more than 1:", positions_above_1_AA/43860, "\n")

#Calculate median count for different variants
Boxplot_1AA_input <- ggplot(allStickles_AA_final_input, aes(y=n,x=Seq))+geom_boxplot(outlier.shape = NA)+ylim(0,65)+ylab("Median Variant Count")+xlab("")+theme_classic()
Boxplot_1AA_input

# Calculate median for each unique sequence
unique_Seq <- unique(allStickles_AA_final_input$Seq)
Medians <- numeric(length(unique_Seq))
for (i in 1:length(unique_Seq)) {
  Medians[i] <- median(allStickles_AA_final_input$n[allStickles_AA_final_input$Seq == unique_Seq[i]])
}

# Create a data frame with unique sequences and their medians
Median_1AA_Seq <- data.frame(seq = unique_Seq, median_count = Medians)
Median_1AA_Seq
ggsave(plot = Boxplot_1AA_input,filename = "Boxplot_1AA_input.pdf",width = 3,height=2)
Saturation_1AA <- allStickles_AA_final_input %>% group_by(Seq) %>% summarize(saturation = sum(n > 1)/2193)
Saturation_1AA



#Supplementary Figure 4. Deep Mutational Scanning of the EV-A71 proteome

#Heatmap for DMS input
# Convert Position to a factor with breaks
breaks_inputheatmap <- seq(0, 2193, by = 250)

# Create a function to assign labels based on position
get_label_input <- function(position) {
  sub_labels <- c("A", "B", "C", "D", "E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V")
  sub_index <- findInterval(position, breaks_inputheatmap)
  return(sub_labels[sub_index])
}

#create a colum sub for the different rows of the heatmap
Fullproteome_input_DMS_Enrich2_long$sub <- sapply(Fullproteome_input_DMS_Enrich2_long$position, get_label_input)

#Use amino acid as factors
Fullproteome_input_DMS_Enrich2_long$AminoAcid=factor(Fullproteome_input_DMS_Enrich2_long$AminoAcid,levels = c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))

#Heatmap of DMS input library
DMS_Heatmap_input <- ggplot(Fullproteome_input_DMS_Enrich2_long) + geom_tile(aes(position,AminoAcid, fill = log10(count))) +
    scale_fill_gradient2(limits = c(0,4),low = muted("white"), high = muted("aquamarine1"),space = "Lab",na.value = "black") +
    theme_classic() +  facet_wrap(sub ~ ., scales = "free_x",ncol=1) + theme(aspect.ratio = 0.1,axis.text.y = element_text(size = 3))+ theme(strip.text.x = element_blank())
DMS_Heatmap_input

#Save figure
ggsave(plot = DMS_Heatmap_input,filename = "DMS_Heatmap_input.pdf",width = 8,height=11)

#Supplementary Figure 5. Reproducibility of the Deep Insertion, Deletion, and Mutational Scanning experiments 

#Comparing the replicates for insertional handle library: Passage 2. Scatter plots of replicates against each other.
Capsid_P2_Scores_shared
Replication_P2_Scores_shared

#Passage 2 Replication

#RepA vs RepB Passage 2 Replication
Scatter_P2_insertion_replication_RepA_RepB <- ggplot(Replication_P2_Scores_shared, aes(x = `RepA_score`, y = `RepB_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_insertion_replication_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_insertion_replication_RepA_RepB,filename = "Scatter_P2_insertion_replication_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepA_RepB_replication <- lm(`RepA_score` ~ `RepB_score`, data = Replication_P2_Scores_shared)
summary_model_P2_insertion_RepA_RepB_replication <- summary(model_P2_insertion_RepA_RepB_replication)
summary_model_P2_insertion_RepA_RepB_replication

#RepA vs RepC Passage 2 Replication
Scatter_P2_insertion_replication_RepA_RepC <- ggplot(Replication_P2_Scores_shared, aes(x = `RepA_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_insertion_replication_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_insertion_replication_RepA_RepC,filename = "Scatter_P2_insertion_replication_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepA_RepC_replication <- lm(`RepA_score` ~ `RepC_score`, data = Replication_P2_Scores_shared)
summary_model_P2_insertion_RepA_RepC_replication <- summary(model_P2_insertion_RepA_RepC_replication)
summary_model_P2_insertion_RepA_RepC_replication

#RepB vs RepC Passage 2 Replication
Scatter_P2_insertion_replication_RepB_RepC <- ggplot(Replication_P2_Scores_shared, aes(x = `RepB_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_insertion_replication_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_insertion_replication_RepB_RepC,filename = "Scatter_P2_insertion_replication_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepB_RepC_replication <- lm(`RepB_score` ~ `RepC_score`, data = Replication_P2_Scores_shared)
summary_model_P2_insertion_RepB_RepC_replication <- summary(model_P2_insertion_RepB_RepC_replication)
summary_model_P2_insertion_RepB_RepC_replication

#Passage 2 Capsid

#RepA vs RepB Passage 2 Capsid
Scatter_P2_insertion_RepA_RepB <- ggplot(Capsid_P2_Scores_shared, aes(x = `RepA_score`, y = `RepB_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_insertion_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_insertion_RepA_RepB,filename = "Scatter_P2_insertion_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepA_RepB <- lm(`RepA_score` ~ `RepB_score`, data = Capsid_P2_Scores_shared)
summary_model_P2_insertion_RepA_RepB <- summary(model_P2_insertion_RepA_RepB)
summary_model_P2_insertion_RepA_RepB

#RepA vs RepC Passage 2 Capsid
Scatter_P2_insertion_RepA_RepC <- ggplot(Capsid_P2_Scores_shared, aes(x = `RepA_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_insertion_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_insertion_RepA_RepC,filename = "Scatter_P2_insertion_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepA_RepC <- lm(`RepA_score` ~ `RepC_score`, data = Capsid_P2_Scores_shared)
summary_model_P2_insertion_RepA_RepC <- summary(model_P2_insertion_RepA_RepC)
summary_model_P2_insertion_RepA_RepC

#RepB vs RepC Passage 2 Capsid
Scatter_P2_insertion_RepB_RepC <- ggplot(Capsid_P2_Scores_shared, aes(x = `RepB_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_insertion_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_insertion_RepB_RepC,filename = "Scatter_P2_insertion_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_insertion_RepB_RepC <- lm(`RepB_score` ~ `RepC_score`, data = Capsid_P2_Scores_shared)
summary_model_P2_insertion_RepB_RepC <- summary(model_P2_insertion_RepB_RepC)
summary_model_P2_insertion_RepB_RepC

#Comparing the replicates for deletion handle library: Passage 2 Scatter plots of replicates against each other.

#Passage 2 Replication shared Deletion

#RepA vs RepB Passage 2 Replication
Scatter_P2_deletion__Replication_RepA_RepB <- ggplot(Replication_P2_Deletion_shared, aes(x = `RepA_score`, y = `RepB_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_deletion__Replication_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_deletion__Replication_RepA_RepB,filename = "Scatter_P2_deletion__Replication_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_replication_RepA_RepB <- lm(`RepA_score` ~ `RepB_score`, data = Replication_P2_Deletion_shared)
summary_model_P2_deletion_replication_RepA_RepB <- summary(model_P2_deletion_replication_RepA_RepB)
summary_model_P2_deletion_replication_RepA_RepB

#RepA vs RepC Passage 2 Replication
Scatter_P2_deletion__Replication_RepA_RepC <- ggplot(Replication_P2_Deletion_shared, aes(x = `RepA_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_deletion__Replication_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_deletion__Replication_RepA_RepC,filename = "Scatter_P2_deletion__Replication_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_replication_RepA_RepC <- lm(`RepA_score` ~ `RepC_score`, data = Replication_P2_Deletion_shared)
summary_model_P2_deletion_replication_RepA_RepC <- summary(model_P2_deletion_replication_RepA_RepC)
summary_model_P2_deletion_replication_RepA_RepC

#RepB vs RepC Passage 2 Replication
Scatter_P2_deletion__Replication_RepB_RepC <- ggplot(Replication_P2_Deletion_shared, aes(x = `RepB_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_deletion__Replication_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_deletion__Replication_RepB_RepC,filename = "Scatter_P2_deletion__Replication_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_replication_RepB_RepC <- lm(`RepB_score` ~ `RepC_score`, data = Replication_P2_Deletion_shared)
summary_model_P2_deletion_replication_RepB_RepC <- summary(model_P2_deletion_replication_RepB_RepC)
summary_model_P2_deletion_replication_RepB_RepC

#Passage 2 Capsid shared Deletion

#RepA vs RepB Passage 2 Capsid
Scatter_P2_deletion_RepA_RepB <- ggplot(Capsid_P2_Deletions_shared, aes(x = `RepA_score`, y = `RepB_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_deletion_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_deletion_RepA_RepB,filename = "Scatter_P2_deletion_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_RepA_RepB <- lm(`RepA_score` ~ `RepB_score`, data = Capsid_P2_Deletions_shared)
summary_model_P2_deletion_RepA_RepB <- summary(model_P2_deletion_RepA_RepB)
summary_model_P2_deletion_RepA_RepB

#RepA vs RepC Passage 2 Capsid
Scatter_P2_deletion_RepA_RepC <- ggplot(Capsid_P2_Deletions_shared, aes(x = `RepA_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_deletion_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_deletion_RepA_RepC,filename = "Scatter_P2_deletion_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_RepA_RepC <- lm(`RepA_score` ~ `RepC_score`, data = Capsid_P2_Deletions_shared)
summary_model_P2_deletion_RepA_RepC <- summary(model_P2_deletion_RepA_RepC)
summary_model_P2_deletion_RepA_RepC

#RepB vs RepC Passage 2 Capsid
Scatter_P2_deletion_RepB_RepC <- ggplot(Capsid_P2_Deletions_shared, aes(x = `RepB_score`, y = `RepC_score`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_deletion_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_deletion_RepB_RepC,filename = "Scatter_P2_deletion_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_deletion_RepB_RepC <- lm(`RepB_score` ~ `RepC_score`, data = Capsid_P2_Deletions_shared)
summary_model_P2_deletion_RepB_RepC <- summary(model_P2_deletion_RepB_RepC)
summary_model_P2_deletion_RepB_RepC


#Passage 2 Replication shared DMS

#RepA vs RepB Passage 2 Replication
Scatter_P2_DMS__Replication_RepA_RepB <- ggplot(Replication_P2_DMS_Enrich2_shared, aes(x = `RepA_DMS_score`, y = `RepB_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_DMS__Replication_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_DMS__Replication_RepA_RepB,filename = "Scatter_P2_DMS__Replication_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_replication_RepA_RepB <- lm(`RepA_DMS_score` ~ `RepB_DMS_score`, data = Replication_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_replication_RepA_RepB <- summary(model_P2_DMS_replication_RepA_RepB)
summary_model_P2_DMS_replication_RepA_RepB

#RepA vs RepC Passage 2 Replication
Scatter_P2_DMS__Replication_RepA_RepC <- ggplot(Replication_P2_DMS_Enrich2_shared, aes(x = `RepA_DMS_score`, y = `RepC_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_DMS__Replication_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_DMS__Replication_RepA_RepC,filename = "Scatter_P2_DMS__Replication_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_replication_RepA_RepC <- lm(`RepA_DMS_score` ~ `RepC_DMS_score`, data = Replication_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_replication_RepA_RepC <- summary(model_P2_DMS_replication_RepA_RepC)
summary_model_P2_DMS_replication_RepA_RepC

#RepB vs RepC Passage 2 Replication
Scatter_P2_DMS__Replication_RepB_RepC <- ggplot(Replication_P2_DMS_Enrich2_shared, aes(x = `RepB_DMS_score`, y = `RepC_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_DMS__Replication_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_DMS__Replication_RepB_RepC,filename = "Scatter_P2_DMS__Replication_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_replication_RepB_RepC <- lm(`RepB_DMS_score` ~ `RepC_DMS_score`, data = Replication_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_replication_RepB_RepC <- summary(model_P2_DMS_replication_RepB_RepC)
summary_model_P2_DMS_replication_RepB_RepC

#Passage 2 Capsid shared DMS

#RepA vs RepB Passage 2 Capsid
Scatter_P2_DMS__Capsid_RepA_RepB <- ggplot(Capsid_P2_DMS_Enrich2_shared, aes(x = `RepA_DMS_score`, y = `RepB_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate A", y = "Replicate B") +theme_classic()
Scatter_P2_DMS__Capsid_RepA_RepB
#Save Figure
ggsave(plot = Scatter_P2_DMS__Capsid_RepA_RepB,filename = "Scatter_P2_DMS__Capsid_RepA_RepB.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_capsid_RepA_RepB <- lm(`RepA_DMS_score` ~ `RepB_DMS_score`, data = Capsid_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_capsid_RepA_RepB <- summary(model_P2_DMS_capsid_RepA_RepB)
summary_model_P2_DMS_capsid_RepA_RepB

#RepA vs RepC Passage 2 Capsid
Scatter_P2_DMS__Capsid_RepA_RepC <- ggplot(Capsid_P2_DMS_Enrich2_shared, aes(x = `RepA_DMS_score`, y = `RepC_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate A", y = "Replicate C") +theme_classic()
Scatter_P2_DMS__Capsid_RepA_RepC
#Save Figure
ggsave(plot = Scatter_P2_DMS__Capsid_RepA_RepC,filename = "Scatter_P2_DMS__Capsid_RepA_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_capsid_RepA_RepC <- lm(`RepA_DMS_score` ~ `RepC_DMS_score`, data = Capsid_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_capsid_RepA_RepC <- summary(model_P2_DMS_capsid_RepA_RepC)
summary_model_P2_DMS_capsid_RepA_RepC

#RepB vs RepC Passage 2 Capsid
Scatter_P2_DMS__Capsid_RepB_RepC <- ggplot(Capsid_P2_DMS_Enrich2_shared, aes(x = `RepB_DMS_score`, y = `RepC_DMS_score`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Replicate B", y = "Replicate C") +theme_classic()
Scatter_P2_DMS__Capsid_RepB_RepC
#Save Figure
ggsave(plot = Scatter_P2_DMS__Capsid_RepB_RepC,filename = "Scatter_P2_DMS__Capsid_RepB_RepC.pdf",device = pdf,width = 2,height=2)
#Creating a linear model
model_P2_DMS_capsid_RepB_RepC <- lm(`RepB_DMS_score` ~ `RepC_DMS_score`, data = Capsid_P2_DMS_Enrich2_shared)
summary_model_P2_DMS_capsid_RepB_RepC <- summary(model_P2_DMS_capsid_RepB_RepC)
summary_model_P2_DMS_capsid_RepB_RepC

#Figure 2. Fitness effects of insertions, deletions, and amino acid changes in EV-A71

#Bar plot of Insertional Handle Enrich2 P2
Scores_Insertional_Handle_Fullproteome_Plot <- ggplot()+stat_summary_bin(data = Scores_Insertional_Handle_Fullproteome, show.legend = F, geom = "bar", aes(indel, score), binwidth =1,fill="#1558a9")+  
geom_rect(data = EV_Features, aes(xmin=(Start-746)/3,xmax=(End-746)/3, ymin = 0, ymax = 8), alpha = rep(times = 1, c(0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0.1, 0)))+
theme_classic()+coord_cartesian(xlim=c(0,2193),ylim=c(-6,8))+xlab("Residue Position")+ylab("Enrich2 Score")+geom_hline(col='darkgrey',  lwd=0.75, yintercept=0)
ggsave(plot = Scores_Insertional_Handle_Fullproteome_Plot,filename = "Scores_Insertional_Handle_Fullproteome.pdf",width = 4,height=1.5)
Scores_Insertional_Handle_Fullproteome_Plot

#Bar plot of Deletion library Enrich2 P2
Scores_Deletions1AA_Fullproteome_plot <- ggplot()+stat_summary_bin(data = Scores_Deletions1AA_Fullproteome_1AA, show.legend = F, geom = "bar", aes(indel, score), binwidth =1,fill="#a91515")+  
geom_rect(data = EV_Features, aes(xmin=(Start-746)/3,xmax=(End-746)/3, ymin = 0, ymax =8), alpha = rep(times = 1, c(0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0.1, 0)))+
theme_classic()+coord_cartesian(xlim=c(0,2193),ylim=c(-6,8))+xlab("Residue Position")+ylab("Enrich2 Score")+geom_hline(col='darkgrey',  lwd=0.75, yintercept=0)
ggsave(plot = Scores_Deletions1AA_Fullproteome_plot,filename = "Scores_Deletions1AA_Fullproteome_plot.pdf",width = 4,height=1.5)
Scores_Deletions1AA_Fullproteome_plot

#Bar plot of DMS library Enrich2 P2
Scores_DMS_Fullproteome_plot <- ggplot()+stat_summary_bin(data = Fullproteome_P2_DMS_Enrich2_long, show.legend = F, geom = "bar", aes(position, score), binwidth =1,fill="#3c7430")+  
geom_rect(data = EV_Features, aes(xmin=(Start-746)/3,xmax=(End-746)/3, ymin = 0, ymax =3), alpha = rep(times = 1, c(0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0.1, 0)))+
theme_classic()+coord_cartesian(xlim=c(0,2193),ylim=c(-5,3))+xlab("Residue Position")+ylab("Enrich2 Score")+geom_hline(col='darkgrey',  lwd=0.75, yintercept=0)
ggsave(plot = Scores_DMS_Fullproteome_plot,filename = "Scores_DMS_Fullproteome_plot.pdf",width = 4,height=1.5)
Scores_DMS_Fullproteome_plot

#Figure 2: DMFE of insertions, deletions, and DMS

#Selecting only the score column and plotting a histogram for insertions
histo_insertions_sel  <- select(Scores_Insertional_Handle_Fullproteome,score)
histo_insertions_sel_plot <- ggplot(histo_insertions_sel, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = ("#1558a9"), color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-6,6))
histo_insertions_sel_plot
#Save Figure
ggsave(plot = histo_insertions_sel_plot,filename = "histo_insertions_sel_plot.pdf",width = 2,height=1.5)
#Filtering df for insertions with a score higher than 0
histo_insertions_sel_nonlethal <- filter(histo_insertions_sel,score>0)
#Plotting a histogram for insertions with a score higher than 0
histo_insertions_sel_nonlethal_plot <- ggplot(histo_insertions_sel_nonlethal, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = ("#1558a9"), color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
histo_insertions_sel_nonlethal_plot
#Save Figure
ggsave(plot = histo_insertions_sel_nonlethal_plot,filename = "histo_insertions_sel_nonlethal_plot.pdf",width = 2,height=1.5)

#Selecting only the score column and plotting a histogram for deletions
histo_del_sel  <- select(Scores_Deletions_Fullproteome,score)
del_histogram <- ggplot(histo_del_sel, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#a91515", color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-6,6))
del_histogram
#Save Figure
ggsave(plot = del_histogram,filename = "del_histogram.pdf",width = 2,height=1.5)
#Filtering df for deletions with a score higher than 0
histo_del_sel_nonlethal <- filter(histo_del_sel,score>0)
#Plotting a histogram for deletions with a score higher than 0
del_histogram_nonlethal <- ggplot(histo_del_sel_nonlethal, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#a91515", color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
del_histogram_nonlethal
#Save Figure
ggsave(plot = del_histogram_nonlethal,filename = "del_histogram_nonlethal.pdf",width = 2,height=1.5)

#Selecting only the score column and plotting a histogram for DMS
histo_DMS_sel  <- select(Fullproteome_P2_DMS_Enrich2_long,score)
histo_DMS_sel_plot <- ggplot(histo_DMS_sel, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = ("#3c7430"), color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-6,6))
histo_DMS_sel_plot
#Save Figure
ggsave(plot = histo_DMS_sel_plot,filename = "histo_DMS_sel_plot.pdf",width = 2,height=1.5)
#Filtering df for DMS with a score higher than 0
histo_DMS_sel_nonlethal <- filter(histo_DMS_sel,score>0)
#Plotting a histogram for DMS with a score higher than 0
histo_DMS_sel_plot_nonlethal_plot <- ggplot(histo_DMS_sel_nonlethal, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = ("#3c7430"), color = "black") +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
histo_DMS_sel_plot_nonlethal_plot
#Save Figure
ggsave(plot = histo_DMS_sel_plot_nonlethal_plot,filename = "histo_DMS_sel_plot_nonlethal.pdf",width = 2,height=1.5)

##enrich2: 
## enrich2<0: low
## 0<enrich2<2 : mid
## enrich2>2: high

#Figure 2: Classification of variants into different enrich2 class scores, and comparing distribution for different viral proteins

# Create a mapping of positions to viral proteins
Viral_protein_Mapping <- data.frame(
  start = c(1, 70, 324,566,863,1013,1112,1441,1527,1549,1732),
  end = c(69, 323, 565,862,1012,1111,1440,1526,1548,1731,2193), 
  viralprotein = c("VP4", "VP2", "VP3","VP1","2A","2B","2C","3A","3B","3C","3D"))

#Setting colors for different enrich2 classes
fitness_class_colors <- c( "#202E38","#7B9CAF","#9DCCED")
fitness_class_colors_del <- c( "#441E1E","#BF7B7B","#FFD9D9")
fitness_class_colors_dms <- c( "#285e1b","#529143","#a6f492")

#Dataframe for insertional handle for comparison of enrich2 classes
Handle_P2_Fitness_Class <- Scores_Insertional_Handle_Fullproteome %>% mutate(viralprotein = sapply(indel, function(pos) {
matching_row <- which(pos >= Viral_protein_Mapping$start & pos <= Viral_protein_Mapping$end)
if (length(matching_row) > 0) {return(Viral_protein_Mapping$viralprotein[matching_row])} else {return("No Protein")}})
)

#Adding scores classes depending on the enrich2 scores and plotting the area plot
Handle_P2_Fitness_Class  <- na.omit(Handle_P2_Fitness_Class)
Handle_P2_Fitness_Class$Class="N"
Handle_P2_Fitness_Class[Handle_P2_Fitness_Class$score<2,]$Class="N"
Handle_P2_Fitness_Class[Handle_P2_Fitness_Class$score<0,]$Class="L"
Handle_P2_Fitness_Class[Handle_P2_Fitness_Class$score>2,]$Class="B"
Handle_P2_Fitness_Class_counts <- table(Handle_P2_Fitness_Class$Class,Handle_P2_Fitness_Class$viralprotein)
chisq.test(Handle_P2_Fitness_Class_counts)
Handle_P2_Fitness_Class$viralprotein <- factor(Handle_P2_Fitness_Class$viralprotein, levels = c("VP4", "VP2", "VP3","VP1","2A","2B","2C","3A","3B","3C","3D"))  
Handle_P2_Fitness_Class$Class <- factor(Handle_P2_Fitness_Class$Class, levels = c("B", "N","L"))  
plot(as.factor(Handle_P2_Fitness_Class$viralprotein),as.factor(Handle_P2_Fitness_Class$Class),ylim=c(0.8,1.0),col = fitness_class_colors,lwd = 0.6)
pdf("Insertion_Fitness_Class.pdf",width = 4, height = 4)
plot(as.factor(Handle_P2_Fitness_Class$viralprotein),as.factor(Handle_P2_Fitness_Class$Class),ylim=c(0.8,1.0),col = fitness_class_colors,lwd = 0.6)
dev.off()

#Dataframe for deletion library for comparison of enrich2 classes
Deletion3bp_P2_Fitness_Class <- Scores_Deletions1AA_Fullproteome_1AA %>% mutate(viralprotein = sapply(indel, function(pos) {
matching_row <- which(pos >= Viral_protein_Mapping$start & pos <= Viral_protein_Mapping$end)
if (length(matching_row) > 0) {return(Viral_protein_Mapping$viralprotein[matching_row])} else {return("No Protein")}})
)

#Adding scores classes depending on the enrich2 scores and plotting the area plot
Deletion3bp_P2_Fitness_Class  <- na.omit(Deletion3bp_P2_Fitness_Class)
Deletion3bp_P2_Fitness_Class$Class="N"
Deletion3bp_P2_Fitness_Class[Deletion3bp_P2_Fitness_Class$score<2,]$Class="N"
Deletion3bp_P2_Fitness_Class[Deletion3bp_P2_Fitness_Class$score<0,]$Class="L"
Deletion3bp_P2_Fitness_Class[Deletion3bp_P2_Fitness_Class$score>2,]$Class="B"
Deletion3bp_P2_Fitness_Class_Counts <- table(Deletion3bp_P2_Fitness_Class$Class,Deletion3bp_P2_Fitness_Class$viralprotein)
chisq.test(Deletion3bp_P2_Fitness_Class_Counts)
Deletion3bp_P2_Fitness_Class$viralprotein <- factor(Deletion3bp_P2_Fitness_Class$viralprotein, levels = c("VP4", "VP2", "VP3","VP1","2A","2B","2C","3A","3B","3C","3D"))  
Deletion3bp_P2_Fitness_Class$Class <- factor(Deletion3bp_P2_Fitness_Class$Class, levels = c("B", "N","L"))  
plot(as.factor(Deletion3bp_P2_Fitness_Class$viralprotein),as.factor(Deletion3bp_P2_Fitness_Class$Class),ylim=c(0.8,1.0),col = fitness_class_colors_del,lwd = 0.6)
pdf("Deletion_Fitness_Class_3bp.pdf",width = 4, height = 4)
plot(as.factor(Deletion3bp_P2_Fitness_Class$viralprotein),as.factor(Deletion3bp_P2_Fitness_Class$Class),ylim=c(0.8,1.0),col = fitness_class_colors_del,lwd = 0.6)
dev.off()

#Dataframe for DMS  for comparison of enrich2 classes
DMS_P2_Fitness_Class <- Fullproteome_P2_DMS_Enrich2_long %>% mutate(viralprotein = sapply(position, function(pos) {
matching_row <- which(pos >= Viral_protein_Mapping$start & pos <= Viral_protein_Mapping$end)
if (length(matching_row) > 0) {return(Viral_protein_Mapping$viralprotein[matching_row])} else {return("No Protein")}})
)

#Adding scores classes depending on the enrich2 scores and plotting the area plot
DMS_P2_Fitness_Class  <- na.omit(DMS_P2_Fitness_Class)
DMS_P2_Fitness_Class$Class="N"
DMS_P2_Fitness_Class[DMS_P2_Fitness_Class$score<2,]$Class="N"
DMS_P2_Fitness_Class[DMS_P2_Fitness_Class$score<0,]$Class="L"
DMS_P2_Fitness_Class[DMS_P2_Fitness_Class$score>2,]$Class="B"
DMS_P2_Fitness_Class_Counts <- table(DMS_P2_Fitness_Class$Class,DMS_P2_Fitness_Class$viralprotein)
chisq.test(DMS_P2_Fitness_Class_Counts)
DMS_P2_Fitness_Class$viralprotein <- factor(DMS_P2_Fitness_Class$viralprotein, levels = c("VP4", "VP2", "VP3","VP1","2A","2B","2C","3A","3B","3C","3D"))  
DMS_P2_Fitness_Class$Class <- factor(DMS_P2_Fitness_Class$Class, levels = c("B", "N","L"))  
plot(as.factor(DMS_P2_Fitness_Class$viralprotein),as.factor(DMS_P2_Fitness_Class$Class),ylim=c(0.5,1.0),col = fitness_class_colors_dms,lwd = 0.6)
pdf("DMS_P2_Fitness_Class.pdf",width = 4, height = 4)
plot(as.factor(DMS_P2_Fitness_Class$viralprotein),as.factor(DMS_P2_Fitness_Class$Class),ylim=c(0.5,1.0),col = fitness_class_colors_dms,lwd = 0.6)
dev.off()


#Supplementary Figure 6.  Enrich2 scores for insertions, deletions, and amino acid changes across the EV-A71 proteome

#Making insertional df compatible for merge
Scores_Insertional_Handle_Fullproteome_formerge <-  Scores_Insertional_Handle_Fullproteome %>% mutate(AminoAcid = "Ins")
Scores_Insertional_Handle_Fullproteome_formerge  <- Scores_Insertional_Handle_Fullproteome_formerge %>% rename(position = colnames(Scores_Insertional_Handle_Fullproteome_formerge)[1])
#Making deletional library df compatible for merge
Scores_Deletions_Fullproteome_formerge  <- Scores_Deletions_Fullproteome %>% rename(AminoAcid = colnames(Scores_Deletions_Fullproteome)[3])
Scores_Deletions_Fullproteome_formerge  <- Scores_Deletions_Fullproteome_formerge %>% rename(position = colnames(Scores_Deletions_Fullproteome_formerge)[1])
#Merging all dataframes
merged_df_indel_DMS <- merge(Scores_Insertional_Handle_Fullproteome_formerge, Scores_Deletions_Fullproteome_formerge, by = c("position","AminoAcid","score"),all = TRUE)
merged_df_indel_DMS <- merge(merged_df_indel_DMS, Fullproteome_P2_DMS_Enrich2_long, by = c("position","AminoAcid","score"),all = TRUE)
merged_df_indel_DMS$AminoAcid=factor(merged_df_indel_DMS$AminoAcid,levels = c("3AAdel","2AAdel","1AAdel","Ins","P","G","Y","W","F","V","L","I","A","T","S","Q","N","M","C","E","D","R","K","H"))
write.table(merged_df_indel_DMS,file="merged_df_indel_DMS.csv",row.names=FALSE,sep=",")
# Convert Position to a factor with breaks
breaks <- seq(0, 2193, by = 250)

# Creating a function to assign labels based on position
get_label <- function(position) {
  sub_labels <- c("A", "B", "C", "D", "E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V")
  sub_index <- findInterval(position, breaks)
  return(sub_labels[sub_index])
}

merged_df_indel_DMS$sub <- sapply(merged_df_indel_DMS$position, get_label)

#Plotting the heatmap of the complete proteome
DMS_indel_Heatmap_fullproteome <- ggplot(merged_df_indel_DMS) +
    geom_tile(aes(position, AminoAcid, fill = score)) +
    scale_fill_gradient2(midpoint = 0,limits = c(-8,8),low = muted("orchid2"), mid = "white",
                         high = muted("aquamarine1"),space = "Lab",na.value = "black") + theme_classic() +
                        facet_wrap(sub ~ ., scales = "free_x",ncol=1) + theme(aspect.ratio = 0.1,axis.text.y = element_text(size = 3))+ 
                        theme(strip.text.x = element_blank())

#Save figure
ggsave(plot = DMS_indel_Heatmap_fullproteome,filename = "DMS_indel_Heatmap_fullproteome.pdf",width = 8,height=11)

DMS_indel_Heatmap_fullproteome

#Figure 3. Impact of altering insertion sequence, deletion length, and amino acid residue substitution on EV-A71 growth

#Creating a factor for Amino acid changes and plotting area plot for different amino acids
DMS_P2_Fitness_Class$AminoAcid=factor(DMS_P2_Fitness_Class$AminoAcid,levels = c("P","G","Y","W","F","V","L","I","A","T","S","Q","N","M","C","E","D","R","K","H"))
DMS_P2_Fitness_Class_Counts_AA <- table(DMS_P2_Fitness_Class$Class,DMS_P2_Fitness_Class$AminoAcid)
chisq.test(DMS_P2_Fitness_Class_Counts_AA)
plot(as.factor(DMS_P2_Fitness_Class$AminoAcid),as.factor(DMS_P2_Fitness_Class$Class),ylim=c(0.5,1.0),col = fitness_class_colors_dms,lwd = 0.6)
pdf("DMS_P2_Fitness_Class_Aminoacids.pdf",width = 4, height = 4)
plot(as.factor(DMS_P2_Fitness_Class$AminoAcid),as.factor(DMS_P2_Fitness_Class$Class),ylim=c(0.5,1.0),col = fitness_class_colors_dms,lwd = 0.6)
dev.off()

#Area plot for deletions of different lengths
Deletion_P2_Fitness_Class <- Scores_Deletions_Fullproteome
Deletion_P2_Fitness_Class  <- na.omit(Deletion_P2_Fitness_Class)
Deletion_P2_Fitness_Class$Class="N"
Deletion_P2_Fitness_Class[Deletion_P2_Fitness_Class$score<2,]$Class="N"
Deletion_P2_Fitness_Class[Deletion_P2_Fitness_Class$score<0,]$Class="L"
Deletion_P2_Fitness_Class[Deletion_P2_Fitness_Class$score>2,]$Class="B"
Deletion_P2_Fitness_Class_Counts <- table(Deletion_P2_Fitness_Class$Class,Deletion_P2_Fitness_Class$dataset)
chisq.test(Deletion_P2_Fitness_Class_Counts)
Deletion_P2_Fitness_Class$Class <- factor(Deletion_P2_Fitness_Class$Class, levels = c("B", "N","L"))  
plot(as.factor(Deletion_P2_Fitness_Class$dataset),as.factor(Deletion_P2_Fitness_Class$Class),ylim=c(0.95,1.0),col = fitness_class_colors_del,lwd = 0.6)
pdf("Deletion_Fitness_Class_3_6_9_bp.pdf",width = 4, height = 4)
plot(as.factor(Deletion_P2_Fitness_Class$dataset),as.factor(Deletion_P2_Fitness_Class$Class),ylim=c(0.95,1.0),col = fitness_class_colors_del,lwd = 0.6)
dev.off()

#Creating a factor for the 5 AA insertion library using inserted amino acid as factor 
Fullproteome_1AA_Enrich2_long$AminoAcid=factor(Fullproteome_1AA_Enrich2_long$AminoAcid,levels = c("P","G","Y","W","F","V","L","I","A","T","S","Q","N","M","C","E","D","R","K","H"))
Fullproteome_1AA_Enrich2_long$sub <- sapply(Fullproteome_1AA_Enrich2_long$insertion, get_label)

#Heat map function for the 5 AA insertion library
create_AA_heatmap <- function(lower_bound, upper_bound) {
  df <- filter(Fullproteome_1AA_Enrich2_long, insertion > lower_bound, insertion < upper_bound)
  
  SingleAA_Heatmap <- ggplot(df) +
    geom_tile(aes(insertion, AminoAcid, fill = score)) +
    scale_fill_gradient2(midpoint = 0,limits = c(-8,8),
                         low = muted("orchid2"),
                         mid = "white",
                         high = muted("aquamarine1"),
                         space = "Lab",
                         na.value = "black") +
    theme_classic() +coord_fixed()
  
  filename <- paste0("SingleAA_Heatmap_", lower_bound, "_to_", upper_bound, ".pdf")
  ggsave(filename, height = 6, width = 6)

  return(SingleAA_Heatmap)
}

#Heatmaps for the sites at N-termini of VP1/2A/3A and VP3
create_AA_heatmap(lower_bound = 399, upper_bound = 451)
create_AA_heatmap(lower_bound = 549, upper_bound = 601)
create_AA_heatmap(lower_bound = 849, upper_bound = 901)
create_AA_heatmap(lower_bound = 1440, upper_bound = 1492)

# Area plot for 5 AA insertion library
SingleAA_fullproteome_Fitness_Class <- Fullproteome_1AA_Enrich2_long
SingleAA_fullproteome_Fitness_Class  <- na.omit(SingleAA_fullproteome_Fitness_Class)
SingleAA_fullproteome_Fitness_Class$Class="N"
SingleAA_fullproteome_Fitness_Class[SingleAA_fullproteome_Fitness_Class$score<2,]$Class="N"
SingleAA_fullproteome_Fitness_Class[SingleAA_fullproteome_Fitness_Class$score<0,]$Class="L"
SingleAA_fullproteome_Fitness_Class[SingleAA_fullproteome_Fitness_Class$score>2,]$Class="B"
SingleAA_fullproteome_Fitness_Class_Counts <- table(SingleAA_fullproteome_Fitness_Class$Class,SingleAA_fullproteome_Fitness_Class$AminoAcid)
chisq.test(SingleAA_fullproteome_Fitness_Class_Counts)
SingleAA_fullproteome_Fitness_Class$Class <- factor(SingleAA_fullproteome_Fitness_Class$Class, levels = c("B", "N","L"))  
plot(as.factor(SingleAA_fullproteome_Fitness_Class$AminoAcid),as.factor(SingleAA_fullproteome_Fitness_Class$Class),ylim=c(0.97,1),col = fitness_class_colors,lwd = 0.6)
pdf("SingleAA_fullproteome_Fitness_Class.pdf",width = 4, height = 4)
plot(as.factor(SingleAA_fullproteome_Fitness_Class$AminoAcid),as.factor(SingleAA_fullproteome_Fitness_Class$Class),ylim=c(0.97,1),col = fitness_class_colors,lwd = 0.6)
dev.off()

#Creating a dataframe with only deletions and then creating bar plots for different deletion lengths at InDel hotspots
merged_df_del_only <- filter(merged_df_indel_DMS,AminoAcid %in% c("1AAdel", "2AAdel","3AAdel"))
merged_df_del_only <- select(merged_df_del_only,position,AminoAcid,score)
merged_df_del_only$AminoAcid=factor(merged_df_del_only$AminoAcid,levels = c("1AAdel", "2AAdel","3AAdel"))

#Bar plot at VP1 N-terminus cleavage site  
Deletion_N_term_VP1  <- filter(merged_df_del_only,position>549,position<601)
vp1_nterminus_del <- ggplot(Deletion_N_term_VP1)+
  geom_bar(aes(factor(position),score,fill=AminoAcid),stat = 'identity',position="dodge")+theme_classic()+scale_fill_manual(values =c("#8E8E8E", "#686868", "#444343"))+
  theme(axis.text.x= element_text(angle=90))+ ylab("Enrich2 score") + xlab("Amino Acid Position")+facet_grid(AminoAcid~.)+theme(legend.position = "none")+coord_fixed(ylim=c(-4,4))
vp1_nterminus_del
#Save figure
ggsave(plot = vp1_nterminus_del,filename = "vp1_nterminus_del_hotspot2.pdf",width = 7,height=4)

#Bar plot at VP1-2A cleavage site  
Deletion_Cleavage_VP1  <- filter(merged_df_del_only,position>849,position<901)
vp1_2a_cleavage <- ggplot(Deletion_Cleavage_VP1)+
  geom_bar(aes(factor(position),score,fill=AminoAcid),stat = 'identity',position="dodge")+theme_classic()+scale_fill_manual(values =c("#8E8E8E", "#686868", "#444343"))+
  theme(axis.text.x= element_text(angle=90))+ ylab("Relative Fitness") + xlab("Amino Acid Position")+facet_grid(AminoAcid~.)+theme(legend.position = "none")+coord_fixed(ylim=c(-4,4))
vp1_2a_cleavage
#Save figure
ggsave(plot = vp1_2a_cleavage,filename = "vp1_2a_cleavage_hotspot3.pdf",width = 7,height=4)

#Bar plot at 3A N-terminus 
Deletion_3A_N_term  <- filter(merged_df_del_only,position>1440,position<1492)
del3Adel <- ggplot(Deletion_3A_N_term)+
  geom_bar(aes(factor(position),score,fill=AminoAcid),stat = 'identity',position="dodge")+theme_classic()+scale_fill_manual(values =c("#8E8E8E", "#686868", "#444343"))+
  theme(axis.text.x= element_text(angle=90))+ ylab("Relative Fitness") + xlab("Amino Acid Position")+facet_grid(AminoAcid~.)+theme(legend.position = "none")+coord_fixed(ylim=c(-4,4))
del3Adel
#Save figure
ggsave(plot = del3Adel,filename = "del3Adel_hotspot4.pdf",width = 7,height=4)


#Figure 4. Structural interpretation of InDel fitness effects for 2A(pro) and 3D(pol) 
#Supplementary Figure 8. Structural interpretation of mutational effects for 3A, 2C, and 3C(pro) 

#Parsing dataframe for exponentiated score as an attribute for chimera for the 2A(pro): Insertional Handle
Insertions_2A <- filter(Scores_Insertional_Handle_Fullproteome,indel>862,indel<1013)
Insertions_2A <- na.omit(Insertions_2A)
Insertions_2A_reindexed <- Insertions_2A %>% mutate(indel=indel-862)
Insertions_2A_reindexed$indel <- paste0(Insertions_2A_reindexed$indel, ".")
Insertions_2A_reindexed$indel <- paste0(Insertions_2A_reindexed$indel, "A")
Insertions_2A_reindexed$indel <- paste(":", Insertions_2A_reindexed$indel, sep = "")
write.csv(Insertions_2A_reindexed,file ="2A_Chimera_Attached_Insertions.csv",)

#Parsing dataframe for exponentiated score as an attribute for chimera for the 2A(pro): Deletion 1 AA 
Deletions_2A <- filter(Scores_Deletions1AA_Fullproteome_1AA,indel>862,indel<1013)
Deletions_2A <- na.omit(Deletions_2A)
Deletions_2A_reindexed <- Deletions_2A %>% mutate(indel=indel-862)
Deletions_2A_reindexed$indel <- paste0(Deletions_2A_reindexed$indel, ".")
Deletions_2A_reindexed$indel <- paste0(Deletions_2A_reindexed$indel, "A")
Deletions_2A_reindexed$indel <- paste(":", Deletions_2A_reindexed$indel, sep = "")
write.csv(Deletions_2A_reindexed,file ="2A_Chimera_Attached_Deletions.csv",)

#Parsing dataframe for mean exponentiated score as an attribute for chimera for the 2A(pro): DMS
DMS_2A <- filter(Fullproteome_P2_DMS_Enrich2_long,position>862,position<1013)
DMS_2A <-  DMS_2A %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
DMS_2A <- DMS_2A %>% replace(is.na(.), "")
DMS_2A_reindexed <- DMS_2A %>% mutate(position=position-862)
DMS_2A_reindexed$position <- paste0(DMS_2A_reindexed$position, ".")
DMS_2A_reindexed$position <- paste0(DMS_2A_reindexed$position, "A")
DMS_2A_reindexed$position <- paste(":", DMS_2A_reindexed$position, sep = "")
write.csv(DMS_2A_reindexed,file ="DMS_2A_reindexed_chimera.csv",)
DMS_2A_reindexed

#STRIDE of 2A(pro)
STRIDEP71_2A <- readr::read_tsv(file.path("STRIDE_3W95.txt"),col_names = TRUE)
STRIDEP71_2A_DMS <- STRIDEP71_2A %>% rename(position = colnames(STRIDEP71_2A)[1])

#Areaplot enrich2 score class and secondary structures for insertions
NS2A_2ndary_insertion <- filter(Scores_Insertional_Handle_Fullproteome,indel>862,indel<1013)
NS2A_2ndary_insertion <- cbind(NS2A_2ndary_insertion, STRIDEP71_2A$Secondary)
colnames(NS2A_2ndary_insertion)[colnames(NS2A_2ndary_insertion) == "STRIDEP71_2A$Secondary"] <- "Secondary"
NS2A_2ndary_insertion <- NS2A_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2A_2ndary_insertion <- NS2A_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS2A_2ndary_insertion_fitness_class  <- na.omit(NS2A_2ndary_insertion)
NS2A_2ndary_insertion_fitness_class$Class="N"
NS2A_2ndary_insertion_fitness_class[NS2A_2ndary_insertion_fitness_class$score<2,]$Class="N"
NS2A_2ndary_insertion_fitness_class[NS2A_2ndary_insertion_fitness_class$score<0,]$Class="L"
NS2A_2ndary_insertion_fitness_class[NS2A_2ndary_insertion_fitness_class$score>2,]$Class="B"
NS2A_2ndary_insertion_fitness_class_Counts <- table(NS2A_2ndary_insertion_fitness_class$Class,NS2A_2ndary_insertion_fitness_class$Secondary)
chisq.test(NS2A_2ndary_insertion_fitness_class_Counts)
NS2A_2ndary_insertion_fitness_class$Class <- factor(NS2A_2ndary_insertion_fitness_class$Class, levels = c("B","N","L"))  
plot(as.factor(NS2A_2ndary_insertion_fitness_class$Secondary),as.factor(NS2A_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.4,1.0),lwd = 0.5)
pdf("insertion_fitness_class_2ndary_2A.pdf",width = 4, height = 4)
plot(as.factor(NS2A_2ndary_insertion_fitness_class$Secondary),as.factor(NS2A_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for deletions
NS2A_2ndary_del <- Scores_Deletions1AA_Fullproteome_1AA
NS2A_2ndary_del <- NS2A_2ndary_del %>% ungroup() 
NS2A_2ndary_del <- filter(NS2A_2ndary_del,indel>862,indel<1013)
NS2A_2ndary_del <- cbind(NS2A_2ndary_del, STRIDEP71_2A$Secondary)
colnames(NS2A_2ndary_del)[colnames(NS2A_2ndary_del) == "STRIDEP71_2A$Secondary"] <- "Secondary"
NS2A_2ndary_del <- NS2A_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2A_2ndary_del <- NS2A_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS2A_2ndary_del_fitness_class  <- na.omit(NS2A_2ndary_del)
NS2A_2ndary_del_fitness_class$Class="N"
NS2A_2ndary_del_fitness_class[NS2A_2ndary_del_fitness_class$score<2,]$Class="N"
NS2A_2ndary_del_fitness_class[NS2A_2ndary_del_fitness_class$score<0,]$Class="L"
NS2A_2ndary_del_fitness_class[NS2A_2ndary_del_fitness_class$score>2,]$Class="B"
NS2A_2ndary_del_fitness_class_Counts <- table(NS2A_2ndary_del_fitness_class$Class,NS2A_2ndary_del_fitness_class$Secondary)
chisq.test(NS2A_2ndary_del_fitness_class_Counts)
NS2A_2ndary_del_fitness_class$Class <- factor(NS2A_2ndary_del_fitness_class$Class, levels = c("B", "N","L"))  
plot(as.factor(NS2A_2ndary_del_fitness_class$Secondary),as.factor(NS2A_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.4,1.0),lwd = 0.5)
pdf("deletions_fitness_class_2ndary_2A.pdf",width = 4, height = 4)
plot(as.factor(NS2A_2ndary_del_fitness_class$Secondary),as.factor(NS2A_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for DMS
NS2A_2ndary_DMS <- Fullproteome_P2_DMS_Enrich2_long
NS2A_2ndary_DMS <- filter(NS2A_2ndary_DMS,position>862,position<1013)
NS2A_2ndary_DMS <- select(NS2A_2ndary_DMS,position,score,AminoAcid)
NS2A_2ndary_DMS <- mutate(NS2A_2ndary_DMS,position=position-862)
NS2A_2ndary_DMS <- merge(NS2A_2ndary_DMS,STRIDEP71_2A_DMS,by="position")
NS2A_2ndary_DMS <- NS2A_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2A_2ndary_DMS <- NS2A_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS2A_2ndary_DMS_fitness_class  <- na.omit(NS2A_2ndary_DMS)
NS2A_2ndary_DMS_fitness_class$Class="N"
NS2A_2ndary_DMS_fitness_class[NS2A_2ndary_DMS_fitness_class$score<2,]$Class="N"
NS2A_2ndary_DMS_fitness_class[NS2A_2ndary_DMS_fitness_class$score<0,]$Class="L"
NS2A_2ndary_DMS_fitness_class[NS2A_2ndary_DMS_fitness_class$score>2,]$Class="B"
NS2A_2ndary_DMS_fitness_class_Counts <- table(NS2A_2ndary_DMS_fitness_class$Class,NS2A_2ndary_DMS_fitness_class$Secondary)
chisq.test(NS2A_2ndary_DMS_fitness_class_Counts)
NS2A_2ndary_DMS_fitness_class$Class <- factor(NS2A_2ndary_DMS_fitness_class$Class, levels = c("B", "N", "L"))  
plot(as.factor(NS2A_2ndary_DMS_fitness_class$Secondary),as.factor(NS2A_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
pdf("DMS_fitness_class_2ndary_2A.pdf",width = 4, height = 4)
plot(as.factor(NS2A_2ndary_DMS_fitness_class$Secondary),as.factor(NS2A_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()


#Function to create areaplots for different amino acids and tolerance at secondary structures in 2A
r2A_mutationtype_function <- function(dfi,AA,plotname) {
df  <- na.omit(dfi)
df  <- filter(df,AminoAcid==AA)
df$Class="N"
df[df$score<2,]$Class="N"
df[df$score<0,]$Class="L"
df[df$score>2,]$Class="B"
df$Class <- factor(df$Class, levels = c("B", "N", "L")) 
df_Counts <- table(df$Class,df$Secondary)
df_Counts
plot(as.factor(df$Secondary),as.factor(df$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
pdf(plotname,width = 4, height = 4)
plot(as.factor(df$Secondary),as.factor(df$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()
    return(df_Counts)
}

#Low alpha helix propensity
r2A_mutationtype_function(NS2A_2ndary_DMS,"G","Glycine_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"P","Proline_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"Y","tyrosine_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"I","IsoLeucine_DMS_2A.pdf")

#High alpha helix propensity
r2A_mutationtype_function(NS2A_2ndary_DMS,"A","Alanine_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"M","met_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"E","Glu_DMS_2A.pdf")
r2A_mutationtype_function(NS2A_2ndary_DMS,"L","Leucine_DMS_2A.pdf")

#Function to Parse dataframe for exponentiated score as an attribute for chimera for 2C
parse2C_function <- function(dataframeoriginal,dataframeindexed,writedataframe,pos) {

dataframeindexed <- filter(dataframeoriginal,indel>1111,indel<1441)
dataframeindexed <- na.omit(dataframeindexed)
dataframeindexed <- mutate(dataframeindexed,score=2^score)
dataframeindexed <- dataframeindexed %>% mutate(indel=indel-1111)
dataframeindexed$indel <- paste0(dataframeindexed$indel, ".")
dataframeindexed$indel <- paste0(dataframeindexed$indel, "A")
dataframeindexed$indel <- paste(":", dataframeindexed$indel, sep = "")
write.csv(dataframeindexed,file = writedataframe,)

}

#Insertion
parse2C_function(Scores_Insertional_Handle_Fullproteome,Insertions_2C_reindexed,"Insertions_2C_reindexed.csv")

#Deletion
parse2C_function(Scores_Deletions1AA_Fullproteome_1AA,Deletions_2C_reindexed,"Deletions_2C_reindexed.csv")

#DMS
Score_DMS_mean_2C <-  Fullproteome_P2_DMS_Enrich2_long %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
Score_DMS_mean_2C <- Score_DMS_mean_2C %>% rename(indel = colnames(Score_DMS_mean_2C)[1])
parse2C_function(Score_DMS_mean_2C,DMS_2C_reindexed,"DMS_2C_reindexed.csv")

#STRIDE of 2C
STRIDEP71_2C <- readr::read_tsv(file.path("STRIDE_5GQ1.txt"),col_names = TRUE)
STRIDEP71_2C <- STRIDEP71_2C %>% complete(Residue = seq(1,329,1))
STRIDEP71_2C_DMS <- STRIDEP71_2C %>% rename(position = colnames(STRIDEP71_2C)[1])

#Areaplot enrich2 score class and secondary structures for insertions
NS2C_2ndary_insertion <- filter(Scores_Insertional_Handle_Fullproteome,indel>1111,indel<1441)
NS2C_2ndary_insertion <- cbind(NS2C_2ndary_insertion, STRIDEP71_2C$Secondary)
colnames(NS2C_2ndary_insertion)[colnames(NS2C_2ndary_insertion) == "STRIDEP71_2C$Secondary"] <- "Secondary"
NS2C_2ndary_insertion <- NS2C_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2C_2ndary_insertion <- NS2C_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS2C_2ndary_insertion_fitness_class  <- na.omit(NS2C_2ndary_insertion)
NS2C_2ndary_insertion_fitness_class$Class="N"
NS2C_2ndary_insertion_fitness_class[NS2C_2ndary_insertion_fitness_class$score<2,]$Class="N"
NS2C_2ndary_insertion_fitness_class[NS2C_2ndary_insertion_fitness_class$score<0,]$Class="L"
#NS2C_2ndary_insertion_fitness_class[NS2C_2ndary_insertion_fitness_class$score>2,]$Class="B"
NS2C_2ndary_insertion_fitness_class_Counts <- table(NS2C_2ndary_insertion_fitness_class$Class,NS2C_2ndary_insertion_fitness_class$Secondary)
chisq.test(NS2C_2ndary_insertion_fitness_class_Counts)
NS2C_2ndary_insertion_fitness_class$Class <- factor(NS2C_2ndary_insertion_fitness_class$Class, levels = c("N","L"))  
plot(as.factor(NS2C_2ndary_insertion_fitness_class$Secondary),as.factor(NS2C_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.8,1.0),lwd = 0.5)
pdf("insertion_fitness_class_2ndary_2C.pdf",width = 4, height = 4)
plot(as.factor(NS2C_2ndary_insertion_fitness_class$Secondary),as.factor(NS2C_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.8,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for deletions
NS2C_2ndary_del <- Scores_Deletions1AA_Fullproteome_1AA
NS2C_2ndary_del <- NS2C_2ndary_del %>% ungroup() 
NS2C_2ndary_del <- filter(NS2C_2ndary_del,indel>1111,indel<1441)
NS2C_2ndary_del <- cbind(NS2C_2ndary_del, STRIDEP71_2C$Secondary)
colnames(NS2C_2ndary_del)[colnames(NS2C_2ndary_del) == "STRIDEP71_2C$Secondary"] <- "Secondary"
NS2C_2ndary_del <- NS2C_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2C_2ndary_del <- NS2C_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS2C_2ndary_del_fitness_class  <- na.omit(NS2C_2ndary_del)
NS2C_2ndary_del_fitness_class$Class="N"
NS2C_2ndary_del_fitness_class[NS2C_2ndary_del_fitness_class$score<2,]$Class="N"
NS2C_2ndary_del_fitness_class[NS2C_2ndary_del_fitness_class$score<0,]$Class="L"
NS2C_2ndary_del_fitness_class[NS2C_2ndary_del_fitness_class$score>2,]$Class="B"
NS2C_2ndary_del_fitness_class_Counts <- table(NS2C_2ndary_del_fitness_class$Class,NS2C_2ndary_del_fitness_class$Secondary)
chisq.test(NS2C_2ndary_del_fitness_class_Counts)
NS2C_2ndary_del_fitness_class$Class <- factor(NS2C_2ndary_del_fitness_class$Class, levels = c("B", "N","L"))  
plot(as.factor(NS2C_2ndary_del_fitness_class$Secondary),as.factor(NS2C_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.8,1.0),lwd = 0.5)
pdf("deletions_fitness_class_2ndary_2C.pdf",width = 4, height = 4)
plot(as.factor(NS2C_2ndary_del_fitness_class$Secondary),as.factor(NS2C_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.8,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for DMS
NS2C_2ndary_DMS <- Fullproteome_P2_DMS_Enrich2_long
NS2C_2ndary_DMS <- filter(NS2C_2ndary_DMS,position>1111,position<1441)
NS2C_2ndary_DMS <- select(NS2C_2ndary_DMS,position,score)
NS2C_2ndary_DMS <- mutate(NS2C_2ndary_DMS,position=position-1111)
NS2C_2ndary_DMS <- merge(NS2C_2ndary_DMS,STRIDEP71_2C_DMS,by="position")
NS2C_2ndary_DMS <- NS2C_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS2C_2ndary_DMS <- NS2C_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
write.csv(NS2C_2ndary_DMS,file ="NS2C_2ndary_DMS.csv",)
NS2C_2ndary_DMS_fitness_class  <- na.omit(NS2C_2ndary_DMS)
NS2C_2ndary_DMS_fitness_class$Class="N"
NS2C_2ndary_DMS_fitness_class[NS2C_2ndary_DMS_fitness_class$score<2,]$Class="N"
NS2C_2ndary_DMS_fitness_class[NS2C_2ndary_DMS_fitness_class$score<0,]$Class="L"
NS2C_2ndary_DMS_fitness_class[NS2C_2ndary_DMS_fitness_class$score>2,]$Class="B"
NS2C_2ndary_DMS_fitness_class_Counts <- table(NS2C_2ndary_DMS_fitness_class$Class,NS2C_2ndary_DMS_fitness_class$Secondary)
chisq.test(NS2C_2ndary_DMS_fitness_class_Counts)
NS2C_2ndary_DMS_fitness_class$Class <- factor(NS2C_2ndary_DMS_fitness_class$Class, levels = c("B", "N", "L"))  
plot(as.factor(NS2C_2ndary_DMS_fitness_class$Secondary),as.factor(NS2C_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.8,1.0),lwd = 0.5)
pdf("DMS_fitness_class_2ndary_2C.pdf",width = 4, height = 4)
plot(as.factor(NS2C_2ndary_DMS_fitness_class$Secondary),as.factor(NS2C_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.8,1.0),lwd = 0.5)
dev.off()

#Function to Parse dataframe for exponentiated score as an attribute for chimera for 3A
parse3A_function <- function(dataframeoriginal,dataframeindexed,writedataframe,pos) {

dataframeindexed <- filter(dataframeoriginal,indel>1440,indel<1527)
dataframeindexed <- na.omit(dataframeindexed)
dataframeindexed <- mutate(dataframeindexed,score=2^score)
dataframeindexed <- dataframeindexed %>% mutate(indel=indel-1440)
dataframeindexed$indel <- paste0(dataframeindexed$indel, ".")
dataframeindexed$indel <- paste0(dataframeindexed$indel, "B")
dataframeindexed$indel <- paste(":", dataframeindexed$indel, sep = "")
write.csv(dataframeindexed,file = writedataframe,)

}

#Insertion
parse3A_function(Scores_Insertional_Handle_Fullproteome,Insertions_3A_reindexed,"Insertions_3A_reindexed.csv")

#Deletion
parse3A_function(Scores_Deletions1AA_Fullproteome_1AA,Deletions_3A_reindexed,"Deletions_3A_reindexed.csv")

#DMS
Score_DMS_mean_3A <-  Fullproteome_P2_DMS_Enrich2_long %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
Score_DMS_mean_3A <- Score_DMS_mean_3A %>% rename(indel = colnames(Score_DMS_mean_3A)[1])
parse3A_function(Score_DMS_mean_3A,DMS_3A_reindexed,"DMS_3A_reindexed.csv")

#Function to Parse dataframe for exponentiated score as an attribute for chimera for 3C

parse3C_function <- function(dataframeoriginal,dataframeindexed,writedataframe,pos) {

dataframeindexed <- filter(dataframeoriginal,indel>1548,indel<1732)
dataframeindexed <- na.omit(dataframeindexed)
dataframeindexed <- mutate(dataframeindexed,score=2^score)
dataframeindexed <- dataframeindexed %>% mutate(indel=indel-1548)
dataframeindexed$indel <- paste0(dataframeindexed$indel, ".")
dataframeindexed$indel <- paste0(dataframeindexed$indel, "A")
dataframeindexed$indel <- paste(":", dataframeindexed$indel, sep = "")
write.csv(dataframeindexed,file = writedataframe,)

}

#Insertion
parse3C_function(Scores_Insertional_Handle_Fullproteome,Insertions_3C_reindexed,"Insertions_3C_reindexed.csv")

#Deletion
parse3C_function(Scores_Deletions1AA_Fullproteome_1AA,Deletions_3C_reindexed,"Deletions_3C_reindexed.csv")

#DMS
Score_DMS_mean_3C <-  Fullproteome_P2_DMS_Enrich2_long %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
Score_DMS_mean_3C <- Score_DMS_mean_3C %>% rename(indel = colnames(Score_DMS_mean_3C)[1])
parse3C_function(Score_DMS_mean_3C,DMS_3C_reindexed,"DMS_3C_reindexed.csv")


#Function to Parse dataframe for exponentiated score as an attribute for chimera for 3D(pol)
parse3D_function <- function(dataframeoriginal,dataframeindexed,writedataframe,pos) {

dataframeindexed <- filter(dataframeoriginal,indel>1731,indel<2194)
dataframeindexed <- na.omit(dataframeindexed)
dataframeindexed <- dataframeindexed %>% mutate(indel=indel-1731)
dataframeindexed <-  mutate(dataframeindexed,score=2^score)
dataframeindexed$indel <- paste0(dataframeindexed$indel, ".")
dataframeindexed$indel <- paste0(dataframeindexed$indel, "A")
dataframeindexed$indel <- paste(":", dataframeindexed$indel, sep = "")
write.csv(dataframeindexed,file = writedataframe,)

}

#Insertion
parse3D_function(Scores_Insertional_Handle_Fullproteome,Insertions_3D_reindexed,"Insertions_3D_reindexed.csv")

#Deletion
parse3D_function(Scores_Deletions1AA_Fullproteome_1AA,Deletions_3D_reindexed,"Deletions_3D_reindexed.csv")

#DMS
Score_DMS_mean_3D <-  Fullproteome_P2_DMS_Enrich2_long %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
Score_DMS_mean_3D <- Score_DMS_mean_3D %>% rename(indel = colnames(Score_DMS_mean_3D)[1])
parse3D_function(Score_DMS_mean_3D,DMS_3D_reindexed,"DMS_3D_reindexed.csv")


#Heatmap showing contact sites of 3D with template RNA

#Reading in txt file with contact sites as identified by chimera
overlap_3D <- read_csv("Overlap_Residues_Selected_filtered.txt")
#Function to create a residue list that contain all sites with the contacts
create_DMS_heatmap_residue_list <- function(residues_df) {
     residues_list <- residues_df[[1]]
        df <- merged_df_indel_DMS
        df <- mutate(df,position=position-1731)
    df <- filter(df, position %in% residues_list)
  df$position <- factor(df$position, levels = residues_list)
       
#Creating heatmap with the contact sites to show mutational robustness
    DMS_Heatmap <- ggplot(df) +
    geom_tile(aes(position, AminoAcid, fill = score)) +
    scale_fill_gradient2(midpoint = 0,limits = c(-8,8),
                         low = muted("hotpink"),
                         mid = "white",
                         high = muted("green"),
                         ,space = "Lab",
                         na.value = "black") +
    theme_classic()+ scale_x_discrete() +
    coord_fixed()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  filename <- paste0("DMS_Heatmap_3D_overlap", paste(collapse = "_"), ".pdf")
  ggsave(filename, height = 4, width = 10)

  return(DMS_Heatmap)
}

#Heatmap for 3D with only the contact sites
create_DMS_heatmap_residue_list(overlap_3D)


#STRIDE of 3D
STRIDEP71_3D <- readr::read_tsv(file.path("STRIDE_3N6L.txt"),col_names = TRUE)
STRIDEP71_3D_DMS <- STRIDEP71_3D %>% rename(position = colnames(STRIDEP71_3D)[1])

#Areaplot enrich2 score class and secondary structures for insertions
NS3D_2ndary_insertion <- filter(Scores_Insertional_Handle_Fullproteome,indel>1731,indel<2194)
NS3D_2ndary_insertion <- cbind(NS3D_2ndary_insertion, STRIDEP71_3D$Secondary)
colnames(NS3D_2ndary_insertion)[colnames(NS3D_2ndary_insertion) == "STRIDEP71_3D$Secondary"] <- "Secondary"
NS3D_2ndary_insertion <- NS3D_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS3D_2ndary_insertion <- NS3D_2ndary_insertion %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS3D_2ndary_insertion_fitness_class  <- na.omit(NS3D_2ndary_insertion)
NS3D_2ndary_insertion_fitness_class$Class="N"
NS3D_2ndary_insertion_fitness_class[NS3D_2ndary_insertion_fitness_class$score<2,]$Class="N"
NS3D_2ndary_insertion_fitness_class[NS3D_2ndary_insertion_fitness_class$score<0,]$Class="L"
#NS3D_2ndary_insertion_fitness_class[NS3D_2ndary_insertion_fitness_class$score>2,]$Class="B"
NS3D_2ndary_insertion_fitness_class_Counts <- table(NS3D_2ndary_insertion_fitness_class$Class,NS3D_2ndary_insertion_fitness_class$Secondary)
chisq.test(NS3D_2ndary_insertion_fitness_class_Counts)
NS3D_2ndary_insertion_fitness_class$Class <- factor(NS3D_2ndary_insertion_fitness_class$Class, levels = c("N","L"))  
plot(as.factor(NS3D_2ndary_insertion_fitness_class$Secondary),as.factor(NS3D_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.4,1.0),lwd = 0.5)
pdf("insertion_fitness_class_2ndary_3D.pdf",width = 4, height = 4)
plot(as.factor(NS3D_2ndary_insertion_fitness_class$Secondary),as.factor(NS3D_2ndary_insertion_fitness_class$Class),col = fitness_class_colors,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for deletions
NS3D_2ndary_del <- Scores_Deletions1AA_Fullproteome_1AA
NS3D_2ndary_del <- NS3D_2ndary_del %>% ungroup() 
NS3D_2ndary_del <- filter(NS3D_2ndary_del,indel>1731,indel<2194)
NS3D_2ndary_del <- cbind(NS3D_2ndary_del, STRIDEP71_3D$Secondary)
colnames(NS3D_2ndary_del)[colnames(NS3D_2ndary_del) == "STRIDEP71_3D$Secondary"] <- "Secondary"
NS3D_2ndary_del <- NS3D_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS3D_2ndary_del <- NS3D_2ndary_del %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
write.csv(NS3D_2ndary_del,file ="NS3D_2ndary_del.csv",)
NS3D_2ndary_del_fitness_class  <- na.omit(NS3D_2ndary_del)
NS3D_2ndary_del_fitness_class$Class="N"
NS3D_2ndary_del_fitness_class[NS3D_2ndary_del_fitness_class$score<2,]$Class="N"
NS3D_2ndary_del_fitness_class[NS3D_2ndary_del_fitness_class$score<0,]$Class="L"
#NS3D_2ndary_del_fitness_class[NS3D_2ndary_del_fitness_class$score>2,]$Class="B"
NS3D_2ndary_del_fitness_class_Counts <- table(NS3D_2ndary_del_fitness_class$Class,NS3D_2ndary_del_fitness_class$Secondary)
chisq.test(NS3D_2ndary_del_fitness_class_Counts)
NS3D_2ndary_del_fitness_class$Class <- factor(NS3D_2ndary_del_fitness_class$Class, levels = c("N","L"))  
plot(as.factor(NS3D_2ndary_del_fitness_class$Secondary),as.factor(NS3D_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.4,1.0),lwd = 0.5)
pdf("deletions_fitness_class_2ndary_3D.pdf",width = 4, height = 4)
plot(as.factor(NS3D_2ndary_del_fitness_class$Secondary),as.factor(NS3D_2ndary_del_fitness_class$Class),col = fitness_class_colors_del,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()

#Areaplot enrich2 score class and secondary structures for DMS
NS3D_2ndary_DMS <- Fullproteome_P2_DMS_Enrich2_long
NS3D_2ndary_DMS <- filter(NS3D_2ndary_DMS,position>1731,position<2194)
NS3D_2ndary_DMS <- select(NS3D_2ndary_DMS,position,score,AminoAcid)
NS3D_2ndary_DMS <- mutate(NS3D_2ndary_DMS,position=position-1731)
NS3D_2ndary_DMS <- merge(NS3D_2ndary_DMS,STRIDEP71_3D_DMS,by="position")
NS3D_2ndary_DMS <- NS3D_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('B', 'T', 'C'), 'L', Secondary))
NS3D_2ndary_DMS <- NS3D_2ndary_DMS %>% mutate(Secondary = ifelse(Secondary %in% c('G', 'H'), 'H', Secondary))
NS3D_2ndary_DMS_fitness_class  <- na.omit(NS3D_2ndary_DMS)
NS3D_2ndary_DMS_fitness_class$Class="N"
NS3D_2ndary_DMS_fitness_class[NS3D_2ndary_DMS_fitness_class$score<2,]$Class="N"
NS3D_2ndary_DMS_fitness_class[NS3D_2ndary_DMS_fitness_class$score<0,]$Class="L"
NS3D_2ndary_DMS_fitness_class[NS3D_2ndary_DMS_fitness_class$score>2,]$Class="B"
NS3D_2ndary_DMS_fitness_class_Counts <- table(NS3D_2ndary_DMS_fitness_class$Class,NS3D_2ndary_DMS_fitness_class$Secondary)
chisq.test(NS3D_2ndary_DMS_fitness_class_Counts)
NS3D_2ndary_DMS_fitness_class$Class <- factor(NS3D_2ndary_DMS_fitness_class$Class, levels = c("B", "N", "L"))  
plot(as.factor(NS3D_2ndary_DMS_fitness_class$Secondary),as.factor(NS3D_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
pdf("DMS_fitness_class_2ndary_3D.pdf",width = 4, height = 4)
plot(as.factor(NS3D_2ndary_DMS_fitness_class$Secondary),as.factor(NS3D_2ndary_DMS_fitness_class$Class),col = fitness_class_colors_dms,ylim=c(0.4,1.0),lwd = 0.5)
dev.off()

#Create areaplots for different amino acids and tolerance at secondary structures in 3D

#Low alphahelix propensity
r2A_mutationtype_function(NS3D_2ndary_DMS,"G","Glycine_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"P","Proline_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"Y","tyrosine_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"I","IsoLeucine_DMS_3D.pdf")

#High alphahelix propensity
r2A_mutationtype_function(NS3D_2ndary_DMS,"A","Alanine_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"M","met_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"E","Glu_DMS_3D.pdf")
r2A_mutationtype_function(NS3D_2ndary_DMS,"L","Leucine_DMS_3D.pdf")

#Figure 5. Structural interpretation of mutational effects of the Capsid

#Creating a dataframe with relative enrichment focused on VP1
Scores_Insertional_Handle_VP1 <- filter(Scores_Insertional_Handle_Fullproteome,indel>565,indel<863)
Scores_Insertional_Handle_VP1 <- mutate(Scores_Insertional_Handle_VP1,score=2^score)
Scores_Deletions_VP1 <- filter(Scores_Deletions_Fullproteome,indel>565,indel<863)
Scores_Deletions_VP1 <- mutate(Scores_Deletions_VP1,score=2^score)
Scores_Deletions1AA_VP1 <- filter(Scores_Deletions_VP1,dataset=="1AAdel")
Scores_Deletions2AA_VP1 <- filter(Scores_Deletions_VP1,dataset=="2AAdel")
Scores_Deletions3AA_VP1 <- filter(Scores_Deletions_VP1,dataset=="3AAdel")
VP1_P2_DMS_Enrich2_long <- filter(Fullproteome_P2_DMS_Enrich2_long,position>565,position<863)
VP1_P2_DMS_Enrich2_long_exponention <- mutate(VP1_P2_DMS_Enrich2_long,score=2^score)


#Function to create geompoint of different mutation types at VP1. VP1 loops are annotated using geomrect
    plot_function_VP1_inserts <- function(data,color,plotname,lim,limm) {
    plot <- ggplot(data)+geom_point(aes(indel,score),color=color,size=0.5)+  
theme_classic() +coord_cartesian(ylim = c(0, 75)) +xlab("Residue Position") +ylab("Relative Enrichment") +  
    geom_rect(aes(xmin = 661, xmax = 670, ymin = 1, ymax = 100),alpha = 0.002,fill = "dimgrey",size = 0) +
    geom_rect(aes(xmin = 706, xmax = 714, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +  
    geom_rect(aes(xmin = 722, xmax = 742, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 748, xmax = 752, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 803, xmax = 811, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0)+
  
  geom_text(data = subset(data, score > 10), 
         aes(indel, score, label = paste(indel)), 
          vjust = -1,hjust=1,size = 1.5)


     ggsave(plot,filename=plotname,width = 4,height=1.5)
        return(plot)

    
    }

# Usage of the plotting function
plot_function_VP1_inserts(Scores_Insertional_Handle_VP1,"#1558a9","Scores_Insertional_Handle_VP1.pdf")
plot_function_VP1_inserts(Scores_Deletions2AA_VP1,"#a91515","Scores_Deletions2AA_VP1.pdf")
plot_function_VP1_inserts(Scores_Deletions3AA_VP1,"#a91515","Scores_Deletions3AA_VP1.pdf")

#Geompoint for 1AA deletion at VP1
Scores_Deletions1AA_VP1_plot <- ggplot(Scores_Deletions1AA_VP1)+geom_point(aes(indel,score),color= "#a91515",size=0.5)+  
theme_classic() +coord_cartesian(ylim = c(0, 120)) +xlab("Residue Position") +ylab("Relative Enrichment") +  
    geom_rect(aes(xmin = 661, xmax = 670, ymin = 1, ymax = 145),alpha = 0.002,fill = "dimgrey",size = 0) +
    geom_rect(aes(xmin = 706, xmax = 714, ymin = 1, ymax = 145), alpha = 0.002, fill = "dimgrey", size = 0) +  
    geom_rect(aes(xmin = 722, xmax = 742, ymin = 1, ymax = 145), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 748, xmax = 752, ymin = 1, ymax = 145), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 803, xmax = 811, ymin = 1, ymax = 145), alpha = 0.002, fill = "dimgrey", size = 0)+
    geom_text(data = subset(Scores_Deletions1AA_VP1, score > 10), 
         aes(indel, score, label = paste(indel)), 
          vjust = -1,hjust=1,size = 1.5)
Scores_Deletions1AA_VP1_plot
#Save figure 
ggsave(Scores_Deletions1AA_VP1_plot,filename="Scores_Deletions1AA_VP1.pdf",width = 4,height=1.5)

#Geompoint for DMS at VP1
DMS_tophit_EV71_VP1 <- ggplot(VP1_P2_DMS_Enrich2_long_exponention)+geom_point(aes(position,score),color="#3c7430",size=0.5)+  
theme_classic() +coord_cartesian(ylim = c(0, 100)) +xlab("Residue Position") +ylab("Relative Enrichment") +  
    geom_rect(aes(xmin = 661, xmax = 670, ymin = 1, ymax = 100),alpha = 0.002,fill = "dimgrey",size = 0) +
    geom_rect(aes(xmin = 706, xmax = 714, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +  
    geom_rect(aes(xmin = 722, xmax = 742, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 748, xmax = 752, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0) +
    geom_rect(aes(xmin = 803, xmax = 811, ymin = 1, ymax = 100), alpha = 0.002, fill = "dimgrey", size = 0)+
  geom_text(data = subset(VP1_P2_DMS_Enrich2_long_exponention, score > 50), 
          aes(position, score, label = paste(position,AminoAcid)), 
           vjust = -2,hjust=1,size = 1.5)
DMS_tophit_EV71_VP1
#Save Figure
ggsave(DMS_tophit_EV71_VP1,filename="DMS_tophit_EV71_VP1.pdf",width = 4,height=1.5)

#Parsing dataframe for exponentiated score as an attribute for chimera for the capsid proteins: Insertional Handle
Insertions_Strucutral_P2 <- filter(Scores_Insertional_Handle_Fullproteome,indel<863)
Insertions_Strucutral_P2 <- Insertions_Strucutral_P2 %>% replace(is.na(.), "")
Insertions_Strucutral_P2 <- mutate(Insertions_Strucutral_P2,score=2^score)
Insertions_Strucutral_P2_reindexed <- Insertions_Strucutral_P2 %>% mutate(indel = ifelse(row_number() %in% c(70:323), indel - 69, indel))
Insertions_Strucutral_P2_reindexed <- Insertions_Strucutral_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(324:565), indel - 323, indel))
Insertions_Strucutral_P2_reindexed <- Insertions_Strucutral_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(566:862), indel - 565, indel))
Insertions_Strucutral_P2_reindexed$indel <- paste0(Insertions_Strucutral_P2_reindexed$indel, ".")
Insertions_Strucutral_P2_reindexed$indel[1:69]  <- paste0(Insertions_Strucutral_P2_reindexed$indel[1:69], "D")
Insertions_Strucutral_P2_reindexed$indel[70:323]  <- paste0(Insertions_Strucutral_P2_reindexed$indel[70:323], "B")
Insertions_Strucutral_P2_reindexed$indel[324:565]  <- paste0(Insertions_Strucutral_P2_reindexed$indel[324:565], "C")
Insertions_Strucutral_P2_reindexed$indel[566:862]  <- paste0(Insertions_Strucutral_P2_reindexed$indel[566:862], "A")
Insertions_Strucutral_P2_reindexed$indel <- paste(":", Insertions_Strucutral_P2_reindexed$indel, sep = "")
write.csv(Insertions_Strucutral_P2_reindexed,file ="Capsid_Chimera_insertions_Attached_p2.csv",)

#Parsing dataframe for exponentiated score as an attribute for chimera for the capsid proteins: 1 AA deletion
Deletions_Strucutral_3bp_P2 <- filter(Scores_Deletions_Fullproteome,dataset=="1AAdel")
Deletions_Strucutral_3bp_P2 <- filter(Deletions_Strucutral_3bp_P2,indel<863)
Deletions_Strucutral_3bp_P2 <- mutate(Deletions_Strucutral_3bp_P2,score=2^score)
Deletions_Strucutral_3bp_P2_reindexed <- Deletions_Strucutral_3bp_P2 %>% mutate(indel = ifelse(row_number() %in% c(70:323), indel - 69, indel))
Deletions_Strucutral_3bp_P2_reindexed <- Deletions_Strucutral_3bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(324:565), indel - 323, indel))
Deletions_Strucutral_3bp_P2_reindexed <- Deletions_Strucutral_3bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(566:862), indel - 565, indel))
Deletions_Strucutral_3bp_P2_reindexed$indel <- paste0(Deletions_Strucutral_3bp_P2_reindexed$indel, ".")
Deletions_Strucutral_3bp_P2_reindexed$indel[1:69]  <- paste0(Deletions_Strucutral_3bp_P2_reindexed$indel[1:69], "D")
Deletions_Strucutral_3bp_P2_reindexed$indel[70:323]  <- paste0(Deletions_Strucutral_3bp_P2_reindexed$indel[70:323], "B")
Deletions_Strucutral_3bp_P2_reindexed$indel[324:565]  <- paste0(Deletions_Strucutral_3bp_P2_reindexed$indel[324:565], "C")
Deletions_Strucutral_3bp_P2_reindexed$indel[566:862]  <- paste0(Deletions_Strucutral_3bp_P2_reindexed$indel[566:862], "A")
Deletions_Strucutral_3bp_P2_reindexed$indel <- paste(":", Deletions_Strucutral_3bp_P2_reindexed$indel, sep = "")
Deletions_Strucutral_3bp_P2_reindexed <- na.omit(Deletions_Strucutral_3bp_P2_reindexed)
write.csv(Deletions_Strucutral_3bp_P2_reindexed,file ="Capsid_Chimera_deletions_Attached_3bp_p2.csv",)

#Parsing dataframe for exponentiated score as an attribute for chimera for the capsid proteins: 2 AA deletion
Deletions_Strucutral_6bp_P2 <- filter(Scores_Deletions_Fullproteome,dataset=="2AAdel")
Deletions_Strucutral_6bp_P2 <- filter(Deletions_Strucutral_6bp_P2,indel<863)
Deletions_Strucutral_6bp_P2 <- mutate(Deletions_Strucutral_6bp_P2,score=2^score)
Deletions_Strucutral_6bp_P2_reindexed <- Deletions_Strucutral_6bp_P2 %>% mutate(indel = ifelse(row_number() %in% c(70:323), indel - 69, indel))
Deletions_Strucutral_6bp_P2_reindexed <- Deletions_Strucutral_6bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(324:565), indel - 323, indel))
Deletions_Strucutral_6bp_P2_reindexed <- Deletions_Strucutral_6bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(566:862), indel - 565, indel))
Deletions_Strucutral_6bp_P2_reindexed$indel <- paste0(Deletions_Strucutral_6bp_P2_reindexed$indel, ".")
Deletions_Strucutral_6bp_P2_reindexed$indel[1:69]  <- paste0(Deletions_Strucutral_6bp_P2_reindexed$indel[1:69], "D")
Deletions_Strucutral_6bp_P2_reindexed$indel[70:323]  <- paste0(Deletions_Strucutral_6bp_P2_reindexed$indel[70:323], "B")
Deletions_Strucutral_6bp_P2_reindexed$indel[324:565]  <- paste0(Deletions_Strucutral_6bp_P2_reindexed$indel[324:565], "C")
Deletions_Strucutral_6bp_P2_reindexed$indel[566:862]  <- paste0(Deletions_Strucutral_6bp_P2_reindexed$indel[566:862], "A")
Deletions_Strucutral_6bp_P2_reindexed$indel <- paste(":", Deletions_Strucutral_6bp_P2_reindexed$indel, sep = "")
Deletions_Strucutral_6bp_P2_reindexed <- na.omit(Deletions_Strucutral_6bp_P2_reindexed)
write.csv(Deletions_Strucutral_6bp_P2_reindexed,file ="Capsid_Chimera_deletions_Attached_6bp_p2.csv",)

#Parsing dataframe for exponentiated score as an attribute for chimera for the capsid proteins: 3 AA deletion
Deletions_Strucutral_9bp_P2 <- filter(Scores_Deletions_Fullproteome,dataset=="3AAdel")
Deletions_Strucutral_9bp_P2 <- filter(Deletions_Strucutral_9bp_P2,indel<863)
Deletions_Strucutral_9bp_P2 <- mutate(Deletions_Strucutral_9bp_P2,score=2^score)
Deletions_Strucutral_9bp_P2_reindexed <- Deletions_Strucutral_9bp_P2 %>% mutate(indel = ifelse(row_number() %in% c(70:323), indel - 69, indel))
Deletions_Strucutral_9bp_P2_reindexed <- Deletions_Strucutral_9bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(324:565), indel - 323, indel))
Deletions_Strucutral_9bp_P2_reindexed <- Deletions_Strucutral_9bp_P2_reindexed %>% mutate(indel = ifelse(row_number() %in% c(566:862), indel - 565, indel))
Deletions_Strucutral_9bp_P2_reindexed$indel <- paste0(Deletions_Strucutral_9bp_P2_reindexed$indel, ".")
Deletions_Strucutral_9bp_P2_reindexed$indel[1:69]  <- paste0(Deletions_Strucutral_9bp_P2_reindexed$indel[1:69], "D")
Deletions_Strucutral_9bp_P2_reindexed$indel[70:323]  <- paste0(Deletions_Strucutral_9bp_P2_reindexed$indel[70:323], "B")
Deletions_Strucutral_9bp_P2_reindexed$indel[324:565]  <- paste0(Deletions_Strucutral_9bp_P2_reindexed$indel[324:565], "C")
Deletions_Strucutral_9bp_P2_reindexed$indel[566:862]  <- paste0(Deletions_Strucutral_9bp_P2_reindexed$indel[566:862], "A")
Deletions_Strucutral_9bp_P2_reindexed$indel <- paste(":", Deletions_Strucutral_9bp_P2_reindexed$indel, sep = "")
Deletions_Strucutral_9bp_P2_reindexed <- na.omit(Deletions_Strucutral_9bp_P2_reindexed)
write.csv(Deletions_Strucutral_9bp_P2_reindexed,file ="Capsid_Chimera_deletions_Attached_9bp_p2.csv",)

#Parsing dataframe for mean exponentiated score as an attribute for chimera for the capsid proteins: DMS
DMS_Strucutral_P2 <- filter(Fullproteome_P2_DMS_Enrich2_long,position<863)
DMS_Strucutral_P2 <- mutate(DMS_Strucutral_P2,score=2^score)
DMS_Strucutral_P2 <-  DMS_Strucutral_P2 %>% group_by(position) %>% summarise(score = mean(score, na.rm = TRUE))
DMS_Strucutral_P2 <- DMS_Strucutral_P2 %>% replace(is.na(.), "")
DMS_Strucutral_P2_reindexed <- DMS_Strucutral_P2 %>% mutate(position = ifelse(row_number() %in% c(70:323), position - 69, position))
DMS_Strucutral_P2_reindexed <- DMS_Strucutral_P2_reindexed %>% mutate(position = ifelse(row_number() %in% c(324:565), position - 323, position))
DMS_Strucutral_P2_reindexed <- DMS_Strucutral_P2_reindexed %>% mutate(position = ifelse(row_number() %in% c(566:862), position - 565, position))
DMS_Strucutral_P2_reindexed$position <- paste0(DMS_Strucutral_P2_reindexed$position, ".")
DMS_Strucutral_P2_reindexed$position[1:69]  <- paste0(DMS_Strucutral_P2_reindexed$position[1:69], "D")
DMS_Strucutral_P2_reindexed$position[70:323]  <- paste0(DMS_Strucutral_P2_reindexed$position[70:323], "B")
DMS_Strucutral_P2_reindexed$position[324:565]  <- paste0(DMS_Strucutral_P2_reindexed$position[324:565], "C")
DMS_Strucutral_P2_reindexed$position[566:862]  <- paste0(DMS_Strucutral_P2_reindexed$position[566:862], "A")
DMS_Strucutral_P2_reindexed$position <- paste(":", DMS_Strucutral_P2_reindexed$position, sep = "")
write.csv(DMS_Strucutral_P2_reindexed,file ="Capsid_Chimera_DMS_Attached_p2.csv",)

#Figure 6. InDels as a contributor of Enterovirus A species evolution

#Function for creating a heatmap
create_DMS_heatmap <- function(lower_bound, upper_bound) {
  df <- filter(merged_df_indel_DMS, position > lower_bound, position < upper_bound)
  
  DMS_Heatmap <- ggplot(df) +
    geom_tile(aes(position, AminoAcid, fill = score)) +
    scale_fill_gradient2(midpoint = 0,limits = c(-8,8),
                         low = muted("hotpink"),
                         mid = "white",
                         high = muted("green"),
                         space = "Lab",
                         na.value = "black") +
    theme_classic() +
    coord_fixed()
  
  filename <- paste0("DMS_Heatmap_", lower_bound, "_to_", upper_bound, ".pdf")
  ggsave(filename, height = 6, width = 6)

  return(DMS_Heatmap)
}

#Create the heatmaps at N-termini of VP1 and 2A
create_DMS_heatmap(lower_bound = 549, upper_bound = 601)
create_DMS_heatmap(lower_bound = 849, upper_bound = 901)

#Supplementary Figure 7. Deep Insertion, Deletion, and Mutational Scanning at InDel hotspot regions 

#Parsing the input competition pooled experiment
Scores_Competition_input_plotting <- Scores_Competition_input
Scores_Competition_input_plotting$position <- gsub("p\\.\\((\\d+).*\\)", "\\1", Scores_Competition_input_plotting$hgvs)
Scores_Competition_input_plotting$Variant <- gsub("p\\.\\(\\d+(.*)\\)", "\\1", Scores_Competition_input_plotting$hgvs)
Scores_Competition_input_plotting$position <- as.numeric(Scores_Competition_input_plotting$position)

#Bar plot for input plasmid library at the N-terminus of VP1
Competition_input_plot_Capsid <-ggplot(Scores_Competition_input_plotting)+
stat_summary_bin(show.legend = F,geom = "bar",aes(position,count),binwidth = 1,fun = function(X) {log10(sum(X))} )+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+ coord_cartesian(xlim = c(577, 612)) 
ggsave("Competition_input_plot_Capsid.pdf",width = 2,height=2)
Competition_input_plot_Capsid

#Bar plot for input plasmid library at the N-terminus of 2A
Competition_input_plot_Replication <-ggplot(Scores_Competition_input_plotting)+
stat_summary_bin(show.legend = F,geom = "bar",aes(position,count),binwidth = 1,fun = function(X) {log10(sum(X))} )+theme_classic()+
ylab(expression(Variant~Count~ (log[10])))+ xlab("Amino acid Position")+ coord_cartesian(xlim = c(863, 910)) 
Competition_input_plot_Replication
ggsave("Competition_input_plot_Replication.pdf",width = 2,height=2)

#Statistics included
#Count positions with count more than 1 
positions_above_1_competition <- sum(Scores_Competition_input_plotting$count > 1)
# Print the result
cat("Number of positions with count more than 1:", positions_above_1_competition/1932, "\n")


#Parsing the enrich2 score dataframes for the competition pooled dataframe at P2
Competition_Pool_Enrich2_Scores$position <- gsub("p\\.\\((\\d+).*\\)", "\\1", Competition_Pool_Enrich2_Scores$hgvs)
Competition_Pool_Enrich2_Scores$Variant <- gsub("p\\.\\(\\d+(.*)\\)", "\\1", Competition_Pool_Enrich2_Scores$hgvs)
Competition_Pool_Enrich2_Scores$position <- as.numeric(Competition_Pool_Enrich2_Scores$position)
Competition_Pool_Enrich2_Scores <- select(Competition_Pool_Enrich2_Scores,position,Variant,score)
Competition_Pool_Enrich2_Scores

#Creating a dataframe with the exponentiated scores
Competition_Pool_Enrich2_Scores$Variant <- factor(Competition_Pool_Enrich2_Scores$Variant, levels = c(levels = "handle","_1AAdel","_2AAdel","_3AAdel","P","G","Y","W","F","V","L","I","A","T","S","Q","N","M","C","E","D","R","K","H"))
Competition_Pool_Enrich2_Scores_exponen <- mutate(Competition_Pool_Enrich2_Scores,score=2^score)
Competition_Pool_Enrich2_Scores_exponen_VP1 <- filter(Competition_Pool_Enrich2_Scores_exponen,position<620)
Competition_Pool_Enrich2_Scores_exponen_2A <- filter(Competition_Pool_Enrich2_Scores_exponen,position>850)

#Boxplot for different variants for the competition pooled experiment
Boxplot_Competition_Variants <- ggplot(Competition_Pool_Enrich2_Scores, aes(y = score, x = Variant))+ylim(-6,6)+geom_boxplot(outlier.shape = NA,width = 0.5)+theme_classic()+labs(y="Median Scores",x="")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Boxplot_Competition_Variants
ggsave("Boxplot_Competition_Variants.pdf", height = 2, width = 4)

#Geompoint for variants at the N-terminus of VP1
DMS_tophit_competitionpool_VP1 <- ggplot(Competition_Pool_Enrich2_Scores_exponen_VP1)+geom_point(aes(position,score),color="#7c7b7b",size=0.5)+  
theme_classic() +coord_cartesian(ylim = c(0, 20)) +xlab("Residue Position") +ylab("Relative Enrichment") +  
    geom_text(data = subset(Competition_Pool_Enrich2_Scores_exponen_VP1, score > 10), 
          aes(position, score, label = paste(position,Variant)), 
           vjust = 0,hjust=0,size = 1) 
#Save figure  
ggsave("DMS_tophit_competitionpool_VP1.pdf", height = 2, width = 2)

#Geompoint for variants at the N-terminus of 2A
DMS_tophit_competitionpool_2A <- ggplot(Competition_Pool_Enrich2_Scores_exponen_2A)+geom_point(aes(position,score),color="#7c7b7b",size=0.5)+  
theme_classic() +coord_cartesian(ylim = c(0, 20)) +xlab("Residue Position") +ylab("Relative Enrichment") +  
    geom_text(data = subset(Competition_Pool_Enrich2_Scores_exponen_2A, score > 10), 
          aes(position, score, label = paste(position,Variant)), 
           vjust = 0,hjust=0,size = 1) 
#Save figure  
ggsave("DMS_tophit_competitionpool_2A.pdf", height = 2, width = 2)

DMS_tophit_competitionpool_VP1
DMS_tophit_competitionpool_2A


#DMFE for deletions
Competition_Pool_Enrich2_Scores_Deletions <- filter(Competition_Pool_Enrich2_Scores,Variant %in% c("_1AAdel", "_2AAdel","_3AAdel"))
histo_competition_pool_Deletions  <- select(Competition_Pool_Enrich2_Scores_Deletions,score)
histo_competition_pool_histo_Deletions <- ggplot(Competition_Pool_Enrich2_Scores_Deletions, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#a91515", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-8,8))
histo_competition_pool_histo_Deletions
ggsave("histo_competition_pool_histo_Deletions.pdf",width = 2,height=1.5)
#Inset for DMFE for deletions
histo_del_sel_nonlethal_competition_pool <- filter(histo_competition_pool_Deletions,score>0)
del_histogram_nonlethal_competition <- ggplot(histo_del_sel_nonlethal_competition_pool, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#a91515", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
del_histogram_nonlethal_competition
#Save figure
ggsave("del_histogram_nonlethal_competition.pdf",width = 2,height=1.5)

#DMFE for Insertions
Competition_Pool_Enrich2_Scores_Insertions <- filter(Competition_Pool_Enrich2_Scores,Variant=="handle")
histo_competition_pool_Insertions  <- select(Competition_Pool_Enrich2_Scores_Insertions,score)
histo_competition_pool_histo_Insertions <- ggplot(histo_competition_pool_Insertions, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#1558a9", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-8,8))
histo_competition_pool_histo_Insertions
ggsave("histo_competition_pool_Insertions.pdf",width = 2,height=1.5)
#Inset for DMFE for Insertions
histo_insertions_sel_nonlethal_competition_pool <- filter(histo_competition_pool_Insertions,score>0)
insertions_histogram_nonlethal_competition <- ggplot(histo_insertions_sel_nonlethal_competition_pool, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#1558a9", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
insertions_histogram_nonlethal_competition
#Save figure
ggsave("insertions_histogram_nonlethal_competition.pdf",width = 2,height=1.5)

#DMFE for DMS
Competition_Pool_Enrich2_Scores_AA <- filter(Competition_Pool_Enrich2_Scores,Variant %in% c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))
histo_competition_pool_AA  <- select(Competition_Pool_Enrich2_Scores_AA,score)
histo_competition_pool_histo_AA <- ggplot(histo_competition_pool_AA, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#3c7430", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(-8,8))
histo_competition_pool_histo_AA
#Save figure
ggsave("histo_competition_pool_histo_AA.pdf",width = 2,height=1.5)
#Inset for DMFE for DMS

histo_AA_sel_nonlethal_competition_pool <- filter(histo_competition_pool_AA,score>0)
AA_histogram_nonlethal_competition <- ggplot(histo_AA_sel_nonlethal_competition_pool, aes(x = score)) +geom_histogram(binwidth = 0.5, fill = "#3c7430", color = "black",size = 0.3) +
labs(x = "Enrich2 Score", y = "Number of Variants") +theme_classic()+coord_cartesian(xlim=c(0,6))
AA_histogram_nonlethal_competition
#Save figure
ggsave("AA_histogram_nonlethal_competition.pdf",width = 2,height=1.5)


#Dataframe for pooled capsid Insertional Handle  
Competition_Pool_Enrich2_Scores_Capsid_Insertion <- filter(Competition_Pool_Enrich2_Scores,position<613)
Competition_Pool_Enrich2_Scores_Capsid_Insertion <- filter(Competition_Pool_Enrich2_Scores_Capsid_Insertion,Variant=="handle")
Competition_Pool_Enrich2_Scores_Capsid_Insertion <- Competition_Pool_Enrich2_Scores_Capsid_Insertion %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Capsid_Insertion)[3])
Competition_Pool_Enrich2_Scores_Capsid_Insertion

#Dataframe for capsid Insertional Handle non-normalized experiment  
Capsid_insertion_N_terminus_VP1 <- filter(merged_df_indel_DMS,position>576 & position<613)
Capsid_insertion_N_terminus_VP1 <- filter(Capsid_insertion_N_terminus_VP1,AminoAcid =="Ins")
Capsid_insertion_N_terminus_VP1$AminoAcid <- gsub("Ins", "handle", Capsid_insertion_N_terminus_VP1$AminoAcid)
Capsid_insertion_N_terminus_VP1 <- Capsid_insertion_N_terminus_VP1 %>% rename(Variant = colnames(Capsid_insertion_N_terminus_VP1)[2])
Capsid_insertion_N_terminus_VP1 <- Capsid_insertion_N_terminus_VP1 %>% rename(score_nonwt = colnames(Capsid_insertion_N_terminus_VP1)[3])
Capsid_insertion_N_terminus_VP1
Merged_Capsid_insertion_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Capsid_Insertion,Capsid_insertion_N_terminus_VP1,all.x=TRUE,all.y=TRUE)
Merged_Capsid_insertion_wt_nonwt

#Scatter plot capsid insertion wildtype vs non wildtype normalized
Scatter_plot_CapsidInsertion_wt_nonwt <- ggplot(Merged_Capsid_insertion_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Enrich2 score", y = "Enrich2 scorewt") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+ geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_CapsidInsertion_wt_nonwt
#Save figure
ggsave(plot = Scatter_plot_CapsidInsertion_wt_nonwt,filename = "Scatter_plot_CapsidInsertion_wt_nonwt.pdf",device = pdf,width = 2,height=2)
#Linear model for wildtype vs nonwildtype experiment
model_Capsid_insertion <- lm(`score_wt`~ `score_nonwt`, data = Merged_Capsid_insertion_wt_nonwt)
summary_model_Capsid_Insertion <- summary(model_Capsid_insertion)
summary_model_Capsid_Insertion

#Dataframe for pooled replication Insertional Handle  
Competition_Pool_Enrich2_Scores_Replication_Insertion <- filter(Competition_Pool_Enrich2_Scores,position>862)
Competition_Pool_Enrich2_Scores_Replication_Insertion <- filter(Competition_Pool_Enrich2_Scores_Replication_Insertion,Variant=="handle")
Competition_Pool_Enrich2_Scores_Replication_Insertion <- Competition_Pool_Enrich2_Scores_Replication_Insertion %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Replication_Insertion)[3])
Competition_Pool_Enrich2_Scores_Replication_Insertion

#Dataframe for replication Insertional Handle non-normalized experiment  
Replication_insertion_N_terminus_2A <- filter(merged_df_indel_DMS,position>862 & position<911)
Replication_insertion_N_terminus_2A <- filter(Replication_insertion_N_terminus_2A,AminoAcid =="Ins")
Replication_insertion_N_terminus_2A$AminoAcid <- gsub("Ins", "handle", Replication_insertion_N_terminus_2A$AminoAcid)
Replication_insertion_N_terminus_2A <- Replication_insertion_N_terminus_2A %>% rename(Variant = colnames(Replication_insertion_N_terminus_2A)[2])
Replication_insertion_N_terminus_2A <- Replication_insertion_N_terminus_2A %>% rename(score_nonwt = colnames(Replication_insertion_N_terminus_2A)[3])
Replication_insertion_N_terminus_2A
Merged_Replication_DMS_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Replication_Insertion,Replication_insertion_N_terminus_2A,all.x=TRUE,all.y=TRUE)

Merged_Replication_DMS_wt_nonwt

#Scatter plot replication insertion wildtype vs non wildtype normalized
Scatter_plot_ReplicationInsertion_wt_nonwt <- ggplot(Merged_Replication_DMS_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#1558a9")+labs(x = "Enrich2 score", y = "Enrich2 scorewt") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_ReplicationInsertion_wt_nonwt
#Save figure
ggsave(plot = Scatter_plot_ReplicationInsertion_wt_nonwt,filename = "Scatter_plot_ReplicationInsertion_wt_nonwt.pdf",device = pdf,width = 2,height=2)
#Linear model for wildtype vs nonwildtype experiment
model_Replication_insertion <- lm(`score_wt`~ `score_nonwt`, data = Merged_Replication_DMS_wt_nonwt)
summary_model_Replication_insertion <- summary(model_Replication_insertion)
summary_model_Replication_insertion


#Dataframe for pooled capsid Deletion   
Competition_Pool_Enrich2_Scores_Capsid_Deletion <- filter(Competition_Pool_Enrich2_Scores,position<613)
Competition_Pool_Enrich2_Scores_Capsid_Deletion <- filter(Competition_Pool_Enrich2_Scores_Capsid_Deletion,Variant %in% c("_1AAdel","_2AAdel","_3AAdel"))
Competition_Pool_Enrich2_Scores_Capsid_Deletion <- Competition_Pool_Enrich2_Scores_Capsid_Deletion %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Capsid_Deletion)[3])
Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant <- gsub("_1AAdel", "1AAdel", Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant <- gsub("_2AAdel", "2AAdel", Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant <- gsub("_3AAdel", "3AAdel", Competition_Pool_Enrich2_Scores_Capsid_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Capsid_Deletion

#Dataframe for capsid Deletion non-normalized experiment  
Capsid_deletion_N_terminus_2A <- filter(merged_df_indel_DMS,position>576 & position<613)
Capsid_deletion_N_terminus_2A <- filter(Capsid_deletion_N_terminus_2A,AminoAcid %in% c("1AAdel","2AAdel","3AAdel"))
Capsid_deletion_N_terminus_2A <- Capsid_deletion_N_terminus_2A %>% rename(Variant = colnames(Capsid_deletion_N_terminus_2A)[2])
Capsid_deletion_N_terminus_2A <- Capsid_deletion_N_terminus_2A %>% rename(score_nonwt = colnames(Capsid_deletion_N_terminus_2A)[3])
Capsid_deletion_N_terminus_2A

Merged_Capsid_deletion_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Capsid_Deletion,Capsid_deletion_N_terminus_2A,all.x=TRUE,all.y=TRUE)
Merged_Capsid_deletion_wt_nonwt

#Scatter plot capsid deletion wildtype vs non wildtype normalized
Scatter_plot_CapsidDeletion_wt_nonwt <- ggplot(Merged_Capsid_deletion_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Enrich2 Score", y = "Enrich2 Scorewt") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+ geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_CapsidDeletion_wt_nonwt
#Save figure
ggsave(plot = Scatter_plot_CapsidDeletion_wt_nonwt,filename = "Scatter_plot_CapsidDeletion_wt_nonwt.pdf",device = pdf,width = 2,height=2)
#Linear model for wildtype vs nonwildtype experiment
model_Capsid_deletion <- lm(`score_wt`~ `score_nonwt`, data = Merged_Capsid_deletion_wt_nonwt)
summary_model_Capsid_deletion <- summary(model_Capsid_deletion)
summary_model_Capsid_deletion

#Dataframe for pooled replication Deletion   
Competition_Pool_Enrich2_Scores_Replication_Deletion <- filter(Competition_Pool_Enrich2_Scores,position>862)
Competition_Pool_Enrich2_Scores_Replication_Deletion <- filter(Competition_Pool_Enrich2_Scores_Replication_Deletion,Variant %in% c("_1AAdel","_2AAdel","_3AAdel"))
Competition_Pool_Enrich2_Scores_Replication_Deletion <- Competition_Pool_Enrich2_Scores_Replication_Deletion %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Replication_Deletion)[3])
Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant <- gsub("_1AAdel", "1AAdel", Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant <- gsub("_2AAdel", "2AAdel", Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant <- gsub("_3AAdel", "3AAdel", Competition_Pool_Enrich2_Scores_Replication_Deletion$Variant)
Competition_Pool_Enrich2_Scores_Replication_Deletion

#Dataframe for replication Deletion non-normalized experiment  
Replication_deletion_N_terminus_2A <- filter(merged_df_indel_DMS,position>862 & position<911)
Replication_deletion_N_terminus_2A <- filter(Replication_deletion_N_terminus_2A,AminoAcid %in% c("1AAdel","2AAdel","3AAdel"))
Replication_deletion_N_terminus_2A <- Replication_deletion_N_terminus_2A %>% rename(Variant = colnames(Replication_deletion_N_terminus_2A)[2])
Replication_deletion_N_terminus_2A <- Replication_deletion_N_terminus_2A %>% rename(score_nonwt = colnames(Replication_deletion_N_terminus_2A)[3])
Replication_deletion_N_terminus_2A

Merged_Replication_deletion_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Replication_Deletion,Replication_deletion_N_terminus_2A,all.x=TRUE,all.y=TRUE)
Merged_Replication_deletion_wt_nonwt

#Scatter plot replication deletion wildtype vs non wildtype normalized
Scatter_plot_ReplicationDeletion_wt_nonwt <- ggplot(Merged_Replication_deletion_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#a91515")+labs(x = "Enrich2 Score", y = "Enrich2 Scorewt") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+ geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_ReplicationDeletion_wt_nonwt
#Save figute
ggsave(plot = Scatter_plot_ReplicationDeletion_wt_nonwt,filename = "Scatter_plot_ReplicationDeletion_wt_nonwt.pdf",device = pdf,width = 2,height=2)
#Linear model for wildtype vs nonwildtype experiment
model_Replication_deletion <- lm(`score_wt`~ `score_nonwt`, data = Merged_Replication_deletion_wt_nonwt)
summary_model_Replication_deletion <- summary(model_Replication_deletion)
summary_model_Replication_deletion


#Dataframe for pooled capsid DMS   
Competition_Pool_Enrich2_Scores_Capsid_DMS <- filter(Competition_Pool_Enrich2_Scores,position<613)
Competition_Pool_Enrich2_Scores_Capsid_DMS <- filter(Competition_Pool_Enrich2_Scores_Capsid_DMS,Variant %in% c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))
Competition_Pool_Enrich2_Scores_Capsid_DMS <- Competition_Pool_Enrich2_Scores_Capsid_DMS %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Capsid_DMS)[3])
Competition_Pool_Enrich2_Scores_Capsid_DMS

#Dataframe for capsid DMS non-normalized experiment  
Capsid_DMS_N_terminus_VP1 <- filter(merged_df_indel_DMS,position>576 & position<613)
Capsid_DMS_N_terminus_VP1 <- filter(Capsid_DMS_N_terminus_VP1,AminoAcid %in% c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))
Capsid_DMS_N_terminus_VP1 <- Capsid_DMS_N_terminus_VP1 %>% rename(Variant = colnames(Capsid_DMS_N_terminus_VP1)[2])
Capsid_DMS_N_terminus_VP1 <- Capsid_DMS_N_terminus_VP1 %>% rename(score_nonwt = colnames(Capsid_DMS_N_terminus_VP1)[3])
Capsid_DMS_N_terminus_VP1
Merged_Capsid_DMS_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Capsid_DMS,Capsid_DMS_N_terminus_VP1,all.x=TRUE,all.y=TRUE)
Merged_Capsid_DMS_wt_nonwt

#Scatter plot capsid DMS wildtype vs non wildtype normalized
Scatter_plot_CapsidDMS_wt_nonwt <- ggplot(Merged_Capsid_DMS_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Enrich2 score", y = "Enrich2 score WT normalized") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+ geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_CapsidDMS_wt_nonwt
#Save figure
ggsave(plot = Scatter_plot_CapsidDMS_wt_nonwt,filename = "Scatter_plot_CapsidDMS_wt_nonwt.pdf",device = pdf,width = 2,height=2)
#Linear model for wildtype vs nonwildtype experiment
model_Capsid_DMS <- lm(`score_wt`~ `score_nonwt`, data = Merged_Capsid_DMS_wt_nonwt)
summary_model_Capsid_DMS <- summary(model_Capsid_DMS)
summary_model_Capsid_DMS


#Dataframe for pooled replication DMS   
Competition_Pool_Enrich2_Scores_Replicationproteins_DMS <- filter(Competition_Pool_Enrich2_Scores,position>862)
Competition_Pool_Enrich2_Scores_Replicationproteins_DMS <- filter(Competition_Pool_Enrich2_Scores_Replicationproteins_DMS,Variant %in% c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))
Competition_Pool_Enrich2_Scores_Replicationproteins_DMS <- Competition_Pool_Enrich2_Scores_Replicationproteins_DMS %>% rename(score_wt = colnames(Competition_Pool_Enrich2_Scores_Replicationproteins_DMS)[3])
Competition_Pool_Enrich2_Scores_Replicationproteins_DMS

#Dataframe for replication DMS non-normalized experiment  
Replication_DMS_N_terminus_2A <- filter(merged_df_indel_DMS,position>862 & position<911)
Replication_DMS_N_terminus_2A <- filter(Replication_DMS_N_terminus_2A,AminoAcid %in% c("Y","W","V","T","S","R","Q","P","N","M","L","K","I","H","G","F","E","D","C","A"))
Replication_DMS_N_terminus_2A <- Replication_DMS_N_terminus_2A %>% rename(Variant = colnames(Replication_DMS_N_terminus_2A)[2])
Replication_DMS_N_terminus_2A <- Replication_DMS_N_terminus_2A %>% rename(score_nonwt = colnames(Replication_DMS_N_terminus_2A)[3])
Replication_DMS_N_terminus_2A
Merged_Replication_DMS_wt_nonwt  <- merge(Competition_Pool_Enrich2_Scores_Replicationproteins_DMS,Replication_DMS_N_terminus_2A,all.x=TRUE,all.y=TRUE)
Merged_Replication_DMS_wt_nonwt

#Scatter plot replication DMS wildtype vs non wildtype normalized
Scatter_plot_ReplicationDMS_wt_nonwt <- ggplot(Merged_Replication_DMS_wt_nonwt, aes(x = `score_nonwt`, y = `score_wt`)) + geom_point(size = 0.05,color="#3c7430")+labs(x = "Enrich2 Score", y = "Enrich2 Scorewt") +theme_classic()+
coord_fixed(xlim=c(-8,8),ylim=c(-8,8))+ geom_smooth(method="lm",color="black",lwd=0.3)
Scatter_plot_ReplicationDMS_wt_nonwt
#Save figure
ggsave(plot = Scatter_plot_ReplicationDMS_wt_nonwt,filename = "Scatter_plot_ReplicationDMS_wt_nonwt.pdf",device = pdf,width = 2,height=2)
model_Replication_DMS <- lm(`score_wt`~ `score_nonwt`, data = Merged_Replication_DMS_wt_nonwt)
#Linear model for wildtype vs nonwildtype experiment
summary_model_Replication_DMS <- summary(model_Replication_DMS)
summary_model_Replication_DMS



