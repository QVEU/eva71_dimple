library(tidyverse)
library(zoo)
library(Biostrings)
library(DescTools)

entropy_aa <- function(i){return(-(i*logb(i, base = 2)))}
compute_shannon_entropy <- function(consensus_matrix){
  frequency_matrix <- consensus_matrix/sum(consensus_matrix[[1]])
  entropy_by_aa <- entropy_aa(frequency_matrix)
  entropy_by_site_int <- colSums(entropy_by_aa, na.rm=TRUE)
  entropy_by_site <- data.frame(entropy_by_site_int) %>% mutate(entropy_by_site_int, aa_pos=row_number(), entropy=entropy_by_site_int) %>% select(aa_pos, entropy)
  return(entropy_by_site)
}

EV71 <- readAAStringSet("Enterovirus_A71_Curated_AA_Aligned.fasta")

mfe <- read.csv("merged_df_indel_DMS.csv") %>% mutate(exp_mfe=2^score)

EV_Features <- read_csv("EV71_4643_Features.csv")

#Identifying which sequence contains the insertion and where it is located, and trimming these residues from the alignment
consensus_matrix_test <- consensusMatrix(EV71) %>% as.data.frame() 

which(letterFrequency(EV71, "-")==0)
matchPattern("---", EV71[[1]])
consensus_matrix_test[1197:1199]

consensus_matrix_EV71 <- consensusMatrix(EV71[c(1:46, 48:482)]) %>% as.data.frame() 

#Trimming the consensus matrix to remove non-AA and gap character options
consensus_matrix_EV71_trimmed <- select(consensus_matrix_EV71, !V1197:V1199)[c(1:20, 28),]
entropy_EV71 <- compute_shannon_entropy(consensus_matrix_EV71_trimmed)

#For plotting: make a sliding scale window for plotting, then plot
entropy_EV71_with_roll <- mutate(entropy_EV71, roll_entropy=rollmean(entropy, 21, fill=NA))

#Processing MFE data
median_mean_mfe <- filter(mfe, AminoAcid!="3AAdel" & AminoAcid!="2AAdel" & AminoAcid!="1AAdel" & AminoAcid !="Ins") %>% group_by(position) %>% summarize(median=median(score, na.rm=TRUE), mean=mean(score, na.rm=TRUE))
median_mean_mfe_EV71_with_roll <- mutate(median_mean_mfe, roll_mfe=(2^rollmean(mean, 21, fill=NA)), aa_pos=position)


#Plotting entropy against MFE
EV71_entropy_plot <- ggplot(entropy_EV71_with_roll) +
  geom_line(aes(x=aa_pos, y=roll_entropy), color="#000000FF") +
  geom_line(data=median_mean_mfe_EV71_with_roll, aes(x=aa_pos, y=roll_mfe), color="#3C7430") +
  
  geom_rect(data = EV_Features, aes(xmin=(Start-746)/3,xmax=(End-746)/3, ymin = 0, ymax = 8), alpha = rep(times = 1, c(0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0.1, 0))) +
  scale_y_continuous(sec.axis = dup_axis(name = "Mean MFE")) +
  theme_classic() +   
  coord_cartesian(xlim=c(0,2193),ylim=c(0,1))+
  #geom_vline(xintercept=862, linetype='dashed', color="#00000080")+
  xlab("Residue Position")+
  ylab("Entropy") +
  theme(axis.line.y.right = element_line(color = "#3C7430"), 
        axis.ticks.y.right = element_line(color = "#3C7430"),
        axis.text.y.right = element_text(color = "#3C7430"), 
        axis.title.y.right = element_text(color = "#3C7430")
  )

EV71_entropy_plot


ggsave("EV71_entropy.pdf", plot=EV71_entropy_plot, width=9, height=1.5)
