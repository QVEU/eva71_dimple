#load dependencies
library(tidyverse)

codonFilter <- function(oligos_programmed_, codon_counts_, wildtype_, start_end_){
  #The arguments specify paths to the dataframes, not dataframes themselves
  
  #These will form the arguments for the overall function
  oligos_programmed <-read_tsv(file(oligos_programmed_))
  codon_counts <- read_tsv(file(codon_counts_))
  wildtype <- readLines(file(wildtype_))
  start_ <- start_end_[1]
  end_ <- start_end_[2]
  proteome_size <- end_ - start_ + 1
  
  #Extract the oligo names (id) and sequence (seq) from oligos_programmed
  oligo_processor <- function(oligos_programmed_){
    oligos_programmed_dirty <- readLines(oligos_programmed_) %>%
      str_c(collapse="") %>% 
      str_split(">") %>% .[[1]] %>% tibble(sequenceid=.) %>% 
      #The first row is empty, which causes problems for separate_wider
      slice(2:nrow(.)) %>% 
      #Separating (id) and (seq)
      separate_wider_regex(
        sequenceid, c(id="^.+_DMS\\-\\d+_[A-z][a-z]{1,2}\\d+[A-z][a-z]{1,2}", seq=".+$"))
    
    oligos_programmed_clean <- oligos_programmed_dirty %>% 
      separate_wider_regex(
        id, c("^.+_DMS\\-\\d+_[A-z][a-z]{1,2}", pos="\\d+", "[A-z][a-z]{1,2}")) %>% 
      mutate(pos=(as.numeric(pos)-1)) %>% 
      separate_wider_regex(
        seq, c("^[A-Z]+[a-z]+", codon="[A-Z]+", "[a-z]+[A-Z]+$"))
    
    return(oligos_programmed_clean)
  }
  
  codons_programmed <- oligo_processor(oligos_programmed_)
  
  #Dealing with and Incorporating the WT sequence
  wildtype_clean <- wildtype[2:length(wildtype)] %>% str_c(collapse="")
  codons_wt <- str_extract_all(wildtype_clean, "[a-z]{3}")[[1]] %>% str_to_upper() %>% tibble(codon=.) %>% mutate(pos=row_number())
  codons_programmed_wt <- codons_programmed %>% add_row(codons_wt)
  
  #Setting up a dataframe to house TRUE/FALSE for each codon
  codon_writer <- function(){
    bases <- c('A', 'T', 'C', 'G')
    codons <- vector(mode="character", length=64)
    counter <- 0
    for(i in 1:4){
      for(j in 1:4){
        for(k in 1:4){
          counter <- counter + 1
          base1 <- bases[i]
          base2 <- bases[j]
          base3 <- bases[k]
          codons[counter] <-
            paste(collapse="", c(base1, base2, base3))
        }
      }
    }
    return(codons)
  }
  
  filler <- function(codons_programmed_wt, proteome_size){
    standard_codons <- codon_writer()
    programmed_empty <- matrix(FALSE, ncol = length(standard_codons),
                               nrow = proteome_size) %>% data.frame()
    colnames(programmed_empty) <- standard_codons
    
    for(i in 1:length(standard_codons)){
      codon_i <- standard_codons[i]
      which_codon <- codons_programmed_wt %>%
        group_by(pos) %>%
        filter(codon==codon_i) %>%
        pull(pos)
      which_codon
      programmed_empty[i][which_codon, ] <- TRUE}
    programmed_full <- programmed_empty
    return(programmed_full)
  }
  
  #By now, a dataframe with TRUE/FALSE for each codon 
  codons_twist_wt <- filler(codons_programmed_wt, proteome_size)
  columns <- c(colnames(codons_twist_wt))
  codons_DMS <- codon_counts %>% select(columns) %>% slice(start_:end_)
  
  filter_codon <- function(codons_twist_wt, codons_DMS){
    
    not_programmed <- apply(codons_twist_wt, 2, function(x) which(x==FALSE)) 
    
    apply(
      codons_DMS, 2,
      function(x) replace(x, not_programmed$x, NA)
    )
    
    for(i in 1:length(columns)){
      codon <- columns[i]
      codons_DMS[[codon]][not_programmed[[codon]]] <- 0
    }
    codons_filtered <- codons_DMS
    return(codons_filtered)
  }
  
  codons_filtered_DMS <- filter_codon(codons_twist_wt, codons_DMS)
  #codons_filtered_AA <- transmute(codons_filtered_DMS, Ala=GCA+GCT+GCG+GCC, Arg=CGT+CGC+CGA+CGG+AGA+AGG, Asn=AAT+AAC, Asp=GAT+GAC, Cys=TGT+TGC, Gln=CAA+CAG, Glu=GAA+GAG, Gly=GGA+GGC+GGT+GGG, His=CAT+CAC, Ile=ATT+ATC+ATA, Lys=AAA+AAG, Met=ATG, Phe=TTT+TTC, Pro=CCT+CCC+CCA+CCG, Leu=TTA+TTG+CTT+CTG+CTA+CTC, Ser=TCT+TCA+TCC+TCG+AGT+AGC, Thr=ACT+ACC+ACA+ACG, Trp=TGG, Tyr=TAT+TAC, Val=GTT+GTC+GTA+GTG)
  codons_filtered_AA <- transmute(codons_filtered_DMS, A=GCA+GCT+GCG+GCC, R=CGT+CGC+CGA+CGG+AGA+AGG, N=AAT+AAC, D=GAT+GAC, C=TGT+TGC, Q=CAA+CAG, E=GAA+GAG, G=GGA+GGC+GGT+GGG, H=CAT+CAC, I=ATT+ATC+ATA, K=AAA+AAG, M=ATG, `F`=TTT+TTC, P=CCT+CCC+CCA+CCG, L=TTA+TTG+CTT+CTG+CTA+CTC, S=TCT+TCA+TCC+TCG+AGT+AGC, `T`=ACT+ACC+ACA+ACG, W=TGG, Y=TAT+TAC, V=GTT+GTC+GTA+GTG)
  invisible(codons_filtered_AA)
}


#Initializing codon files from Twist, Wildtype nucleotide sequences, and start/end positions (amino acid)
oligos_programmed_capsid <- "input_files/All_Oligos_Capsid.fasta"
wildtype_capsid <- "input_files/bbfree_2-746-3331.fasta"

start_end_capsid <- c(1,862)

oligos_programmed_replication <- ("input_files/All_Oligos_Replication.fasta")
wildtype_replication <- ("input_files/bbfree-3332-7324.fasta")
start_end_replication <- c(863, 2193)

#Running codonFilter for all of the data.
capsid_P0_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_input_DMS.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P1_A_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P1_RepA_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P1_B_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P1_RepB_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P1_C_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P1_RepC_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P2_A_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P2_RepA_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P2_B_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P2_RepB_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)
capsid_P2_C_filtered <- codonFilter(oligos_programmed_capsid, "input_files/codon_counts/Capsid_P2_RepC_DMS_Nextseq.codonCounts", wildtype_capsid, start_end_capsid)

#Runs codonFilter on all the replication proteins codon read files
replication_P0_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_input_DMS_DMS.codonCounts", wildtype_replication, start_end_replication)
replication_P1_A_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P1_RepA_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)
replication_P1_B_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P1_RepB_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)
replication_P1_C_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P1_RepC_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)
replication_P2_A_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P2_RepA_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)
replication_P2_B_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P2_RepB_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)
replication_P2_C_filtered <- codonFilter(oligos_programmed_replication, "input_files/codon_counts/Replication_P2_RepC_DMS_Nextseq.codonCounts", wildtype_replication, start_end_replication)

#Removing wildtype codons (because Enrich2 analysis will be performed without a wild type reference)
remove_wt <- function(codon_counts){
  max_aa <- apply(codon_counts, MARGIN=1, which.max)
  for(i in 1:length(max_aa)){
    codon_counts[i, max_aa[i]] <- 0
  }
  return(codon_counts)
}

#Removing WT codons from the Capsid libraries
capsid_P0_filtered_wt_removed <- remove_wt(codon_counts=capsid_P0_filtered)
capsid_P1_A_filtered_wt_removed <- remove_wt(codon_counts=capsid_P1_A_filtered)
capsid_P1_B_filtered_wt_removed <- remove_wt(codon_counts=capsid_P1_B_filtered)
capsid_P1_C_filtered_wt_removed <- remove_wt(codon_counts=capsid_P1_C_filtered)
capsid_P2_A_filtered_wt_removed <- remove_wt(codon_counts=capsid_P2_A_filtered)
capsid_P2_B_filtered_wt_removed <- remove_wt(codon_counts=capsid_P2_B_filtered)
capsid_P2_C_filtered_wt_removed <- remove_wt(codon_counts=capsid_P2_C_filtered)

#removing WT codons from the Replication libraries
replication_P0_filtered_wt_removed <- remove_wt(codon_counts=replication_P0_filtered)
replication_P1_A_filtered_wt_removed <- remove_wt(codon_counts=replication_P1_A_filtered)
replication_P1_B_filtered_wt_removed <- remove_wt(codon_counts=replication_P1_B_filtered)
replication_P1_C_filtered_wt_removed <- remove_wt(codon_counts=replication_P1_C_filtered)
replication_P2_A_filtered_wt_removed <- remove_wt(codon_counts=replication_P2_A_filtered)
replication_P2_B_filtered_wt_removed <- remove_wt(codon_counts=replication_P2_B_filtered)
replication_P2_C_filtered_wt_removed <- remove_wt(codon_counts=replication_P2_C_filtered)

#Converts filtered capsid counts to hvgs format for Enrich2
capsid_saver <- function(filtered_counts, write_path){
  capsid_aa <- read_csv("input_files/capsid_aa.fasta", skip=1, col_names=FALSE)$X1 %>% str_split("")
  capsid_aa_one <- capsid_aa[[1]]
  filtered_long <- mutate(filtered_counts, pos=row_number()) %>% select(pos, everything()) %>%
    pivot_longer(cols=(A:V))
  filtered_hgvs <- transmute(filtered_long, hgvs=str_c("p.(", capsid_aa_one[pos], pos, name, ")"), count=value)
  write_tsv(filtered_hgvs, write_path)
}

replication_saver <- function(filtered_counts, write_path){
  replication_aa <- read_csv("input_files/capsid_aa.fasta", skip=1, col_names=FALSE)$X1 %>% str_split("")
  replication_aa_one <- replication_aa[[1]]
  filtered_long <- mutate(filtered_counts, pos=row_number()) %>% select(pos, everything()) %>%
    pivot_longer(cols=(A:V))
  filtered_hgvs <- transmute(filtered_long, hgvs=str_c("p.(", replication_aa_one[pos], pos, name, ")"), count=value)
  write_tsv(filtered_hgvs, write_path)
}

capsid_saver(capsid_P0_filtered_wt_removed, "output_files/capsid_P0_hgvs.tsv")
capsid_saver(capsid_P1_A_filtered_wt_removed, "output_files/capsid_P1_A_hgvs.tsv")
capsid_saver(capsid_P1_B_filtered_wt_removed, "output_files/capsid_P1_B_hgvs.tsv")
capsid_saver(capsid_P1_C_filtered_wt_removed, "output_files/capsid_P1_C_hgvs.tsv")
capsid_saver(capsid_P2_A_filtered_wt_removed, "output_files/capsid_P2_A_hgvs.tsv")
capsid_saver(capsid_P2_B_filtered_wt_removed, "output_files/capsid_P2_B_hgvs.tsv")
capsid_saver(capsid_P2_C_filtered_wt_removed, "output_files/capsid_P2_C_hgvs.tsv")

replication_saver(replication_P0_filtered_wt_removed, "output_files/replication_P0_hgvs.tsv")
replication_saver(replication_P1_A_filtered_wt_removed, "output_files/replication_P1_A_hgvs.tsv")
replication_saver(replication_P1_B_filtered_wt_removed, "output_files/replication_P1_B_hgvs.tsv")
replication_saver(replication_P1_C_filtered_wt_removed, "output_files/replication_P1_C_hgvs.tsv")
replication_saver(replication_P2_A_filtered_wt_removed, "output_files/replication_P2_A_hgvs.tsv")
replication_saver(replication_P2_B_filtered_wt_removed, "output_files/replication_P2_B_hgvs.tsv")
replication_saver(replication_P2_C_filtered_wt_removed, "output_files/replication_P2_C_hgvs.tsv")
