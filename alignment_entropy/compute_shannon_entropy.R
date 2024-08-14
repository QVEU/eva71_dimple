compute_shannon_entropy <- function(consensus_matrix){
  
  entropy_aa <- function(i){return(-(i*logb(i, base = 2)))}
  
  frequency_matrix <- consensus_matrix/sum(consensus_matrix[[1]])
  entropy_by_aa <- entropy_aa(frequency_matrix)
  entropy_by_site_int <- colSums(entropy_by_aa, na.rm=TRUE)
  entropy_by_site <- data.frame(entropy_by_site_int) %>% mutate(entropy_by_site_int, aa_pos=row_number(), entropy=entropy_by_site_int) %>% select(aa_pos, entropy)
  return(entropy_by_site)
}