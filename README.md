# EVA71 DIMPLE Alignment Analysis and Codon Count Processing Scripts

## Overview

This repository contains two scripts used to process data and generate figures in Bakhache et. al. 2024:

- dms_processing: raw codon counts from mutational scanning experiments and scripts for filtering codon counts based on codons introduced via synthetic biology. 
- alignment_entropy: sequence files and scripts for generating Shannon entropy from curated EV-A71 sequence set and comparing this to mean mutational fitness effect as measured via mutational scanning.

### Prerequisites

- R 4.4.3+
- R packages: tidyverse, zoo, Biostrings, DescTools

### Instructions 
#### alignment_entropy:

Running instructions: Set the working directory to this folder before running the script.

input_files descriptions:
- EV71_4643_Features.csv: EV71_4643_Features.csv: Contains list of features of EV-A71 genome used to plot boundaries of viral proteins.
- Enterovirus_A71_Curated_AA_Aligned.fasta: 482 complete amino sequences of EV-A71 genomes.
- merged_df_indel_DMS.csv: Dataframe including enrich2 scores for insertion, deletion, and amino acid variants. Position references the position of the variant. AminoAcid references the type of variant introduced. Score is the enrich2 score for this (position, variant).

EVA71_EntropyAnalysis.R: R script to calculate Shannon entropy and to plot entropy calculations along with the mean fitness effects derived from the deep mutational scanning experiments. 
compute_shannon_entropy.R: Function to calculate Shannon entropy from consensus matrix

output_files: pdf image of entropy/DMS MFE trace for EV-A71.

#### dms_processing

Running instructions: Set the working directory to this folder before running the script.

input_files descriptions:
- All_Oligos_Capsid.fasta: Oligopools containing every possible amino acid change in the capsid proteins.
- All_Oligos_Replication.fasta: Oligopools containing every possible amino acid change in the replication proteins.
- bbfree_2-746-3331.fasta: Nucleotide sequence of the EV-A71 capsid proteins region.
- bbfree-3332-7324.fasta: Nucleotide sequence of the EV-A71 replication proteins region.
- capsid_aa.fasta: Coding sequence (amino acid) of the EV-A71 capsid proteins region.
- replication_aa.fasta: Coding sequence (amino acid) of the EV-A71 replication proteins region.
- codon_counts: Files in this directory GATK/Analyze Saturation Mutagenesis codon output for all amino acid scanning experiments. These files contain input (P0) and passaged (P1 and P2) condon count data for capsid and replication protein experiments. Each column represents a codon, each row is a position in the viral genome, and the measurements represent the sequencing counts at each position for a certain codon.

DMS_Processing_Workflow.R: R script that will filter for only designed codon variants from GATK/Analyze Saturation Mutagenesis codon output. 
output_files: Contains all filtered amino acid counts for all amino acid scanning experiments. hgvs references the nomenclature used to name the variants, and count represents the sequencing count for each variant. 

### License

This project is licensed under the MIT License - see the LICENSE.md file for details.

### Corresponding Author Contact

Email: patrick.dolan@nih.gov

### Acknowledgements
We wish to acknowledge Willow Coyote-Maestas, originator of the DIMPLE pipeline; see Macdonald, C.B., Nedrud, D., Grimes, P.R. et al. DIMPLE: deep insertion, deletion, and missense mutation libraries for exploring protein variation in evolution, disease, and biology. Genome Biol 24, 36 (2023). https://doi.org/10.1186/s13059-023-02880-6.
