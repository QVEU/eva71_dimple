# EVA71 DIMPLE Alignment Analysis and Codon Count Processing Scripts

## Overview

This repository contains three scripts used to process data and generate figures in Bakhache et. al. 2024:

- dms_processing: raw codon counts from mutational scanning experiments and scripts for filtering codon counts based on codons introduced via synthetic biology. 
- alignment_entropy: sequence files and scripts for generating Shannon entropy from curated EV-A71 sequence set and comparing this to mean mutational fitness effect as measured via mutational scanning.
- figure_generation_script: script generates all the figures displayed in the paper.
- *Note:* Other python scripts used for read mapping of the sequenced libraries can be found in [InDel_Toolkit](https://github.com/QVEU/InDel_Toolkit) 
### Prerequisites

- R 4.4.3+
- R packages: tidyverse, zoo, Biostrings, DescTools, ggplot2, tidyr, ggpubr, dplyr, ggridges, ineq, RColorBrewer, stringr, gglorenz, readr, scales

### Instructions 
#### alignment_entropy:

Running instructions: Run R script in r or rstudio from the `alignment_entropy` directory.

e.g.:
```r
setwd alignment_entropy/
```

input_files descriptions:
- EV71_4643_Features.csv: EV71_4643_Features.csv: Contains list of features of EV-A71 genome used to plot boundaries of viral proteins.
- Enterovirus_A71_Curated_AA_Aligned.fasta: 482 complete amino sequences of EV-A71 genomes.
- merged_df_indel_DMS.csv: Dataframe including enrich2 scores for insertion, deletion, and amino acid variants. Position references the position of the variant. AminoAcid references the type of variant introduced. Score is the enrich2 score for this (position, variant).

EVA71_EntropyAnalysis.R: R script to calculate Shannon entropy and to plot entropy calculations along with the mean fitness effects derived from the deep mutational scanning experiments. 
compute_shannon_entropy.R: Function to calculate Shannon entropy from consensus matrix

output_files: pdf image of entropy/DMS MFE trace for EV-A71.

#### dms_processing:

##### Running instructions: Run R script in r or rstudio from the `dms_processing` directory.
##### Note: we recommend running interactively or in a [conda](http://anaconda.org) environment with tidyverse and r installed. 

```{bash}
conda create -n newenv
conda activate newenv
conda install -c conda-forge r --solver=classic
conda install -c conda-forge r-tidyverse --solver=classic
Rscript DMS_Processing_Workflow.R
```


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

#### figure_generation_script:

Running instructions: Run Rscript in r or rstudio from the `figure_generation_script` directory.

Contains R scripts required to generate figures for deep mutational scanning experiments. 

Executable_R_script_and_dataframes: “Data_Generation_Markdown_InDel_Manuscript.r” is the R script used for data analysis and generation of figures in the paper. R script can be launched directly from the directory where all the files required are present: /Dryad_Repository_InDel_Paper_v3/3_R_analysis_scripts/Figure_Generationscript/Executable_R_script_and_dataframes/. 
The input data required for running the R script are: Dataframes corresponding to enrich2 outputs, STRIDE analysis on viral protein structures, and formatted data for chimera analysis. The variables in these dataframes include: indel, insertion, position, or residue: position of variant, score: relative enrichment values or enrich2 scores, dataset: dataset where the values were measured, n: count, Seq, or Amino Acid: type of variant, hgvs: nomenclature used to name the variants, SE: standard error, epsilon: change in the standard error after the last iteration of the random-effects model, Secondary: secondary structure assignment using STRIDE, Residue_Overlap_3D: contains residues in 3D(pol) that interact with the template RNA.

HTML_version: HTML version of the executable R script where the script can be seen along with the figures generated.

Output_Figures: Contains all output figures generated by the R script.

### License

This project is licensed under the MIT License - see the LICENSE.md file for details.

### Corresponding Author Contact

Email: patrick.dolan@nih.gov

### Acknowledgements
We wish to acknowledge Willow Coyote-Maestas, originator of the DIMPLE pipeline; see Macdonald, C.B., Nedrud, D., Grimes, P.R. et al. DIMPLE: deep insertion, deletion, and missense mutation libraries for exploring protein variation in evolution, disease, and biology. Genome Biol 24, 36 (2023). https://doi.org/10.1186/s13059-023-02880-6.
