# Developmental and housekeeping transcriptional programs display distinct modes of enhancer-enhancer cooperativity in Drosophila

This repository contains all custom scripts supporting the conclusions from "Developmental and housekeeping transcriptional programs display distinct modes of enhancer-enhancer cooperativity in Drosophila" (V. Loubiere, Bernardo P. De Almeida, M. Pagani, A. Stark). doi: https://doi.org/10.1101/2023.10.10.561770

The "main_analyses_paper.R" file lists all the scripts that are needed to reproduce the results presented in the article.

1/ The "pSTARR-Seq PROCESSING" section contains the code that was used to:
  - Align sequencing reads
  - Compute the activity of individual enhancers and enhancer pairs
  - Fit linear and LASSO regression

2/ The "Figures" section lists, for each main and supplementary figure, the different scripts that were used to generate the corresponding panels. Concerning raw and processed data files, they can be retrieved as follows: 
  - Raw sequencing data have been deposited on GEO and will be made publicly available upon publication (GSE245033).
  - Processed data can be found in supplementary tables 4,7,9,10 (https://doi.org/10.1101/2023.10.10.561770). 
      Supplementary table 4: WT oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib002"" in the different scripts.
      Supplementary table 7: Mutated oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib029"" in the different scripts.
      Supplementary table 9: Focused oligo pool STARR-Seq activities, housekeeping Core Promoter (Rps12). Referred to as "vllib016"" in the different scripts.
      Supplementary table 10: Focused oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib015"" in the different scripts.

For any reasonable further request, please contact vincent.loubiere@imp.ac.at.
