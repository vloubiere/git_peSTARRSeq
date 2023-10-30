# Developmental and housekeeping transcriptional programs display distinct modes of enhancer-enhancer cooperativity in Drosophila

This repository contains all custom scripts supporting the conclusions from "Developmental and housekeeping transcriptional programs display distinct modes of enhancer-enhancer cooperativity in Drosophila" (V. Loubiere, Bernardo P. De Almeida, M. Pagani, A. Stark). doi: https://doi.org/10.1101/2023.10.10.561770

System requirements:
  - All the custom scripts generated for this study were written in R (version 4.2.0) using the R studio IDE (https://www.R-project.org/).
  - No special hardware should be required.

Installation guide:
  - R and RStudio can be downloaded at https://posit.co/download/rstudio-desktop/. Installation time is approximately 20min.

Instructions for use:

The "main_analyses_paper.R" file lists all the scripts that are needed to process the raw sequencing data, analyse the results and plot the panels presented in the manuscript.

1/ The "pSTARR-Seq PROCESSING" section contains the code that was used to:
  - Align raw sequencing data (which were deposited on GEO (GSE245033) and will be made publicly available upon publication)
  - Compute the activity of individual enhancers and enhancer pairs
  - Fit linear and LASSO regressions
  - Total processing time can be up to 1h on a "normal" desktop computer

2/ The "Figures" section lists, for each main and supplementary figure, the scripts that were used for the downstream analyses of processed data and to generate the corresponding figure panels. Input processed data files containing the STARR-Seq activities of all tested pairs can be retrieved in supplementary tables 4,7,9,10 (https://www.biorxiv.org/content/biorxiv/early/2023/10/12/2023.10.10.561770/DC1/embed/media-1.xlsx?download=true).Detailed description:
  - Supplementary table 4: WT oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib002" in the different scripts.
  - Supplementary table 7: Mutated oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib029" in the different scripts.
  - Supplementary table 9: Focused oligo pool STARR-Seq activities, housekeeping Core Promoter (Rps12). Referred to as "vllib016" in the different scripts.
  - Supplementary table 10: Focused oligo pool STARR-Seq activities, developmental Core Promoter (DSCP). Referred to as "vllib015" in the different scripts.
Running times should be less than one minute for each individual panel.

Further details about the different oligo pools and their processing can be found in the "Methods" section of the manuscript (https://doi.org/10.1101/2023.10.10.561770).

For any reasonable further request, please contact vincent.loubiere@imp.ac.at.
