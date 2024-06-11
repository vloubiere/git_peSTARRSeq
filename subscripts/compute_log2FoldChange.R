setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)
require(data.table)

# Import metadata ----
meta <- readRDS("Rdata/metadata_processed.rds")

# Cat umi counts files ----
meta <- meta[, .(umi_counts= paste0(unique(umi_counts), collapse= ",")), .(screen, vllib, cdition)]
meta <- merge(meta[cdition=="screen"],
              meta[cdition=="input", .(vllib, umi_counts)],
              by= "vllib",
              suffixes= c("_screen", "_input"))

# Retrieve ID patterns and counts file ----
meta[, sublibRegexprL:= fcase(vllib=="vllib002", "_.*_", # DSCP large WT library
                              vllib=="vllib006", "_A_", # DSCP long spacer
                              vllib=="vllib015", "_A_", # DSCP focused WT library
                              vllib=="vllib016", "_A_", # RpS12 focused WT library
                              vllib=="vllib029", "_A_", # DSCP mutant library
                              vllib=="vllib030", "_B_")] # DSCP DHS sites library
meta[, sublibRegexprR:= fcase(vllib=="vllib002", "_.*_",
                              vllib=="vllib006", "_B_",
                              vllib=="vllib015", "_A_",
                              vllib=="vllib016", "_A_",
                              vllib=="vllib029", "_A_",
                              vllib=="vllib030", "_B_")]
meta[, sublibRegexprControlPairs:= fcase(vllib=="vllib002", "^control.*__control",
                                         vllib=="vllib006", "^control.*__control",
                                         vllib=="vllib015", "^control.*__control",
                                         vllib=="vllib016", "^control.*__control",
                                         vllib=="vllib029", "^control_WT.*__control_WT",
                                         vllib=="vllib030", "^control.*__control")]

# Input and output file names ----
meta[, dds_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/dds/", screen, "_DESeq2.dds")]
meta[, FC_DESeq2_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/", screen, "_FC_DESeq2.rds")]
meta[, FC_ratio_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/", screen, "_FC_ratio.rds")]

# Compute log2FoldChange using DESeq2 ----
meta[, DESeq2_cmd:= {
  if(!file.exists(FC_DESeq2_file))
  {
    # [required] 1/ A list of comma-separated list of input counts (each file considered as a separate replicate) \n
    # [required] 2/ A list of comma-separated list of screen counts (each file considered as a separate replicate) \n
    # [required] 3/ The method to be used. One of 'DESeq2' or 'ratio' (the later one tolerates one replicate) \n
    # [required] 4/ Regexpr to be applied to 5' (Left) IDs \n
    # [required] 5/ Regexpr to be applied to 3' (Right) IDs \n
    # [required] 6/ Output folder\n
    # [required] 7/ Output file prefix\n"
    # paste("module load build-env/2020; module load r/3.6.2-foss-2018b;
    #       Rscript /groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R",
    paste("/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
          "/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R",
          umi_counts_input,
          umi_counts_screen,
          "DESeq2",
          sublibRegexprL,
          sublibRegexprR,
          sublibRegexprControlPairs,
          dirname(FC_DESeq2_file),
          gsub("_FC_DESeq2.rds", "", basename(FC_DESeq2_file)))
  }else
    as.character(NA)
}, .(umi_counts_screen, umi_counts_input, sublibRegexprL, sublibRegexprR, sublibRegexprControlPairs, FC_DESeq2_file)]

# Compute log2FoldChange using ratio ----
meta[, ratio_cmd:= {
  if(!file.exists(FC_ratio_file))
  {
    # [required] 1/ A list of comma-separated list of input counts (each file considered as a separate replicate) \n
    # [required] 2/ A list of comma-separated list of screen counts (each file considered as a separate replicate) \n
    # [required] 3/ The method to be used. One of 'DESeq2' or 'ratio' (the later one tolerates one replicate) \n
    # [required] 4/ Regexpr to be applied to 5' (Left) IDs \n
    # [required] 5/ Regexpr to be applied to 3' (Right) IDs \n
    # [required] 6/ Output folder\n
    # [required] 7/ Output file prefix\n"
    # paste("module load build-env/2020; module load r/3.6.2-foss-2018b;
    #       Rscript /groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R",
    paste("/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
          "/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R",
          umi_counts_input,
          umi_counts_screen,
          "ratio",
          sublibRegexprL,
          sublibRegexprR,
          sublibRegexprControlPairs,
          dirname(FC_DESeq2_file),
          gsub("_FC_ratio.rds", "", basename(FC_ratio_file)))
  }else
    as.character(NA)
}, .(umi_counts_screen, umi_counts_input, sublibRegexprL, sublibRegexprR, FC_ratio_file)]

# Submit ----
cores <- 4
mem <- 16
run <- na.omit(melt(meta, "screen", patterns("_cmd$"), value.name = "cmd"))
if(nrow(run))
{
  run[, {
    .c <- vl_bsub(cmd, 
            cores= cores, 
            m = mem, 
            name = screen, 
            t = '1-00:00:00',
            o= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            execute = T)
  }, .(cmd, screen)]
}
