setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)
require(data.table)

# Update exp data (dropbox folder) and import metadata ----
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- names(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

# Select usable libraries ----
sel <- c("vllib002", # Large dev library
         "vllib006", # Reverse pSTARR-Seq
         "vllib015", # Restricted dev library
         "vllib016", # Restricted hk library
         "vllib025", # highHk CP
         "vllib026", # lowHk CP
         "vllib027", # lowDev CP
         "vllib028", # highDev CP
         "vllib029", # mutant lib
         "vllib030") # DHS lib
meta <- as.data.table(meta)[(DESeq2) & vllib %in% sel]


# Retrieve ID patterns and counts file ----
meta[, sublibRegexprL:= fcase(vllib=="vllib006", "_A_",
                              vllib=="vllib015", "_A_",
                              vllib=="vllib016", "_A_",
                              vllib=="vllib029", "_A_",
                              vllib=="vllib030", "_B_",
                              default = "_.*_")]
meta[, sublibRegexprR:= fcase(vllib=="vllib006", "_B_",
                              vllib=="vllib015", "_A_",
                              vllib=="vllib016", "_A_",
                              vllib=="vllib029", "_A_",
                              vllib=="vllib030", "_B_",
                              default = "_.*_")]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", vllib, "_", cdition, "_", DESeq2_pseudo_rep, ".txt")]

# Input and output file names ----
meta[, dds_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/dds/", vllib, ".dds"), vllib]
meta[, FC_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/", vllib, "_DESeq2.rds"), vllib]

# For vllib029 and vllib030, reads from scren rep 3 were randomly split and added to rep 1 and 2
# script -> "git_peSTARRSeq/subscripts/vllib029_vllib030_pseudo_reps.R"
rm <- which(meta[, vllib %in% c("vllib029", "vllib030") & rep=="rep3"])
meta <- meta[-rm]
meta[vllib %in% c("vllib029", "vllib030") & cdition=="screen", umi_counts:= gsub("rep", "pseudoRep", umi_counts)]

# Compute log2FoldChange using DESeq2 ----
meta[, cmd:= {
  if(!file.exists(FC_file))
  {
    # [required] 1/ A list of comma-separated list of input counts (each file considered as a separate replicate) \n
    # [required] 2/ A list of comma-separated list of screen counts (each file considered as a separate replicate) \n
    # [required] 3/ The method to be used. One of 'DESeq2' or 'ratio' (the later one tolerates one replicate) \n
    # [required] 4/ Regexpr to be applied to 5' (Left) IDs \n
    # [required] 5/ Regexpr to be applied to 3' (Right) IDs \n
    # [required] 6/ Output DESeq2 dds file (.dds)\n
    # [required] 7/ FoldChange output (.rds)\n")
    paste("module load build-env/2020; module load r/3.6.2-foss-2018b; Rscript /groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R",
          paste0(unique(umi_counts[cdition=="input"]), collapse = ","),
          paste0(unique(umi_counts[cdition=="screen"]), collapse = ","),
          "DESeq2",
          sublibRegexprL,
          sublibRegexprR,
          dds_file,
          FC_file)
  }else
    as.character(NA)
}, FC_file]

# Submit ----
cores <- 8
mem <- 32
run <- meta[!is.na(cmd)]
if(nrow(run))
{
  run[, {
    vl_bsub(cmd, 
            cores= cores, 
            m = mem, 
            name = "vlloub", 
            t = '1-00:00:00',
            o= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/")
  }, cmd]
}