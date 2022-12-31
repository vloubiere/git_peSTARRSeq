setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

#--------------------------------------------------------------#
# Update exp data
# Fetch dropbox folder containing my metadata and update local files
#--------------------------------------------------------------#
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)[(DESeq2)]
cols <- colnames(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

#-------------------------------------------------------------#
# vllib002
#-------------------------------------------------------------#
lib <- "vllib002"
type <-  "pe-STARR-Seq"
index <-  "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8"
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)[(DESeq2) & vllib=="vllib002", .(BAM_path, i5, cdition, rep= DESeq2_pseudo_rep)]
tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
fwrite(meta, tmp)
sublibRegexprL <- ".*"
sublibRegexprR <- ".*"
Rcmd <- paste("module load r/4.1.2-foss-2021b; Rscript", 
              normalizePath("git_peSTARRSeq/subscripts/pipeline_peSTARRSeq.R"),
              lib, type, index, tmp, sublibRegexprL, sublibRegexprR)
vl_bsub(Rcmd, 
        cores = 12, 
        name = lib, 
        m= 100, 
        o= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
        t = '2-00:00:00',
        e= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/")

#-------------------------------------------------------------#
# vllib006
#-------------------------------------------------------------#
lib <- "vllib006"
type <-  "rev-pe-STARR-Seq"
index <-  "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8"
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)[(DESeq2) & vllib=="vllib006", .(BAM_path, i5, cdition, rep= DESeq2_pseudo_rep)]
tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
fwrite(meta, tmp)
sublibRegexprL <- "_A_"
sublibRegexprR <- "_B_"
Rcmd <- paste("module load r/4.1.2-foss-2021b; Rscript", 
              normalizePath("git_peSTARRSeq/subscripts/pipeline_peSTARRSeq.R"),
              lib, type, index, tmp, sublibRegexprL, sublibRegexprR)
vl_bsub(Rcmd, 
        cores = 8, 
        name = lib, 
        m= 40, 
        o= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
        t = '2-00:00:00',
        e= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/")

#-------------------------------------------------------------#
# Twist12
#-------------------------------------------------------------#
for(lib in c("vllib015", "vllib016"))
{
  type <-  "pe-STARR-Seq"
  index <-  "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12"
  meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
  meta <- as.data.table(meta)[(DESeq2) & vllib==lib, .(BAM_path, i5, cdition, rep= DESeq2_pseudo_rep)]
  tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
  fwrite(meta, tmp)
  sublibRegexprL <- ".*"
  sublibRegexprR <- ".*"
  Rcmd <- paste("module load r/4.1.2-foss-2021b; Rscript", 
                normalizePath("git_peSTARRSeq/subscripts/pipeline_peSTARRSeq.R"),
                lib, type, index, tmp, sublibRegexprL, sublibRegexprR)
  vl_bsub(Rcmd, 
          cores = 8, 
          name = lib, 
          m= 40, 
          o= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
          t = '2-00:00:00',
          e= "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/")
}
