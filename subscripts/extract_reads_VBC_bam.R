setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# Update exp data (dropbox folder) and import metadata ----
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- names(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

# Select used libraries ----
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
meta <- meta[, .(vllib, type, CP, library, cdition, rep= DESeq2_pseudo_rep, i5, BAM_path, Read_Type, Platform)]

# generate fq names and command lines ----
meta[, fq1:= paste0("/scratch/stark/vloubiere/fastq/", gsub(".bam$", "", basename(BAM_path)), "_", i5, "_1.fq.gz")]
meta[, fq2:= paste0("/scratch/stark/vloubiere/fastq/", gsub(".bam$", "", basename(BAM_path)), "_", i5, "_2.fq.gz")]
meta[, cmd:= {
  if(any(!file.exists(c(fq1, fq2))))
  {
    revComp <- Biostrings::DNAStringSet(gsub(".*-(.*)", "\\1", i5))
    revComp <- Biostrings::reverseComplement(revComp)
    revComp <- as.character(revComp)
    paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; /groups/stark/software-all/shell/demultiplexPE.sh -i ", 
           BAM_path,
           " -o ", gsub("_1.fq.gz", "", fq1), # Output prefix
           " -b ", "\"", revComp, "\"", " -u TRUE")
  }else
    as.character(NA)
}, .(BAM_path, i5, fq1, fq2)]

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

# Clean and SAVE ----
meta$cmd <- NULL
# Add vllib029 screen rep3 to rep1 ----
meta[vllib=="vllib029" & rep=="rep3", rep:= "rep1"]
# Selected libraries ---
meta[, paper:= vllib %in% c("vllib002", "vllib029", "vllib015", "vllib016")]
# Add sublib ID patterns ----
saveRDS(meta,
        "Rdata/metadata.rds")

