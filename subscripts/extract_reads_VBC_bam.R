setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)
require(data.table)


# Import metadata ----
if(F) # Update exp data (dropbox folder) 
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)[(DESeq2)]
cols <- colnames(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

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
            o= "/groups/stark/vloubiere/projects/ORFTRAP_1/logs/",
            e= "/groups/stark/vloubiere/projects/ORFTRAP_1/logs/")
  }, cmd]
}