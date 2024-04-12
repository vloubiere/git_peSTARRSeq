setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# Update exp data (dropbox folder) and import metadata ----
meta <- read_xlsx("Rdata/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)

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
cores <- 2
mem <- 8
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
meta$cmd <- NULL

# SAVE ----
saveRDS(meta,
        "Rdata/metadata.rds")