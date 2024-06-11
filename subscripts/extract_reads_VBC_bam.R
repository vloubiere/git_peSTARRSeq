setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# Update exp data (dropbox folder) and import metadata ----
meta <- read_xlsx("Rdata/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)

# generate fq names ----
meta[, fq1:= paste0("/scratch/stark/vloubiere/fastq/", gsub(".bam$", "", basename(BAM_path)), "_", i5, "_1.fq.gz")]
meta[, fq2:= paste0("/scratch/stark/vloubiere/fastq/", gsub(".bam$", "", basename(BAM_path)), "_", i5, "_2.fq.gz")]
meta[, fq1_trimmed:= gsub("_1.fq.gz", "_trimmed_1.fq.gz", fq1)]
meta[, fq2_trimmed:= gsub("_2.fq.gz", "_trimmed_2.fq.gz", fq2)]

# SAVE ----
saveRDS(meta,
        "Rdata/metadata.rds")

# Extract reads ----
meta[, extract_cmd:= {
  if(any(!file.exists(c(fq1, fq2))))
  {
    revComp <- Biostrings::DNAStringSet(gsub(".*-(.*)", "\\1", i5))
    revComp <- Biostrings::reverseComplement(revComp)
    revComp <- as.character(revComp)
    paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; /groups/stark/software-all/shell/demultiplexPE.sh -i ", 
           BAM_path,
           " -o ", gsub("_1.fq.gz", "", fq1), # Output prefix
           " -b ", "\"", revComp, "\"", " -u TRUE")
  }
}, .(BAM_path, i5, fq1, fq2)]

# Trim reads ----
meta[, trim_cmd:= {
  if(any(!file.exists(c(fq1_trimmed, fq2_trimmed))))
  {
    if(type=="pe-STARR-Seq")
    {
      fw <- paste("module load build-env/2020;",
                  "module load cutadapt/1.18-foss-2018b-python-2.7.15;",
                  "cutadapt -g GCGTCTCTCACCG",
                  "-m", 10, # Remove trimmed reads shorter than 10bp
                  "-o", fq1_trimmed, fq1)
      rev <- paste("cutadapt -g TCGATCCTAGG -g GCTTCAAAAGC -g TAAGGCACAGG",
                   "-m", 10, # Remove trimmed reads shorter than 10bp
                   "-o", fq2_trimmed, fq2)
    }else if (type=="rev-pe-STARR-Seq")
    {
      fw <- paste("module load build-env/2020;",
                  "module load cutadapt/1.18-foss-2018b-python-2.7.15;",
                  "cutadapt -g TCGATCCTAGG -g GCTTCAAAAGC -g TAAGGCACAGG",
                  "-m", 10, # Remove trimmed reads shorter than 10bp
                  "-o", fq1_trimmed, fq1)
      rev <- paste("cutadapt -g GCGTCTCTCACCG",
                   "-m", 10, # Remove trimmed reads shorter than 10bp
                   "-o", fq2_trimmed, fq2)
    }
    paste(c(fw, rev), collapse = "; ")
  }
}, .(fq1, fq1_trimmed, fq2, fq2_trimmed, type)]

# Submit ----
cores <- 2
mem <- 8
run <- meta[, .(cmd= paste0(unique(na.omit(unlist(.SD))), collapse = "; ")), .(fq1, fq2, fq1_trimmed, fq2_trimmed), .SDcols= patterns("_cmd$")]
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
