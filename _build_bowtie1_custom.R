setwd("/groups/stark/vloubiere/genomes/")
source("/groups/stark/vloubiere/scripts/shell/R_wrappers/my_wrappers.R")
require(data.table)
require(Biostrings)

# Load lib 
lib <- as.data.table(readRDS("/groups/stark/vloubiere/projects/0002_library_design/Rdata/TWIST008_library_design_112019/vl_library_112019.rds"))

# Load templates for template switching control
ts <- fread("/groups/stark/vloubiere/projects/0002_library_design/Rdata/TWIST008_library_design_112019/template_switching_templates.txt")
# Add extra nt for consistency with the TWIST library
ts[, oligo_full_sequence:= paste0("G", oligo_full_sequence, "CACCG")]

# Bind library and TS
lib <- as.data.table(rbind(lib, ts))

# Save fasta
sequences <- DNAStringSet(lib$oligo_full_sequence)
names(sequences) <- lib$ID_vl
writeXStringSet(sequences, "/groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/fasta/combinations_peSTARRSeq.fa") 

# generate bowtie index
cmd <- paste0("module load bowtie/1.2.2-foss-2018b; 
              /software/2020/software/bowtie/1.2.2-foss-2018b/bin/bowtie-build --threads 8 ", #bowtie
              "-f /groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/fasta/combinations_peSTARRSeq.fa ",  # fasta
              "/groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/index_bowtie1/custom_peSTARR1") # output

bsub(cmd, cores= 8)
