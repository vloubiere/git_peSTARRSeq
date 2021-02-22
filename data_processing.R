setwd("/groups/stark/vloubiere/projects/pe_STARRSeq_2/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
require(data.table)
require(Rsubread)
require(Biostrings)
require(parallel)
require(gridExtra)
require(pheatmap)

# Extract fastq from VBC bam ####
fq <- fread("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.txt")
fq <- fq[grepl("002|012|013|014", Sample_ID)]
fq[, fq_prefix:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq_2/db/fastq/", `Flowcell ID`, "_", Sample_ID)]
fq[, i5:= as.character(reverseComplement(DNAStringSet(unlist(tstrsplit(i5, "-", keep=2)))))]
fq <- fq[, .(fq_file= paste0(fq_prefix, c("_1.fq.gz", "_2.fq.gz"))), fq]
fq <- fq[!file.exists(fq_file), .(cmd= extract_reads(BAM_path, fq_prefix, i5)), .(BAM_path, fq_prefix, i5)]
mclapply(fq$cmd, system, mc.preschedule = F, mc.cores = getDTthreads()-1)
####

# Create Rsubread index ####
check <- list.files("/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/index/")
if(!any(grepl("^peSTARR", check))){
  # Import lib
  lib <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/library/vl_library_112019.rds"))
  lib <- lib[, .(ID_vl, oligo_full_sequence)]
  # add template switch sequences
  add <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt", key= "name")
  lib <- rbind(lib, data.table(ID_vl= c("ts_SCR2_01002", "ts_HAM1_01003", "ts_SUP1_01004"), 
                               oligo_full_sequence= c(paste0(add[c("Flink_lib", "SCR2", "R1link_lib"), sequence], collapse= ""),
                                                      paste0(add[c("Flink_lib", "HAM1", "R1link_lib"), sequence], collapse= ""),
                                                      paste0(add[c("Flink_lib", "SUP1", "R1link_lib"), sequence], collapse= ""))))
  # Save fasta
  sequences <- DNAStringSet(lib$oligo_full_sequence)
  names(sequences) <- lib$ID_vl
  writeXStringSet(sequences, "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/peSTARRSeq_sequences.fa")
  
  buildindex(basename= "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/index/peSTARR_idx", 
             reference = "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/peSTARRSeq_sequences.fa")
}
####

# alignment ####
aln <- data.table(file= list.files("db/fastq/", "fq.gz", full.names = T))
aln[, bam:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq_2/db/bam/", gsub(".fq.gz", ".bam", basename(file)))]
aln[, {
  if(!file.exists(bam)){
    print(paste("START", file))
    stats <- capture.output(align(index = "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/index/peSTARR_idx",
                                  readfile1 = file, type= "dna", output_file = bam, maxMismatches= 0, unique= T, nTrim5 = 15, 
                                  nthreads= getDTthreads()-1))
    writeLines(stats, gsub(".bam$", "_stats.txt", bam))
    print(paste(bam, "DONE!"))
  }
}, .(file, bam)]
####

# Convert to text files ####
convert <- data.table(bam= list.files("db/bam", ".bam$", full.names = T))
convert[, txt:= gsub(".bam$", ".txt", bam)]
convert[, cmd:= paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; samtools view ", bam)]
convert <- convert[!file.exists(txt), .(cmd= paste0(cmd, " > ", txt))]
mclapply(convert$cmd, system, mc.preschedule = F, mc.cores = getDTthreads()-1)

# Alignment statistics ####
if(!file.exists("pdf/alignment_statistics.pdf")){
  stats <- data.table(file= list.files("db/bam", "stats.txt", full.names = T))
  stats <- stats[, .(Mapped_reads= readLines(file)), file]
  stats <- stats[grepl("Uniquely_mapped_reads|Unmapped_reads", Mapped_reads)]
  stats[, counts:= {current <- tstrsplit(Mapped_reads, " "); current[length(current)]}, Mapped_reads]
  stats[, variable:= ifelse(grepl("Uniquely_mapped_reads", Mapped_reads), "Uniquely_mapped_reads", "Unmapped_reads"), Mapped_reads]
  stats[, counts:= round(as.numeric(counts)/1e6, 2)]
  stats <- dcast(stats, file~variable, value.var = "counts")
  stats[, perc:= round(Uniquely_mapped_reads/(Uniquely_mapped_reads+Unmapped_reads)*100, 1)] 
  
  pdf("pdf/alignment_statistics.pdf", width = 15, height = 30)
  grid.table(stats)
  dev.off()
}
####

# Map fine alignment of sub-libs ####
if(!file.exists("pdf/sub_libraries_fine_alignment.pdf"))
{
  refine <- data.table(file= list.files("db/bam/", ".txt$", full.names = T))
  refine[, lib:= tstrsplit(basename(file), "_", keep= 2)]
  refine <- refine[!grepl("stats", file)]
  refine <- refine[, fread(file, nrows = 10000, fill= T), refine]
  refine <- refine[V3!="*" & V5>30]
  refine[grep("_A_", V3), class:= "R1"]
  refine[grep("_B_", V3), class:= "R2"]
  refine[grep("_C_", V3), class:= "R3"]
  refine[grep("_HAM1_", V3), class:= "HAM1"]
  refine[grep("_SCR2_", V3), class:= "SCR2"]
  refine[grep("_SUP1_", V3), class:= "SUP1"]
  heat <- dcast(refine, lib+class~V4)
  gaps <- cumsum(heat[, .N, lib][, N])
  heat[, lib:= paste0(lib, "__", class)]
  heat <- heat[, !"class"]
  mat <- as.matrix(heat, 1)
  pheatmap(log2(mat+1), filename = "pdf/sub_libraries_fine_alignment.pdf",
           cluster_rows = F, cluster_cols = F, width = 20, height = 3, gaps_row = gaps)
}
####
# read 1= 16 to 19
# read 2= 256 to 303

# retrieve pairs and compute counts ####
pa <- data.table(file= list.files("db/bam", ".txt$", full.names = T))
pa <- pa[!grepl("stats", file)]
pa[, cdition:= gsub("_1.txt|_2.txt", "", basename(file))]
pa[, read:= ifelse(grepl("_1.txt$", file), "read_1", "read_2")]
pa <- dcast(pa, cdition~read, value.var = "file")
pa[, cdition:= substr(cdition, 11, nchar(cdition)), cdition]
pa[, umi_counts:= paste0("db/count/", cdition, "_umi_count.txt"), cdition]
pa[!file.exists(umi_counts), {
  .L <- rbindlist(lapply(read_1, function(x) fread(x, fill = T, select = c(1,3,4,5))))
  .L <- .L[V3!="*" & between(V4, 16, 19, incbounds = T) & V5>30]
  .R <- rbindlist(lapply(read_2, function(x) fread(x, fill = T, select = c(1,3,4,5))))
  .R <- .R[V3!="*" & between(V4, 256, 303, incbounds = T) & V5>30]
  .pair <- merge(.L, .R, by.x= "V1", by.y= "V1", suffixes= c("_L", "_R"))
  .pair[, UMI:= tstrsplit(V1, "_", keep= 2)]
  .pair <- unique(.pair[, .(L= V3_L, R= V3_R, UMI)])
  .pair <- .pair[, .(count= .N), .(L, R)]
  fwrite(.pair, umi_counts, col.names = T, row.names = F, sep= "\t", quote= F)
  print(paste0(umi_counts, " DONE!"))
}, umi_counts]
####

# Check saturation ####