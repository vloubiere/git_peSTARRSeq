setwd("/groups/stark/vloubiere/projects/pe_STARRSeq_2/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
# Packs ####
require(data.table)
require(Rsubread)
require(Biostrings)
require(parallel)
require(gridExtra)
require(pheatmap)
require(DESeq2)
####

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
pa[, umi_counts:= paste0("db/count/", cdition, "_all_umi_counts.txt"), cdition]
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

# Clean pairs and keep only expected pairs ####
clean <- data.table(file= list.files("db/count", "_all_umi_counts.txt$", full.names = T))
clean[, lib:= tstrsplit(basename(file), "_", keep= 1) ]
clean[grepl("DSCP", file), cdition:= "DSCP"]
clean[grepl("input", file), cdition:= "input"]
clean[, rep:= substr(file, nchar(file)-22, nchar(file)-19)]
clean[, clean:= paste0("db/count/", lib, "_", cdition, "_", rep, "_clean_counts.txt")]
clean[!file.exists(clean), {
  .c <- fread(file)
  if(lib=="vllib012"){
    .c <- .c[(grepl("_C_", L) & grepl("_SCR2_", R)) 
             | (grepl("_SCR2_", L) & grepl("_C_", R))]
  }
  if(lib=="vllib013"){
    .c <- .c[(grepl("_A_", L) & grepl("_A_", R)) 
             | (grepl("_C_", L) & grepl("_SCR2_", R)) 
             | (grepl("_SCR2_", L) & grepl("_C_", R)) ]
  }
  if(lib=="vllib014"){
    .c <- .c[(grepl("_B_", L) & grepl("_B_", R)) 
             | (grepl("_C_", L) & grepl("_SCR2_", R)) 
             | (grepl("_SCR2_", L) & grepl("_C_", R)) ]
  }
  fwrite(.c, paste0("db/count/", lib, "_", cdition, "_", rep, "_clean_umi_counts.txt"))
}, .(file, clean)]
####

# Check saturation ####
sat <- data.table(file= list.files("db/count", "clean_umi_counts.txt$", full.names = T))
sat[, sample:= tstrsplit(basename(file), "_clean", keep= 1)]
sat[, lib:= tstrsplit(sample, "_", keep= 1)]
sat <- sat[, fread(file), sat]
pdf("pdf/sequencing_saturation.pdf")
sat[, {
  .c <- .SD[, .(dens= list(density(log2(count)+1))), sample]
  cc1 <- colorRampPalette(c("tomato", "darkred"))(.c[grepl("input", sample), .N])
  cc2 <- colorRampPalette(c("cornflowerblue", "royalblue2"))(.c[grepl("DSCP", sample), .N])
  .c[grepl("input", sample), Cc:= cc1[.GRP], sample]
  .c[grepl("DSCP", sample), Cc:= cc2[.GRP], sample]
  .x <- .c[, range(dens[[1]]$x), sample][, range(V1)]
  .y <- .c[, range(dens[[1]]$y), sample][, range(V1)]
  plot(NA, xlim= .x, ylim= .y, las= 1, xlab= "log2(counts+1)", ylab= "density", main= lib)
  legend("topright", fill= .c$Cc, legend = .c$sample, bty= "n")
  .c[, lines(dens[[1]], col= Cc), .(sample, Cc)]
  print("")
}, lib]
dev.off()
####




# Differential analysis ####
# 1- Make master counts table
if(!file.exists("Rdata/master_counts_table.txt"))
{
  raw <- data.table(file= list.files("db/count/", "clean", full.names = T))
  raw[, c("lib", "cdition", "rep"):= tstrsplit(basename(file), "_", keep= 1:3)]
  raw <- raw[, fread(file), raw]
  raw <- dcast(raw, L+R~lib+cdition+rep, value.var = "count")
  
  # Add pseudo replicates when missing
  vl13_long <- raw[!is.na(vllib013_input_rep1)][rep(seq(.N), vllib013_input_rep1), .(L, R)]
  raw[vl13_long[sample(.N, 6e6)][, .N, .(L, R)], vllib013_pinput_rep1:= i.N, on= c("L", "R")]
  raw[vl13_long[sample(.N, 6e6)][, .N, .(L, R)], vllib013_pinput_rep2:= i.N, on= c("L", "R")]
  vl14_long <- raw[!is.na(vllib014_input_rep1)][rep(seq(.N), vllib014_input_rep1), .(L, R)]
  raw[vl14_long[sample(.N, 7.5e6)][, .N, .(L, R)], vllib014_pinput_rep1:= i.N, on= c("L", "R")]
  raw[vl14_long[sample(.N, 7.5e6)][, .N, .(L, R)], vllib014_pinput_rep2:= i.N, on= c("L", "R")]
  
  # Clean counts table
  counts <- data.table(raw[, .(L, R)],
                       libvl002_DSCP_rep1= raw$libvl002_DSCP_rep1+raw$libvl002_DSCP_rep2,
                       libvl002_DSCP_rep2= raw$libvl002_DSCP_rep3+raw$libvl002_DSCP_rep4,
                       libvl002_input_rep1= raw$libvl002_input_rep1+raw$libvl002_input_rep3+raw$libvl002_input_rep4,
                       libvl002_input_rep2= raw$libvl002_input_rep2+raw$libvl002_input_rep5+raw$libvl002_input_rep6,
                       libvl013_input_rep1= raw$vllib013_pinput_rep1,
                       libvl013_input_rep2= raw$vllib013_pinput_rep2,
                       libvl013_DSCP_rep1= raw$vllib013_DSCP_rep1,
                       libvl013_DSCP_rep2= raw$vllib013_DSCP_rep2,
                       libvl014_input_rep1= raw$vllib014_pinput_rep1,
                       libvl014_input_rep2= raw$vllib014_pinput_rep2,
                       libvl014_DSCP_rep1= raw$vllib014_DSCP_rep1,
                       libvl014_DSCP_rep2= raw$vllib014_DSCP_rep2)
  fwrite(counts, "Rdata/master_counts_table.txt", col.names = T, row.names = F, sep= "\t", quote= F)
}
# 2- run DESeq
DE <- data.table(lib= colnames(fread("Rdata/master_counts_table.txt", nrows = 0))[-c(1,2)])
DE[, c("lib", "cdition", "rep"):= tstrsplit(lib, "_", keep= 1:3)]
DE[, dds_file:= paste0("db/DE_analysis/", lib, "_dds.rds")]
if(any(!file.exists(DE$dds_file)))
{
  if(!exists("master_counts")){
    master_counts <- fread("Rdata/master_counts_table.txt")
    master_counts[, ID:= paste0(L, "_vs_", R)]
    master_counts$L <- NULL
    master_counts$R <- NULL
    setcolorder(master_counts, "ID")
  }
  DE[!file.exists(dds_file), {
    # counts
    sel <- c("ID", grep(lib, colnames(master_counts), value = T))
    .c <- na.omit(master_counts[, which(colnames(master_counts) %in% sel), with= F])
    .c <- data.frame(.c[rowSums(.c[, -1])>10], row.names = 1)
    # sampleTable
    .s <- data.frame(.SD[, .(cdition, rep)], row.names = paste0(lib, "_", .SD$cdition, "_", .SD$rep))
    # DESeq2
    dds <- DESeqDataSetFromMatrix(countData= .c, colData= .s, design= ~rep+cdition)
    # SizeFactors
    sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(.c[grep("control.*vs.*control", rownames(.c)),]))
    # dds result
    res <- DESeq(dds)
    saveRDS(res, dds_file)
  }, .(lib, dds_file)]
}
# make diff table
DE[, diff:= gsub("_dds.rds$", "_DE.txt", dds_file)]
DE[!file.exists(diff), {
  .c <- readRDS(dds_file)
  # Differential expression
  .c <- as.data.table(as.data.frame(results(.c, contrast= c("cdition", "DSCP", "input"))), keep.rownames= T)
  .c[, c("L", "R"):= tstrsplit(rn, "_vs_")]
  .c <- .c[, !"rn"]
  setcolorder(.c, c("L", "R"))
  fwrite(.c, diff, col.names = T, row.names = F, sep= "\t", quote= F)
  print(paste0(diff, " DONE!"))
}, .(dds_file, diff)]
####

# Compute expected score ####
FC <- data.table(file= list.files("db/DE_analysis/", "DE.txt$", full.names = T))
FC[,  lib:= gsub("_DE.txt$", "", basename(file))]
FC[,  FC:= paste0("db/DE_analysis/", lib, "_add_scores.txt")]
FC[!file.exists(FC), {
  .c <- fread(file)
  .c[, median_L:= .SD[grep("control", R), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], L]
  .c[grepl("_C_", L) & R=="ts_SCR2_01002", median_L:= log2FoldChange, L]
  .c[, median_R:= .SD[grep("control", L), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], R]
  .c[L=="ts_SCR2_01002" & grepl("_C_", R), median_R:= log2FoldChange, R]
  .c <- na.omit(.c[!grepl("control", L) & !grepl("control", R)])
  .c[, add:= log2(sum(2^median_L+2^median_R)), .(L, R)]
  fwrite(.c, FC, col.names = T, row.names = F, sep= "\t", quote= F)
}, .(file, FC)]
####

# Comparison Bernardo ####
if(!file.exists("pdf/comparison_BA.pdf"))
{
  lib <- readRDS("../pe_STARRSeq/Rdata/library/uniq_library_final.rds")
  res <- data.table(file= list.files("db/DE_analysis/", "add_scores.txt", full.names = T))
  res[, lib:= tstrsplit(basename(file), "_", keep= 1)]
  res <- res[, fread(file), res]
  res[lib, single_L:= i.dev_log2FoldChange, on= "L==ID"]
  res[lib, single_R:= i.dev_log2FoldChange, on= "R==ID"]
  pdf("pdf/comparison_BA.pdf", width = 6, height = 10)
  par(mfrow= c(3,2), pch= 16)
  res[, {
    xl <- "twist-STARR-Seq"
    yl <- "pe-STARR-Seq"
    .c <- unique(.SD[, .(L, median_L, single_L)])
    .p <- seq(-2, 12, 0.1)
    plot(.c$median_L~.c$single_L, .c, main= paste0(lib, " L"), xlab= xl, ylab= yl)
    abline(0,1)
    lines(.p, predict(loess(.c$median_L~.c$single_L), .p), col= "red")
    .c <- unique(.SD[, .(R, median_R, single_R)])
    plot(median_R~single_R, .c, main= paste0(lib, " R"), xlab= xl, ylab= yl)
    lines(.p, predict(loess(.c$median_R~.c$single_R), .p), col= "red")
    abline(0,1)
  }, lib]
  dev.off()
}
####
# Comparison 