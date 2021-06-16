setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
exp_data_dropbox <- "https://www.dropbox.com/sh/rp60dty3vpudmci/AADU9O_QpnJ-smjB46RW7tfxa?dl=0"
dir_exp_data <- "/groups/stark/vloubiere/exp_data/"
dir_fq <- normalizePath("db/fastq/")
dir_sam <- normalizePath("db/sam/")
dir_count <- normalizePath("db/counts/")
dir_dds <- normalizePath("db/dds/")
dir_FC <- normalizePath("db/FC_tables/")
dir_allCounts <- normalizePath("db/umi_counts/")
dir_mergedCounts <- normalizePath("db/merged_counts/")
dir_annotatedCounts <- normalizePath("db/annotated_counts/")
subread_index <- paste0(normalizePath("db/subread_index"), "/vllib001-014")

#--------------------------------------------------------------#
# Update exp data
# Fetch dropbox folder containing my metadata and update local files
#--------------------------------------------------------------#
if(F)
{
  tmp <- tempfile()
  system(paste("wget", exp_data_dropbox, "-O", tmp))
  system(paste("unzip -o", tmp, "-d", dir_exp_data))
}
constructs <- fread(paste0(dir_exp_data, "vl_constructs_sequences.txt"), key = "name")
libs <- fread(paste0(dir_exp_data, "vl_libraries.txt"))
twist8 <- as.data.table(readRDS(paste0("Rdata/vl_library_twist008_112019.rds")))
meta <- fread(paste0(dir_exp_data, "vl_sequencing_metadata.txt"))
meta[, output_prefix:= paste0("/", my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]

#--------------------------------------------------------------#
# Extract from VBC bam file
# For each sequencing run, extract my reads from the bam containing the full lane
#--------------------------------------------------------------#
dir.create(dir_fq, showWarnings = F)
mcmapply(function(b, o, i){
  fq_files <- paste0(o, c("_1.fq.gz", "_2.fq.gz"))
  if(any(!file.exists(fq_files)))
  {
    cmd <- vlfunctions::vl_extract_reads_VBC(bam= b, 
                                             rev_comp_i5 = i)
    system(cmd)
  }
  print(paste0(fq_files, " -->>DONE!"))
},
b= meta[, BAM_path],
o= meta[, paste0(dir_fq, output_prefix)],
i= meta[, as.character(reverseComplement(DNAStringSet(gsub(".*-(.*)", "\\1", i5))))], 
mc.preschedule = F, 
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Create Rsubread index
# Do only once 
#--------------------------------------------------------------#
if(!file.exists("/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_index/vllib001-014.log"))
{
  seqs <- c(twist8$oligo_full_sequence,
            paste0(constructs["Flink_lib", sequence], 
                   constructs[c("SCR1", "SCR2", "HAM1", "SUP1"), sequence], 
                   constructs["R1link_lib", sequence]))
  names <- c(twist8$ID_vl,
             "ts_SCR1_01001", 
             "ts_SCR2_01002",
             "ts_HAM1_01003",
             "ts_SUP1_01004")
  write.fasta(as.list(seqs), 
              names, 
              file.out = "db/fasta/vllib001-014.fasta",
              as.string = T)
  Rsubread::buildindex(basename = subread_index, 
                       reference = "db/fasta/vllib001-014.fasta")
}

#--------------------------------------------------------------#
# Alignment
# Align esach fastq file and produces SAM ouptut
#--------------------------------------------------------------#
dir.create(dir_sam, showWarnings = F)
meta[, {
  sam <- paste0(dir_sam, output_prefix, ".sam")
  if(!file.exists(sam))
  {
   align(index = subread_index,
         readfile1 = paste0(dir_fq, "/", output_prefix, "_1.fq.gz"),
         readfile2 = paste0(dir_fq, "/", output_prefix, "_2.fq.gz"), 
         output_format = "SAM", 
         output_file = sam, 
         maxMismatches = 0, 
         unique = T, 
         nTrim5 = 3, 
         nthreads = getDTthreads()) 
  }
  print(paste0(sam, " -->>DONE"))
}, output_prefix]

#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and collapsed reads using UMIs (0 differences)
#--------------------------------------------------------------#
dir.create(dir_allCounts, showWarnings = F)
meta[, {
  counts <- paste0(dir_allCounts, output_prefix, ".txt")
  if(!file.exists(counts))
  {
    .c <- fread(paste0(dir_sam, output_prefix, ".sam"), 
                skip = 1006,
                header= F, 
                fill= T, 
                select = c(1,3,4,7,8))
    # Keep pairs with two mates are aligned
    .c <- .c[V3!="*" & V7!="*"]
    # Select firt read
    .c <- .c[, .SD[1], V1]
    # Keep pairs with reads on opposite enh ends
    if(type=="peSTARRSeq")
      .c <- .c[(V8-V4)>230]
    if(type=="revpeSTARRSeq")
      .c <- .c[(V4-V8)>230]
    # Extract UMI
    .c <- .c[, .(L= V3, R= V7, UMI= gsub(".*_([A-Z]{10}).*", "\\1", V1))]
    # ompute statistics and collapse
    stat <- data.table(total_reads= nrow(.c))
    .c <- unique(.c)
    stat[, umi_collapsed_reads:= nrow(.c)]
    # SAVE
    fwrite(.c, counts)
    fwrite(stat, gsub(".txt$", "_summary.txt", counts))
  }
  print(paste0(counts, " -->>DONE"))
}, .(output_prefix, type)]

#--------------------------------------------------------------#
# Basic sequencing statistics
# Outputs a barplot showing total and collapsed reads for quality check
#--------------------------------------------------------------#
stats <- data.table(file= list.files("db/umi_counts/", "summary.txt$", full.names = T))
stats <- stats[, fread(file), file]
stats[, output_prefix:= paste0("/", gsub("_summary.txt", "", basename(file)))]
stats[, output_prefix:= paste0(output_prefix, "_", meta[.BY, used, on= "output_prefix"]), output_prefix]
stats[, total_reads:= total_reads-umi_collapsed_reads]
pdf("pdf/alignment_statistics.pdf", 
    height = nrow(stats)/5, 
    width = 10)
par(mar= c(7,30,2,2))
barplot(t(stats[, .(umi_collapsed_reads, total_reads)]), 
        col= c("black", "white"), 
        names.arg = basename(stats$output_prefix),
        cex.names= 0.5, 
        horiz = T)
abline(v= 1e6, 
       lty= 2)
mtext(text = "N reads", 
      side = 1, 
      line = 5)
legend("bottomright", 
       bty= "n",
       fill= c("black", "white"), 
       legend = c("UMI collapsed", "all"))
dev.off()

#--------------------------------------------------------------#
# Merged counts
# Takes counts from separated runs and merge them per condition
# +UMI collapsing with 1nt difference
#--------------------------------------------------------------#
dir.create(dir_mergedCounts, showWarnings = F)
meta[(used), {
  counts <- paste0(dir_mergedCounts, "/", my_ID, "_merged.txt")
  if(!file.exists(counts))
  {
    files <- paste0(dir_allCounts, output_prefix, ".txt")
    .c <- lapply(files, fread)
    .c <- unique(rbindlist(.c))
    # Remove problematic UMIs
    .c <- .c[!grepl("GGGGGGGGG", UMI)]
    # Advanced UMI collapsing
    for(i in 0:9)
    {
      .c <- .c[, .(UMI= UMI[1]), .(L, R, sub(paste0("(.{", i, "})."), "\\1", UMI))]
    }
    .c <- .c[, .(umi_counts= .N), .(L, R)]
    fwrite(.c, counts)
  }
  print(paste0(counts, " -->>DONE"))
}, my_ID]

#-------------------------------#
# Annotated merged counts 
# identify spikein, pairs, artifacts based on library design
#-------------------------------#
dir.create(dir_annotatedCounts, showWarnings = F)
meta[(used), {
  counts <- paste0(dir_annotatedCounts, "/", my_ID, "_annotated.txt")
  if(!file.exists(counts))
  {
    .c <- fread("db/merged_counts/vllib002_DSCP_vllib002_spike_NA_peSTARRSeq_rep1_merged.txt")
    # When red 1 and 2 are frome the same enhancer, read 2 seqnems is "="
    .c[R=="=", R:= L]
    # Check if exists and is spike in
    .ex <- libs[lib_ID %in% c(vllib, Spike_in)]
    .ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]], 
                   R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
    .ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]], 
                    R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
    .ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
                                          paste0("_", R_pattern, "_"))]
    .ex$L <- NULL
    .ex$R <- NULL
    .ex[, Spike:= ifelse(lib_ID==Spike_in, T, F)]
    .ex[, {
      .c[grepl(L_pattern, L) & grepl(R_pattern, R), type:= ifelse(Spike, "spike-in", "pair")]
    }, .(L_pattern, R_pattern, Spike)]
    .c[is.na(type), type:= "switched"]
    fwrite(.c, counts)
  }
  print(paste0(counts, " -->>DONE"))
}, .(my_ID, vllib, Spike_in)]


