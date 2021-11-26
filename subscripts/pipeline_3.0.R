setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
# require(DESeq2)
require(readxl)

# PATHS
dir_fq <- normalizePath("db/fastq/", mustWork = F)
dir.create(dir_fq, showWarnings = F)
dir_sam <- normalizePath("db/sam/", mustWork = F)
dir.create(dir_sam, showWarnings = F)
dir_umiCounts <- normalizePath("db/umi_counts/", mustWork = F)
dir.create(dir_umiCounts, showWarnings = F)
dir_mergedCounts <- normalizePath("db/merged_counts/", mustWork = F)
dir.create(dir_mergedCounts, showWarnings = F)
dir_dds <- normalizePath("db/dds/", mustWork = F)
dir.create(dir_dds, showWarnings = F)
dir_FC <- normalizePath("db/FC_tables/", mustWork = F)
dir.create(dir_FC, showWarnings = F)
dir_final <- normalizePath("db/final_tables_exp_model/", mustWork = F)
dir.create(dir_final, showWarnings = F)
dir.create("pdf/alignment", showWarnings = F)
dir.create("db/final_tables_exp_model/counts_norm/", showWarnings = F)
dir.create("db/final_tables_exp_model/replicates_counts_norm/", showWarnings = F)

#--------------------------------------------------------------#
# Update exp data
# Fetch dropbox folder containing my metadata and update local files
#--------------------------------------------------------------#
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", key = "name")
libs <- as.data.table(read_excel("/groups/stark/vloubiere/exp_data/vl_libraries.xlsx"))
twist8 <- as.data.table(readRDS(paste0("Rdata/vl_library_twist008_112019.rds")))
twist12 <- as.data.table(readRDS(paste0("Rdata/vl_library_twist12_210610.rds")))
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- colnames(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]
meta[, output_prefix:= paste0("/", my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]
# IMPORTANT!!
if(any(meta[, .N, output_prefix]$N>1))
  stop(paste0("Some output prefixes are not unique and cannot be used. Check metadata table (replicates?)!\n"))

#--------------------------------------------------------------#
# Extract from VBC bam file
# For each sequencing run, extract my reads from the bam containing the full lane
#--------------------------------------------------------------#
mcmapply(function(b, o, i){
  fq_files <- paste0(o, c("_1.fq.gz", "_2.fq.gz"))
  if(any(!file.exists(fq_files)))
  {
    cmd <- vlfunctions::vl_extract_reads_VBC(bam= b,
                                             output_prefix = o,
                                             rev_comp_i5 = i)
    system(cmd)
    print(paste0(fq_files, " -->> DONE!"))
  }else
    print(paste0(fq_files, " -->> ALREADY EXISTS!"))
},
b= meta[, BAM_path],
o= meta[, paste0(dir_fq, output_prefix)],
i= meta[, as.character(reverseComplement(DNAStringSet(gsub(".*-(.*)", "\\1", i5))))], 
mc.preschedule = F, 
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
# Align each fastq file and produces SAM ouptut
#--------------------------------------------------------------#
meta[, {
  sam <- paste0(dir_sam, output_prefix, ".sam")
  if(file.exists(sam))
    print(paste0(sam, " -->> ALREADY EXISTS"))
  else
  {
    if(library=="T8")
      subread_index <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8"
    else if(library=="T12")
      subread_index <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12"
    align(index = subread_index,
          readfile1 = paste0(dir_fq, "/", output_prefix, "_1.fq.gz"),
          readfile2 = paste0(dir_fq, "/", output_prefix, "_2.fq.gz"),
          output_format = "SAM",
          output_file = sam,
          maxMismatches = 3,
          unique = T,
          nTrim5 = 3,
          nthreads = getDTthreads())
    print(paste0(sam, " -->> DONE"))
  }
}, .(output_prefix, library)]

#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and collapse reads using UMIs (0 mistmatches)
#--------------------------------------------------------------#
meta[, {
  counts <- paste0(dir_umiCounts, output_prefix, ".txt")
  if(file.exists(counts))
    print(paste0(counts, " -->> ALREADY EXISTS"))
  else
  {
    .c <- fread(paste0(dir_sam, output_prefix, ".sam"), 
                skip = switch(library, "T8"= 1006, "T12"= 395),
                header= F, 
                fill= T, 
                select = c(1,3,4,7,8))
    # Keep pairs with two mates are aligned
    .c <- .c[V3!="*" & V7!="*"]
    # Select firt read
    .c <- .c[, .SD[1], V1]
    # Keep pairs with reads on opposite enh ends
    if(type=="pe-STARR-Seq")
      .c <- .c[(V8-V4)>230]
    if(type=="rev-pe-STARR-Seq")
      .c <- .c[(V4-V8)>230]
    # Extract UMI
    .c <- .c[, .(L= V3, R= V7, UMI= gsub(".*_([A-Z]{10}).*", "\\1", V1))]
    # Compute total reads
    stat <- data.table(total_reads= nrow(.c))
    # UMI collapsing
    .c <- unique(.c)
    # When red 1 and 2 are frome the same enhancer, read 2 seqnems is "="
    .c[R=="=", R:= L]
    # Compute collapsed reads
    stat[, umi_collapsed_reads:= nrow(.c)]
    # SAVE
    fwrite(.c, counts)
    fwrite(stat, gsub(".txt$", "_summary.txt", counts))
    print(paste0(counts, " -->>DONE"))
  }
}, .(output_prefix, type, library)]

#--------------------------------------------------------------#
# Merged counts
# Takes counts from separated runs and merge them per condition
# +UMI collapsing with 1nt difference
# +annotation based on library
#--------------------------------------------------------------#
meta[!is.na(DESeq2_group), {
  if(.N>0)
  {
    counts <- paste0(dir_mergedCounts, "/", DESeq2_group, "_", cdition, "_rep", DESeq2_pseudo_rep, "_merged.txt")
    if(file.exists(counts))
      print(paste0(counts, " -->> ALREADY EXISTS")) else
      {
        print(paste0("Start ", counts))
        files <- paste0(dir_umiCounts, output_prefix, ".txt")
        .c <- lapply(files, fread)
        .c <- unique(rbindlist(.c))
        # Remove problematic UMIs
        .c <- .c[!grepl("GGGGGGGGG", UMI)]
        # Advanced UMI collapsing
        for(i in 0:9) 
          .c <- .c[, .(UMI= UMI[1]), .(L, R, sub(paste0("(.{", i, "})."), "\\1", UMI))]
        # UMI collapsing
        .c <- .c[, .(umi_counts= .N), .(L, R)]
        # Compute patterns to identify library pairs
        .ex <- libs[lib_ID %in% c(vllib, Spike_in)]
        .ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]], 
                       R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
        .ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]], 
                        R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
        .ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
                                              paste0("_", R_pattern, "_"))]
        .ex$L <- NULL
        .ex$R <- NULL
        .ex[, Spike:= ifelse(lib_ID==Spike_in & !is.na(Spike_in), T, F)]
        # Check if pair exists and is spike in
        .ex[, {
          .c[grepl(L_pattern, L) & grepl(R_pattern, R), type:= ifelse(Spike, "spike-in", "pair")]
        }, .(L_pattern, R_pattern, Spike)]
        .c[is.na(type), type:= "switched"]
        # Generate read summary
        sum_files <- data.table(pattern= gsub(".txt", "", basename(files)))
        sum_files[, file:= list.files(dir_sam, paste0(pattern, ".*.summary$"), full.names = T), pattern]
        sum_files[, mapped:= fread(file)[V1=="Mapped_fragments", V2], file]
        summary <- sum_files[, .(mapped= sum(mapped), collapsed= sum(.c$umi_counts))]
        # SAVE
        fwrite(.c, counts)
        fwrite(summary, gsub(".txt$", ".summary.txt", counts))
        print(paste0(counts, " -->> DONE"))
      }
  }
}, .(DESeq2_group, cdition, CP, DESeq2_pseudo_rep, vllib, Spike_in)]

#--------------------------------------------------------------#
# Method #1 using counts norm (Tolerates one rep)
#--------------------------------------------------------------#
meta[!is.na(DESeq2_group), {
  if(.N>0)
  {
    FC <- paste0(dir_final, "/counts_norm/", DESeq2_group, "_counts_norm_final_oe.txt")
    if(file.exists(FC))
      print(paste0(FC, "  -->> ALREADY EXISTS")) else
      {
        # Import counts
        dat <- data.table(file= list.files(dir_mergedCounts, paste0(DESeq2_group, ".*_merged.txt$"), full.names = T))
        dat <- dat[, fread(file), file][type=="pair"]
        dat[, rep:= paste0("log2FC_rep", gsub(".*_rep([1-9]{1})_merged.txt$", "\\1", file)), file]
        dat[grepl("input", basename(file)), cdition:= "input"]
        dat[grepl("screen", basename(file)), cdition:= "screen"]
        # Filter
        dat[, check_counts:= sum(umi_counts)>=10, .(L, R)]
        dat <- dat[(check_counts), !c("file", "check_counts", "type")]
        # dat <- dat[L!=R] # COULD BE USEFUL (SAFER?)!!!!!!!!!!!!
        # Add merged counts
        dat <- rbind(dat, 
                     dat[, .(umi_counts= sum(umi_counts), 
                             rep= "log2FC_merge"), .(L, R, cdition)])
        # Normalize
        dat[, norm:= (umi_counts+1)/sum(umi_counts)*1e6, .(rep, cdition)]
        # Compute FC 
        dat <- dat[, .(FC= log2(norm[cdition=="screen"]/norm[cdition=="input"])), .(L, R, rep)]
        dat <- rbind(dat,
                     dat[!grepl("merge$", rep), data.table(rep= "log2FoldChange",
                                                           FC= mean(FC)), .(L, R)])
        # Cast
        res <- dcast(dat, 
                     L+R~rep, 
                     value.var= "FC")
        # Subtract basal activity (center controls on 0)
        cols <- grep("^log2", names(res), value = T)
        res[, (cols):= lapply(.SD, function(x) x-median(x[grepl("^control", L) & grepl("^control", R)], na.rm = T)), .SDcols= cols]
        # Compute expected
        med_L <- gsub("log2FoldChange", "median", paste0(gsub("log2FC", "median", cols), "_L"))
        res[, (med_L):= lapply(.SD, function(x) {
          x <- na.omit(x[grepl("control", R)])
          if(length(x)>=5)
            median(x) else
              as.numeric(NA)
        }), L, .SDcols= cols]
        med_R <- gsub("log2FoldChange", "median", paste0(gsub("log2FC", "median", cols), "_R"))
        res[, (med_R):= lapply(.SD, function(x) {
          x <- na.omit(x[grepl("control", L)])
          if(length(x)>=5)
            median(x) else
              as.numeric(NA)
        }), R, .SDcols= cols]
        # SAVE
        fwrite(res, FC)
        print(paste0(FC, "  -->> DONE"))
      }
  }
}, DESeq2_group]

# #--------------------------------------------------------------#
# # Method # 2 using DESeq2 ==> DEPRECATED
# #--------------------------------------------------------------#
# # dds --------------#
# meta[!is.na(DESeq2_group), {
#   for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
#   {
#     dds_file <- paste0(dir_dds, "/", DESeq2_group, "_", stringency, "_dds.rds")
#     if(file.exists(dds_file))
#       print(paste0(dds_file, "  -->> ALREADY EXISTS"))
#     else
#     {
#       # Import counts
#       file <- list.files(dir_mergedCounts, DESeq2_group, full.names = T)
#       .c <- lapply(file, fread)
#       names(.c) <- gsub(".*_(.*)_(.*)_merged.txt", "\\1_\\2", basename(file))
#       .c <- rbindlist(.c, idcol = T)
#       # Cast and make DF
#       mat <- dcast(.c[type!="switched"], L+R~.id, value.var = "umi_counts", fill= 0)
#       DF <- data.frame(mat[, -c(1,2)])
#       rownames(DF) <- paste0(mat$L, "_vs_", mat$R)
#       if(stringency=="low_cutoff")
#         DF <- DF[rowSums(DF)>5,]
#       if(stringency=="medium_cutoff")
#         DF <- DF[rowSums(DF)>25,]
#       if(stringency=="high_cutoff")
#         DF <- DF[rowSums(DF)>5 & apply(DF, 1, function(x) all(x)>0),]
#       if(stringency=="max_cutoff")
#         DF <- DF[rowSums(DF)>25 & apply(DF, 1, function(x) all(x)>0),]
#       # SampleTable
#       sampleTable <- data.frame(cdition= unlist(tstrsplit(colnames(DF), "_rep", keep= 1)),
#                                 rep= unlist(tstrsplit(colnames(DF), "_rep", keep= 2)))
#       rownames(sampleTable) <- colnames(DF)
#       # DESeq2
#       dds <- try(DESeq2::DESeqDataSetFromMatrix(countData= DF, 
#                                                 colData= sampleTable, 
#                                                 design= ~rep+cdition), 
#                  silent = T)
#       if(class(dds)=="try-error")
#         print(paste0(DESeq2_group, " could not pass DESeqDataSetFromMatrix command and will be skipped -->> Check N replicates?"))
#       else
#       {
#         # SizeFactors
#         ctls_idx <- grep("control.*vs.*control", rownames(DF))
#         sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[ctls_idx,]))
#         # dds result
#         res <- try(DESeq2::DESeq(dds), silent = T)
#         if(class(res)=="try-error")
#           print(paste0("No replicates --> ", DESeq2_group, "SKIPED"))
#         else
#         {
#           saveRDS(res, dds_file) 
#           print(paste0(dds_file, "  -->> DONE"))
#         }
#       }
#     }
#   }
# }, DESeq2_group]
# 
# # Differential expression----------#
# meta[!is.na(DESeq2_group), {
#   for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
#   {
#     FC <- paste0(dir_FC, "/", DESeq2_group, "_", stringency, "_log2FC.txt")
#     if(file.exists(FC))
#       paste0(FC, "  -->> ALREADY EXISTS")
#     else
#     {
#       # Check if models exists
#       dds_file <- paste0(dir_dds, "/", DESeq2_group, "_", stringency, "_dds.rds")
#       if(!file.exists(dds_file))
#         print(paste0(FC, "  -->> No model --> ", DESeq2_group, " -->> SKIPPED"))
#       else
#       {
#         # Import dds
#         dds <- readRDS(dds_file)
#         # Compoute FC
#         .c <- DESeq2::results(dds, contrast= c("cdition", "DSCP", "input"))
#         # Reformat
#         .c <- as.data.table(as.data.frame(.c), keep.rownames= T)
#         .c[, c("L", "R"):= tstrsplit(rn, "_vs_")]
#         .c <- .c[, !"rn"]
#         setcolorder(.c, c("L", "R"))
#         # save
#         fwrite(.c, FC, col.names = T, row.names = F, sep= "\t", quote= F)
#         print(paste0(FC, "  -->> DONE!"))
#       }
#     }
#   }
# }, DESeq2_group]
# 
# # Compute Left and right activities -----------------------#
# meta[!is.na(DESeq2_group), {
#   for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
#   {
#     final <- paste0(dir_final, "/DESeq2/", DESeq2_group, "_", stringency, "_final_oe.txt")
#     if(file.exists(final))
#       paste0(final, "  -->> ALREADY EXISTS")
#     else
#     {
#       FC_file <- paste0(dir_FC, "/", DESeq2_group, "_", stringency, "_log2FC.txt")
#       if(!file.exists(FC_file))
#         print(paste0(FC_file, " does not exist. -->> SKIPPED"))
#       else
#       {
#         dat <- fread(FC_file)
#         if(type=="peSTARRSeq")
#         {
#           #---------- Select robust controls -----------#
#           cor <- data.table(ID= unique(grep("^control", c(dat$L, dat$R), value = T)))
#           # Compute L/R PCC for each control fragment
#           cor[, PCC:= {
#             .c <- merge(dat[.BY, .(enh= L, L= log2FoldChange), on= "R==ID"],
#                         dat[.BY, .(enh= R, R= log2FoldChange), on= "L==ID"])
#             .c <- na.omit(.c)
#             if(nrow(.c)>=5)
#               cor.test(.c$L, .c$R)$estimate
#             else
#               as.numeric(NA)
#           }, ID]
#           # Compute outliers frequency
#           out <- copy(dat)
#           out[grepl("^control", L), out_R:= log2FoldChange %in% boxplot(log2FoldChange, plot= F)$out, R]
#           out[grepl("^control", L), out_fL:= length(which(out_R))/.N, L]
#           out[grepl("^control", R), out_L:= log2FoldChange %in% boxplot(log2FoldChange, plot= F)$out, L]
#           out[grepl("^control", R), out_fR:= length(which(out_L))/.N, R]
#           # Robust controls list
#           dat[, ctl_L:= L %in% cor[PCC>0.75, ID] & L %in% out[out_fL<0.025, L]]
#           dat[, ctl_R:= R %in% cor[PCC>0.75, ID] & R %in% out[out_fR<0.025, R]]
#           #---------- Compute expected values -----------#
#           dat[, length(which(ctl_R)), L]
#         }else if(type=="revpeSTARRSeq")
#         {
#           dat[, ctl_L:= grepl("^control", L)]
#           dat[, ctl_R:= grepl("^control", R)]
#         }
#         dat[, median_L:= ifelse(length(which(ctl_R))>=5, median(log2FoldChange[ctl_R], na.rm = T), as.numeric(NA)), L]
#         dat[, median_R:= ifelse(length(which(ctl_L))>=5, median(log2FoldChange[ctl_L], na.rm = T), as.numeric(NA)), R]
#         fwrite(dat, final)
#         print(final, "  -->> DONE")
#       }
#     }
#   }
# }, .(DESeq2_group, type)]
