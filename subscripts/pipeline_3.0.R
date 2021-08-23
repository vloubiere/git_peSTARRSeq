setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
require(DESeq2)
require(readxl)
dir_exp_data <- "/groups/stark/vloubiere/exp_data/"
dir_fq <- normalizePath("db/fastq/", mustWork = F)
dir_sam <- normalizePath("db/sam/", mustWork = F)
dir_allCounts <- normalizePath("db/umi_counts/", mustWork = F)
dir_mergedCounts <- normalizePath("db/merged_counts/", mustWork = F)
dir_dds <- normalizePath("db/dds/", mustWork = F)
dir_FC <- normalizePath("db/FC_tables/", mustWork = F)
dir_final <- normalizePath("db/final_tables_exp_model/", mustWork = F)

#--------------------------------------------------------------#
# Update exp data
# Fetch dropbox folder containing my metadata and update local files
#--------------------------------------------------------------#
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
constructs <- fread(paste0(dir_exp_data, "vl_constructs_sequences.txt"), key = "name")
libs <- fread(paste0(dir_exp_data, "vl_libraries.txt"))
twist8 <- as.data.table(readRDS(paste0("Rdata/vl_library_twist008_112019.rds")))
twist12 <- as.data.table(readRDS(paste0("Rdata/vl_library_twist12_210610.rds")))
meta <- read_xlsx(paste0(dir_exp_data, "vl_sequencing_metadata.xlsx"))
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
dir.create(dir_fq, showWarnings = F)
mcmapply(function(b, o, i){
  fq_files <- paste0(o, c("_1.fq.gz", "_2.fq.gz"))
  if(any(!file.exists(fq_files)))
  {
    cmd <- vlfunctions::vl_extract_reads_VBC(bam= b,
                                             output_prefix = o,
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
if(!file.exists("/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8.log"))
  source("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist8_subread_index.R")
if(!file.exists("/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12.log"))
  source("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/create_twist12_subread_index.R")

#--------------------------------------------------------------#
# Alignment
# Align each fastq file and produces SAM ouptut
#--------------------------------------------------------------#
dir.create(dir_sam, showWarnings = F)
meta[, {
  sam <- paste0(dir_sam, output_prefix, ".sam")
  if(file.exists(sam))
    print(paste0(sam, " -->> ALREADY EXISTS")) 
  else
  {
    if(twist_lib==8)
      subread_index <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8"
    else if(twist_lib==12)
      subread_index <- "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12"
    align(index = subread_index,
          readfile1 = paste0(dir_fq, "/", output_prefix, "_1.fq.gz"),
          readfile2 = paste0(dir_fq, "/", output_prefix, "_2.fq.gz"), 
          output_format = "SAM", 
          output_file = sam, 
          maxMismatches = 0, 
          unique = T, 
          nTrim5 = 3, 
          nthreads = getDTthreads()) 
    print(paste0(sam, " -->> DONE")) 
  }
}, .(output_prefix, twist_lib)]

#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and collapsed reads using UMIs (0 mistmatches)
#--------------------------------------------------------------#
dir.create(dir_allCounts, showWarnings = F)
meta[, {
  counts <- paste0(dir_allCounts, output_prefix, ".txt")
  if(file.exists(counts))
    print(paste0(counts, " -->> ALREADY EXISTS"))
  else
  {
    .c <- fread(paste0(dir_sam, output_prefix, ".sam"), 
                skip = ifelse(twist_lib==8, 1006, 395),
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
}, .(output_prefix, type, twist_lib)]

#--------------------------------------------------------------#
# Basic sequencing statistics
# Outputs a barplot showing total and collapsed reads for quality check
#--------------------------------------------------------------#
stats <- data.table(file= list.files("db/umi_counts/", "summary.txt$", full.names = T))
stats <- stats[, fread(file), file]
stats[, name:= gsub("_summary.txt", "", basename(file))]
stats[, Cc:= ifelse(is.na(meta[grep(name, output_prefix), DESeq2_group]), "red", "green"), name]

dir.create("pdf/alignment", 
           showWarnings = F)

pdf("pdf/alignment/alignment_statistics.pdf", 
    height = nrow(stats)/5, 
    width = 10)
par(mar= c(7,30,2,2))
barplot(t(stats[, .(umi_collapsed_reads, total_reads)]), 
        beside = T,
        col= unlist(lapply(stats$Cc, function(x) c(x, "white"))), 
        names.arg = basename(stats$name),
        cex.names= 0.5, 
        # col.names= stats$Cc,
        horiz = T, 
        las= 1)
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
# +annotation based on library
#--------------------------------------------------------------#
dir.create(dir_mergedCounts, showWarnings = F)
meta[!is.na(DESeq2_group), {
  counts <- paste0(dir_mergedCounts, "/", DESeq2_group, "_", cdition, "_rep", DESeq2_pseudo_rep, "_merged.txt")
  if(file.exists(counts))
    print(paste0(counts, " -->> ALREADY EXISTS"))
  else
  {
    print(paste0("Start ", counts))
    files <- paste0(dir_allCounts, output_prefix, ".txt")
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
    # SAVE
    fwrite(.c, counts)
    print(paste0(counts, " -->> DONE"))
  }
}, .(DESeq2_group, cdition, DESeq2_pseudo_rep, vllib, Spike_in)]

#--------------------------------------------------------------#
# Method #1 using counts norm
#--------------------------------------------------------------#
dir.create("db/final_tables_exp_model/counts_norm/", 
           showWarnings = F)
meta[!is.na(DESeq2_group), {
  FC <- paste0(dir_final, "/counts_norm/", DESeq2_group, "_counts_norm_final_oe.txt")
  if(file.exists(FC))
    print(paste0(FC, "  -->> ALREADY EXISTS"))
  else
  {
    # Import counts
    file <- list.files(dir_mergedCounts, DESeq2_group, full.names = T)
    .c <- lapply(file, fread)
    names(.c) <- gsub(".*_(.*)_(.*)_merged.txt", "\\1_\\2", basename(file))
    .c <- rbindlist(.c, idcol = T)
    # Cast and make matrix
    mat <- dcast(.c[type!="switched"], L+R~.id, value.var = "umi_counts", fill= 0)
    mat <- mat[rowSums(mat[, -c(1,2)])>=10]
    # melt
    dat <- melt(mat, id.vars = c("L", "R"))
    # dat <- dat[L!=R] # COULD BE USEFUL (SAFER?)!!!!!!!!!!!!
    dat[, c("cdition", "rep"):= tstrsplit(variable, "_rep")]
    # Check that input and screen both exist for each rep
    if(any(dat[, length(unique(cdition)), rep]$V1 != 2))
      print(paste0(DESeq2_group, " has missing replicates --> SKIPED")) else
        dat[, norm:= (value+1)/sum(value)*1e6, .(rep, cdition)]
    # Compute FC and scale
    res <- unique(dat[, .(L, R, rep)])
    res$FC <- dat[cdition!="input", norm]/dat[cdition=="input", norm]
    res <- res[, .(log2FoldChange= log2(mean(FC))), .(L, R)]
    # Subtract basal activity
    res[, log2FoldChange:= log2FoldChange-mean(res[grepl("^control", L) & grepl("^control", R), log2FoldChange])]
    # Compute expected
    res[, ctl_L:= grepl("^control", L)]
    res[, ctl_R:= grepl("^control", R)]
    res[, median_L:= ifelse(length(which(ctl_R))>=5, median(log2FoldChange[ctl_R], na.rm = T), as.numeric(NA)), L]
    res[, median_R:= ifelse(length(which(ctl_L))>=5, median(log2FoldChange[ctl_L], na.rm = T), as.numeric(NA)), R]
    # SAVE
    fwrite(res, FC)
    print(paste0(FC, "  -->> DONE"))
  }
}, DESeq2_group]

#--------------------------------------------------------------#
# Method # 2 using DESeq2
#--------------------------------------------------------------#
# dds --------------#
dir.create(dir_dds, showWarnings = F)
meta[!is.na(DESeq2_group), {
  for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
  {
    dds_file <- paste0(dir_dds, "/", DESeq2_group, "_", stringency, "_dds.rds")
    if(file.exists(dds_file))
      print(paste0(dds_file, "  -->> ALREADY EXISTS"))
    else
    {
      # Import counts
      file <- list.files(dir_mergedCounts, DESeq2_group, full.names = T)
      .c <- lapply(file, fread)
      names(.c) <- gsub(".*_(.*)_(.*)_merged.txt", "\\1_\\2", basename(file))
      .c <- rbindlist(.c, idcol = T)
      # Cast and make DF
      mat <- dcast(.c[type!="switched"], L+R~.id, value.var = "umi_counts", fill= 0)
      DF <- data.frame(mat[, -c(1,2)])
      rownames(DF) <- paste0(mat$L, "_vs_", mat$R)
      if(stringency=="low_cutoff")
        DF <- DF[rowSums(DF)>5,]
      if(stringency=="medium_cutoff")
        DF <- DF[rowSums(DF)>25,]
      if(stringency=="high_cutoff")
        DF <- DF[rowSums(DF)>5 & apply(DF, 1, function(x) all(x)>0),]
      if(stringency=="max_cutoff")
        DF <- DF[rowSums(DF)>25 & apply(DF, 1, function(x) all(x)>0),]
      # SampleTable
      sampleTable <- data.frame(cdition= unlist(tstrsplit(colnames(DF), "_rep", keep= 1)),
                                rep= unlist(tstrsplit(colnames(DF), "_rep", keep= 2)))
      rownames(sampleTable) <- colnames(DF)
      # DESeq2
      dds <- try(DESeq2::DESeqDataSetFromMatrix(countData= DF, 
                                                colData= sampleTable, 
                                                design= ~rep+cdition), 
                 silent = T)
      if(class(dds)=="try-error")
        print(paste0(DESeq2_group, " could not pass DESeqDataSetFromMatrix command and will be skipped -->> Check N replicates?"))
      else
      {
        # SizeFactors
        ctls_idx <- grep("control.*vs.*control", rownames(DF))
        sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[ctls_idx,]))
        # dds result
        res <- try(DESeq2::DESeq(dds), silent = T)
        if(class(res)=="try-error")
          print(paste0("No replicates --> ", DESeq2_group, "SKIPED"))
        else
        {
          saveRDS(res, dds_file) 
          print(paste0(dds_file, "  -->> DONE"))
        }
      }
    }
  }
}, DESeq2_group]

# Differential expression----------#
dir.create(dir_FC, showWarnings = F)
meta[!is.na(DESeq2_group), {
  for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
  {
    FC <- paste0(dir_FC, "/", DESeq2_group, "_", stringency, "_log2FC.txt")
    if(file.exists(FC))
      paste0(FC, "  -->> ALREADY EXISTS")
    else
    {
      # Check if models exists
      dds_file <- paste0(dir_dds, "/", DESeq2_group, "_", stringency, "_dds.rds")
      if(!file.exists(dds_file))
        print(paste0(FC, "  -->> No model --> ", DESeq2_group, " -->> SKIPPED"))
      else
      {
        # Import dds
        dds <- readRDS(dds_file)
        # Compoute FC
        .c <- DESeq2::results(dds, contrast= c("cdition", "DSCP", "input"))
        # Reformat
        .c <- as.data.table(as.data.frame(.c), keep.rownames= T)
        .c[, c("L", "R"):= tstrsplit(rn, "_vs_")]
        .c <- .c[, !"rn"]
        setcolorder(.c, c("L", "R"))
        # save
        fwrite(.c, FC, col.names = T, row.names = F, sep= "\t", quote= F)
        print(paste0(FC, "  -->> DONE!"))
      }
    }
  }
}, DESeq2_group]

# Compute Left and right activities -----------------------#
dir.create(dir_final, showWarnings = F)
meta[!is.na(DESeq2_group), {
  for(stringency in c("low_cutoff", "medium_cutoff", "high_cutoff", "max_cutoff"))
  {
    final <- paste0(dir_final, "/DESeq2/", DESeq2_group, "_", stringency, "_final_oe.txt")
    if(file.exists(final))
      paste0(final, "  -->> ALREADY EXISTS")
    else
    {
      FC_file <- paste0(dir_FC, "/", DESeq2_group, "_", stringency, "_log2FC.txt")
      if(!file.exists(FC_file))
        print(paste0(FC_file, " does not exist. -->> SKIPPED"))
      else
      {
        dat <- fread(FC_file)
        if(type=="peSTARRSeq")
        {
          #---------- Select robust controls -----------#
          cor <- data.table(ID= unique(grep("^control", c(dat$L, dat$R), value = T)))
          # Compute L/R PCC for each control fragment
          cor[, PCC:= {
            .c <- merge(dat[.BY, .(enh= L, L= log2FoldChange), on= "R==ID"],
                        dat[.BY, .(enh= R, R= log2FoldChange), on= "L==ID"])
            .c <- na.omit(.c)
            if(nrow(.c)>=5)
              cor.test(.c$L, .c$R)$estimate
            else
              as.numeric(NA)
          }, ID]
          # Compute outliers frequency
          out <- copy(dat)
          out[grepl("^control", L), out_R:= log2FoldChange %in% boxplot(log2FoldChange, plot= F)$out, R]
          out[grepl("^control", L), out_fL:= length(which(out_R))/.N, L]
          out[grepl("^control", R), out_L:= log2FoldChange %in% boxplot(log2FoldChange, plot= F)$out, L]
          out[grepl("^control", R), out_fR:= length(which(out_L))/.N, R]
          # Robust controls list
          dat[, ctl_L:= L %in% cor[PCC>0.75, ID] & L %in% out[out_fL<0.025, L]]
          dat[, ctl_R:= R %in% cor[PCC>0.75, ID] & R %in% out[out_fR<0.025, R]]
          #---------- Compute expected values -----------#
          dat[, length(which(ctl_R)), L]
        }else if(type=="revpeSTARRSeq")
        {
          dat[, ctl_L:= grepl("^control", L)]
          dat[, ctl_R:= grepl("^control", R)]
        }
        dat[, median_L:= ifelse(length(which(ctl_R))>=5, median(log2FoldChange[ctl_R], na.rm = T), as.numeric(NA)), L]
        dat[, median_R:= ifelse(length(which(ctl_L))>=5, median(log2FoldChange[ctl_L], na.rm = T), as.numeric(NA)), R]
        fwrite(dat, final)
        print(final, "  -->> DONE")
      }
    }
  }
}, .(DESeq2_group, type)]
