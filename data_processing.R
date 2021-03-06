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
require(seqinr)
require(rtracklayer)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(motifmatchr)
require(PWMEnrich)
require(TFBSTools)
require(seqLogo)
require(motifmatchr)
require(kohonen)
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

# Merge vl BA libraries ####
if(!file.exists("Rdata/uniq_library_final.rds"))
{
  # Merge BA and my TWIST into a single non-redundant clean object
  vl <- as.data.table(readRDS("Rdata/vl_library_112019.rds"))
  BA <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
  lib <- merge(BA[, .(ID_BA= ID, BA_group= Group, BA_enhancer_group= enhancer_group, BA_enh_group_detail= enhancer_group_detail, 
                      dev_log2FoldChange, hk_log2FoldChange, seqnames, start, end, width, strand)], 
               vl[, .(ID_BA= BA_ID, ID_vl, group, detail, linker_ID, seqnames, start, end, width, strand)], all.x= T, all.y= T)
  # Uniq IDs
  lib[, vl:= ifelse(!is.na(ID_vl), T, F)]
  lib[, uniq_id:= ifelse(is.na(ID_vl), as.character(ID_BA), as.character(ID_vl))]
  # Uniq groups
  lib[BA_enhancer_group=="dev", uniq_group:= "dev"]
  lib[BA_enhancer_group=="hk", uniq_group:= "hk"]
  lib[BA_enhancer_group=="Inducible", uniq_group:= "inducible"]
  lib[BA_enhancer_group=="shared", uniq_group:= "shared"]
  lib[BA_enhancer_group=="Controls", uniq_group:= "control"]
  lib[group=="OSC", uniq_group:= "OSC"]
  lib[group=="control", uniq_group:= "control"]
  # Uniq details
  lib[uniq_group== "dev", uniq_detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,10), labels = c("inactive", "weak", "medium", "strong"))]
  lib[uniq_group== "shared", uniq_detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,10), labels = c("inactive", "weak", "medium", "strong"))]
  lib[uniq_group== "hk", uniq_detail:= cut(hk_log2FoldChange, quantile(hk_log2FoldChange, seq(0, 1, length.out= 5)), include.lowest = T, labels = c("inactive", "weak", "medium", "strong"))]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "Ecdysone", uniq_detail:= "ecdysone"]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "HeatShock", uniq_detail:= "heatshock"]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "Cadmium", uniq_detail:= "cadmium"]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "Paraquat", uniq_detail:= "paraquat"]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "PGN", uniq_detail:= "PGN"]
  lib[uniq_group== "inducible" & BA_enh_group_detail== "Wnt", uniq_detail:= "wnt"]
  lib[uniq_group== "OSC", uniq_detail:= "OSC_specific"]
  lib[uniq_group== "control" & detail== "exon", uniq_detail:= "exon"]
  lib[uniq_group== "control" & detail== "Ecoli", uniq_detail:= "ecoli"]
  lib[uniq_group== "control" & is.na(uniq_detail), uniq_detail:= "flat_genomic_region"]
  # colors
  clean <- lib[, .(ID= uniq_id, group= uniq_group, detail= uniq_detail, vl, linker_ID, seqnames, start, end, strand, dev_log2FoldChange, hk_log2FoldChange)]
  class_Cc <- data.table(group= c("hk", "shared", "dev", "OSC", "inducible", "control"),
                         col= c("tomato", "royalblue2", "#74C27A", "black", "gold", "lightgrey"))
  clean <- clean[class_Cc, , on= "group"]
  # SAVE
  saveRDS(clean, "Rdata/uniq_library_final.rds")
}
####
# Motifs clustering  SOM ####
lib <- readRDS("Rdata/uniq_library_final.rds")
lib <- lib[!detail=="ecoli"]
# Count hits low cutoff
if(!file.exists("Rdata/master_motifs_counts.rds"))
{
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds, GRanges(lib$seqnames, IRanges(lib$start, lib$end)), 
                     genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds)
  rownames(counts) <- lib$ID
  counts <- as.data.table(counts, keep.rownames = T)
  saveRDS(counts, "Rdata/master_motifs_counts.rds")
}
# Select motifs enriched in at least one of the lib groups
if(!file.exists("Rdata/enriched_motifs_counts.rds"))
{
  counts <- readRDS("Rdata/master_motifs_counts.rds")
  
  dat <- melt(counts, id.vars= "rn")
  dat[lib, group:= i.group, on= "rn==ID"]
  sel <- dat[group!="control", any(value>0), .(variable, group)][(V1), .(variable, group)]
  setkey(sel, group)
  
  mot <- list()
  for(c_gr in c("dev", "hk", "inducible", "OSC"))
  {
    sub <- dat[group %in% c(c_gr, "control") & variable %in% sel[c_gr, variable]]
    sub <- sub[, fisher.test(table(group!="control", value>0))[c("estimate", "p.value")], variable]
    mot[[c_gr]] <- sub[estimate>2 & p.value<0.00001, variable]
  }
  enriched_motifs_counts <- dat[variable %in% unique(unlist(mot))]
  saveRDS(enriched_motifs_counts, "Rdata/enriched_motifs_counts.rds")
}
# SOM clustering
if(!file.exists("Rdata/som_enriched_motifs.rds"))
{
  dat <- readRDS("Rdata/enriched_motifs_counts.rds")
  dmat <- dcast(dat, variable~rn, value.var = "value")
  mat <- log2(as.matrix(dmat, 1)+1)
  mat <- scale(mat) 
  mat <- apply(mat, 2, function(x)
  {
    lim <- quantile(x, 0.025, 0.975)
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    return(x)
  })
  mygrid <- somgrid(xdim= 6, ydim= 6, topo = 'hexagonal', toroidal = T)
  set.seed(1234)
  som.model <- supersom(mat, grid = mygrid, rlen = 500)
  
  # Diag plots
  pdf("pdf/SOM_motifs_clustering_diag_plots.pdf")
  par(mfrow= c(4, 4))
  for(type in c("changes", "dist.neighbours", "counts", "mapping", "quality"))
  {
    plot(som.model, type= type, palette.name = colorRampPalette(c("blue", "yellow")), shape= "straight", border= "black")
  }
  dev.off()
  
  # Heatmap
  pdf("pdf/SOM_motifs_clustering_heatmap.pdf")
  pl <- mat[order(som.model$unit.classif),]
  my_heatmap(pl, cluster_rows = F, col = colorRampPalette(c("blue", "yellow", "white"))(100), plot_dendro_col = F)
  abline(h= 1-which(diff(som.model$unit.classif[order(som.model$unit.classif)])!=0)/nrow(mat), lwd= 1, col= "white")
  dev.off()
  
  # Identify representative motifs
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  cl <- data.table(cl= som.model$unit.classif, 
                   dist= som.model$distances, 
                   motif= rownames(som.model$data[[1]]))
  cl[, name:= TF_clusters_PWMs$metadata$Dmel[match(motif, TF_clusters_PWMs$metadata$motif_name)]]
  cl[, best_match:= motif[which.min(dist)], cl]
  cl[, Dmel:= paste0(unique(na.omit(name)), collapse= "__"), cl]
  cl[, check:= gsub("A|C|G|T|R|Y|S|W|K|M|B|D|H|V|N", "", name)=="", name]
  cl[, Dmel_s:= paste0(unique(na.omit(.SD[!(check), name])), collapse= "_"), keyby= cl]
  cl[, Dmel_s:= paste0(.GRP, ".", Dmel_s), keyby= .(cl, Dmel_s)]
  
  # SAVE
  som.model$info <- cl
  saveRDS(som.model, "Rdata/som_enriched_motifs.rds")
}
####

# Compute enhancer features ####
if(!file.exists("Rdata/master_lib_features.rds"))
{
  lib <- as.data.table(readRDS("Rdata/uniq_library_final.rds"))
  # Closest promoter
  gene <- import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf")
  seqlevelsStyle(gene) <- "UCSC"
  gene <- as.data.table(gene)[type=="gene", .(seqnames, start, end, strand, gene_id, symbol= gene_name)]
  tss <- copy(gene)
  tss[strand=="+", end:= start]
  tss[strand=="-", start:= end]
  ctss <- tss[lib[, .(seqnames, start, end, ID)], .(ID, gene_id= gene_id[which.min(abs(start-((i.start+i.end)/2)))]), .EACHI, on= "seqnames"]
  ctss <- ctss[gene, , on= "gene_id", nomatch= NULL]
  ctss <- ctss[, .(ID, gene_id, symbol, gene_coor= paste0(i.seqnames, ":", start, "-", end, ":", strand))]
  lib <- ctss[lib, , on= "ID", nomatch= NA]
  
  # Multiple enhancer promoters
  selg <- gene[c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo"), , on= "symbol"]
  selg <- selg[, .(seqnames, start= start-5000, end= end+5000, strand, gene_id, symbol)]
  setkeyv(selg, c("seqnames", "start", "end"))
  ov <- foverlaps(lib, selg, nomatch = NULL)[, .(ID, me_gene_id= gene_id, me_symbol= symbol)]
  lib <- merge(lib, ov, by= "ID", all.x= T)
  
  # STARR-Seq enrichment
  if(!file.exists("Rdata/master_STARR-Seq_counts_table.rds"))
  {
    STARR <- data.table(file= c("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_DSCP_200bp_cut.bed",
                                "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed", 
                                "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed",
                                "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/input_RPS12_200bp_cut.bed",
                                "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed", 
                                "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/raw_data/RpS12_200bp_gw_Rep1.UMI_cut.bed",
                                "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_input.merged.uniq.bed",
                                "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_1.UMI.bed", 
                                "/groups/stark/vloubiere/projects/rep_STARRSeq_Lorena/db/sgl_gw_200bp/sgl_DSCP.libLOH009_enh10-2_DSCP_gw200_rep_2.UMI.merged.bed"), 
                        CP= c("DSCP", "DSCP", "DSCP", "RPS12", "RPS12", "RPS12", "DSCP", "DSCP", "DSCP"),
                        cdition= c("input", "screen", "screen", "input", "screen", "screen", "input", "screen", "screen"),
                        assay= c("STARR-Seq","STARR-Seq","STARR-Seq","STARR-Seq","STARR-Seq","STARR-Seq","rep-STARR-Seq","rep-STARR-Seq","rep-STARR-Seq"),
                        rep=  c("rep1","rep1","rep2","rep1","rep1","rep2","rep1","rep1","rep2"),
                        name= rep(c("devSTARRSeq_log2FC", "hkSTARRSeq_log2FC", "sglrepSTARRSeq_log2FC"), each= 3))
    setkeyv(lib, c("seqnames", "start", "end"))
    STARR <- STARR[, {
      .c <- fread(file)
      colnames(.c)[1:3] <- c("seqnames", "start", "end")
      .c <- foverlaps(.c, lib, nomatch= NULL)
      .c[, .N, ID]
    }, STARR]
    STARR <- STARR[, .(count= sum(N)), .(name, cdition, ID)]
    STARR <- STARR[, {
      .c <- dcast(.SD, ID~cdition, value.var = "count")
      .c <- .c[input+screen>100, .(ID, log2FC= log2(screen)-log2(input))]
    }, name]
    STARR <- dcast(STARR, ID~name, value.var= "log2FC")
    saveRDS(STARR, "Rdata/master_STARR-Seq_counts_table.rds")
  }else{STARR <- readRDS("Rdata/master_STARR-Seq_counts_table.rds")}
  lib <- STARR[lib, , on= "ID"]
  
  # ChIP-Seq enrichment
  if(!file.exists("Rdata/master_ChIP_counts_table.rds"))
  {
    # Compute ChIP-Seq ATAC-Seq enrichment
    avail <- fread("/groups/stark/vloubiere/projects/available_data_dm3/db/metadata/available_data_metadata.txt")
    avail[, file:= list.files("/groups/stark/vloubiere/projects/available_data_dm3/db/bed/", paste0(uniq_name, ".*.bed"), full.names = T), uniq_name]
    avail <- avail[sample!="input", .(file, cdition= sample)]
    setkeyv(lib, c("seqnames", "start", "end"))
    avail <- avail[, {
      .c <- fread(file)
      colnames(.c)[1:3] <- c("seqnames", "start", "end")
      .c <- foverlaps(.c, lib, nomatch= NULL)
      .c[, .N, ID]
    }, avail]
    avail <- avail[, .(count= sum(N)), .(cdition, ID)]
    avail[, cdition:= paste0(cdition, "_log2FC"), cdition]
    avail[, log2FC:= log2(count)-log2(mean(.SD[grepl("control|_NegativeRegions", ID), count])), cdition]
    avail <- dcast(avail, ID~cdition, value.var = "log2FC")
    saveRDS(avail, "Rdata/master_ChIP_counts_table.rds")
  }else{avail <- readRDS("Rdata/master_ChIP_counts_table.rds")}
  lib <- avail[lib, , on= "ID"]
  
  # Counts most informative motifs
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  som <- readRDS("Rdata/som_enriched_motifs.rds")
  sel <- match(unique(som$info$best_match), name(TF_clusters_PWMs$All_pwms_log_odds))
  # counts
  sublib <- copy(lib)
  sublib <- sublib[!detail=="ecoli"]
  hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], GRanges(sublib$seqnames, IRanges(sublib$start, sublib$end)), 
                     genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds[sel])
  rownames(counts) <- sublib$ID
  counts <- as.data.table(counts, keep.rownames = T)
  colnames(counts)[-1] <- paste0("motif__", colnames(counts)[-1])
  # add to lib and add ecoli sequences
  sublib <- cbind(sublib, counts[, !"rn"])
  lib <- rbind(sublib, lib[detail=="ecoli"], fill= T)
  
  # SAVE
  setcolorder(lib, c("ID", "group", "detail", "vl", "linker_ID", "col", "seqnames", "start", "end", "strand", 
                     "gene_id", "symbol", "gene_coor",  "me_gene_id", "me_symbol", 
                     "dev_log2FoldChange", "hk_log2FoldChange"))
  saveRDS(lib, "Rdata/master_lib_features.rds")
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
pa[, cdition:= gsub("^[^_]*_(.*)", "\\1", cdition), cdition]
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
clean[, clean:= paste0("db/count/", lib, "_", cdition, "_", rep, "_clean_umi_counts.txt")]
clean[!file.exists(clean), {
  .c <- fread(file)
  if(lib=="libvl002"){
    .c <- .c[(grepl("^ts", L) & L==R) | (!grepl("^ts", L) & !grepl("^ts", R))]
  }
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
  fwrite(.c, clean)
}, .(file, clean)]
####

# Check saturation ####
if(!file.exists("pdf/sequencing_saturation.pdf"))
{
  sat <- data.table(file= list.files("db/count", "clean_umi_counts.txt$", full.names = T))
  sat[, sample:= tstrsplit(basename(file), "_clean", keep= 1)]
  sat[, lib:= tstrsplit(sample, "_", keep= 1)]
  sat <- sat[, fread(file), sat]
  pdf("pdf/sequencing_saturation.pdf", 5, 5)
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
}
####

# Spike-ins percentage ####
if(!file.exists("pdf/spike_ins_percentage_barplot.pdf"))
{
  spike <- data.table(file= list.files("db/count/", "_clean_umi_counts.txt", full.names = T))
  spike[, lib:= tstrsplit(basename(file), "_", keep=1), file]
  spike <- spike[lib != "vllib012"]
  spike <- spike[, fread(file), spike]
  spike[grepl("^libvl002", lib) & grepl("^ts_", L), spike:= T]
  spike[!grepl("^libvl002", lib) & (grepl("^ts_|_C_", L) | grepl("^ts_|_C_", R)), spike:= T]
  spike[is.na(spike), spike:= F]
  spike <- spike[, .(perc= sum(.SD[(spike), count])/sum(.SD[, count])*100), lib]
  
  pdf("pdf/spike_ins_percentage_barplot.pdf", width = 2.5, height = 3)
  par(mar= c(4,4,1,1), las= 2)
  barplot(spike$perc, names.arg = spike$lib, ylim = c(0, 10), ylab= "% Spike-ins")
  dev.off()
}
####
# Template switching ####
# Compute
if(!file.exists("Rdata/master_ts_table.rds"))
{
  # SampleTable
  sampleTable <- data.table(file= list.files("db/count", "_all_", full.names = T))
  sampleTable[, lib:= tstrsplit(basename(file), "_", keep= 1)]
  sampleTable[, cdition:= ifelse(grepl("input", file), "input", "DSCP")]
  sampleTable[, rep:= tstrsplit(basename(file), "DSCP_|input_|_all", keep = 2)]
  # build tables for counting
  seq <- names(read.fasta("/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq/peSTARRSeq_sequences.fa"))
  seq <- CJ(L= seq, R= seq)
  ts <- c("ts_SCR2_01002", "ts_HAM1_01003", "ts_SUP1_01004")
  tab <- sampleTable[, {
    if(lib=="libvl002"){seq[L %in% ts | R %in% ts, .(L, R, ts= ifelse(L %in% ts & L==R, "ok", "switched"))]}
    else{seq[L %in% ts | R %in% ts | grepl("_C_", L) | grepl("_C_", R), 
             .(L, R, ts= ifelse((L=="ts_SCR2_01002" & grepl("_C_", R)) | (R=="ts_SCR2_01002" & grepl("_C_", L)), "ok", "switched"))]}
  }, sampleTable]
  # Percentage barplot
  raw <- sampleTable[, fread(file), sampleTable]
  tab[raw, count:= i.count, on= c("lib", "cdition", "rep", "L", "R")]
  tab <- na.omit(tab)
  saveRDS(tab, "Rdata/master_ts_table.rds")
}

# barplot percentage template switching
if(!file.exists("pdf/template_switching_barplot.pdf"))
{
  tab <- readRDS("Rdata/master_ts_table.rds")
  res <- tab[, .(count= sum(count, na.rm= T)), .(lib, cdition, rep, ts)]
  res <- dcast(res, lib+cdition+rep~ts, value.var = "count")
  res[, perc:= switched/(ok+switched)*100]
  
  pdf("pdf/template_switching_barplot.pdf", height = 6)
  par(mar= c(10,4,2,2), las= 2)
  bar <- barplot(res$perc, ylim= c(0, 20), ylab= "% template switching", 
                 names.arg = res[, paste0(lib, "_", cdition, "_" , rep)])
  text(bar, res$perc, round(res$perc, 1), pos= 3)
  dev.off()
}

# Scatterplot vllib 014 bad and "expected" reads
if(!file.exists("pdf/template_switching_scatterplot_vllib014.pdf"))
{
  tab <- readRDS("Rdata/master_ts_table.rds")
  sub <- tab[lib=="vllib014" & cdition=="DSCP" & rep=="rep1" & !is.na(count)]
  setorderv(sub, "count")
  switched <- sub[ts=="switched"]
  ok <- sub[ts=="ok"]
  switched[ok, correct_L:= i.count, on= "L"]
  switched[ok, correct_R:= i.count, on= "R"]
  
  x <- switched[, log2(count)]
  
  pdf("pdf/template_switching_scatterplot_vllib014.pdf", 4.5, 4.5)
  layout(matrix(c(2,4,1,3), ncol=2, byrow = T), heights = c(0.25,1), widths = c(1, 0.25))
  for(dir in c("correct_L", "correct_R"))
  {
    y <- log2(switched[[dir]])
    par(mar=c(4,4,1,1), pch= 16, las= 1)
    plot(x, y, xlim= c(0, 12), ylim= c(0, 12), col= adjustcolor("lightgrey", 0.5),
         xlab= "Switched counts", ylab= "Putative correct counts")
    abline(0, 1)
    abline(1, 1, lty= 2)
    abline(-1, 1, lty= 2)
    par(mar=c(1,4,1,1))
    plot(density(x), xlim= c(0, 12), xaxt= "n", ylab= "density", main= "")
    .dens <- density(na.omit(y))
    par(mar=c(4,1,1,1))
    plot(.dens$y, .dens$x, type= "l", ylim= c(0, 12), yaxt= "n", xlab= "density", main= "")
    plot.new()
  }
  dev.off()
}
####

# Differential analysis ####
# 1- Make master counts table
if(!file.exists("Rdata/master_counts_table.txt"))
{
  raw <- data.table(file= list.files("db/count/", "clean_umi_counts", full.names = T))
  raw[, c("lib", "cdition", "rep"):= tstrsplit(basename(file), "_", keep= 1:3)]
  raw <- raw[, fread(file), raw]
  raw <- dcast(raw, L+R~lib+cdition+rep, value.var = "count")
  
  # Clean counts table
  counts <- data.table(raw[, .(L, R)],
                       libvl002_DSCP_rep1= raw$libvl002_DSCP_rep1+raw$libvl002_DSCP_rep2,
                       libvl002_DSCP_rep2= raw$libvl002_DSCP_rep3+raw$libvl002_DSCP_rep4,
                       libvl002_input_rep1= raw$libvl002_input_rep1+raw$libvl002_input_rep3+raw$libvl002_input_rep4,
                       libvl002_input_rep2= raw$libvl002_input_rep2+raw$libvl002_input_rep5+raw$libvl002_input_rep6,
                       libvl013_input_rep1= raw$vllib013_input_rep1,
                       libvl013_input_rep2= raw$vllib013_input_rep2,
                       libvl013_DSCP_rep1= raw$vllib013_DSCP_rep1,
                       libvl013_DSCP_rep2= raw$vllib013_DSCP_rep2,
                       libvl014_input_rep1= raw$vllib014_input_rep1,
                       libvl014_input_rep2= raw$vllib014_input_rep2,
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
  DE[, {
    if(!file.exists(dds_file))
    {
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
    }
  }, .(lib, dds_file)]
}
# make diff table
DE[, diff:= gsub("_dds.rds$", "_DE.txt", dds_file)]
DE[, {
  if(!file.exists(diff))
  {
    .c <- readRDS(dds_file)
    # Differential expression
    .c <- as.data.table(as.data.frame(results(.c, contrast= c("cdition", "DSCP", "input"))), keep.rownames= T)
    .c[, c("L", "R"):= tstrsplit(rn, "_vs_")]
    .c <- .c[, !"rn"]
    setcolorder(.c, c("L", "R"))
    fwrite(.c, diff, col.names = T, row.names = F, sep= "\t", quote= F)
    print(paste0(diff, " DONE!"))
  }
}, .(dds_file, diff)]
####

# Compute expected score ####
FC <- data.table(file= list.files("db/DE_analysis/", "DE.txt$", full.names = T))
FC[,  lib:= gsub("_DE.txt$", "", basename(file))]
FC[,  FC_file:= paste0("db/DE_analysis/", lib, "_add_scores.txt")]
FC[, {
  if(!file.exists(FC_file))
  {
    .c <- fread(file)
    .c[, spike_in:= ifelse(grepl("^ts", L) | grepl("^ts", R), T, F)]
    .c[, median_L:= .SD[grep("control", R), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], L]
    .c[grepl("_C_", L) & R=="ts_SCR2_01002", median_L:= log2FoldChange, L]
    .c[, median_R:= .SD[grep("control", L), ifelse(.N<5, as.numeric(NA), median(log2FoldChange, na.rm = T))], R]
    .c[L=="ts_SCR2_01002" & grepl("_C_", R), median_R:= log2FoldChange, R]
    .c <- na.omit(.c[!grepl("control", L) & !grepl("control", R)])
    .c[, add:= log2(sum(2^median_L+2^median_R)), .(L, R)]
    fwrite(.c, FC_file, col.names = T, row.names = F, sep= "\t", quote= F)
  }
}, .(file, FC_file)]
if(!file.exists("Rdata/master_results_peSTARRSeq.rds"))
{
  .c <- FC[, fread(FC_file), FC]
  .c$file <- NULL
  .c$FC_file <- NULL
  saveRDS(.c, "Rdata/master_results_peSTARRSeq.rds")
}
####

# Comparison Bernardo ####
if(!file.exists("pdf/comparison_BA.pdf"))
{
  lib <- readRDS("Rdata/uniq_library_final.rds")
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
# Smoothscatter models ####
if(!file.exists("Rdata/linear_models_prediction.rds"))
{
  feat <- readRDS("Rdata/master_lib_features.rds")
  cols <- c("ID", grep("^motif", colnames(feat), value= T))
  feat <- feat[, ..cols]
  tab <- readRDS("Rdata/master_results_peSTARRSeq.rds")
  tab <- tab[!(spike_in)]
  tab <- feat[tab, , on= "ID==L"]
  tab <- merge(tab, feat, by.x= "R", by.y= "ID", all.x= T, suffixes= c("___L", "___R"))
  cols <- unique(gsub("___L|___R", "", grep("^motif", colnames(tab), value = T)))
  tab[, (paste0(cols, "___sum")):= lapply(cols, function(x) rowSums(tab[, grep(x, colnames(tab)), with= F], na.rm= T))]
  mots1mots2 <- paste0(c(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T)), collapse= "+")
  motssum <- paste0(grep("___sum", colnames(tab), value = T), collapse= "+")
  motspw <- CJ(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T))[, paste0(V1, "*", V2)]
  motspw <- paste0(motspw, collapse= "+")
  models <- data.table(lib= unique(tab$lib))
  models <- models[, .(title= c("lm_L", 
                                "lm_R", 
                                "lm_add", 
                                "lm_L+R", 
                                "lm_L*R", 
                                "lm_L*R+add", 
                                "lm_L*R*add", 
                                "lm_mots1+mot2", 
                                "lm_motsums", 
                                "lm_motspw", 
                                "lm_L+R+mots1+mot2", 
                                "lm_L+R+motsums", 
                                "lm_L+R+motspw")), lib]
  models[, file:= paste0("db/models/", lib, "_", title, ".rds"), models]
  if(any(!file.exists(models$file)))
  {
    models[, formula:= c("log2FoldChange~median_L",
                         "log2FoldChange~median_R",
                         "log2FoldChange~add",
                         "log2FoldChange~median_L+median_R",
                         "log2FoldChange~median_L*median_R",
                         "log2FoldChange~median_L*median_R+add",
                         "log2FoldChange~median_L*median_R*add",
                         paste0("log2FoldChange~", mots1mots2),
                         paste0("log2FoldChange~", motssum),
                         paste0("log2FoldChange~", motspw),
                         paste0("log2FoldChange~median_L+median_R+", mots1mots2),
                         paste0("log2FoldChange~median_L+median_R+", motssum),
                         paste0("log2FoldChange~median_L+median_R+", motspw)), lib]
    setkeyv(tab, "lib")
    models[, {if(!file.exists(file)){saveRDS(lm(as.formula(formula), tab[lib]), file)}}, models]
    models$formula <- NULL
  }
  models[, rsq:= summary(readRDS(file))$r.squared, file]
  models[, obs:= .(list(readRDS(file)$model$log2FoldChange)), file]
  models[, exp:= .(list(predict(readRDS(file)))), file]
  models[, PCC:= cor.test(unlist(obs), unlist(exp))$estimate, .(lib, title)]
  saveRDS(models, "Rdata/linear_models_prediction.rds")
}else{
  models <- readRDS("Rdata/linear_models_prediction.rds")
}

pdf("pdf/linear_models_peSTARRSeq.pdf", width = 52, height = 4)
par(mfrow= c(1,13), las= 1)
models[, {
  smoothScatter(unlist(exp), unlist(obs), main= paste(lib, title))
  legend("topleft", c(paste0("R2= ", round(rsq, 2)), paste0("PCC= ", round(PCC, 2))), bty= "n")
  print("")
}, .(lib, title, PCC, rsq)]
dev.off()
####

# Smoothscatter models residuals ####
if(!file.exists("Rdata/linear_models_residuals_prediction.rds"))
{
  feat <- readRDS("Rdata/master_lib_features.rds")
  cols <- c("ID", grep("^motif", colnames(feat), value= T))
  feat <- feat[, ..cols]
  tab <- readRDS("Rdata/master_results_peSTARRSeq.rds")
  tab <- tab[!(spike_in)]
  tab <- feat[tab, , on= "ID==L"]
  tab <- merge(tab, feat, by.x= "R", by.y= "ID", all.x= T, suffixes= c("___L", "___R"))
  cols <- unique(gsub("___L|___R", "", grep("^motif", colnames(tab), value = T)))
  tab[, (paste0(cols, "___sum")):= lapply(cols, function(x) rowSums(tab[, grep(x, colnames(tab)), with= F], na.rm= T))]
  tab[, diff:= log2FoldChange-add]
  mots1mots2 <- paste0(c(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T)), collapse= "+")
  motssum <- paste0(grep("___sum", colnames(tab), value = T), collapse= "+")
  motspw <- CJ(grep("___L", colnames(tab), value = T), grep("___R", colnames(tab), value = T))[, paste0(V1, "*", V2)]
  motspw <- paste0(motspw, collapse= "+")
  models <- data.table(lib= unique(tab$lib))
  models <- models[, .(title= c("lm_L", 
                                "lm_R", 
                                "lm_add", 
                                "lm_L+R", 
                                "lm_L*R", 
                                "lm_L*R+add", 
                                "lm_L*R*add", 
                                "lm_mots1+mot2", 
                                "lm_motsums", 
                                "lm_motspw", 
                                "lm_L+R+mots1+mot2", 
                                "lm_L+R+motsums", 
                                "lm_L+R+motspw")), lib]
  models[, file:= paste0("db/models/", lib, "_", title, "_residuals.rds"), models]
  if(any(!file.exists(models$file)))
  {
    models[, formula:= c("diff~median_L",
                         "diff~median_R",
                         "diff~add",
                         "diff~median_L+median_R",
                         "diff~median_L*median_R",
                         "diff~median_L*median_R+add",
                         "diff~median_L*median_R*add",
                         paste0("diff~", mots1mots2),
                         paste0("diff~", motssum),
                         paste0("diff~", motspw),
                         paste0("diff~median_L+median_R+", mots1mots2),
                         paste0("diff~median_L+median_R+", motssum),
                         paste0("diff~median_L+median_R+", motspw)), lib]
    setkeyv(tab, "lib")
    models[, {if(!file.exists(file)){saveRDS(lm(as.formula(formula), tab[lib]), file)}}, models]
    models$formula <- NULL
  }
  models[, rsq:= summary(readRDS(file))$r.squared, file]
  models[, obs:= .(list(readRDS(file)$model$diff)), file]
  models[, exp:= .(list(predict(readRDS(file)))), file]
  models[, PCC:= cor.test(unlist(obs), unlist(exp))$estimate, .(lib, title)]
  saveRDS(models, "Rdata/linear_models_residuals_prediction.rds")
}else{
  models <- readRDS("Rdata/linear_models_residuals_prediction.rds")
}

pdf("pdf/linear_models_residuals_peSTARRSeq.pdf", width = 52, height = 4)
par(mfrow= c(1,13), las= 1)
models[, {
  smoothScatter(unlist(exp), unlist(obs), main= paste(lib, title))
  legend("topleft", c(paste0("R2= ", round(rsq, 2)), paste0("PCC= ", round(PCC, 2))), bty= "n")
  print("")
}, .(lib, title, PCC, rsq)]
dev.off()

####








# Close versus distant ####
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in)]
feat <- readRDS("Rdata/master_lib_features.rds")
dat[feat, c("seqnames_L", "pos_L", "gene_L"):= .(i.seqnames, round((i.end+i.start)/2), symbol), on= "L==ID"]
dat[feat, c("seqnames_R", "pos_R", "gene_R"):= .(i.seqnames, round((i.end+i.start)/2), symbol), on= "R==ID"]
dat[seqnames_L==seqnames_R, dist:= abs(pos_L-pos_R)]
dat[seqnames_L!=seqnames_R, dist:= Inf]
pdf("pdf/close_vs_distant_pairs.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
cutoff <- 30000
dat[, {
  close <- .SD[L!=R & dist<cutoff, .(L, R, median_L, median_R, diff= log2FoldChange-add)]
  rest <- .SD[dist>=cutoff]
  far <- data.table()
  for(i in seq(nrow(close)))
  {
    sel <- which.min(abs(close[i, median_L]-rest$median_L)+abs(close[i, median_R]-rest$median_R))
    far <- rbind(far, rest[sel])
    rest <- rest[-sel]
  }
  close <- close$diff
  far <- far[, log2FoldChange-add]
  boxplot(close, far, notch= T, names= c("close", "far"), main= lib)
  text(0.5, max(close), labels = length(close), pos= 4, offset = 0)
  text(1.5, max(far), labels = length(far), pos= 4, offset = 0)
  print("")
}, lib]
dev.off()
####
# multiple enhancer genes vs few enhancer genes ####
dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
dat <- dat[!(spike_in)]
feat <- readRDS("Rdata/master_lib_features.rds")
dat[feat, gene_L:= me_symbol, on= "L==ID"]
dat[feat, gene_R:= me_symbol, on= "R==ID"]
pdf("pdf/multiple_enhancer_genes_vs_few_enhancer.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
dat[, {
  close <- .SD[!is.na(gene_L) & gene_L==gene_R, log2FoldChange-add]
  far <- .SD[is.na(gene_L) | is.na(gene_R) | gene_L!=gene_R, log2FoldChange-add]
  boxplot(close, far, notch= T, names= c("mult. enh.", "few enh."), main= lib)
  text(0.5, max(close), labels = length(close), pos= 4, offset = 0)
  text(1.5, max(far), labels = length(far), pos= 4, offset = 0)
  print("")
}, lib]
dev.off()
####
# Homotypic vs heterotypic pairs ####
dat <- data.table(file= list.files("db/DE_analysis/", "DE.txt", full.names = T))
dat <- dat[, fread(file), file]
pdf("pdf/homotypic_pairs.pdf", height = 3.5)
par(mfrow= c(1, 3), las= 2, mar= c(7,4,2,2))
dat[, {
  controls <- .SD[grepl("control", L) & grepl("control", R)]
  enh <- .SD[!grepl("control", L) & !grepl("control", R)]
  boxplot(controls[L==R, log2FoldChange], 
          controls[L!=R, log2FoldChange], 
          enh[L==R, log2FoldChange], 
          enh[L!=R, log2FoldChange], 
          main= basename(file), notch= T,
          names= c("homo_ctl", "hetero_ctl", "homo_enh", "hetero_enh"))
  pval <- wilcox.test(enh[L==R, log2FoldChange], enh[L!=R, log2FoldChange])$p.value
  text(3.5, 10, round(pval, 3))
  print("")
}, file]
dev.off()
####
# Heatmap log2FC ####
if(!file.exists("pdf/Heatmaps.pdf"))
{
  dat <- readRDS("Rdata/master_results_peSTARRSeq.rds")
  dat <- dat[!(spike_in) & lib=="libvl013" & median_L>0.5 & median_R>0.5]
  dat[, c("L", "R"):= .(substr(L, nchar(L)-4, nchar(L)), substr(R, nchar(R)-4, nchar(R)))]
  dat[, diff:= log2FoldChange-add]
  mat1 <- dcast(dat, -median_L+L~median_R+R, value.var = "log2FoldChange", fill = NA)
  ind_L <- -(rev(mat1$median_L))
  ind_R <- as.numeric(unlist(tstrsplit(colnames(mat1), "_", keep=1))[-c(1,2)])
  colnames(mat1) <- unlist(tstrsplit(colnames(mat1), "_", keep = 2))
  mat1 <- as.matrix(mat1[, -1], 1)
  mat2 <- dcast(dat, -median_L+L~median_R+R, value.var = "diff", fill = NA)
  mat2 <- as.matrix(mat2[, -c(1, 2)])
  rownames(mat2) <- rownames(mat1)
  colnames(mat2) <- colnames(mat1)
  pdf("pdf/Heatmaps.pdf", 20, 20)
  layout(matrix(1:4, ncol= 2), heights = c(0.25,1), widths = c(1,0.25))
  #FC
  par(mar= c(1,4.1,1,5.1), xaxs= "i", yaxs= "i")
  barplot(ind_L)
  par(mar= c(4.1,4.1,1,5.1))
  my_heatmap(mat1, cluster_cols = F, cluster_rows = F, breaks = c(-2, 0, 10))
  plot.new()
  par(mar= c(4.1,1,1,1))
  barplot(ind_R, horiz = T)
  #diff
  par(mar= c(1,4.1,1,5.1), xaxs= "i", yaxs= "i")
  barplot(ind_L)
  par(mar= c(4.1,4.1,1,5.1))
  my_heatmap(mat2, cluster_cols = F, cluster_rows = F, breaks = c(-4, 0, 4))
  plot.new()
  par(mar= c(4.1,1,1,1))
  barplot(ind_R, horiz = T)
  dev.off()
}
####
