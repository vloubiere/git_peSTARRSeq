# Merge vl BA libraries ####
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

#--------------------------------------------------------------------------####
# Motifs clustering  SOM ####
lib <- readRDS("Rdata/uniq_library_final.rds")
lib <- lib[!detail=="ecoli"]

# Count hits low cutoff
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds, GRanges(lib$seqnames, IRanges(lib$start, lib$end)), 
                   genome= "dm3", p.cutoff= 5e-4, bg="even", out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds)
rownames(counts) <- lib$ID
counts <- as.data.table(counts, keep.rownames = T)
saveRDS(counts, "Rdata/master_motifs_counts.rds")

# Select motifs enriched in at least one of the lib groups
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

# SOM clustering
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

#--------------------------------------------------------------------------####
# Compute enhancer features ####
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

lib <- STARR[lib, , on= "ID"]

# ChIP-Seq enrichment
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

#--------------------------------------------------------------------------####