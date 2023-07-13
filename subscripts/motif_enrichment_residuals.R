require(vlfunctions)
require(data.table)

dat <- readRDS("db/linear_models/FC_vllib002_lm_predictions.rds")
# dat <- readRDS("db/linear_models/FC_vllib002_actPairs_lm_predictions.rds")
mot <- readRDS("db/motif_counts/twist008_motif_counts.rds")
rdm <- readRDS("db/motif_counts/twist008_motif_counts.rds")

enh <- unique(data.table(enh= dat[, c(L, R)]))
enh <- merge(enh, 
             dat[, .(meanL= mean(residuals)), .(enh= L, indL)],
             by= "enh")
enh <- merge(enh, 
             dat[, .(meanR= mean(residuals)), .(enh= R, indR)],
             by= "enh")
enh <- merge(enh, 
             mot[, .(enh= ID, cisbp__M6212)],
             by= "enh")
enh[, group:= fcase(meanL>quantile(meanL, 0.8) & meanR>quantile(meanR, 0.8), "stronger",
                    meanL<quantile(meanL, 0.2) & meanR<quantile(meanR, 0.2), "weaker",
                    default = "bulk")]
counts <- rbind(mot[enh$enh, on="ID"], rdm)
enr <- vl_motif_cl_enrich(counts_list = split(counts[, !"ID"], c(as.character(enh$group), rep("ctl", nrow(rdm)))),
                          control_cl = "ctl")
par(las= 1)
plot(enr,
     padj_cutoff= 0.5)


#-----------------------------------------------#
# Select motifs of interest and make counts matrix
#-----------------------------------------------#
# Import motif counts
enh <- readRDS("db/motif_counts/twist008_motif_counts_low_stringency_no_collapsing.rds")
rdm <- readRDS("db/motif_counts/random_controls_1000_low_stringency_no_collapsing.rds")
# For each motif cluster enriched in enhancers vs controls, select top enriched
counts <- rbindlist(list(enh= enh, rdm= rdm), idcol = "group")
counts <- melt(counts, 
               id.vars = c("ID", "group"),
               variable.factor = F)
enr <- counts[, {
  fisher.test(table(group=="enh", value>0), 
              alternative = "greater")[c("estimate", "p.value")]
}, variable]
enr[vl_Dmel_motifs_DB_full, cluster:= i.motif_cluster, on= "variable==motif_ID"]
sel <- enr[p.adjust(p.value)<0.05, .SD[which.max(estimate), variable], cluster]$V1
# Counts matrix
motL <- enh[dat$L, sel, on= "ID", with= F]
setnames(motL, gsub("$", "__L", names(motL)))
motR <- enh[dat$R, sel, on= "ID", with= F]
setnames(motR, gsub("$", "__R", names(motR)))
mat <- as.matrix(cbind(motL, motR))
rm(list= c("motL", "motR"))
gc()

