setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(data.table)

#---------------------------------------------#
# Metadata
#---------------------------------------------#
meta <- fread("Rdata/metadata_processed.txt")
unique(meta[, .(vllib, spacer, CP, library)])

#---------------------------------------------#
# Longer spacer
#---------------------------------------------#
# Import counts
dat <- meta[vllib %in% c("vllib015", "vllib018", "vllib023",
                         "vllib016", "vllib020", "vllib024") & (DESeq2),
            .(file= pairs_counts,
              rep= DESeq2_pseudo_rep,
              cdition,
              spacer,
              CP)]
dat[, spacer:= switch(spacer,
                      "SCR1"= "no_intron",
                      "shortened-intron4"= "short_intron",
                      "intron4"= "long_intron"), spacer]
dat[cdition=="screen", cdition:= paste0(cdition, "_rep", rep)]
dat <- unique(dat)
dat <- dat[, fread(file), (dat)]
# Collapse input counts and screen reps
dat <- dat[, .(umi_counts= sum(umi_counts)), .(cdition, CP, spacer, L, R)]
# Cast counts
res <- dcast(dat,
             L+R+CP~spacer+cdition,
             value.var= "umi_counts",
             fill= 0,
             sep = "__")
cols <- setdiff(names(res), c("L", "R", "CP"))
res[, check:= rowSums(.SD), .SDcols= cols]
res <- res[check>20, !"check"]
# Normalize
res[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
res[, (cols):= lapply(.SD, function(x) (x+1)/sum(x)*1e6), CP, .SDcols= cols]
# Format and compute FC
res <- melt(res, id.vars = c("L", "R", "CP"))
res[, c("spacer", "cdition"):= tstrsplit(variable, "__")]
res <- dcast(res, L+R+CP+spacer~cdition, value.var = "value")
res[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
# Check if individual enhancer is active
control_pairs_log2FC <- res[grepl("control", L) & grepl("control", R), log2FoldChange, .(CP, spacer)]
res[, act_wilcox_L:= {
  .c <- log2FoldChange[grepl("control", R)]
  if(length(.c)>5)
    wilcox.test(.c,
                control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "spacer")],
                alternative = "greater")$p.value else
                  as.numeric(NA)
}, .(CP, spacer, L)]
res[, act_wilcox_R:= {
  .c <- log2FoldChange[grepl("control", L)]
  if(length(.c)>5)
    wilcox.test(.c,
                control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "spacer")],
                alternative = "greater")$p.value else
                  as.numeric(NA)
}, .(CP, spacer, R)]
# Subtract basal activity (center controls on 0)
res[, log2FoldChange:= log2FoldChange-median(control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "spacer")]), .(CP, spacer)]
# Compute expected
median_L <- res[grepl("control", R) , .(check= .N>5, median_L= median(log2FoldChange)), .(L, CP, spacer)][(check)]
res[median_L, median_L:= i.median_L, on= c("CP", "spacer", "L")]
median_R <- res[grepl("control", L) , .(check= .N>5, median_R= median(log2FoldChange)), .(R, CP, spacer)][(check)]
res[median_R, median_R:= i.median_R, on= c("CP", "spacer", "R")]
# Expected
res[, additive:= log2(2^median_L+2^median_R)]
res[, multiplicative:= median_L+median_R]
# SAVE
res <- na.omit(res[, .(CP, spacer, L, R, log2FoldChange, median_L, median_R, act_wilcox_L, act_wilcox_R, additive, multiplicative)])
saveRDS(res,
        "Rdata/final_results_table_spacer_size.rds")

#---------------------------------------------#
# Core promoters
#---------------------------------------------#
# Import counts
dat <- meta[vllib %in% c("vllib015", "vllib027", "vllib028",
                         "vllib016", "vllib026", "vllib025") & (DESeq2), 
            .(file= pairs_counts, 
              rep= DESeq2_pseudo_rep,
              cdition,
              CP)]
dat[, type:= fcase(CP %in% c("DSCP", "devLow", "devHigh"), "dev",
                   CP %in% c( "RpS12", "hkLow", "hkHigh"), "hk"), CP]
dat[, CP:= switch(CP,
                 "DSCP"= "ref",
                 "devLow"= "low",
                 "devHigh"= "high",
                 "RpS12"= "ref",
                 "hkLow"= "low",
                 "hkHigh"= "high"), CP]
setnames(dat, 
         c("CP", "type"), 
         c("basalAct", "CP"))
dat[cdition=="screen", cdition:= paste0(cdition, "_rep", rep)]
dat <- unique(dat)
dat <- dat[, fread(file), (dat)]
# Collapse input counts and screen reps
dat <- dat[, .(umi_counts= sum(umi_counts)), .(cdition, CP, basalAct, L, R)]
# Cast counts
res <- dcast(dat,
             L+R+CP~basalAct+cdition,
             value.var= "umi_counts",
             fill= 0, 
             sep = "__")
cols <- setdiff(names(res), c("L", "R", "CP"))
res[, check:= rowSums(.SD), .SDcols= cols]
res <- res[check>20, !"check"]
# Normalize
res[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
res[, (cols):= lapply(.SD, function(x) (x+1)/sum(x)*1e6), CP, .SDcols= cols]
# Format and compute FC
res <- melt(res, id.vars = c("L", "R", "CP"))
res[, c("basalAct", "cdition"):= tstrsplit(variable, "__")]
res <- dcast(res, L+R+CP+basalAct~cdition, value.var = "value") 
res[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
# Check if individual enhancer is active
control_pairs_log2FC <- res[grepl("control", L) & grepl("control", R), log2FoldChange, .(CP, basalAct)]
res[, act_wilcox_L:= {
  .c <- log2FoldChange[grepl("control", R)]
  if(length(.c)>5)
    wilcox.test(.c, 
                control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "basalAct")], 
                alternative = "greater")$p.value else
                  as.numeric(NA)
}, .(CP, basalAct, L)]
res[, act_wilcox_R:= {
  .c <- log2FoldChange[grepl("control", L)]
  if(length(.c)>5)
    wilcox.test(.c, 
                control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "basalAct")], 
                alternative = "greater")$p.value else
                  as.numeric(NA)
}, .(CP, basalAct, R)]
# Subtract basal activity (center controls on 0)
res[, log2FoldChange:= log2FoldChange-median(control_pairs_log2FC[.BY, log2FoldChange, on= c("CP", "basalAct")]), .(CP, basalAct)]
# Compute expected
median_L <- res[grepl("control", R) , .(check= .N>5, median_L= median(log2FoldChange)), .(L, CP, basalAct)][(check)]
res[median_L, median_L:= i.median_L, on= c("CP", "basalAct", "L")]
median_R <- res[grepl("control", L) , .(check= .N>5, median_R= median(log2FoldChange)), .(R, CP, basalAct)][(check)]
res[median_R, median_R:= i.median_R, on= c("CP", "basalAct", "R")]
# Expected
res[, additive:= log2(2^median_L+2^median_R)]
res[, multiplicative:= median_L+median_R]
# SAVE
res <- na.omit(res[, .(CP, basalAct, L, R, log2FoldChange, median_L, median_R, act_wilcox_L, act_wilcox_R, additive, multiplicative)])
saveRDS(res, 
        "Rdata/final_results_table_CP_basalAct.rds")
