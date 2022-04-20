setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)


# Metadata
meta <- fread("Rdata/metadata_processed.txt")
unique(meta[, .(vllib, spacer, CP, library)])

# Import counts
dat <- meta[vllib %in% c("vllib015", "vllib016") & (DESeq2),
            .(file= pairs_counts,
              rep= DESeq2_pseudo_rep,
              cdition,
              CP)]
dat <- unique(dat)
# Collapse input counts and screen reps
dat <- dat[, fread(file), (dat)]
dat <- dat[, .(umi_counts= sum(umi_counts)), .(rep, cdition, CP, L, R)]

# hk rougly 3 times less sequenced than dev
dat[, sum(umi_counts), .(cdition, CP, rep)]

# Down sample dev 3 times, 10 times
dat[, down_3:= {
  if(cdition=="input")
    umi_counts else
    {
      set.seed(1)
      res <- .SD[rep(seq(nrow(.SD)), umi_counts)][sample(.N, sum(umi_counts)/3)][, .(sub= .N), .(L, R)]
      res[.SD, sub, on= c("L", "R"), nomatch= NA]
    }
}, .(CP, cdition, rep), .SDcols= c("L", "R")]
dat[, down_10:= {
  if(cdition=="input")
    umi_counts else
    {
      set.seed(1)
      res <- .SD[rep(seq(nrow(.SD)), umi_counts)][sample(.N, sum(umi_counts)/10)][, .(sub= .N), .(L, R)]
      res[.SD, sub, on= c("L", "R"), nomatch= NA]
    }
}, .(CP, cdition, rep), .SDcols= c("L", "R")]
dat <- melt(dat, id.vars = c("cdition", "CP", "rep", "L", "R"))
dat <- na.omit(dat)

# Compute activity and expected score
res <- dat[CP=="DSCP" & variable=="umi_counts", {
  res <- dcast(.SD,
               L+R~cdition+rep,
               value.var= "value",
               fill= 0,
               sep = "__")
  
  sampleTable <- data.table(rn= setdiff(names(res), c("L", "R")))
  sampleTable[, c("cdition", "rep"):= tstrsplit(rn, "__")]
  sampleTable <- data.frame(sampleTable, 
                            row.names = "rn")
  
  # DF
  res[, rn:= paste0(L, "__", R)]
  DF <- data.frame(res[, !c("L", "R")], 
                   row.names = "rn")
  DF <- DF[rowSums(DF)>20,]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                        colData= sampleTable,
                                        design= ~rep+cdition)
  sizeFactors(dds)= DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*__control.*", rownames(DF)),]))
  
  res <- DESeq2::DESeq(dds)
  
  # Differential expression
  FC <- as.data.frame(DESeq2::results(res, contrast= c("cdition", "screen", "input")), keep.rownames= "rn")
  FC <- as.data.table(FC, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]

  # Compute expected
  median_L <- FC[grepl("control", R), .(check= .N>5, median_L= median(log2FoldChange, na.rm = T)), L][(check)]
  FC[median_L, median_L:= i.median_L, on= "L"]
  median_R <- FC[grepl("control", L) , .(check= .N>5, median_R= median(log2FoldChange, na.rm = T)), R][(check)]
  FC[median_R, median_R:= i.median_R, on= "R"]
  FC[, additive:= log2(2^median_L+2^median_R)]
  FC[, multiplicative:= median_L+median_R]
  
  # RETURN
  print(paste0(CP, ".../...", variable))
  na.omit(FC[, .(L, R, log2FoldChange, padj, median_L, median_R, additive, multiplicative)])
}, .(CP, variable)]

pl <- dcast(res, CP+L+R~variable, value.var = list("log2FoldChange", "padj", "median_L", "median_R", "additive", "multiplicative"))
par(mfrow= c(2,2))
smoothScatter(pl[CP=="DSCP", .(log2FoldChange_umi_counts, log2FoldChange_down_3)])
smoothScatter(pl[CP=="DSCP", .(log2FoldChange_umi_counts, log2FoldChange_down_10)])
points(pl[CP=="DSCP" & grepl("control", L) & grepl("control", R) & L==R, .(log2FoldChange_umi_counts, log2FoldChange_down_10)], pch= 19, cex= 0.1, col= "red")
points(pl[CP=="DSCP" & grepl("control", L) & grepl("control", R) & L!=R, .(log2FoldChange_umi_counts, log2FoldChange_down_10)], pch= 19, cex= 0.1, col= "green")
boxplot(pl[CP=="DSCP", .(log2(2^median_L_umi_counts+2^median_R_umi_counts),
                         log2(2^median_L_down_3+2^median_R_down_3),
                         log2(2^median_L_down_10+2^median_R_down_10))])
boxplot(pl[CP=="RpS12", .(log2(2^median_L_umi_counts+2^median_R_umi_counts),
                          log2(2^median_L_down_3+2^median_R_down_3),
                          log2(2^median_L_down_10+2^median_R_down_10))])
