setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
source("git_peSTARRSeq/functions/plot_transgene.R")

#-----------------------------------------------#
# Count matched
#-----------------------------------------------#
dat <- data.table(file= list.files("db/merged_counts/", c(""), full.names = T))
dat <- SJ(file= pairs_counts, 
          cdition= cdition, 
          rep= DESeq2_pseudo_rep)
setkeyv(dat, c("cdition", "rep"))
dat <- dat[, fread(file), (dat)]
counts <- dcast(dat, 
                L+R~cdition+rep, 
                value.var = "umi_counts", 
                fun.aggregate = sum, 
                sep= "_rep")
# Remove homotypic pairs and cutoff low counts
counts <- counts[L!=R & rowSums(counts[, !c("L", "R")])>20]

####### DESeq2 ########
if(!file.exists(FC_file_DESeq) & 
   all(c("input_rep1", "input_rep2", "screen_rep1", "screen_rep2") %in% names(counts)))
{
  # Format sampleTable
  sampleTable <- SJ(name= setdiff(names(counts), c("L","R")))
  sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_rep")], 
                            row.names = "name")
  # Format DF
  DF <- counts[, name:= paste0(L, "__", R)][, !c("L", "R")]
  DF <- data.frame(DF, 
                   row.names = "name")
  # DESeq
  dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                        colData= sampleTable,
                                        design= ~rep+cdition)
  sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*__control.*", rownames(DF)),]))
  dds <- DESeq2::DESeq(dds)
  
  # Differential expression
  FC <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
  FC <- as.data.table(FC, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
  
  # Compute expected
  FC[, median_L:= .SD[grepl("^control", R), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], L]
  FC[, median_R:= .SD[grepl("^control", L), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], R]
  FC[, additive:= log2(2^median_L+2^median_R)]
  FC[, multiplicative:= median_L+median_R]
  
  # SAVE
  saveRDS(dds, dds_file)
  fwrite(FC, FC_file_DESeq, sep= "\t", na = NA)
  print(paste0(FC_file_DESeq, "  -->> DONE"))
}

####### Ratios ########
if(!file.exists(FC_file_ratio))
{
  # Format counts
  norm <- copy(counts)
  inputs <- grep("^input", names(norm), value = T)
  norm[, input:= rowSums(.SD), .SDcols= inputs]
  norm <- norm[, !(inputs), with= F]
  # Pseudocount normalized counts
  cols <- grep("^screen|^input$", names(norm), value= T)
  norm[, (cols):= lapply(.SD, function(x) (x+0.5)/sum(x+0.5)*1e6), .SDcols= cols]
  # FoldChange
  norm[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
  # Check if individual enhancer is active
  control_pairs_log2FC <- norm[grepl("control", L) & grepl("control", R), log2FoldChange]
  norm[, act_wilcox_L:= {
    .c <- log2FoldChange[grepl("control", R)]
    if(length(.c)>5)
      wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
        as.numeric(NA)
  }, L]
  norm[, act_wilcox_R:= {
    .c <- log2FoldChange[grepl("control", L)]
    if(length(.c)>5)
      wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
        as.numeric(NA)
  }, R]
  # Subtract basal activity (center controls on 0)
  norm[, log2FoldChange:= log2FoldChange-median(control_pairs_log2FC)]
  # Compute expected
  norm[, median_L:= .SD[grepl("^control", R), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], L]
  norm[, median_R:= .SD[grepl("^control", L), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], R]
  norm[, additive:= log2(2^median_L+2^median_R)]
  norm[, multiplicative:= median_L+median_R]
  # SAVE
  norm <- na.omit(norm[, .(L, R, log2FoldChange, median_L, median_R, act_wilcox_L, act_wilcox_R, additive, multiplicative)])
  fwrite(norm, FC_file_ratio, na = NA, sep= "\t")
  print(paste0(FC_file_ratio, "  -->> DONE"))

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
# dat <- readRDS("Rdata/final_results_table_spacer_size.rds")
# dat <- dat[act_wilcox_L<0.001 & median_L>log2(1.5) 
#            & act_wilcox_R<0.001 
#            & median_R>log2(1.5) 
#            & grepl("^dev|^hk", L) 
#            & grepl("^dev|^hk", R)]
dat <- dat[grepl("^dev|^hk", L) 
           & grepl("^dev|^hk", R)]
dat[, diff:= log2FoldChange-additive]
set.seed(1)
dat <- dat[sample(nrow(dat), nrow(dat))]
dat[, class:= factor(
  fcase(grepl("dev", L) & grepl("dev", R), "dev/dev",
        grepl("dev", L) & grepl("hk", R), "dev/hk",
        grepl("hk", L) & grepl("dev", R), "hk/dev",
        grepl("hk", L) & grepl("hk", R), "hk/hk"),
  levels= c("hk/hk",
            "hk/dev",
            "dev/hk",
            "dev/dev"))]
dat[class=="dev/dev", Cc:= "#74C27A"]
dat[class=="hk/dev", Cc:= "royalblue2"]
dat[class=="dev/hk", Cc:= "cyan"]
dat[class=="hk/hk", Cc:= "tomato"]
dat[, CP:= factor(CP, c("DSCP", "RpS12"))]
dat[, spacer:= factor(spacer, c("no_intron", "short_intron", "long_intron"))]

pdf("pdf/draft/Figure_4A.pdf",
    height= 3,
    width = 3.75)
layout(matrix(1:8, 
              nrow= 2, 
              byrow = T),
       widths = c(1.3,1,1,1),
       heights = c(1, 1.375))
for(var in c("log2FoldChange", "diff"))
{
  ylab <- NA
  if(var=="log2FoldChange")
  {
    ylim <- c(-4, 12)
    ylab <- "Activity (log2)"
  }else if(var=="diff")
  {
    ylim <- c(-5, 5.5)
    ylab <- "Obs./Exp. Add. (log2)"
  }
  dat[,
      {
        margins <- c(1,2,1.5,0.1)
        xaxt <- "n"
        if(class=="hk/hk")
          margins[2] <- 4 else
            ylab <- NA
        if(CP=="RpS12")
        {
          margins[c(1, 3)] <- c(6, 0.05)
          xaxt <- "o"
        }
        par(mar= margins,
            lwd= 0.5,
            mgp= c(1.5, 0.5, 0),
            tcl= -0.2,
            las= 1)
        vl_boxplot(split(get(var), spacer),
                   ylim= ylim,
                   boxcol= Cc,
                   xaxt= xaxt,
                   compute_pval= list(c(1,2), c(2,3)),
                   tilt.names= T,
                   notch= T,
                   ylab= ylab)
        transgene(grconvertX(1, "lines", "user"),
                  -3.6, 
                  CP= ifelse(CP=="DSCP", "dev", "hk"), 
                  rotate = T,
                  cex= 0.5)
        abline(h= 0,
               lty= 2)
        if(CP=="DSCP")
          title(class)
        print("")
      }, keyby= .(CP, class, Cc)]
}
dev.off()

