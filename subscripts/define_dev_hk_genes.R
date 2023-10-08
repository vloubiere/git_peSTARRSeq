setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

dat <- fread("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm3/gene_rpkm_report_fb_2014_03.tsv.gz",
             fill= T,
             skip = 5)
colnames(dat) <- gsub("\\#", "", colnames(dat))
dat <- dcast(dat[Parent_library_name == "modENCODE_mRNA-Seq_U"], FBgn+GeneSymbol~RNASource_name, value.var = "RPKM_value")
dat <- melt(dat,
            id.vars = c("FBgn", "GeneSymbol"),
            value.name = "RPKM",
            variable.name = "condition")
# Compute quantile
dat[, quant:= quantile(RPKM, 0.4), condition]
res <- dat[, .(class= ifelse(all(RPKM>quant), "housekeeping", "developmental")), FBgn]

fwrite(res,
       "db/public/hk_dev_FB_r5.57_2014_03.txt",
       sep= "\t",
       quote=F,
       col.names = T,
       row.names = F,
       na= NA)
