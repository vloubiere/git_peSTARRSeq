setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Metadata ----
dat <- readRDS("Rdata/metadata_processed.rds")
dat <- dat[(paper)]
dat[, name:= paste0(switch(vllib,
                           "vllib002"= "WT_oligo_pool_",
                           "vllib029"= "Mutated_oligo_pool_",
                           "vllib015"= "Focused_oligo_pool_",
                           "vllib016"= "Focused_oligo_pool_"),
                    switch(CP,
                           "DSCP"= "dCP_",
                           "RpS12"= "hkCP_"),
                    cdition, "_",
                    rep), .(CP, vllib)]
# Fq files ----
dat[, fq_geo_1:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/GEO/",
                        name,
                        "_1.fq.gz"), .(CP, vllib)]
dat[, fq_geo_2:= gsub("_1.fq.gz$", "_2.fq.gz", fq_geo_1)]
dat[, {
  if(!file.exists(fq_geo_1))
  {
    cmd <- paste(c("cat", unique(fq1), ">", fq_geo_1), collapse= " ")
    system(cmd)
  }
  print("done")
  if(!file.exists(fq_geo_2))
  {
    cmd <- paste(c("cat", unique(fq2), ">", fq_geo_2), collapse= " ")
    system(cmd)
  }
  print("done")
  .SD
}, .(fq_geo_1, fq_geo_2)]

# Umi counts ----
dat[, umi_counts_geo:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/GEO/",
                              name,
                              "_umi_counts.txt")]
dat[, {
  if(!file.exists(umi_counts_geo))
    file.copy(umi_counts, umi_counts_geo)
}, .(umi_counts, umi_counts_geo)]

# FC tables ----
dat[, FC_table:= list.files("db/FC_tables/", paste0(vllib, ".*.rds$"), full.names = T), vllib]
dat[, FC_table_geo:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/GEO/",
                            sub("(_input|_screen).*", "", name),
                            "_FC_table.txt"), .(name, FC_table)]
dat[, {
  if(!file.exists(FC_table_geo))
  {
    .c <- readRDS(FC_table)
    .c$groupL <- .c$groupR <- .c$actL <- .c$actR <- NULL
    setcolorder(.c,
                c("L", "R", "ctlL", "ctlR", "indL", "padjL", "indR", "padjR"))
    fwrite(.c,
           FC_table_geo,
           col.names = T,
           row.names = F,
           sep= "\t",
           na=NA,
           quote=F)
  }
}, .(FC_table, FC_table_geo)]

# GEO submission metadata
samples <- dat[, .(name,
                   fq1= basename(fq_geo_1),
                   fq2= basename(fq_geo_2),
                   umi_counts= basename(umi_counts_geo),
                   FC_table= basename(FC_table_geo))]
samples <- unique(samples)
fwrite(samples,
       "db/GEO/samples.txt",
       col.names = T,
       row.names = F,
       sep= "\t")

# GEO checksums
md5 <- melt(dat,
            id.vars = "name",
            measure.vars = patterns("geo"))
md5[, name:= sub("(_input|_screen).*", "", name)]
md5 <- unique(md5)
md5[, checksum:= tools::md5sum(value), value]
md5[, value:= basename(value)]
fwrite(md5,
       "db/GEO/checksums.txt",
       col.names = T,
       row.names = F,
       sep= "\t")
