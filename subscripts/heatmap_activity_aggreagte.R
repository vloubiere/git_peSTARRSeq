setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

dat <- data.table(file= list.files("db/final_tables_exp_model/counts_norm/", full.names = T))
dat[, cdition:= tstrsplit(basename(file), "_", keep= 1)]
dat <- dat[, fread(file), (dat)]
lib8 <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
dat[cdition %in% c("vllib002", "vllib006", "vllib013", "vllib014"), ]


dat[, {
  vl_heatmap()
}, cdition]