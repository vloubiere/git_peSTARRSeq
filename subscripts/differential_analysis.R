# Differential analysis ####
# 1- Make master counts table
# Template switching 
raw <- data.table(file= list.files("db/count/", "_umi_counts.txt", full.names = T))
raw[, lib:= tstrsplit(basename(file), "_", keep=1), file]
raw[, cdition:= ifelse(grepl("input", file), "input", "DSCP")]
raw[, rep:= tstrsplit(basename(file), "DSCP_|input_|_umi", keep = 2)]
raw <- raw[, fread(file), raw]
raw <- raw[(exist)] # Only retain correct pairs
raw <- dcast(raw, L+R~lib+cdition+rep, value.var = "count")

# Clean counts table
counts <- data.table(raw[, .(L, R)],
                     libvl002_DSCP_rep1= raw$vllib002_DSCP_rep1+raw$vllib002_DSCP_rep2,
                     libvl002_DSCP_rep2= raw$vllib002_DSCP_rep3+raw$vllib002_DSCP_rep4,
                     libvl002_input_rep1= raw$vllib002_input_rep1+raw$vllib002_input_rep3+raw$vllib002_input_rep4,
                     libvl002_input_rep2= raw$vllib002_input_rep2+raw$vllib002_input_rep5+raw$vllib002_input_rep6,
                     libvl013_input_rep1= raw$vllib013_input_rep1,
                     libvl013_input_rep2= raw$vllib013_input_rep2,
                     libvl013_DSCP_rep1= raw$vllib013_DSCP_rep1,
                     libvl013_DSCP_rep2= raw$vllib013_DSCP_rep2,
                     libvl014_input_rep1= raw$vllib014_input_rep1,
                     libvl014_input_rep2= raw$vllib014_input_rep2,
                     libvl014_DSCP_rep1= raw$vllib014_DSCP_rep1,
                     libvl014_DSCP_rep2= raw$vllib014_DSCP_rep2)
fwrite(counts, "Rdata/master_counts_table.txt", col.names = T, row.names = F, sep= "\t", quote= F)

# 2- run DESeq
DE <- data.table(lib= colnames(fread("Rdata/master_counts_table.txt", nrows = 0))[-c(1,2)])
DE[, c("lib", "cdition", "rep"):= tstrsplit(lib, "_", keep= 1:3)]
DE[, dds_file:= paste0("db/DE_analysis/", lib, "_dds.rds")]
DE[, diff:= gsub("_dds.rds$", "_DE.txt", dds_file)]

master_counts <- fread("Rdata/master_counts_table.txt")
master_counts[, ID:= paste0(L, "_vs_", R)]
master_counts$L <- NULL
master_counts$R <- NULL
setcolorder(master_counts, "ID")

DE[, {
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
  # Differential expression
  .c <- as.data.table(as.data.frame(results(res, contrast= c("cdition", "DSCP", "input"))), keep.rownames= T)
  .c[, c("L", "R"):= tstrsplit(rn, "_vs_")]
  .c <- .c[, !"rn"]
  setcolorder(.c, c("L", "R"))
  fwrite(.c, diff, col.names = T, row.names = F, sep= "\t", quote= F)
  print(paste0(diff, " DONE!"))
}, .(lib, dds_file, diff)]

