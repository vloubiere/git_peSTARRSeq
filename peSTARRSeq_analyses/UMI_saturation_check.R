setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

.c <- data.table(file= c("db/read_counts/libvl002_SCR1_DSCP_rep1.all.rds", "db/read_counts/libvl002_SCR1_DSCP_rep1.uniq.UMI.rds"))
.c <- .c[, readRDS(file), .c]

agg <- .c[, .N, .(file, enh_L, enh_R)]
agg <- dcast(agg, enh_L+enh_R~file, value.var = "N")

agg <- unique(agg[, `db/read_counts/libvl002_SCR1_DSCP_rep1.all.rds`:`db/read_counts/libvl002_SCR1_DSCP_rep1.uniq.UMI.rds`])
plot(agg)
