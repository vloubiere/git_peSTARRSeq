setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(gridExtra)

if(!exists("temp"))
{
  dat <- readRDS("db/read_counts/all_uniq_counts.rds")
  temp <- dat[grepl("temp_switch", enh_L) | grepl("temp_switch", enh_R)]
}

res <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/library/template_switching_templates.txt", 
             select = c(9,15), col.names = c("name", "template"))
res[, no_switching_input:= sum(temp[enh_L==template & enh_R==template & grepl("input", file), counts]), res]
res[, no_switching_DSCP:= sum(temp[enh_L==template & enh_R==template & grepl("DSCP", file), counts]), res]
res[, switching_input:= sum(temp[xor(enh_L==template, enh_R==template) & grepl("input", file), counts]), res]
res[, switching_DSCP:= sum(temp[xor(enh_L==template, enh_R==template) & grepl("DSCP", file), counts]), res]

.m <- melt(res, measure.vars = patterns(no_switching= "^no_switching", switching= "^switching"))
.m[, variable:= c("input", "DSCP")[variable]]
.m[, switching_percentage:= round(switching/sum(no_switching+switching)*100, 1), (.m)]
.m <- .m[, !"template"]

pdf("pdf/ts_and_sanger_checks/template_switching_stats.pdf")
grid.table(.m)
dev.off()