setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

# Import
dat <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")
lib <- readRDS("Rdata/library/lib_features.rds")
dat[, group:= paste0(lib[.BY, group, on= "ID==enh_L"], "~", lib[.BY, group, on= "ID==enh_R"]), .(enh_L, enh_R)]
dat <- dat[group=="control~control" | (!grepl("^control", enh_L) & !grepl("^control", enh_R))]
setorderv(dat, "group")
dat[, x:= .GRP*2-1, group]
dat[group=="control~control", x:= 1]
dat[, xjit:= jitter(x, 1, 0.15), .(Sample_ID, group, x)]


pdf("pdf/luciferase_validations/connected_scatterplot_luc_validations.pdf", 15, 5)
par(las= 2, mar= c(5,4,2,1))

plot(NA, xlim = c(1, 15), ylim= range(c(dat$luc_norm, dat$luc_add), na.rm= T), ylab= "Normalized luciferase activity", las= 1, xaxt= "n", xlab= "")
mtext(unique(dat$group), side= 1, line= 3, at= unique(dat$x)+c(0, rep(-0.5, length(unique(dat$group))-1)), las= 1)
axis(1, at= unique(dat$x), labels= rep("obs.", length(unique(dat$x))))
axis(1, at= unique(dat$x)[-1]-1, labels= rep("add.", length(unique(dat$x))-1))

dat[, {
  points(xjit, mean(luc_norm, na.rm = T))
  if(group!="control~control")
  {
    points(xjit-1, mean(luc_add, na.rm = T))
    segments(xjit, mean(luc_norm, na.rm = T), xjit-1, mean(luc_add, na.rm = T))
  }
}, .(Sample_ID, group, xjit)]

dev.off()

