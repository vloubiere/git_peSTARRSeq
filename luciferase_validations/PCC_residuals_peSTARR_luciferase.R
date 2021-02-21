setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

# Import
dat <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")
lib <- readRDS("Rdata/library/lib_features.rds")
dat[, group:= paste0(lib[.BY, group, on= "ID==enh_L"], "~", lib[.BY, group, on= "ID==enh_R"]), .(enh_L, enh_R)]
dat <- dat[group=="control~control" | (!grepl("^control", enh_L) & !grepl("^control", enh_R))]
setorderv(dat, "group")
dat[, Cc= , group]
dat <- dat[!grepl("control", enh_L) & !grepl("control", enh_R), 
           lapply(.SD, mean, na.rm= T), 
           .(Sample_ID, enh_L, enh_R, group), 
           .SDcols= patterns("^luc_")]
dat[, diff_luc:=  log2(luc_norm)-log2(luc_add)]
dat[, Cc:= adjustcolor(rainbow(length(unique(dat$group))), 0.6)[.GRP], group]


# Add peSTARR-Seq
peS <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat[peS, diff_peSTARR:= i.diff, on= c("enh_L", "enh_R")]

pdf("pdf/luciferase_validations/comparison_residuals_luciferase_peSTARR.pdf")
plot(diff_peSTARR~diff_luc, dat, col= dat$Cc, pch= 19, las= 1, 
     xlab= "luciferase o/e (log2)", ylab= "pe-STARR-Seq o/e (log2)")
legend("topleft", legend = unique(dat$group), col= unique(dat$Cc), bty= "n", pch= 19)
dev.off()

# correlation / group
pdf("pdf/luciferase_validations/PCC_residuals_luciferase_peSTARR.pdf")
dat[, {
  plot(diff_peSTARR~diff_luc, .SD, col= Cc, pch= 19, las= 1, 
       xlab= "luciferase o/e (log2)", ylab= "pe-STARR-Seq o/e (log2)", main= group)
  if(.N>2)
  {
    .lm <- lm(diff_peSTARR~diff_luc, .SD)
    abline(.lm, lty= 2)
    leg <- c(paste0("R2= ", round(summary(.lm)$r.squared, 2)), 
             paste0("PCC= ", round(cor.test(.SD$diff_luc, .SD$diff_peSTARR)$estimate, 2)))
    legend("topleft", legend = leg, col= Cc, bty= "n", pch= 19)
  }
  print(group)
}, .(group, Cc)]
dev.off()





