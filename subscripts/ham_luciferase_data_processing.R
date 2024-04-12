setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(readxl)
require(plater)

# 1- Import and format data ----
dat <- data.table(scheme= list.files("db/luciferase/HAM_pilot/", "scheme", full.names = T))
dat <- dat[, read_plate(scheme), scheme]
dat[, date:= tstrsplit(basename(scheme), "_", keep= 1)]
dat[, c("L", "R", "prep", "tech_replicate"):= tstrsplit(rownames, "_")]
dat[, rep:= .GRP, .(date, prep)]
var <- data.table(file= list.files("db/luciferase/HAM_pilot/", "lum.*.csv", full.names = T))
var <- var[, read_plate(file), file]
var[, date:= tstrsplit(basename(file), "_", keep= 1)]
var[, lum:= gsub(".csv", "", unlist(tstrsplit(basename(file), "_", keep= 3)))]
var <- dcast(var, date+Wells~lum, value.var = "values")
dat[var, c("luc", "ren"):= .(i.lum1, i.lum2), on= c("date", "Wells")]
dat[, L:= switch(L, "HAM1"= "A", "SUP1"= "B", "SCR2"= "control"), L]
dat[, R:= switch(R, "HAM1"= "A", "SUP1"= "B", "SCR2"= "control"), R]
comb <- CJ(L= c("control", "A", "B"),
           R= c("control", "A", "B"),
           sorted = F)
dat <- merge(comb, dat[, .(L, R, rep, luc, ren)])

# 2- Compute enrichment ----
dat <- dat[ren>2500]
dat <- dat[, .(act= mean(log2(luc/ren))), .(L, R, rep)]
dat <- dat[, .(L, R, log2FoldChange= act-act[L=="control" & R=="control"]), rep]
dat[, c("indL", "sdL"):= .(mean(log2FoldChange[R=="control"]), sd(log2FoldChange[R=="control"])), L]
dat[, c("indR", "sdR"):= .(mean(log2FoldChange[L=="control"]), sd(log2FoldChange[L=="control"])), R]
dat <- unique(dat[, .(log2FoldChange= mean(log2FoldChange),
                      sd= sd(log2FoldChange)), .(L, R, indL, indR, sdL, sdR)])
dat[, additive:= log2(2^indL+2^indR-1)]
dat[, multiplicative:= indL+indR]
dat[, best:= ifelse(abs(log2FoldChange-additive)<abs(log2FoldChange-multiplicative), "Additive", "Multiplicative")]

saveRDS(dat, "Rdata/ham_luciferase_final_table.rds")