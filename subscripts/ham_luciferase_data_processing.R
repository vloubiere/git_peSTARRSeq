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
dat[, L:= switch(L, "HAM1"= "B", "SUP1"= "A", "SCR2"= "control"), L]
dat[, R:= switch(R, "HAM1"= "B", "SUP1"= "A", "SCR2"= "control"), R]
comb <- CJ(L= c("control", "A", "B"),
           R= c("control", "A", "B"),
           sorted = F)
dat <- merge(comb, dat[, .(L, R, rep, luc, ren)])

# 2- Compute enrichment ----
dat <- dat[ren>2500]
dat <- dat[, .(act= mean(luc/ren)), .(L, R, rep)]
dat <- dat[, .(L, R, act= act/act[L=="control" & R=="control"]), rep]
# dat[, act:= log2(act)]
dat[, indL:= act[R=="control"], .(rep, L)]
dat[, indR:= act[L=="control"], .(rep, R)]
dat[, additive:= indL+indR]
# dat[, additive:= log2(2^indL+2^indR-1)]

# Melt
.m <- melt(dat[(L %in% c("A", "B") & R %in% c("A", "B"))], c("rep", "L", "R"))
.m[, variable:= factor(variable, c("indL", "indR", "additive", "act"))]
setorderv(.m, "variable")
.m <- .m[, .(mean= mean(value), sd= sd(value), all= .(value)), .(L, R, variable)]
.m[, name:= factor(paste0(L, R), c("AB", "BA", "AA", "BB"))]
setorderv(.m, "name")

saveRDS(.m, "Rdata/ham_luciferase_final_table.rds")
