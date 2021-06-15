setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(plater)

# import plates
dat <- data.table(scheme= list.files("db/luciferase/HAM_pilot/", "scheme", full.names = T))
dat <- dat[, read_plate(scheme), scheme]
dat[, date:= tstrsplit(basename(scheme), "_", keep= 1)]
dat[, c("enh_L", "enh_R", "prep", "tech_replicate"):= tstrsplit(rownames, "_")]
dat[, replicate:= .GRP, .(date, prep)]
var <- data.table(file= list.files("db/luciferase/HAM_pilot/", "lum.*.csv", full.names = T))
var <- var[, read_plate(file), file]
var[, date:= tstrsplit(basename(file), "_", keep= 1)]
var[, lum:= gsub(".csv", "", unlist(tstrsplit(basename(file), "_", keep= 3)))]
var <- dcast(var, date+Wells~lum, value.var = "values")
dat[var, c("luc", "ren"):= .(i.lum1, i.lum2), on= c("date", "Wells")]
dat[, name:= paste0(enh_L, "_", enh_R), .(enh_L, enh_R)]
dat <- dat[ren>5000 & name!="Ubi_GFP", .(enh_L, enh_R, name, replicate, value= luc/ren)]
dat <- dat[, .(value= mean(value)), .(enh_L, enh_R, name, replicate)]

# Add additive scores
dat[, mean_L:= .(mean(.SD[enh_R=="SCR2", value])), .(replicate, enh_L)]
dat[, mean_R:= .(mean(.SD[enh_L=="SCR2", value])), .(replicate, enh_R)]
dat[, add:= mean_L+mean_R]
dat[!(enh_L %in% c("HAM1", "SUP1") | enh_R %in% c("HAM1", "SUP1")), c("add", "mean_L", "mean_R") := .(NA, NA, NA)]
setkeyv(dat, "name")

# PLOT
ord <- c("DSCP_ZFH1", "p002_empty", "SCR2_SCR2", "HAM1_HAM1", "HAM1_SUP1", "SUP1_HAM1", "SUP1_SUP1")
pl <- dat[ord, .(mean= mean(value),
                 mean_all= list(value),
                 sd= sd(value),
                 max= max(value),
                 add= mean(add),
                 add_all= list(add),
                 add_sd= sd(add),
                 left= mean(mean_L),
                 pval= ifelse(length(which(!is.na(add)))>1, t.test(value, add, paired = T)$p.value, as.numeric(NA)),
                 FC= round(mean(value)/mean(add), 1)), name]
yl <- "Activity (luciferase/renilla)"
Cc <- c(rep("grey", 3), "mediumturquoise", "gold", "limegreen", "tomato")
spa1 <- c(0.5,0.5,0.5,2,2,2,2)
spa2 <- c(0.5,0.5,0.5,1,2,2,2)

pdf("pdf/barplot_luciferase.pdf", height = 5)
par(las= 1)
bar1 <- barplot(pl[, mean], space = spa1, col= Cc, ylab= yl, ylim= c(0, 14))

pl[, b1:= bar1[,1]]
pl[, points(jitter(rep(b1, length(unlist(mean_all)))), unlist(mean_all), pch= 16, cex= 0.6), b1]
arrows(pl$b1, pl$mean, pl$b1, pl[, mean+sd], length = 0.075, angle = 90)
arrows(pl$b1, pl$mean, pl$b1, pl[, mean-sd], length = 0.075, angle = 90)

text(pl$b1, pl$max, my_pval_format(pl$pval)$p.value, offset = 1, pos = 3)
text(pl$b1, pl$max, ifelse(is.na(pl$FC), NA, paste0("x", pl$FC)), offset = 0.7, pos = 3, cex= 0.8)

bar2 <- barplot(pl[, add], space = spa2, col= Cc, ylim= c(0, 14), add= T)
pl[, b2:= bar2[,1]]
pl[, points(jitter(rep(b2, length(unlist(add_all)))), unlist(add_all), pch= 16, cex= 0.6), b2]
arrows(pl$b2, pl$add, pl$b2, pl[, add+add_sd], length = 0.075, angle = 90)
arrows(pl$b2, pl$add, pl$b2, pl[, add-add_sd], length = 0.075, angle = 90)

barplot(pl[, left], space = spa2, col= "black", ylim= c(0, 14), add= T, density = 10, angle = 45)

at <- rowMeans(cbind(bar1, bar2))
text(at, -0.2, pl$name, srt= 45, pos= 2, xpd= T, offset = -0.25)
dev.off()

# Left and righ activities plotted next each other
pdf("pdf/barplot_left_right_activity.pdf", height = 4, width = 5)
sel <- copy(dat)
setkeyv(sel, "name")
sel <- sel[c("DSCP_ZFH1", "p002_empty", "SCR2_alone", "SCR2_SCR2", "HAM1_alone", "HAM1_SCR2", "SCR2_HAM1", "SUP1_alone", "SUP1_SCR2", "SCR2_SUP1")]
Cc <- c("limegreen", "black", "lightgrey", "lightgrey", rep("gold", 3), rep("tomato", 3))
bar <- barplot(sel[, mean(value), name][, V1], ylim= c(0,12), ylab = yl, col= Cc, las=1)
text(bar, -0.5, unique(sel$name), srt= 45, xpd= T, offset = -0.5, pos= 2)
dev.off()