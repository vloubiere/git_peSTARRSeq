setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)
require(plater)
dir_pdf <- normalizePath("pdf/luciferase", mustWork = F)
dir.create(dir_pdf, showWarnings = F)
constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", key= "name")
pl <- fread("/groups/stark/vloubiere/exp_data/vl_plasmids.txt")
pl <- pl[Experiment=="ham_pilot_luc"]

#----------------------------#
# sanger sequencing
#----------------------------#
p1 <- "db/sanger_sequencing/ham_pilot/100220/"
pl[, file_1:= .(list(list.files(p1, x, full.names= T))), .(x= paste0("pluc", gsub(".* (.*$)", "\\1_", labbook)))]
p2 <- "db/sanger_sequencing/ham_pilot/240220_recloning/"
pl[, file_2:= .(list(list.files(p2, x, full.names = T))), .(x= paste0("_p", formatC(as.numeric(gsub(".* (.*$)", "\\1", labbook)), width = 2, flag = "0")))]
p3 <- "db/sanger_sequencing/ham_pilot/240220_reseq_last_sample/"
pl[, file_3:= .(list(list.files(p3, x, full.names= T))), .(x= paste0("pluc", gsub(".* (.*$)", "\\1_", labbook)))]
pl <- melt(pl, measure.vars = patterns("file"))
pl <- pl[lengths(value)>0, .(value= unlist(value)), setdiff(colnames(pl), "value")]
upstream <- vl_digest(constructs["DSCP_pluc002", sequence], "KpnI")[1]
downstream <- vl_digest(constructs["DSCP_pluc002", sequence], "KpnI")[2]
pl[, seq:= paste0(substr(upstream, 
                         start = nchar(upstream)-200,
                         stop = nchar(upstream)),
                  constructs["illumina_F", sequence],
                  constructs["Flink_+0", sequence],
                  constructs[gsub("(^.*)-(.*)-(.*$)", "\\1", contains), sequence],
                  constructs["R1link+0", sequence],
                  constructs["CGCov_F", sequence],
                  constructs[gsub("(^.*)-(.*)-(.*$)", "\\2", contains), sequence],
                  constructs["CGCov_R", sequence],
                  constructs["Flink_+0", sequence],
                  constructs[gsub("(^.*)-(.*)-(.*$)", "\\3", contains), sequence],
                  constructs["R1link+0", sequence],
                  substr(downstream, 1, 200)), contains]
pl[, rev:= ifelse(grepl("CASeq001", value), F, T)]

pdf(paste0(dir_pdf, "/sanger_sequencing.pdf"), height = 5)
par(mar= c(2,20,5,2))
pl[, 
   {
     vl_sanger_align(refseq = seq, 
                     abfiles = value, 
                     revcomp = rev, 
                     feat_sequences = constructs[c("HAM1", "SCR1", "SUP1"), sequence], 
                     feat_names = c("HAM1", "SCR1", "SUP1"),
                     feat_cols = c("green", "black", "red"))
     mtext(contains)
   }, .(seq, dirname(value))]
dev.off()

#----------------------------#
# Luciferase assays
#----------------------------#
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

pdf(paste0(dir_pdf, "/barplot_luciferase.pdf"), height = 5)
par(las= 1)
bar1 <- barplot(pl[, mean], space = spa1, col= Cc, ylab= yl, ylim= c(0, 14))

pl[, b1:= bar1[,1]]
pl[, points(jitter(rep(b1, length(unlist(mean_all)))), unlist(mean_all), pch= 16, cex= 0.6), b1]
arrows(pl$b1, pl$mean, pl$b1, pl[, mean+sd], length = 0.075, angle = 90)
arrows(pl$b1, pl$mean, pl$b1, pl[, mean-sd], length = 0.075, angle = 90)

vl_plot_pval(x = pl$b1, y = pl$max, pl$pval, offset = 1, pos = 3, stars_only = T)
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
pdf(paste0(dir_pdf, "/barplot_left_right_activity.pdf"), height = 4, width = 5)
sel <- copy(dat)
setkeyv(sel, "name")
sel <- sel[c("DSCP_ZFH1", "p002_empty", "SCR2_alone", "SCR2_SCR2", "HAM1_alone", "HAM1_SCR2", "SCR2_HAM1", "SUP1_alone", "SUP1_SCR2", "SCR2_SUP1")]
Cc <- c("limegreen", "black", "lightgrey", "lightgrey", rep("gold", 3), rep("tomato", 3))
bar <- barplot(sel[, mean(value), name][, V1], ylim= c(0,12), ylab = yl, col= Cc, las=1)
text(bar, -0.5, unique(sel$name), srt= 45, xpd= T, offset = -0.5, pos= 2)
dev.off()