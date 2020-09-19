setwd("/groups/stark/vloubiere/projects/0003_luciferase_HAM1_SUP1_SCR2/")
require(data.table)

dat <- readRDS("Rdata/luciferase_validation.rds")
collapse <- unique(dat[, .(name, mean_activity, sd_activity, order)])

# Barplot
Cc <- c(rep("grey", 5), rep("yellow", 3), rep("green", 3), rep("darkgoldenrod1", 2), rep("cornflowerblue", 2), rep("tomato", 2), rep("darkolivegreen2", 2))

pdf("pdf/C_luciferase_barplot.pdf", height = 4.5)
par(mar= c(8,5,2,2))
bar <- barplot(collapse$mean_activity, ylim= c(0, 18), col= Cc, ylab = "Normalized luciferase activity", las= 1)
text(bar, -0.5, collapse$name, srt=45, xpd= T, pos= 2, offset= -0.25)

# sd bars
collapse[, xpos:= bar[order]]
arrows(collapse$xpos, collapse$mean_activity, collapse$xpos, collapse$mean_activity+collapse$sd_activity, length = 0.075, angle = 90)
arrows(collapse$xpos, collapse$mean_activity, collapse$xpos, collapse$mean_activity-collapse$sd_activity, length = 0.075, angle = 90)

# exp_L
exp_L <-  unique(dat[, .(name, mean_exp_L= mean(exp_L)), order])
barplot(exp_L$mean_exp_L, add=T, density = 30, angle = 45, col= "black", border= NA, yaxt= "n")

# add points
dat[, xpos:= bar[order]]
set.seed(1)
dat[, xjit:= jitter(bar[order])]
points(dat$xjit, dat$activity, pch= 19, cex= 0.7, col= adjustcolor("black", 0.6))
dat <- dat[!is.na(order)]
setorderv(dat, "order")

# Stat
stat <- copy(dat[order %in% 12:19])
stat <- stat[, cdition := ifelse(grepl("exp.", name), "exp", "obs")]
stat <- dcast(stat, enh_L+enh_R+idx~cdition, value.var= c("activity", "xpos"))
stat[, pval:= cut(t.test(.SD[,activity_obs], .SD[,activity_exp])$p.value, c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), c("****", "***", "**", "*", "N.S")), .(enh_L, enh_R)]
stat[, diff:= round(mean(activity_obs)/mean(activity_exp), 1), .(enh_L, enh_R)]
offset <- 1.5
cols <- c("activity_exp", "activity_obs")
stat[, segments(xpos_exp[1], max(.SD)+offset, xpos_obs[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, segments(xpos_exp[1], max(activity_exp)+0.2, xpos_exp[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, segments(xpos_obs[1], max(activity_obs)+0.2, xpos_obs[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, text(mean(c(xpos_exp[1], xpos_obs[1])), max(.SD)+offset, paste0("x", diff[1]), xpd= T, pos= 3), .(enh_L, enh_R), .SDcols= cols]
stat[, text(mean(c(xpos_exp[1], xpos_obs[1])), max(.SD)+offset+0.4, pval[1], xpd= T, pos= 3, cex= 1.2), .(enh_L, enh_R), .SDcols= cols]
leg <- c("* pval<5e-2", "** pval<1e-2", "*** pval<1e-3", "**** pval<1e-5")
legend("topleft", paste0("            ", leg), bty= "n")
dev.off()

# Simplified barplot

pdf("pdf/C_luciferase_barplot_simplified.pdf", width = 5.5, height = 4.5)
par(mar= c(8,5,2,2))

simp <- collapse[c(3,5,12,13,14,15,16,17,18,19)]
dat_simp <- dat[simp, , on= "name"]
simp[c(1, 2), name:= c("empty plasmid", "control x 2")]

Cc1 <- c("grey", "grey", "darkgoldenrod2", "gold1", "midnightblue", "cornflowerblue", "seagreen", "limegreen", "tomato", "lightcoral")
Cc2 <- c(NA, NA, "yellow", NA, "mediumturquoise", NA, "greenyellow", NA, "plum1")

bar <- barplot(simp$mean_activity, ylim= c(0, 18), col= Cc1, ylab = "Normalized luciferase activity", las= 1, space = rep(c(1, 0.1), 5))
text(bar, -0.5, simp$name, srt=45, xpd= T, pos= 2, offset= -0.25)

# sd bars
simp[, xpos:= bar]
arrows(simp$xpos, simp$mean_activity, simp$xpos, simp$mean_activity+simp$sd_activity, length = 0.075, angle = 90)
arrows(simp$xpos, simp$mean_activity, simp$xpos, simp$mean_activity-simp$sd_activity, length = 0.075, angle = 90)

# exp_L
exp_L <-  c(0, 0, collapse[name=="hamlet x control", mean_activity], 0,
            collapse[name=="hamlet x control", mean_activity], 0,
            collapse[name=="sup-add x control", mean_activity], 0,
            collapse[name=="sup-add x control", mean_activity])
barplot(exp_L, add=T, col= Cc2, border= NA, yaxt= "n", space = rep(c(1, 0.1), 5))
barplot(exp_L, add=T, density = 10, angle = 45, col= "black", yaxt= "n", space = rep(c(1, 0.1), 5))

# add points
dat_simp[, order:= .GRP, name]
dat_simp[, xpos:= bar[order], name]
set.seed(1)
dat_simp[, xjit:= jitter(xpos)]
points(dat_simp$xjit, dat_simp$activity, pch= 19, cex= 0.7, col= adjustcolor("black", 0.6))
dat_simp <- dat_simp[!is.na(order)]
setorderv(dat_simp, "order")

# Stat
stat <- copy(dat_simp[order %in% 3:10])
stat <- stat[, cdition := ifelse(grepl("exp.", name), "exp", "obs")]
stat <- dcast(stat, enh_L+enh_R+idx~cdition, value.var= c("activity", "xpos"))
stat[, pval:= cut(t.test(.SD[,activity_obs], .SD[,activity_exp])$p.value, c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), c("****", "***", "**", "*", "N.S")), .(enh_L, enh_R)]
stat[, diff:= round(mean(activity_obs)/mean(activity_exp), 1), .(enh_L, enh_R)]
offset <- 1.5
cols <- c("activity_exp", "activity_obs")
stat[, segments(xpos_exp[1], max(.SD)+offset, xpos_obs[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, segments(xpos_exp[1], max(activity_exp)+0.2, xpos_exp[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, segments(xpos_obs[1], max(activity_obs)+0.2, xpos_obs[1], max(.SD)+offset, xpd= T), .(enh_L, enh_R), .SDcols= cols]
stat[, text(mean(c(xpos_exp[1], xpos_obs[1])), max(.SD)+offset, paste0("x", diff[1]), xpd= T, pos= 3), .(enh_L, enh_R), .SDcols= cols]
stat[, text(mean(c(xpos_exp[1], xpos_obs[1])), max(.SD)+offset+0.4, pval[1], xpd= T, pos= 3, cex= 1.2), .(enh_L, enh_R), .SDcols= cols]
leg <- c("* pval<5e-2", "** pval<1e-2", "*** pval<1e-3", "**** pval<1e-5")
legend("topleft", paste0(leg), bty= "n", cex= 0.6)
dev.off()

# ### Scatterplot
# collapse_scatter <- stat[, .(mean_exp= mean(activity_exp), mean_obs= mean(activity_obs)), .(enh_L, enh_R)]
# 
# pdf("pdf/C_luciferase_scatterplot.pdf", 4, 4.5)
# Cc <- c("darkgoldenrod1", "cornflowerblue", "tomato", "darkolivegreen2")
# plot(collapse_scatter[, mean_exp:mean_obs], xlim= c(0, 10), ylim= c(0, 10), las= 1, col= Cc, cex= 1.5, pch= 19, xlab= "Expected additive activity", "Activity")
# legend("bottomright", c("hamlet x 2", "hamlet x sup-add.", "sup-add. x hamlet", "sup-add. x 2"), text.col = Cc, bty= "n")
# abline(0, 1, lty= 2)
# dev.off()



