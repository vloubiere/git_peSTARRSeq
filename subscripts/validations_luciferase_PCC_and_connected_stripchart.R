dat <- readRDS("Rdata/validations_luciferase_final_table.rds")

# Add features
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
dat[feat, group_L:= i.group, on= "L==ID"]
dat[feat, group_R:= i.group, on= "R==ID"]

# Collapse per enhancer pair
coll <- dat[, .(mean= mean(luc_norm, na.rm= T), 
                sd= sd(luc_norm, na.rm= T),
                mean_add= mean(luc_mean_L+luc_mean_R, na.rm= T),
                sd_add= sd(luc_mean_L+luc_mean_R, na.rm= T)), .(L, R, group_L, group_R)]

# Add STARR-Seq data
# STARR <- fread("db/final_tables_exp_model/vllib002_short_spacer_peSTARRSeq_max_cutoff_final_oe.txt")
STARR <- fread("db/final_tables_exp_model/vllib002_short_spacer_peSTARRSeq_low_cutoff_final_oe.txt") # Use low cutoff to keep control pairs
coll[STARR, log2FoldChange:= i.log2FoldChange, on= c("L", "R")]

# Compute colors and pch
Cc <- data.table(group= c("OSC", "dev", "hk", "inducible", "shared", "control"),
                 col= c("black", "#74C27A", "tomato", "gold", "royalblue2", "lightgrey"))
coll[Cc, col_L:= i.col, on= "group_L==group"]
coll[Cc, col_R:= i.col, on= "group_R==group"]
coll[xor(group_L=="control", group_R=="control"), col:= ifelse(group_L=="control", col_R, col_L)]
coll[xor(group_L=="control", group_R=="control"), pch:= ifelse(group_L=="control", 17, 15)]
coll[!xor(group_L=="control", group_R=="control") & group_L!=group_R, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
coll[!xor(group_L=="control", group_R=="control") & group_L==group_R, col:= col_L, .(col_L, col_R)]
coll[!xor(group_L=="control", group_R=="control"), pch:= 19]

# Add connected strpchart classes
coll[group_L=="control" & group_R=="control", class:= "A_control"]
coll[is.na(class) & (group_L=="hk" | group_R=="hk") & group_L!="control" & group_R!="control", class:= "B_hk containing"]
coll[is.na(class) & (group_L=="shared" | group_R=="shared") & group_L!="control" & group_R!="control", class:= "C_shared containing"]
coll[is.na(class) & group_L=="dev" & group_R=="dev", class:= "D_developmental pairs"]
setorderv(coll, "class")
coll[!is.na(class), stripchart_x:= c(1,2,4,6)[.GRP], class]
coll[, class:= gsub(".*_(.*)", "\\1", class), class]

.lm <- lm(log2(mean)~log2FoldChange, coll)
.lo <- loess(log2(mean)~log2FoldChange, coll)
leg <- unique(coll[, .(group=paste0(group_L,  "+", group_R), col, pch)])
leg <- leg[order(col)]

#----------------------------------------------#
# PLOT
#----------------------------------------------#
# PCC ---------------#
pdf("pdf/luciferase/PCC_luciferase_validations_peSTARRSeq.pdf", 4.5, 5)
par(las= 1)
plot(NA, 
     xlim= range(coll$log2FoldChange, na.rm = T),
     ylim= range(log2(c(coll[, mean+sd], coll[, mean-sd]))),
     xlab= "pe-STARR-Seq activity (log2)",
     ylab= "Normalized luciferase activity (log2)")
coll[, {
  points(log2FoldChange, 
         log2(mean), 
         pch= pch,
         col= col)
  segments(log2FoldChange, 
           log2(mean-sd), 
           log2FoldChange, 
           log2(mean+sd), 
           col= col)
  print("")
}, .(col, pch)]
abline(.lm, lty=2)
lines(seq(-5, 15, 0.1), 
      predict(.lo, newdata = seq(-5, 15, 0.1)), 
      lty=2, 
      lwd=2, 
      col= "red")
legend("bottomright", 
       legend = leg$group, 
       col= leg$col, 
       bty= "n", 
       pch= leg$pch,
       cex= 0.6)
legend("topleft", 
       legend = paste("R²=", round(summary(.lm)$r.squared, 2)), 
       lty= 2, 
       bty= "n")
dev.off()

# STRIPCHART ---------------#
pdf("pdf/luciferase/Stripchart_luciferase_validations_peSTARRSeq.pdf", 4.5, 5)
par(mar= c(7.1, 4.1, 2.1, 2.1),
    las= 2)
plot(NA,
     xlim= c(0.6,7.4),
     xaxt= "n",
     ylim= c(0, 100),
     xlab= NA,
     ylab= "Normalized luciferase activity (log2)")
axis(1, 
     at= 1:7, 
     labels = c("Negative pairs",
                "Exp. additive",
                "Obecsrved",
                "Exp. additive",
                "Obecsrved",
                "Exp. additive",
                "Obecsrved"))
coll[!is.na(class), {
  jit_amount <- 0.2
  jit <- jitter(rep(0, nrow(.SD)),
                amount = jit_amount)
  if(class=="control")
  {
    points(stripchart_x+jit,
           mean,
           col= adjustcolor(col, 0.6),
           pch= pch)
    segments(stripchart_x+jit, 
             mean-sd, 
             stripchart_x+jit,
             mean+sd,
             col= col)
  }else
  {
    points(stripchart_x+jit+1,
           mean,
           col= adjustcolor(col, 0.6),
           pch= pch)
    segments(stripchart_x+jit+1, 
             mean-sd, 
             stripchart_x+jit+1,
             mean+sd,
             col= col)
    points(stripchart_x+jit,
           mean_add,
           col= adjustcolor(col, 0.6),
           pch= pch)
    segments(stripchart_x+jit, 
             mean_add-sd_add, 
             stripchart_x+jit,
             mean_add+sd_add,
             col= col)
    segments(stripchart_x+jit, 
             mean_add, 
             stripchart_x+1+jit, 
             mean, 
             col= adjustcolor(col, 0.6))
    
  }
  print("")
}, .(class, stripchart_x)]

legend(0.75, 
       105,
       bty= "n",
       col= leg[group=="control+control", col],
       pch= leg[group=="control+control", pch],
       legend = "ctl+ctl",
       cex= 0.6)
legend(2.25, 
       105,
       bty= "n",
       col= leg[grepl("hk", group) & !grepl("control", group), col],
       pch= leg[grepl("hk", group) & !grepl("control", group), pch],
       legend = leg[grepl("hk", group) & !grepl("control", group), group],
       cex= 0.6)
legend(4.25, 
       105,
       bty= "n",
       col= leg[grepl("shared", group) & !grepl("hk|control", group), col],
       pch= leg[grepl("shared", group) & !grepl("hk|control", group), pch],
       legend = leg[grepl("shared", group) & !grepl("hk|control", group), group],
       cex= 0.6)
legend(6.25,
       105,
       bty= "n",
       col= leg[group=="dev+dev", col],
       pch= leg[group=="dev+dev", pch],
       legend = leg[group=="dev+dev", group],
       cex= 0.6)
dev.off()
