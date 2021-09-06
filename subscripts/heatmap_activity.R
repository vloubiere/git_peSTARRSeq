setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import
if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
if("plot_group" %in% names(dat))
  dat$plot_group <- NULL
dat[cdition=="vllib002" & ((group_L=="OSC" & group_R=="OSC" & median_L>1 & median_R>1)), plot_group:= "vllib002: OSC x OSC"]
dat[cdition=="vllib002" & ((group_L=="OSC" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: OSC x dev"]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x OSC"]
dat[cdition=="vllib002" & ((group_L=="hk" & group_R=="hk" & median_L>1 & median_R>1)), plot_group:= "vllib002: hk x hk"]
dat[cdition=="vllib002" & ((group_L=="hk" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: hk x dev"]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="hk" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x hk"]
dat[cdition=="vllib002" & ((group_L=="dev" & group_R=="dev" & median_L>1 & median_R>1)), plot_group:= "vllib002: dev x dev"]

dat[cdition=="vllib015" & group_L=="hk" & group_R=="hk", plot_group:= "vllib015: hk x hk"]
dat[cdition=="vllib015" & group_L=="hk" & group_R=="dev", plot_group:= "vllib015: hk x dev"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="hk", plot_group:= "vllib015: dev x hk"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="dev", plot_group:= "vllib015: dev x dev"]
dat[cdition=="vllib015" & group_L=="DHS" & group_R=="DHS", plot_group:= "vllib015: DHS x DHS"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="DHS", plot_group:= "vllib015: dev x DHS"]
dat[cdition=="vllib015" & group_L=="DHS" & group_R=="dev", plot_group:= "vllib015: DHS x dev"]
dat[cdition=="vllib015" & group_L=="Silencer" & group_R=="Silencer", plot_group:= "vllib015: Sil. x Sil."]
dat[cdition=="vllib015" & group_L=="Silencer" & group_R=="dev", plot_group:= "vllib015: Sil. x dev"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="Silencer", plot_group:= "vllib015: dev x Sil."]
dat[cdition=="vllib015" & group_L=="CP" & group_R=="CP", plot_group:= "vllib015: CP x CP"]
dat[cdition=="vllib015" & group_L=="CP" & group_R=="dev", plot_group:= "vllib015: CP x dev"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="CP", plot_group:= "vllib015: dev x CP"]
dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib015: Ins. x Ins."]
dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="dev", plot_group:= "vllib015: Ins. x dev"]
dat[cdition=="vllib015" & group_L=="dev" & group_R=="SUHW_peak", plot_group:= "vllib015: dev x Ins."]
dat[cdition=="vllib015" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib015: Ins. x Ins."]

dat[cdition=="vllib016" & group_L=="hk" & group_R=="hk", plot_group:= "vllib016: hk x hk"]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="dev", plot_group:= "vllib016: hk x dev"]
dat[cdition=="vllib016" & group_L=="dev" & group_R=="hk", plot_group:= "vllib016: dev hk"]
dat[cdition=="vllib016" & group_L=="dev" & group_R=="dev", plot_group:= "vllib016: dev x dev"]
dat[cdition=="vllib016" & group_L=="DHS" & group_R=="DHS", plot_group:= "vllib016: DHS x DHS"]
dat[cdition=="vllib016" & group_L=="DHS" & group_R=="hk", plot_group:= "vllib016: DHS x hk"]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="DHS", plot_group:= "vllib016: hk x DHS"]
dat[cdition=="vllib016" & group_L=="Silencer" & group_R=="Silencer", plot_group:= "vllib016: Sil. x Sil."]
dat[cdition=="vllib016" & group_L=="Silencer" & group_R=="hk", plot_group:= "vllib016: Sil. x hk"]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="Silencer", plot_group:= "vllib016: hk x Sil."]
dat[cdition=="vllib016" & group_L=="SUHW_peak" & group_R=="SUHW_peak", plot_group:= "vllib016: Ins. x Ins."]
dat[cdition=="vllib016" & group_L=="SUHW_peak" & group_R=="hk", plot_group:= "vllib016: Ins. x hk"]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="SUHW_peak", plot_group:= "vllib016: hk x Ins."]
dat[cdition=="vllib016" & group_L=="CP" & group_R=="CP", plot_group:= "vllib016: CP x CP"]
dat[cdition=="vllib016" & group_L=="CP" & group_R=="hk", plot_group:= "vllib016: CP x hk"]
dat[cdition=="vllib016" & group_L=="hk" & group_R=="CP", plot_group:= "vllib016: hk x CP"]


# PLOT
pdf("pdf/all_pairs_activity_additivity.pdf", height = 4.75, width = 5*3)
par(mfrow= c(1, 3))
dat[!is.na(plot_group) & L!=R, {
  # Use expected additive values to order the heatmap
  .add <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "additive")
  .add <- as.matrix(.add, 1)
  row_ord <- order(-rowMeans(.add, na.rm = T))
  col_ord <- order(colMeans(.add, na.rm = T))
  .add <- .add[row_ord,]
  .add <- .add[, col_ord]
  vl_heatmap(mat = .add, 
             main= plot_group, 
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             show_colnames = F,
             col = viridis::plasma(10),
             legend_title = "Exp. Add. (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # Activity
  .act <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "log2FoldChange")
  .act <- as.matrix(.act, 1)
  .act <- .act[row_ord,]
  .act <- .act[, col_ord]
  vl_heatmap(mat = .act, 
             main= plot_group, 
             cluster_rows = F,
             cluster_cols = F,
             col = viridis::plasma(10),
             show_rownames = F,
             show_colnames = F,
             legend_title = "Activity (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # Additivity
  .diff <- dcast(data = .SD, 
                formula= L~R, 
                value.var= "diff")
  .diff <- as.matrix(.diff, 1)
  .diff <- .diff[row_ord,]
  .diff <- .diff[, col_ord]
  .breaks <- max(abs(quantile(.diff, c(0.05, 0.975), na.rm= T)))
  vl_heatmap(mat = .diff, 
             main= plot_group, 
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             show_colnames = F,
             breaks = c(-.breaks, -0.5, 0.5, .breaks),
             col= c("cornflowerblue", "white", "white", "tomato"),
             legend_title = "Additivity (log2)")
  rect(0,0,1,1,lwd= 0.5)
  # # Motif
  # .mot <- dcast(data = .SD[, .(L, R, value= log2(motif__flyfactorsurvey__suHw_FlyReg_FBgn0003567_R+1))], 
  #                formula= L~R, 
  #                value.var= "value")
  # .mot <- as.matrix(.mot, 1)
  # .mot <- .mot[row_ord,]
  # .mot <- .mot[, col_ord]
  # .breaks <- max(abs(quantile(.mot, c(0.05, 0.975), na.rm= T)))
  # vl_heatmap(mat = .mot, 
  #            main= plot_group, 
  #            cluster_rows = F,
  #            cluster_cols = F,
  #            show_rownames = F,
  #            show_colnames = F,
  #            # col= c("cornflowerblue", "white", "white", "tomato"),
  #            col = viridis::plasma(10),
  #            # breaks = c(-.breaks, -0.5, 0.5, .breaks),
  #            legend_title = "Motif counts")
  # rect(0,0,1,1)
  print("DONE")
}, plot_group]
dev.off()


