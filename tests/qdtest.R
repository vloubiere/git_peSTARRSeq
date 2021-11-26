setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")
setorderv(dat, c("active_group", "cdition"))

# SmoothScatter observed/expected
par(mfrow= c(3,4))
dat[!is.na(active_group), {
  if(nrow(.SD)>10)
  {
    smoothScatter(additive_merge, 
                  log2FC_merge, 
                  main= paste(active_group[1], cdition[1]))
    abline(0, 1, lty= 2)
    legend("topleft", 
           legend = round(cor.test(additive_merge, log2FC_merge)$estimate, 2), 
           bty= "n")
  }
  print(cdition)
}, .(active_group, cdition)]

# SmoothScatter observed/expected
par(mfrow= c(1,2), 
    las= 2,
    mar= c(16,4.1,4,1))
dat[!is.na(active_group), {
  boxplot(diff_merge~cdition, 
          .SD,
          outline= F, 
          xlab= NA,
          main= active_group[1],
          ylim= c(-6, 6), 
          notch= T)
  abline(h= 0, lty= 2)
  print("")
}, active_group]

# Heatmap left/right activity
pl <- merge(unique(dat[!is.na(active_group_L), .(enh= L, cdition, median_L= median_merge_L)]),
            unique(dat[!is.na(active_group_R), .(enh= R, cdition, median_R= median_merge_R)]))
pl[, diff:= median_R-median_L]
pl[, c("lib", "CP", "spacer_name", "spacer_length"):= tstrsplit(cdition, "_", keep= c(1, 2, 4, 5))]
pl <- pl[, .(value= median(diff, na.rm= T)), .(lib, CP, spacer_name, spacer_length)]
mat <- dcast(pl, 
             CP~paste0(lib, "_", spacer_name, "_", spacer_length), 
             value.var = "value")
vl_heatmap(mat, 
           display_numbers=T,
           cluster_rows= F,
           cluster_cols= F,
           breaks= c(-1,0,1))

boxplot(dat[cdition=="vllib021_DSCP_T12_enh-intron1_300" 
            & grepl("^control", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib021_DSCP_T12_enh-intron1_300" 
            & grepl("^SU", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib021_DSCP_T12_enh-intron1_300" 
            & grepl("^control", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib021_DSCP_T12_enh-intron1_300" 
            & grepl("^SU", L)
            & grepl("^SU", R), log2FC_merge])

boxplot(dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^control", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^control", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^SU", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^SU", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^control", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^control", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^SU", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^SU", L)
            & grepl("^SU", R), log2FC_merge], 
        notch= T)


test <- merge(dat[cdition=="vllib018_DSCP_T12_intron4_2000", .(L, R, log2FoldChange)],
              dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000", .(L, R, log2FoldChange)],
              by= c("L", "R"))
smoothScatter(test[[3]], test[[4]], notch= T)
abline(0,1)

test <- merge(dat[cdition=="vllib018_DSCP_T12_intron4_2000", .(L, R, diff)],
              dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000", .(L, R, diff)],
              by= c("L", "R"))
smoothScatter(test[[3]], test[[4]], notch= T)
abline(0,1)

boxplot(dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^control", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^control", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^SU", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^SU", L)
            & grepl("^control", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^control", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^control", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib018_DSCP_T12_intron4_2000" 
            & grepl("^SU", L)
            & grepl("^SU", R), log2FC_merge],
        dat[cdition=="vllib022_DSCP_T12_enh-intron3_2000" 
            & grepl("^SU", L)
            & grepl("^SU", R), log2FC_merge], 
        notch= T)



