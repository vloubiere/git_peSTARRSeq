setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

if(!exists("dat"))
  dat <- readRDS("Rdata/final_results_table.rds")

sub <- merge(unique(dat[group_L %in% c("dev","hk"), .(cdition, coor= coor_L, group= group_L, enhancer= L, median_merge_L)]),
             unique(dat[group_R %in% c("dev","hk"), .(cdition, coor= coor_R, group= group_R, enhancer= R, median_merge_R)]))
sub[, group:= factor(group, levels= c("hk", "dev"))]
sub[, c("lib", "CP", "spacer"):= tstrsplit(cdition, "_", keep= c(1,2,5))]
sub[, diff:= median_merge_R-median_merge_L]

pl <- sub[, .(diff= list(na.omit(diff))), .(lib, CP, spacer, group)]

par(mar= c(15,5,2,2),
    mfrow= c(1, 2)) 
pl[, {
  boxplot(.SD$diff, 
          names= paste0(.SD$lib, "_", .SD$CP, "_", .SD$spacer, " (", lengths(.SD$diff), ")"),
          las= 2, 
          main= group, 
          notch= T)
  abline(h= 0)
  print("")
}, group]

par(mfrow= c(3,3),
    mar= c(5.1,4.1,4.1,2.1))
for(expr in c("vllib015|vllib018", "vllib021|vllib022", "vllib018|vllib022"))
{
  pcc <- sub[grepl(expr, cdition) & group== "dev",
             {
               plot(rep(c(1, 2), each= .N),
                    c(median_merge_L, median_merge_R), 
                    xlim= c(0, 3),
                    pch= 21,
                    main= cdition[1])
               
               segments(1,
                        c(median_merge_L),
                        2,
                        c(median_merge_R))
               .(enhancer, 
                 median_merge_R-median_merge_L, 
                 rowMeans(.SD[, .(median_merge_L, median_merge_R)]))
             }, cdition]
  pcc <- dcast(pcc, enhancer~cdition, value.var = "V2")
  names(pcc) <- gsub("-", "_", names(pcc))
  plot(pcc[[2]], 
       pcc[[3]],
       )
  .f <- paste0(names(pcc)[3], "~", names(pcc)[2])
  abline(h= 0, lty=2)
  abline(v= 0, lty=2)
  abline(lm(as.formula(.f), pcc))
}

sub[, pref:= ifelse(median_merge_L>median_merge_R, "L", "R")]
mot <- vl_motif_counts(unique(as.data.table(GRanges(sub$coor, name= sub$enhancer))), "dm3")
mot <- dcast(mot, seqnames+start+end+strand+name~motif+motif_name, value.var= "motif_counts")
mot[, coor:= paste0(seqnames, ":", start, "-", end, ":", strand)]
sub <- merge(sub, 
             mot[, `bergman__Abd-B_Abd-B`:coor],
             by= "coor",
             allow.cartesian=TRUE)
sub[dat, , cdition]
