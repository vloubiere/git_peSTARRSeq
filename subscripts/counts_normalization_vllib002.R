setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- data.table(file= list.files("db/merged_counts", "vllib002", full.names = T))
dat <- dat[, fread(file), file]
dat[, cdition:= tstrsplit(basename(file), "peSTARRSeq_|_merged", keep= 2), file]
dat <- dat[type=="pair"]

mat <- dcast(dat,
             L+R~cdition,
             value.var = "umi_counts",
             fill = 0)
check <- mat[, apply(.SD, 1, sum), .SDcols= patterns("rep")]
mat <- mat[(check>20)]
cols <- grep("rep", names(mat))
mat[, (cols):= lapply(.SD, function(x) (x+1)/sum(x)*1e6), .SDcols= cols]
mat[, FC_rep1:= DSCP_rep1/input_rep1]
mat[, FC_rep2:= DSCP_rep2/input_rep2]

res <- merge(mat[!xor(grepl("^control", L), grepl("^control", R)), .(L, R, FC_rep1, FC_rep2)],
             mat[grepl("^control", R), .(median_rep1_L= median(FC_rep1),
                                         median_rep2_L= median(FC_rep2)), L],
             by= "L")
res <- merge(res,
             mat[grepl("^control", L), .(median_rep1_R= median(FC_rep1),
                                         median_rep2_R= median(FC_rep2)), R],
             by= "R")

sel_act <- res[, !grepl("control", L) & 
                 !grepl("control", R) &
                 log2(median_rep1_L) > quantile(log2(res[grepl("^control", L), median_rep1_L]), 0.95) &
                 log2(median_rep2_L) > quantile(log2(res[grepl("^control", L), median_rep2_L]), 0.95) &
                 log2(median_rep1_R) > quantile(log2(res[grepl("^control", R), median_rep1_R]), 0.95) &
                 log2(median_rep2_R) > quantile(log2(res[grepl("^control", R), median_rep2_R]), 0.95)]
sel_ctl <- res[,grepl("control", L) & grepl("control", R) & L!=R]

par(mfrow= c(2,2))
smoothScatter(log2(res[sel_act, median_rep1_L+median_rep1_R]),
              log2(res[sel_act, FC_rep1]))
abline(0, 1)
smoothScatter(log2(res[sel_act, median_rep2_L+median_rep2_R]),
              log2(res[sel_act, FC_rep2]))
abline(0, 1)
smoothScatter(log2(res[sel_act, median_rep1_L+median_rep1_R]-median(res[sel_ctl,FC_rep1])),
              res[sel_act, log2(FC_rep1)])
abline(0, 1)
smoothScatter(log2(res[sel_act, median_rep2_L+median_rep2_R]-median(res[sel_ctl,FC_rep2])),
              res[sel_act, log2(FC_rep2)])
abline(0, 1)

mod <- data.table(y= res[sel_act, log2(FC_rep1)],
                  x= log2(res[sel_act, median_rep1_L+median_rep1_R]-median(res[sel_ctl,FC_rep1])))
summary(lm(y~x, mod))

mod <- data.table(y= res[sel_act, log2(FC_rep1)],
                  x= log2(res[sel_act, median_rep1_L+median_rep1_R]))
summary(lm(y~x, mod))
