setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)
require(vlfunctions)
require(data.table)

######################################################
# BERNCHMAK (see original code from Bernardo lower down)
######################################################
# Import and sort based on counts
test <- fread("db/umi_counts/vllib015_pe-STARR-Seq-2.0_DSCP_T12_SCR1_300_1000_T_screen_NA_NA_final_rep1__H2KMLBGXK_1_20210831B_20210901.txt.gz")
set.seed(1)
test <- test[sample(nrow(test), 2e6)]
test[, ID:= paste0(L, "__", R)]

######------ Stark way -------#######
t1 <- Sys.time()
stark <- test[, .(umi_N= .N), .(L, R, UMI)]
setorderv(stark, "umi_N", order = -1)
stark <- mclapply(split(stark$UMI, stark[, .(L,R)]),
                  function(x)
                  {
                    names(x) <- x
                    keep <- vector("numeric", length(x))
                    names(keep) <- x
                    while (length(x)>0) {
                      rmv = which(stringdist(x[1], x, method="hamming", nthread= 1) <= 1)
                      keep[names(x)[1]] <- length(rmv)
                      x <- x[-rmv]
                    }
                    return(as.data.table(keep, keep.rownames = "UMI"))
                  })
stark <- rbindlist(stark, idcol= T)
stark[, c("L", "R"):= tstrsplit(.id, "[.]")]
stark <- stark[keep>0][order(L, R, UMI), .(L, R, UMI)]
t2 <- Sys.time()

######------ My way -------#######
t3 <- Sys.time()
.c <- test[, .(umi_N= .N), .(L, R, UMI)]
setorderv(.c, "umi_N", order = -1)
.c[, collapsed:= .N==1, .(L, R)]
while(any(!.c$collapsed))
{
  .c[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist(UMI[1],
                       UMI,
                       method="hamming",
                       nthread= getDTthreads()-1)<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, .(L, R)]
}
.c <- unique(.c[, .(L, R, UMI)])
t4 <- Sys.time()

# Compare
paste0("Stark:", t2-t1)
paste0("VL:", t4-t3)

res <- merge(test[, .(before= .N), .(L, R)],
             stark[, .(stark= .N), .(L, R)],
             by= c("L", "R"))
res <- merge(res,
             .c[, .(vl= .N), .(L, R)],
             by= c("L", "R"))
identical(res$stark, res$vl)
# par(mfrow=c(2,2))
# res[, {
#   plot(before, stark)
#   abline(0,1)
#   plot(before, vl)
#   abline(0,1)
#   plot(stark, vl)
#   abline(0,1)
# }]

######################################################
# BERNARDO
######################################################
# #!/usr/bin/env Rscript
# 
# options(stringsAsFactors=FALSE)
# options("scipen"=100, "digits"=4)
# 
# #test options
# opt=list()
# #opt$input <- "/scratch/stark/almeida/20190315_starrseq/data/tmp.InfiNOD0Tr/collapsed_frags_sub.bed"
# opt$input <- "/Volumes/clustertmp/stark/almeida/20190315_starrseq/data/tmp.InfiNOD0Tr/collapsed_frags_sub.bed"
# opt$MM <- 1
# opt$core <- 1
# opt$out <- "filtered.bed"
# 
# # OPTION PARSING
# suppressPackageStartupMessages(library("optparse"))
# 
# option_list <- list(
#   make_option(c("-i", "--input"),  action = "store", type="character", default=NULL,
#               help="$TMP/collapsed_frags.bed", metavar="character"),
#   make_option(c("-m", "--MM"),  action = "store", type="numeric", default=1,
#               help="$UMI_MM", metavar="character"),
#   make_option(c("-c", "--core"),  action = "store", type="numeric", default=1,
#               help="number of cores", metavar="character"),
#   make_option(c("-o", "--out"),  action = "store", type="character", default="average.txt",
#               help="$TMP/reads.filtered.3.bed", metavar="character")
# )
# 
# opt_parser <- OptionParser(
#   usage = "%prog [options]",
#   option_list=option_list,
#   description = "UMI collapsing"
# )
# arguments <- parse_args(opt_parser, positional_arguments = TRUE)
# opt <- arguments$options
# 
# # LIBRARIES
# suppressPackageStartupMessages(library("rtracklayer"))
# suppressPackageStartupMessages(library("parallel"))
# suppressPackageStartupMessages(library("stringdist"))
# 
# # print options
# cat("\nRunning UMI collapsing\n")
# 
# opt
# 
# # Prepare data
# test_big <- import.bed(opt$input)
# test_big$ID <- paste(test_big@seqnames, test_big@ranges@start, test_big@ranges@start+test_big@ranges@width-1, test_big@strand, sep="_")
# lvl <- names(sort(table(test_big$ID)))
# test_big$ID <- factor(test_big$ID, levels = rev(lvl))
# test_big_sorted <- test_big[order(test_big$ID)]
# 
# f3 = function(bar,c=1){
#   keep <- vector("numeric",length(bar))
#   names(keep) = names(bar)
#   while (length(bar)>0) {
#     rmv = which(stringdist(bar[1],bar,method="hamming",nthread =c)<=opt$MM)
#     keep[names(bar)[1]] = length(rmv)
#     bar = bar[-rmv]
#   }
#   return(keep)
# }
# 
# # Run
# # Calculate the number of cores
# no_cores <- opt$core
# # Initiate cluster
# cl <- makeCluster(no_cores)
# 
# clusterExport(cl, "stringdist")
# clusterExport(cl, "f3")
# clusterExport(cl, "opt")
# 
# results <- unlist(parLapply(cl, split(test_big_sorted$name,test_big_sorted$ID),
#                             function(x){
#                               names(x) <- 1:length(x)
#                               f3(x)
#                             })
# )
# 
# # Finish
# stopCluster(cl)
# 
# test_big_sorted$counts <- as.numeric(results)
# # print to see top
# test_big_sorted
# test_big_sorted <- test_big_sorted[test_big_sorted$counts>0]
# test_big_sorted
# 
# out <- as.data.frame(test_big_sorted)
# out$start <- out$start-1
# out$name <- paste(out$name, out$score, sep="_")
# 
# write.table(out[,c(1:3,6,9,5)], opt$out, sep="\t", row.names = F, col.names = F, quote = F)
# 
# sessionInfo()
#
