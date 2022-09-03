setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)
require(vlfunctions)
require(data.table)

######################################################
# BERNCHMAK (see original code from Bernardo lower down)
######################################################
# Import and sort based on counts
test <- fread("db/umi_counts/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_1000_F_input_vllib002-spike_0.01_final_rep1__HJFH3BGXF_1_20200403B_20200404.txt.gz")
set.seed(1)
test <- test[sample(.N, 1e6)]
test[, UMI_counts := .N, .(L, R, UMI)]
test[, pair_counts := .N, .(L, R)]
setorderv(test, 
          c("pair_counts", "UMI_counts"), 
          order = c(-1, -1))
test[, ID:= paste0(L, "__", R)]
test[, ID:= factor(ID, unique(ID))]

######------ Stark way -------#######
t1 <- Sys.time()
results <- mclapply(split(test$UMI, test$ID),
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
results <- rbindlist(results, idcol= T)
setnames(results, c("ID", "UMI", "counts"))
t2 <- Sys.time()

######------ My way -------#######
t3 <- Sys.time()
test[, collapsed:= .N==1, .(L, R)]
while(any(!test$collapsed))
{
  test[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist(UMI[1],
                       UMI,
                       method="hamming",
                       nthread= getDTthreads()-1)<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, .(L, R)]
}
t4 <- Sys.time()


# Compare
paste0("Stark:", t2-t1)
paste0("VL:", t4-t3)

results[, ID:= factor(ID, levels= levels(test$ID))]
stark <- results[counts>0][order(ID, UMI), .(ID, UMI)]
vl <- unique(test[order(ID, UMI), .(ID, UMI)])
identical(stark, vl)

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
