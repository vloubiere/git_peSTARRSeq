#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############--------------------------------------------------------------------##############
############  Compute individual and paired activities from pSTARR-Seq counts   ##############
############--------------------------------------------------------------------##############
# Takes as input umi collapsed counts and computes log2FoldChange, using either DESeq2 approach or log2 ratios

# test if there is at least 2 args: if not, return an error
if (length(args)!=7) {
  stop(paste0(length(args), " / 7 required arguments were provided:\n",
              paste0(args, "\n"),
              "Please specify:\n
              [required] 1/ A list of comma-separated list of input counts (each file considered as a separate replicate) \n
              [required] 2/ A list of comma-separated list of screen counts (each file considered as a separate replicate) \n
              [required] 3/ The method to be used. One of 'DESeq2' or 'ratio' (the later one tolerates one replicate) \n
              [required] 4/ Regexpr to be applied to 5' (Left) IDs \n
              [required] 5/ Regexpr to be applied to 3' (Right) IDs \n
              [required] 6/ Output DESeq2 dds file (.dds)\n
              [required] 7/ FoldChange output (.rds)\n"))
}

.libPaths(c("/users/vincent.loubiere/R/x86_64-pc-linux-gnu-library/3.6",  "/software/2020/software/r/3.6.2-foss-2018b/lib64/R/library"))
require(DESeq2)
require(data.table)

# Variables and checks ----
inputs <- unique(unlist(tstrsplit(args[1], ",")))
screens <- unique(unlist(tstrsplit(args[2], ",")))
method <- args[3]
sublibRegexprL <- args[4]
sublibRegexprR <- args[5]
dds_file <- args[6]
FC_file <- args[7]
if(!method %in% c("DESeq2", "ratio"))
  stop("Third argument 'method' should be one of 'DESeq2' or 'ratio'.")
if(!grepl(".dds$", dds_file))
  stop("Sixth argument 'dds_file' should finish with .dds extension.")
if(!grepl(".rds$", FC_file))
  stop("Seventh argument 'FC_file' should finish with .rds extension.")

# # Tests
# inputs <- c("db/umi_counts/vllib002_input_rep1.txt", "db/umi_counts/vllib002_input_rep2.txt")
# screens <- c("db/umi_counts/vllib002_screen_rep1.txt", "db/umi_counts/vllib002_screen_rep2.txt")
# sublibRegexprL <- ".*"
# sublibRegexprR <- ".*"
# dds_file <- 

# Import counts ---
dat <- list(input= inputs, screen= screens)
dat <- lapply(dat, function(x) data.table(file= x))
dat <- rbindlist(dat, idcol = "cdition")
dat[, rep:= rowid(cdition)]
dat <- dat[, fread(file, colClasses = c("character", "character", "numeric", "numeric")), .(cdition, rep, file)]
# Initiate report
report <- dat[, .(total_counts= sum(total_counts),
                  umi_counts= sum(umi_counts)), .(cdition, rep, file)]
# Select reads with correct combinations
dat <- dat[grepl(sublibRegexprL, L) & grepl(sublibRegexprR, R)]
report$umi_counts_correct_combinations <- dat[, sum(umi_counts), .(cdition, rep)]$V1
# Cast counts
counts <- dcast(dat, 
                L+R~cdition+rep, 
                value.var = "umi_counts", 
                fun.aggregate = sum)
# Remove low counts and homotypic pairs
check <- apply(counts[, .SD, .SDcols= patterns("input")], 1, function(x) all(x>=5))
Nthreshold <- paste0(sum(check), " pairs out of ", length(check), " had enough reads to pass DESeq2 threshold")
# Remove homotypic pairs
Nhomo <- paste0(nrow(counts[L==R]), " homotypic pairs removed")
counts <- counts[L!=R & (check)]
cols <- setdiff(names(counts), c("L", "R"))
# Add pseudocount
counts[, (cols):= lapply(.SD, function(x) x+1), .SDcols= cols]

# Compute log2FoldChanges ----
if(method=="DESeq2")
{
  # DESeq2  ----
  DF <- data.frame(counts[, !c("L", "R")], row.names = counts[, paste0(L, "__", R)])
  sampleTable <- SJ(name= colnames(DF))
  sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_")], 
                            row.names = "name")
  dds <- DESeqDataSetFromMatrix(countData= DF,
                                colData= sampleTable,
                                design= ~rep+cdition)
  sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[grepl("^control.*__control.*", rownames(DF)),]))
  dds <- DESeq(dds)
  # Differential expression
  norm <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
  norm <- as.data.table(norm, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
}else if(method=="ratio"){
  ## Ratio method (tolerates one replicate) ----
  norm <- copy(counts)
  cols <- grep("^input|^screen", names(norm), value= T)
  norm[, norm_input:= rowSums(.SD), .SDcols= patterns("^input")]
  norm[, norm_input:= norm_input/sum(norm_input)*1e6]
  norm[, norm_screen:= rowSums(.SD), .SDcols= patterns("^screen")]
  norm[, norm_screen:= norm_screen/sum(norm_screen)*1e6]
  # FoldChange
  norm[, log2FoldChange:= log2(norm_screen/norm_input)]
}

# Center and compute inidividual activities ----
# Select usable, inactive control pairs
ctlPairs <- norm[grepl("^control", L) & grepl("^control", R)]
inactCtlL <- ctlPairs[, mean(log2FoldChange), L][between(scale(V1), -1, 1), L]
inactCtlR <- ctlPairs[, mean(log2FoldChange), R][between(scale(V1), -1, 1), R]
norm[, ctlL:= L %in% inactCtlL]
norm[, ctlR:= R %in% inactCtlR]
NctlL <- paste0("Out of ", length(unique(ctlPairs[, L])), " 5' control sequences, ", length(unique(norm[(ctlL), L])), " passed quality thresholds and were used")
NctlR <- paste0("Out of ", length(unique(ctlPairs[, R])), " 3' control sequences, ", length(unique(norm[(ctlR), R])), " passed quality thresholds and were used")
if(method=="ratio") # Ratio method -> center on inactive pairs (already done earlier for DESeq2)
  norm[, log2FoldChange:= log2FoldChange-median(norm[ctlL & ctlR, log2FoldChange])]
# Compute individual act and pval (vs ctlPairs)
ctlPairs <- norm[ctlL & ctlR, log2FoldChange]
Left <- norm[(ctlR), 
             .(.N>=10, 
               wilcox.test(log2FoldChange, ctlPairs, alternative = "greater")$p.value,
               mean(log2FoldChange)), L][(V1)]
Left[, padj:= p.adjust(V2, "fdr")]
norm[Left, c("indL", "padjL"):= .(i.V3, i.padj), on= "L"]
Right <- norm[(ctlL), 
              .(.N>=10, 
                wilcox.test(log2FoldChange, ctlPairs, alternative = "greater")$p.value,
                mean(log2FoldChange)), R][(V1)]
Right[, padj:= p.adjust(V2, "fdr")]
norm[Right, c("indR", "padjR"):= .(i.V3, i.padj), on= "R"]
# Remove pairs for which combined or ind act could not be computed accurately
Npairs <- paste0("Individual + combined activity of ", nrow(norm[!is.na(indL) & !is.na(indR)]), " pairs out of ", nrow(norm), " could be properly inferred")
norm <- norm[!is.na(indL) & !is.na(indR)]
# Add activity classes
norm[, groupL:= fcase(grepl("^control", L), "Random seq.", default = "Candidate seq.")]
norm[, groupR:= fcase(grepl("^control", R), "Random seq.", default = "Candidate seq.")]
# Split Active enhancers based on their strength
norm[, actL:= fcase(between(indL, 1, 2) & padjL<0.05, "Low",
                    between(indL, 2, 4) & padjL<0.05, "Medium",
                    between(indL, 4, Inf) & padjL<0.05, "Strong",
                    default = "Inactive")]
norm[, actR:= fcase(between(indR, 1, 2) & padjR<0.05, "Low",
                    between(indR, 2, 4) & padjR<0.05, "Medium",
                    between(indR, 4, Inf) & padjR<0.05, "Strong",
                    default = "Inactive")]
cols <- c("actL", "actR")
norm[, (cols):= lapply(.SD, function(x) factor(x, c("Inactive", "Low", "Medium", "Strong"))), .SDcols= cols]
# Handle missing columns
cols <- c("L", "R",
          "groupL", "groupR", "actL", "actR",
          "indL", "padjL", "indR", "padjR", 
          "log2FoldChange", "padj",
          "ctlL", "ctlR")
cols <- cols[cols %in% names(norm)]
norm <- norm[, cols, with= F]

## SAVE ----
saveRDS(dds, dds_file)
saveRDS(norm, FC_file)
report <- rbind(report, 
                data.table(cdition= c("", Nthreshold, Nhomo, NctlL, NctlR, Npairs)),
                fill= T)
fwrite(report,
       gsub(".rds$", "_report.txt", FC_file))

