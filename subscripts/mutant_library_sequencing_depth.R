old <- rbindlist(list(screen_rep1= fread("db/umi_counts/vllib002_screen_rep1.txt"),
                      screen_rep2= fread("db/umi_counts/vllib002_screen_rep2.txt"),
                      input_rep1= fread("db/umi_counts/vllib002_input_rep1.txt"),
                      input_rep2= fread("db/umi_counts/vllib002_input_rep2.txt")), idcol = T)
old <- old[, sum(umi_counts), .id]
new <- rbindlist(list(screen_rep1= fread("db/umi_counts/vllib029_screen_rep1.txt"),
                      screen_rep2= fread("db/umi_counts/vllib029_screen_rep2.txt"),
                      input_rep1= fread("db/umi_counts/vllib029_input_rep1.txt"),
                      input_rep2= fread("db/umi_counts/vllib029_input_rep2.txt")), idcol = T)
new <- new[, sum(umi_counts), .id]
new[, c("cdition", "rep"):= tstrsplit(.id, "_")]

# Show impact of sequencing on individual activity
dat <- data.table(file= list.files("db/umi_counts", "vllib002_", full.names = T))
dat[, c("cdition", "rep"):= tstrsplit(basename(file), "_|[.]", keep= 2:3)]
setkeyv(dat, c("cdition", "rep"))
dat <- dat[, fread(file), (dat)]
dat[new, size:= i.V1, on= c("cdition","rep")]
# subsampling
sub <- dat[, {
  set.seed(1)
  .s <- sample(.N, size= size, replace = T, prob = umi_counts)
  .SD[.s, .(umi_counts= .N), .(L, R)]
}, .(cdition, rep, size)]
counts <- dcast(sub, 
                L+R~cdition+rep, 
                value.var = "umi_counts", 
                fun.aggregate = sum)
# Remove homotypic pairs and cutoff low input counts
check <- apply(counts[, .SD, .SDcols= patterns("input")], 1, function(x) all(x>=5))
counts <- counts[L!=R & (check)]
cols <- grep("rep", names(counts), value = T)
# Add pseudocount
counts[, (cols):= lapply(.SD, function(x) x+1), .SDcols= cols]
# DF and sampleTable
DF <- data.frame(counts[, !c("L", "R")], row.names = counts[, paste0(L, "__", R)])
sampleTable <- SJ(name= colnames(DF))
sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_rep")], 
                          row.names = "name")
# DESeq
dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                      colData= sampleTable,
                                      design= ~rep+cdition)
sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grepl("^control.*__control.*", rownames(DF)),]))
dds <- DESeq2::DESeq(dds) # Check if missing replicates -> skip
# Differential expression
norm <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
norm <- as.data.table(norm, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
# Select usable, inactive control pairs
ctlPairs <- norm[grepl("^control", L) & grepl("^control", R)]
inactCtlL <- ctlPairs[, mean(log2FoldChange), L][between(scale(V1), -1, 1), L]
inactCtlR <- ctlPairs[, mean(log2FoldChange), R][between(scale(V1), -1, 1), R]
norm[, ctlL:= L %in% inactCtlL]
norm[, ctlR:= R %in% inactCtlR]
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

full <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
compare <- merge(full, norm, by= c("L", "R"))
par(mfrow=c(2,1))
plot(unique(compare[, .(L, `Full seq. depth`= indL.x, `Sub-sampled`= indL.y)])[, -1])
abline(0,1)
plot(unique(compare[, .(R, `Full seq. depth`= indR.x, `Sub-sampled`= indR.y)])[, -1])
abline(0,1)