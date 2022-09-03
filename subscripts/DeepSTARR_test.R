setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

###########################
# Import
###########################
# Import feat
feat <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
# Add chromosome sizes and compute sets
chr_sizes <- as.data.table(GenomicFeatures::getChromInfoFromUCSC("dm3"))
feat[chr_sizes, length:= i.length, on= "seqnames==chrom"]
feat[, set:= fcase(seqnames=="chr2R" & between(start, 1, length/2), "validation",
                   seqnames=="chr2R" & between(start, length/2, length), "test", 
                   default= "train"), .(seqnames, length)]

# Import vllib002 and merge with feat
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002"]
lib[feat, c("set_L", "seq_L"):= .(i.set, i.enh_sequence), on= "L==ID_vl"]
lib[feat, c("set_R", "seq_R"):= .(i.set, i.enh_sequence), on= "R==ID_vl"]
# Add residuals from CV lm
model <- readRDS("Rdata/CV_linear_model_vllib002.rds")$pred
lib[model, residuals:= log2FoldChange-i.predicted, on= c("L", "R")]
# Keep only the lines for which both enhancers are in the same set
dat <- lib[set_L==set_R, .(L, R, log2FoldChange, residuals, set= set_L, seq_L= seq_L, seq_R= seq_R, median_L, median_R, norm_input, norm_screen_rep1, norm_screen_rep2)]

###########################
# Save fasta files
###########################
if(F)
{
  dat[, seqinr::write.fasta(sequences = as.list(seq_L), 
                            names = L, 
                            file.out = paste0("db/for_Bernado/deepSTARR/", set, "_left.fasta"), 
                            as.string = T), set]
  dat[, seqinr::write.fasta(sequences = as.list(seq_R), 
                            names = R, 
                            file.out = paste0("db/for_Bernado/deepSTARR/", set, "_right.fasta"), 
                            as.string = T), set]
  dat[, fwrite(.(log2FoldChange),
               paste0("db/for_Bernado/deepSTARR/", set, "_activity.txt"),
               quote= F,
               col.names = F), set]
  dat[, fwrite(.(residuals),
               paste0("db/for_Bernado/deepSTARR/", set, "_residuals.txt"),
               quote= F,
               col.names = F), set]
  dat[, fwrite(.(median_L),
               paste0("db/for_Bernado/deepSTARR/", set, "_activity_L.txt"),
               quote= F,
               col.names = F), set]
  dat[, fwrite(.(median_R),
               paste0("db/for_Bernado/deepSTARR/", set, "_activity_R.txt"),
               quote= F,
               col.names = F), set]
}

###########################
# Train linear model and CV
###########################
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= dat[set=="train"])
dat[, predicted:= predict(model, 
                          newdata = dat)]
par(mfrow=c(2,2))
dat[set!="train", {
  # model
  smoothScatter(log2FoldChange, 
                predicted, 
                main= paste0("model ", set))
  legend("topleft",
         paste0("PCC= ", round(cor.test(log2FoldChange, predicted)$estimate, 2)))
  abline(0,1)
  # replicates
  x <- log2(norm_screen_rep1/norm_input)
  y <- log2(norm_screen_rep2/norm_input)
  smoothScatter(x, 
                y, 
                main= paste0("rep ", set))
  legend("topleft",
         paste0("PCC= ", round(cor.test(x, y)$estimate, 2)))
  abline(0,1)
  print("")
}, set]