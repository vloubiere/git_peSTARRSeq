setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

downsize <- 3
for(file in c("db/umi_counts/DSCP_large_WT_input_rep1.txt",
              "db/umi_counts/DSCP_large_WT_input_rep2.txt",
              "db/umi_counts/DSCP_large_WT_screen_rep1.txt",
              "db/umi_counts/DSCP_large_WT_screen_rep2.txt"))
{
  # Output file
  output <- gsub(".txt$", "_downsampled.txt", file)
  if(!file.exists(output))
  {
    umi <- fread(file)
    idx <- sample(seq(nrow(umi)),
                  round(sum(umi$umi_counts)/downsize),
                  replace = T,
                  prob= umi$umi_counts)
    umi <- umi[idx, .(umi_counts= .N), .(L, R, total_counts)]
    fwrite(umi,
           output,
           col.names = T,
           row.names = F,
           sep= "\t",
           quote= F,
           na = NA)
  }
  print(output)
}

# Compute log2FC
FC_file <- "db/FC_tables/large_WT_downsampled_FC_DESeq2.rds"
if(!file.exists(FC_file))
{
  cmd <- "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript /groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_compute_activities.R /groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/DSCP_large_WT_input_rep1_downsampled.txt,/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/DSCP_large_WT_input_rep2_downsampled.txt /groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/DSCP_large_WT_screen_rep1_downsampled.txt,/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/DSCP_large_WT_screen_rep2_downsampled.txt DESeq2 _.*_ _.*_ /groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/ large_WT_downsampled"
  vl_bsub(cmd)
}
down <- readRDS(FC_file)[!grepl("^control", L) & !grepl("^control", R)]
down[, exp:= log2(2^indL+2^indR-1)]
full <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")[!grepl("^control", L) & !grepl("^control", R)]
full[, exp:= log2(2^indL+2^indR-1)]

# Plot 
pdf("pdf/draft/review_impact_seq_depth.pdf", width = 8, height = 3.5)
layout(matrix(1:4, ncol= 4),
       widths= c(2,1,2,1))
full[, {
  vl_par(mai= c(.9, .6, .9, .4))
  smoothScatter(exp,
                log2FoldChange,
                xlim= c(-4, 7),
                main= "Full")
  abline(0, 1)
  vl_par(mai= c(.9, .2, .9, .8))
  .SD[actL!="Inactive" & actR!="Inactive", {
    vl_boxplot(log2FoldChange-exp, violin= T)
    med <- median(log2FoldChange-exp, na.rm= T)
    text(1,
         med,
         round(med, 3),
         pos= 4,
         xpd= T)
  }]
  .SD
}]
down[, {
  vl_par(mai= c(.9, .6, .9, .4))
  smoothScatter(exp,
                log2FoldChange,
                xlim= c(-4, 7),
                main= "Sub-sampled")
  abline(0, 1)
  vl_par(mai= c(.9, .2, .9, .8))
  .SD[actL!="Inactive" & actR!="Inactive", {
    vl_boxplot(log2FoldChange-exp, violin= T)
    med <- median(log2FoldChange-exp, na.rm= T)
    text(1,
         med,
         round(med, 3),
         pos= 4,
         xpd= T)
  }]
  .SD
}]
dev.off()