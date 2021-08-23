setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggrepel)

dat <- fread("Rdata/vllib002_modelling.txt")

#------------------------#
# Residuals lm3
#------------------------#
dat[, residuals3:= log2FoldChange-predict(.lm3)]
res <- melt(dat,
            id.vars = c("L", "R", "residuals3"),
            measure.vars = patterns("motif"= "^motif"))
res[, variable:= gsub("__L|__R", "", variable)]
res[, value:= 2^value-1]
res <- res[, .(value= log2(sum(value)+1)), L:variable]
# high pairs
res[, c("p.value_high", "OR_high"):= {
  .c <- fisher.test(value>0, residuals3>2)
  .(.c$p.value, .c$estimate)
}, variable]
# low pairs
res[, c("p.value_low", "OR_low"):= {
  .c <- fisher.test(value>0, residuals3<(-2))
  .(.c$p.value, .c$estimate)
}, variable]
# Melt table, add colors and names
res <- melt(res,
            id.vars = "variable",
            variable.name = "trend",
            measure.vars = patterns("pval"= "^p.val", "OR"= "^OR"))
res <- unique(res)
res[trend== 1, Cc:= "tomato"]
res[trend== 2, Cc:= "cornflowerblue"]
res[, padj:= p.adjust(pval, method= "fdr"), trend]
res[, name:= variable]
for(i in seq(nrow(som)))
  res[, name:= gsub(som[i, motif], som[i, BA_cluster], name), name]


#------ Volcano plot motifs ------#
pdf("pdf/modeling/vllib002_modeling_volcano_motifs.pdf", 6, 6)
ggplot(res, aes(log2(OR), -log10(padj), label = trend)) +
  geom_point(color = res$Cc, cex= 2) +
  geom_label_repel(data= res[padj<0.05], 
                   aes(label= gsub("motif__", "", name)),
                   fontface = 'bold', 
                   segment.color = 'black',
                   color = res[padj<0.05, Cc])
dev.off()


