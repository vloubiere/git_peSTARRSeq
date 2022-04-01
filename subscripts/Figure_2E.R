setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

# Import
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")
screen <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class== "enh./enh."]
screen[, diff:= log2FoldChange-additive]
screen <- feat$add_feature(screen, feat$lib)
cl <- readRDS("Rdata/vllib002_clustering_additive_scores_draft_figure.rds")
motifs <- rbindlist(list(L= cl$mot_enr_L$enr,
                         R= cl$mot_enr_R$enr), 
                    idcol = T)
motifs[vl_Dmel_motifs_DB_full, motif_cluster:= i.Motif_cluster_name, on= "variable==uniqName_noSpecialChar"]
motifs <- motifs[order(abs(log2OR), decreasing = T)]
motifs <- motifs[padj<0.001, .SD[1], .(.id, motif_cluster)]
motifs <- unique(motifs[, .(.id, motif_cluster, variable)])

# Build data
counts_L <- as.data.table(cl$mot_enr_L$counts[, motifs[.id=="L", variable]])
setnames(counts_L, function(x) paste0("motif_", x, "_L"))
counts_R <- as.data.table(cl$mot_enr_R$counts[, motifs[.id=="R", variable]])
setnames(counts_R, function(x) paste0("motif_", x, "_R"))
dat <- cbind(screen,
             counts_L[match(screen$L, cl$rows$name),],
             counts_R[match(screen$R, cl$cols$name),])

# Sample for CV
set.seed(1)
sel_L <- dat[, 1-.N/dat[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- dat[, 1-.N/dat[,.N], R][, sample(R, round(.N/10), prob = V1)]
dat[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]

# Train linear model
mot <- grep("^motif_", names(dat), value= T)
form <- paste0("log2FoldChange~median_L*median_R+", paste0(mot, collapse= "+"))
model <- lm(formula = as.formula(form), 
            data= dat[set=="train"])

# Predict and CV
rsq <- vl_model_eval(observed = dat[set=="test", log2FoldChange], 
                     predicted = predict(model, new= dat[set=="test"]))
pred <- predict(model, new= dat)

# Compute percentage of explained variance
af <- anova(model)
af$PctExp <- af$"Sum Sq"/sum(af$"Sum Sq")*100
af <- as.data.table(af, keep.rownames = "variable")
setorderv(af, "PctExp", 1)
af[, name:= ifelse(grepl("^motif", variable),
                   gsub("^motif_|_L$|_R$", "", variable),
                   variable)]
af[motifs, name:= paste0(i.motif_cluster, 
                         fcase(grepl("median|Residuals", variable), "",
                               grepl("_L$", variable), "_5'",
                               grepl("_R$", variable), "_3'",
                               default= "")), on= "name==variable"]
af[, name:= gsub("median_L", "median_5'", name)]                               
af[, name:= gsub("median_R", "median_3'", name)]        

#---------------------------------------#
# PLOT
#---------------------------------------#
pdf("pdf/draft/Figure_2E.pdf", 
    height = 4.5, 
    width = 5.6)
layout(matrix(c(1,2), nrow= 1), 
       widths = c(1,0.4))
smoothScatter(pred, 
              dat$log2FoldChange,
              xlab= "Linear model prediction (log2)",
              ylab= "Observed (log2)",
              las= 1)
legend("topleft", 
       paste0("CV R²= ", round(rsq$Rsquare, 2)), 
       bty= "n")
abline(0, 1, lty= 2)
par(mar= c(5,-0.1,4,4)+0.1)
bar <- barplot(af$PctExp, 
        beside = T, 
        las= 1,
        xlab= "% explained\nvariance",
        border= NA,
        horiz= T,
        axes= F)
axis(1, at= c(0, 30), labels = c(0, 30))
text(af$PctExp, 
     bar[,1],
     af$name,
     pos= 4,
     xpd= T,
     cex= 0.6, 
     offset= 0.25)
dev.off()