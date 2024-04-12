setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import counts ----
dat <- data.table(file= list.files("db/dds/", ".dds", full.names = T))
dat[, library:= gsub(".dds$", "", basename(file))]
dat <- dat[, {
  as.data.table(DESeq2::counts(readRDS(file)))
}, library]

# Format ----
dat <- melt(dat, id.vars = "library")
dat[, value:= log2(value)]
dat[, c("cdition", "rep"):= tstrsplit(variable, "_"), variable]
dat <- dcast(dat, 
             library+cdition~rep, 
             value.var = "value", 
             fun.aggregate = list)

# Compute PCC and lm ----
dat[, c("PCC", "lm", "rsq"):= {
  x <- unlist(`1`)
  y <- unlist(`2`)
  .lm <- lm(y~x)
  .(round(cor.test(x, y)$estimate, 2), 
    .(.lm),
    round(summary(.lm)$r.squared, 2))
}, .(library, cdition)]

pdf("pdf/draft/PCC_peSTARRSeq_replicates.pdf", 
    height = 3, 
    width = 6)
par(mfrow= c(1,2),
    mai= rep(.9, 4), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    pty= "s",
    lend= 2,
    font.main= 1)
dat[, {
  smoothScatter(unlist(`1`),
                unlist(`2`),
                xlab= "Replicate 1",
                ylab= "Replicate 2",
                main= paste(cdition, library),
                col= adjustcolor(blues9[9], .3),
                xaxt= "n")
  axis(1, padj = -1.25)
  vl_plot_coeff(value= PCC,
                type= "pcc",
                cex= 7/12)
  print(.GRP)
}, .(library, cdition)]
dev.off()