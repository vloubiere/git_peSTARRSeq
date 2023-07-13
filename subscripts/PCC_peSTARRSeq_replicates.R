setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import counts
dat <- data.table(file= list.files("db/dds/", full.names = T))
dat[, library:= gsub(".dds$", "", basename(file))]
dat <- dat[, {
  as.data.table(DESeq2::counts(readRDS(file)))
}, library]
# Format
dat <- melt(dat, id.vars = "library")
dat[, value:= log2(value)]
dat[, c("cdition", "rep"):= tstrsplit(variable, "_"), variable]
dat <- dcast(dat, 
             library+cdition~rep, 
             value.var = "value", 
             fun.aggregate = list)
# Compute PCC and lm
dat[, c("PCC", "lm", "rsq"):= {
  x <- unlist(rep1)
  y <- unlist(rep2)
  .lm <- lm(y~x)
  .(round(cor.test(x, y)$estimate, 2), 
    .(.lm),
    round(summary(.lm)$r.squared, 2))
}, .(library, cdition)]

pdf("pdf/draft/PCC_peSTARRSeq_replicates.pdf", 
    height = 5.5, 
    width = 2.5*length(unique(dat$library)))
vl_par(mfcol= c(2, length(unique(dat$library))),
       cex= 1,
       mar= c(3,3,3,1),
       bty= "n",
       mgp= c(1.5,0.5,0))
dat[, {
  smoothScatter(unlist(rep1),
                unlist(rep2),
                xlab= "Replicate 1",
                ylab= "Replicate 2",
                main= paste(cdition, library))
  legend("topleft",
         c(paste("PCC=", PCC)),
           # paste("R2=", rsq)),
         bty= "n",
         cex= 0.8)
  # abline(lm[[1]], lty= 2)
  print(.GRP)
}, .(library, cdition)]
dev.off()

file.show("pdf/draft/PCC_peSTARRSeq_replicates.pdf")