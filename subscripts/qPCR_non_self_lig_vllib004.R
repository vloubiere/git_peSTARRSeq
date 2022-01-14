setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

#------------------------------------------------#
# PCR 1, 12-01-2021, too high quantities 
#------------------------------------------------#
# Import
dat <- fread("db/qPCR/20210112_non_selflig_vllib004/20210112_qPCR_vllib004_nonself_ligation -  Quantification Cq Results_0.txt", 
             select = c(2, 8))
id <- fread("db/qPCR/20210112_non_selflig_vllib004/20210112_non_selflig_vllib004_plate_template.txt")
melt <- fread("db/qPCR/20210112_non_selflig_vllib004/20210112_qPCR_vllib004_nonself_ligation -  Melt Curve Derivative Results_SYBR.txt")
melt <- melt(melt[, -1], id.vars = "Temperature")
melt[, variable:= paste0(substr(variable, 1,1), 
                         formatC(as.numeric(substr(variable, 2, nchar(as.character(variable)))), 
                                 width = 2, 
                                 flag = "0"))]

# Master dat
dat <- dat[id, , on= "Well==well"]
dat[, melt_temp:= .(list(melt[.BY, Temperature, on= "variable==Well"])), by= Well]
dat[, melt_value:= .(list(melt[.BY, value, on= "variable==Well"])), by= Well]
dat[, Cc:= rainbow(4)[.GRP], primer]

# Prepare PCR plot obj
pl1 <- dat[, .(mean_cq= mean(Cq), sd= sd(Cq)), .(primer, dna, Cc)]
pl1[, dna_quant:= as.numeric(gsub("ng", "", unlist(tstrsplit(dna, "_",  keep=4))))]
pl1[dna_quant>0, lm:= .(list(lm(mean_cq~log(dna_quant), .SD))), primer]
pl1[dna_quant>0, efficiency:= sapply(lm, function(x) round((-1/(10^summary(x)$coefficients[2,1]-1))*100))]

# Prepare melting curve plot obj
pl2 <- dat[!grepl("0ng$", dna), .(temp= unlist(melt_temp[[1]]), value= unlist(melt_value[[1]])), .(Well, primer, Cc)]

# PLOT
pdf("pdf/design/qPCR_non_self_ligation_vllib004_20210112.pdf", width = 11, height = 4)
layout(matrix(1:3, ncol=3), widths = c(1,0.5,1))
plot(NA, xlim= c(-3.5, 2), ylim= c(5, 20), xlab= "log(DNA)", ylab= "Cq")
pl1[dna_quant>0, 
    {
      points(log(dna_quant), mean_cq, col= Cc, pch= 19)
      abline(lm[[1]], col= Cc, 0.5)
      leg.space <- paste0(rep("\n", .GRP-1), collapse= "")
      legend("topleft", legend = paste(leg.space, "Efficiency=",  efficiency), bty= "n", text.col= Cc[1])
      legend("topright", legend = paste(leg.space, primer), bty= "n", text.col= Cc[1])
      print(legend)
    }, .(primer, Cc, efficiency)]

diff <- dcast(pl1[dna_quant>0], dna_quant~primer, value.var = "mean_cq")
boxplot(2^(diff$non_self_ligation_1-diff$self_ligation_1), 2^(diff$non_self_ligation_2-diff$self_ligation_2), ylab= "FC self-lig/non-self-lig")

plot(NA, xlim = range(pl2$temp), ylim = range(pl2$value), xlab= "temperature", ylab= "Derivative", main= "melting curve")
pl2[, lines(temp, value, col= Cc), .(Well, Cc)]

dev.off()

#------------------------------------------------#
# PCR 2, 13-01-2021, too high quantities 
# illF-SCR2 and SCR2-illR melting curves are bad -> several amplicons
# -> illF-SCR2 amplifies both  illF-SCR2 and illF-SCR2-SCR2 frags?
#------------------------------------------------#
# Import
dat <- fread("db/qPCR/20210113_non_selflig_vllib004/20210113_qPCR_vllib004_nonself_ligation -  Quantification Cq Results_0.txt", 
             select = c(2, 8))
id <- fread("db/qPCR/20210113_non_selflig_vllib004/20210113_non_selflig_vllib004_plate_template.txt")
melt <- fread("db/qPCR/20210113_non_selflig_vllib004/20210113_qPCR_vllib004_nonself_ligation -  Melt Curve Derivative Results_SYBR.txt")
melt <- melt(melt[, -1], id.vars = "Temperature")
melt[, variable:= paste0(substr(variable, 1,1), 
                         formatC(as.numeric(substr(variable, 2, nchar(as.character(variable)))), 
                                 width = 2, 
                                 flag = "0"))]

# Master dat
dat <- dat[id, , on= "Well==well"]
dat[, melt_temp:= .(list(melt[.BY, Temperature, on= "variable==Well"])), by= Well]
dat[, melt_value:= .(list(melt[.BY, value, on= "variable==Well"])), by= Well]
dat[, Cc:= rainbow(4)[.GRP], primer]

# Prepare PCR plot obj
pl1 <- dat[, .(mean_cq= mean(Cq), sd= sd(Cq)), .(primer, dna, Cc)]
pl1[, dna_quant:= 1/as.numeric(unlist(tstrsplit(dna, ":",  keep=2)))]
pl1[dna_quant>0, lm:= .(list(lm(mean_cq~log(dna_quant), .SD))), primer]
pl1[dna_quant>0, efficiency:= sapply(lm, function(x) round((-1/(10^summary(x)$coefficients[2,1]-1))*100))]

# Prepare melting curve plot obj
pl2 <- dat[!grepl("H2O$", dna), .(temp= unlist(melt_temp[[1]]), value= unlist(melt_value[[1]])), .(Well, primer, Cc)]

# PLOT
pdf("pdf/design/qPCR_non_self_ligation_vllib004_20210113.pdf", width = 11, height = 4)
layout(matrix(1:4, ncol=4), widths = c(1,0.5,0.5,1))
plot(NA, xlim= c(-10, -5), ylim= c(8, 30), xlab= "log(DNA)", ylab= "Cq")
pl1[dna_quant>0, 
    {
      points(log(dna_quant), mean_cq, col= Cc, pch= 19)
      abline(lm[[1]], col= Cc, 0.5)
      leg.space <- paste0(rep("\n", .GRP-1), collapse= "")
      legend("topleft", legend = paste(leg.space, "Efficiency=",  efficiency), bty= "n", text.col= Cc[1])
      legend("topright", legend = paste(leg.space, primer), bty= "n", text.col= Cc[1])
      print(legend)
    }, .(primer, Cc, efficiency)]

diff <- dcast(pl1[dna_quant>0], dna_quant~primer, value.var = "mean_cq")
boxplot(diff$non_self_ligation_1-diff$self_ligation_1, diff$self_ligation_1-diff$illF_SCR2, ylab= "FC self-lig/non-self-lig (log2)")
boxplot(2^(diff$non_self_ligation_1-diff$self_ligation_1), 2^(diff$self_ligation_1-diff$illF_SCR2), ylab= "FC self-lig/non-self-lig")

plot(NA, xlim = range(pl2$temp), ylim = range(pl2$value), xlab= "temperature", ylab= "Derivative", main= "melting curve")
pl2[, lines(temp, value, col= Cc), .(Well, Cc)]

dev.off()
