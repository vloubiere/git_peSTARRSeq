setwd("/groups/stark/vloubiere/projects/0003_luciferase_HAM1_SUP1_SCR2/")
require(data.table)
require(plater)

# import plate schemes
scheme <- data.table(scheme= list.files("data/", "scheme", full.names = T))
scheme <- scheme[, read_plate(scheme), scheme]
scheme[, sample:= substr(basename(scheme), 1, 13), scheme]
scheme[, cdition:= rownames]
scheme$rownames <- NULL
scheme[, cdition:= gsub("_r1?|_r2?", "", cdition)]
scheme[cdition %in% c("DSCP_ZFH1", "Ubi_GFP"), cdition := gsub("_", "", cdition)]
scheme[, c("enh_L", "enh_R"):= tstrsplit(cdition, "_")]

# Import data
lum <- data.table(file= list.files("data/", "lum", full.names = T))
lum <- lum[, read_plate(file), file]
lum[, variable:= ifelse(grepl("lum1.csv?", file), "luciferase", "renilla")]
lum[, sample:= gsub("_lum1.csv|_lum2.csv", "", basename(file)), file]
lum <- dcast(lum, sample+Wells~variable, value.var = "values")

# Create object
dat <- scheme[lum, , on= c("sample", "Wells")]
dat[, idx:= .SD[,.I], .(enh_L, enh_R)]
dat[, activity:= luciferase/renilla]

# Annotate
ord <- unique(dat[, .(enh_L, enh_R)])
ord[, order:= c(13,5,15,6,7,9,17,4,19,1,10,3,8,2,11)]
ord[, name:= c("hamlet x 2 obs.", "control x 2 obs.", "hamlet x sup-add obs.", "hamlet", "hamlet x control", "sup-add",
               "sup-add x hamlet obs.", "control", "sup-add x2 obs.", "ZFH1", "sup-add x control", "p002", "control x hamlet", "UbiGFP", "control x sup-add")]
dat <- ord[dat, , on= c("enh_L", "enh_R")]

# Add additive scores
exp_L <- na.omit(dat[, .SD[enh_R=="SCR2", .(exp_L= activity)], .(enh_L, sample, idx)])
exp_R <- na.omit(dat[, .SD[enh_L=="SCR2", .(exp_R= activity)], .(enh_R, sample, idx)])
exp <- merge(exp_L, exp_R, by= c("sample", "idx"), allow.cartesian = T)
exp <- exp[enh_L != "SCR2" & enh_R != "SCR2"]

ord <- unique(exp[, .(enh_L, enh_R)])
ord[, order:= c(12, 14, 16, 18)]
ord[, name := c("hamlet x 2 exp.", "hamlet x sup-add. exp.", "sup-add. x hamlet exp.", "sup-add x 2 exp.")]
exp <- ord[exp, , on= c("enh_L", "enh_R")]
exp[, activity := exp_L+exp_R]

dat <- rbind(dat, exp, fill= T)
setkeyv(dat, "order")

# Compute informative values and save
dat[, mean_activity:= mean(activity), order]
dat[, sd_activity:= sd(activity), order]

saveRDS(dat, "Rdata/luciferase_validation.rds")




