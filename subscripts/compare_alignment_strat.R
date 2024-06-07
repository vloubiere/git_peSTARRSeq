setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

old <- data.table(old= list.files("db/umi_counts/", full.names = T))
old[, name:= basename(old)]
new <- data.table(new= list.files("db/umi_counts/test/", full.names = T))
new[, name:= basename(new)]

dat <- merge(old, new)
dat[, c("total_old", "umi_old"):= fread(old)[, .(sum(umi_counts), sum(total_counts))], old]
dat[, c("total_new", "umi_new"):= fread(new)[, .(sum(umi_counts), sum(total_counts))], new]

pl <- melt(dat, id.vars = "name", measure.vars = c("umi_old", "umi_new"))
pl[, perc:= round(value/value[variable=="umi_old"], 3), name]
setorderv(pl, c("name", "variable"))

par(las= 2,
    mai= c(3,1,0.9,0.9),
    cex.lab= 9/12,
    cex.axis= 8/12)
bar <- vl_barplot(pl$value,
                  names.arg= gsub(".txt$", "", pl$name))
text(bar, pl$value, ifelse(pl$perc==1, NA, pl$perc), pos= 3, xpd= T)
