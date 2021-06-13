# Spike-ins percentage ####
spike <- data.table(file= list.files("db/count/", "_umi_counts.txt", full.names = T))
spike[, lib:= tstrsplit(basename(file), "_", keep=1), file]
spike <- spike[lib!="vllib012"]
spike <- spike[, fread(file), spike]
spike <- spike[, .(perc= round(sum(.SD[spike=="ok", count])/sum(.SD[(exist), count])*100, 2)), lib]

pdf("pdf/spike_ins_percentage_barplot.pdf", width = 2.5, height = 3)
par(mar= c(4,4,1,1), las= 2)
barplot(spike$perc, names.arg = spike$lib, ylim = c(0, 10), ylab= "% Spike-ins")
dev.off()
#####
