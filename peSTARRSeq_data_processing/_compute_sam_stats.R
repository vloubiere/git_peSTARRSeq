setwd("/groups/stark/vloubiere/data/pe_STARRSeq/")
source("/groups/stark/vloubiere/scripts/R_functions/R_shell_cmd_wrap.R")
require(data.table)

reads <- data.table(file= list.files("/groups/stark/vloubiere/data/pe_STARRSeq/sam", full.names = T))
reads[, sam_stats:= paste0("/groups/stark/vloubiere/data/pe_STARRSeq/sam_stats/", gsub(".sam", ".txt", basename(file)))]
reads[, bsub(paste0("module load samtools/1.9-foss-2018b; /software/2020/software/samtools/1.9-foss-2018b/bin/samtools stats ", file, " | grep ^SN | cut -f 2- > ", sam_stats)), .(file, sam_stats)]
