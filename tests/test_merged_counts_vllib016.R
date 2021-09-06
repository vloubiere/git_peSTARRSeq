setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

files <- list.files("db/umi_counts/", "vllib016_screen.*rep2", full.names= T)
files <- files[!grepl("summary", files)]
.c <- lapply(files, fread)
.c <- unique(rbindlist(.c))
# Remove problematic UMIs
.c <- .c[!grepl("GGGGGGGGG", UMI)]
# Advanced UMI collapsing
for(i in 0:9) 
  .c <- .c[, .(UMI= UMI[1]), .(L, R, sub(paste0("(.{", i, "})."), "\\1", UMI))]
# UMI collapsing
.c <- .c[, .(umi_counts= .N), .(L, R)]
# Compute patterns to identify library pairs
libs <- fread("/groups/stark/vloubiere/exp_data/vl_libraries.txt")
.ex <- libs[lib_ID == "vllib016"]
.ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]], 
               R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
.ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]], 
                R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
.ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
                                      paste0("_", R_pattern, "_"))]
.ex$L <- NULL
.ex$R <- NULL
.ex[, Spike:= F]
# Check if pair exists and is spike in
.ex[, {
  .c[grepl(L_pattern, L) & grepl(R_pattern, R), type:= ifelse(Spike, "spike-in", "pair")]
}, .(L_pattern, R_pattern, Spike)]
.c[is.na(type), type:= "switched"]
# Generate read summary
sum_files <- data.table(pattern= gsub(".txt", "", basename(files)))
sum_files[, file:= list.files("db/sam/", paste0(pattern, ".*.summary$"), full.names = T), pattern]
sum_files[, mapped:= fread(file)[V1=="Mapped_fragments", V2], file]
summary <- sum_files[, .(mapped= sum(mapped), collapsed= sum(.c$umi_counts))]

real <- fread("db/merged_counts/vllib016_short_spacer_peSTARRSeq_screen_RpS12_rep2_merged.txt")
barplot(matrix(c(.c[, sum(umi_counts)], 
                 .c[L!=R, sum(umi_counts)],
                 real[, sum(umi_counts)],
                 real[L!=R, sum(umi_counts)]), 
               ncol= 2), beside= T)

