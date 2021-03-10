# Map fine alignment of sub-libs ####
refine <- data.table(file= list.files("db/bam/", ".txt$", full.names = T))
refine[, lib:= tstrsplit(basename(file), "_", keep= 2)]
refine <- refine[!grepl("stats", file)]
refine <- refine[, fread(file, nrows = 10000, fill= T), refine]
refine <- refine[V3!="*" & V5>30]
refine[grep("_A_", V3), class:= "R1"]
refine[grep("_B_", V3), class:= "R2"]
refine[grep("_C_", V3), class:= "R3"]
refine[grep("_HAM1_", V3), class:= "HAM1"]
refine[grep("_SCR2_", V3), class:= "SCR2"]
refine[grep("_SUP1_", V3), class:= "SUP1"]
heat <- dcast(refine, lib+class~V4)
gaps <- cumsum(heat[, .N, lib][, N])
heat[, lib:= paste0(lib, "__", class)]
heat <- heat[, !"class"]
mat <- as.matrix(heat, 1)
pheatmap(log2(mat+1), filename = "pdf/sub_libraries_fine_alignment.pdf",
         cluster_rows = F, cluster_cols = F, width = 20, height = 3, gaps_row = gaps)
# read 1 libvl002= 17
# read 2 HAM1 libvl002= 303
# read 2 SCR2 libvl002= 257
# read 2 SUP1 libvl002= 298
# read 2 R1   libvl002= 260
# read 2 R2   libvl002= 260
# read 2 R3   libvl002= 264

# read 1 R1/2/3 libvl012/13/14= 16-19
# read 1 SCR2   libvl012/13/14= 17
# read 2 SCR2   libvl012/13/14= 257
# read 2 R1     libvl0012/13/14= 260-263
# read 2 R2     libvl0012/13/14= 260-263
# read 2 R3     libvl0012/13/14= 261-264
#----------------------------------------------------------------------------####
# retrieve pairs and compute counts ####
pa <- data.table(file= list.files("db/bam", ".txt$", full.names = T))
pa <- pa[!grepl("stats", file)]
pa[, cdition:= gsub("_1.txt|_2.txt", "", basename(file))]
pa[, read:= ifelse(grepl("_1.txt$", file), "read_1", "read_2")]
pa <- dcast(pa, cdition~read, value.var = "file")
pa[, cdition:= gsub("^[^_]*_(.*)", "\\1", cdition), cdition]
pa[, umi_counts:= paste0("db/count/", cdition, "_umi_counts.txt"), cdition]
pa[, lib:= tstrsplit(cdition, "_", keep= 1), cdition]
pa[, {
  .L <- rbindlist(lapply(read_1, function(x) fread(x, fill = T, select = c(1,2,3,4,5))))
  .L <- .L[V2==0 & V5>30]
  .R <- rbindlist(lapply(read_2, function(x) fread(x, fill = T, select = c(1,2,3,4,5))))
  .R <- .R[V2==16 & V5>30]
  if(lib=="vllib002")
  {
    .L <- .L[between(V4, 16, 18, incbounds = T)]
    .R <- .R[(grepl("_HAM1_", V3) & between(V4, 302, 304, incbounds = T))
             | (grepl("_SCR2_", V3) & between(V4, 256, 258, incbounds = T))
             | (grepl("_SUP1_", V3) & between(V4, 297, 299, incbounds = T))
             | (grepl("_A_|_B_", V3) & between(V4, 259, 261, incbounds = T))
             | (grepl("_C_", V3) & between(V4, 263, 265, incbounds = T))]
  }else
  {
    .L <- .L[(grepl("_SCR2_", V3) & between(V4, 16, 18, incbounds = T)) 
             | (!grepl("_SCR2_", V3) & between(V4, 15, 20, incbounds = T))]
    .R <- .R[(grepl("_SCR2_", V3) & between(V4, 256, 258, incbounds = T))
             | (grepl("_A_|_B_", V3) & between(V4, 259, 264, incbounds = T))
             | (grepl("_C_", V3) & between(V4, 260, 265, incbounds = T))]
  }
  .pair <- merge(.L, .R, by.x= "V1", by.y= "V1", suffixes= c("_L", "_R"))
  .pair[, UMI:= gsub("^[^_]*_(.*)", "\\1", V1)]
  .pair <- unique(.pair[, .(L= V3_L, R= V3_R, UMI)])
  .pair <- .pair[, .(count= .N), .(L, R)]
  # Check if pair suppositely exists and if spike-in / switching
  if(lib=="vllib002")
  {
    .pair[, exist:= (grepl("^ts", L) & L==R) | (!grepl("^ts", L) & !grepl("^ts", R))]
    .pair[grepl("^ts", L) | grepl("^ts", R), spike:= ifelse(L==R, "ok", "switched")]
  }
  if(lib=="vllib012")
  {
    .pair[, exist:= (grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))]
    .pair[grepl("_SCR2_|_C_", L) | grepl("_SCR2_|_C_", R), 
          spike:= ifelse(((grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))), "ok", "switched")]
  }
  if(lib=="vllib013")
  {
    .pair[, exist:= (grepl("_A_", L) & grepl("_A_", R)) | (grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))]
    .pair[grepl("_SCR2_|_C_", L) | grepl("_SCR2_|_C_", R), 
          spike:= ifelse(((grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))), "ok", "switched")]
  }
  if(lib=="vllib014")
  {
    .pair[, exist:= (grepl("_B_", L) & grepl("_B_", R)) | (grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))]
    .pair[grepl("_SCR2_|_C_", L) | grepl("_SCR2_|_C_", R), 
          spike:= ifelse(((grepl("_C_", L) & grepl("_SCR2_", R)) | (grepl("_SCR2_", L) & grepl("_C_", R))), "ok", "switched")]
  }
  fwrite(.pair, umi_counts, col.names = T, row.names = F, sep= "\t", quote= F, na = NA)
  print(paste0(umi_counts, " DONE!"))
}, .(umi_counts, lib)]

