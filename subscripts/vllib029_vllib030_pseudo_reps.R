setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# For vllib029 and vllib030, reads from scren rep 3 were randomly split and added to rep 1 and 2
counts <-  meta[vllib %in% c("vllib029", "vllib030") & cdition=="screen"]
counts <- counts[, {
  reps <- rbindlist(list(`1`= fread(umi_counts[rep=="rep1"]),
                         `2`= fread(umi_counts[rep=="rep2"])),
                    idcol= "rep")
  # Split rep 3 in 2
  .c <- fread(umi_counts[rep=="rep3"])
  .c <- .c[, .(umi_counts= rep(T, umi_counts)), .(L, R, total_counts)]
  set.seed(1)
  idx <- sample(1:2, nrow(.c), replace = T)
  .c <- split(.c, idx)
  .c <- rbindlist(.c, idcol = "rep")
  .c <- .c[, .(umi_counts= .N), .(L, R, rep)]
  # Merge into two pseudo rep
  res <- rbind(reps, .c, fill= T)
  res <- res[, .(total_counts= as.numeric(NA), umi_counts= sum(umi_counts)), .(L, R, rep)]
  # Save
  basen <- gsub("rep1.txt$", "", umi_counts[rep=="rep1"])
  res[, {
    file <- paste0(basen, "pseudoRep", rep, ".txt")
    fwrite(.SD, file)
    print(file)
  }, rep]
  .SD
}, .(vllib, cdition)]