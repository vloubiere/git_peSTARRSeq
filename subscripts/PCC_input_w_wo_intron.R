setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

# Compute merged counts for input with intron
libs <- as.data.table(read_excel("/groups/stark/vloubiere/exp_data/vl_libraries.xlsx"))
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- colnames(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]
meta[, output_prefix:= paste0("/", my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]
meta <- meta[grepl("vllib017.*input|vllib019.*input|vllib021.*input", my_ID)]
meta[suffix=="with-intron"]
meta[, {
  if(.N>0)
  {
    counts <- paste0("db/merged_counts_with_intron/", vllib, "_", cdition, "_rep", rep, "_merged.txt")
    if(file.exists(counts))
      print(paste0(counts, " -->> ALREADY EXISTS")) else
      {
        print(paste0("Start ", counts))
        files <- paste0("db/umi_counts/", output_prefix, ".txt")
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
        .ex <- libs[lib_ID %in% c(vllib, Spike_in)]
        .ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]], 
                       R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
        .ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]], 
                        R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
        .ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
                                              paste0("_", R_pattern, "_"))]
        .ex$L <- NULL
        .ex$R <- NULL
        .ex[, Spike:= ifelse(lib_ID==Spike_in & !is.na(Spike_in), T, F)]
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
        # SAVE
        fwrite(.c, counts)
        fwrite(summary, gsub(".txt$", ".summary.txt", counts))
        print(paste0(counts, " -->> DONE"))
      }
  }
}, .(cdition, CP, rep, vllib, Spike_in)]

# Compute inputs with and without introns
dat <- data.table(file= list.files(c("db/merged_counts/", "db/merged_counts_with_intron/"), "vllib017.*input.*merged.txt|vllib019.*input.*merged.txt|vllib021.*input.*merged.txt", full.names = T))
dat[, cdition:= tstrsplit(basename(file), "_", keep= 1), file]
dat[grepl("with_intron", file), intron:= "intron"]
dat[!grepl("with_intron", file), intron:= "no_intron"]
dat <- dat[, fread(file), (dat)]

pdf("pdf/alignment/PCC_with_without_introns.pdf", width = 5)
par(pty= "s", las= 1)
dat[, {
  .c <- dcast(.SD,
               L+R~intron, 
               value.var = "umi_counts")
  x <- log2(.c$intron)
  y <- log2(.c$no_intron)
  smoothScatter(x, 
                y, 
                xlab= "With intron (log2 counts)",
                ylab= "Without intron (log2 counts)",
                main= cdition)
  legend("topleft", 
         legend = paste0("PCC= ", round(cor.test(x, y)$estimate, 2)),
         bty= "n")
  print("")
}, cdition]
dev.off()
