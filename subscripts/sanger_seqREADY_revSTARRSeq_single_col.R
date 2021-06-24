setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(sangerseqR)

# Extract
dat <- data.table(file= list.files("db/sanger_sequencing/revSTARRSeq_single_colonies/", full.names = T))
dat[, cdition:= tstrsplit(basename(file), "_", keep= 2)]
dat[, abdir:= ifelse(grepl("_T3_", file), "F", "R")]
dat[, abseq:= as.character(primarySeq(readsangerseq(file))), (dat)]

# Find illumina primers
seq <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", 
             key= "name")
dat[, c("start_F", "end_F"):= 
      as.data.table(matchPattern(DNAString(seq["illumina_F", sequence]), subject = abseq)@ranges), abseq]
dat[, c("start_F", "end_F"):= 
      as.data.table(matchPattern(reverseComplement(DNAString(seq["illumina_F", sequence])), subject = abseq)@ranges), abseq]
dat[, c("start_R", "end_R"):= 
      as.data.table(matchPattern(DNAString(seq["illumina_R", sequence]), subject = abseq)@ranges), abseq]
dat[, c("start_R", "end_R"):=
      as.data.table(matchPattern(reverseComplement(DNAString(seq["illumina_R", sequence])), subject = abseq)@ranges), abseq]
dat[, dir:= ifelse(end_F<start_R, "F", "R")]

# Compute sequence pattern 
enh <- readDNAStringSet("/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_revSTARRSeq/revSTARRSeq_sequences.fa")
dat[dir=="F", pL:= as.character(reverseComplement(DNAString(substr(abseq, end_F+5, end_F+35)))), abseq]
dat[dir=="F", pR:= as.character(reverseComplement(DNAString(substr(abseq, start_R-35, start_R-5)))), abseq]
dat[dir=="R", pL:= as.character(DNAString(substr(abseq, start_F-35, start_F-5))), abseq]
dat[dir=="R", pR:= as.character(DNAString(substr(abseq, end_R+5, end_R+35))), abseq]
# Identify enhancers
dat[!is.na(pL), enhF:= names(enh)[which(lengths(vmatchPattern(DNAString(pL), enh))==1)], pL]
dat[!is.na(pR), enhR:= names(enh)[which(lengths(vmatchPattern(DNAString(pR), enh))==1)], pR]
# Compute refseq
dat[!is.na(enhF) & !is.na(enhR), refseq:= paste0(seq["illumina_F", sequence], 
                                                 as.character(reverseComplement(enh[match(enhF, names(enh))])),
                                                 "CCGCGCGCGCTCATCAATGTATCTTATCATGTCTGCTCGAAGCGGCCGGCCGAATTCG", # Ring linker sequence
                                                 as.character(reverseComplement(enh[match(enhR, names(enh))])),
                                                 seq["illumina_R", sequence]), .(enhF, enhR)]

# Features table
feat <- rbind(seq[c("illumina_F",
                    "Flink_-2"), .(name, sequence)],
              seq[c("illumina_R",
                    "R1link+0",
                    "R2link+0",
                    "R3link-3"), .(name, sequence= as.character(reverseComplement(DNAStringSet(sequence))))])
feat <- rbind(feat, data.table(name= "spacer", sequence= "CCGCGCGCGCTCATCAATGTATCTTATCATGTCTGCTCGAAGCGGCCGGCCGAATTCG"))
feat[, Cc:= c("cornflowerblue",
              "tomato",
              "pink",
              "orange",
              "gold",
              "purple",
              "green")]

# Plot result
pdf("pdf/sanger_sequencing/sanger_single_colonies_STARRSeq_vllib006.pdf", height = 2.5, width = 10)
par(xaxs= "i", yaxs= "i", mar= c(4,12,4,12))
dat[!is.na(refseq), {
  .c <- c(.SD$refseq[1])
  names <- c("refseq")
  
  # Create sequence vector
  if(nrow(.SD[dir=="F"])>0)
  {
    .c <- c(.c, .SD[dir=="F", abseq])
    names <- c(names, gsub(".ab1", "", basename(.SD[dir=="F", file])))
  }
  if(nrow(.SD[dir=="R"])>0)
  {
    .c <- c(.c, as.character(reverseComplement(DNAString(.SD[dir=="R", abseq]))))
    names <- c(names, gsub(".ab1", "", basename(.SD[dir=="R", file])))
  }
  
  # Align
  align <- msa(DNAStringSet(.c), order = "input", method= "ClustalOmega")
  align <- rbindlist(lapply(as.character(align@unmasked), function(x) strsplit(x, "")), idcol = T)
  align[, idx:= .SD[, .I], .id]
  align <- dcast(align, idx~.id, value.var = "V1")
  colnames(align)[-1] <- c("refseq", gsub(".ab1", "", basename(file)))
  
  # make matrix
  mat <- as.matrix(align, 1)
  mat[mat=="-"] <- "white"
  for(i in 2:ncol(mat))
  {
    mat[mat[, i] != "white", i] <- ifelse(mat[mat[, i] != "white", i]==mat[mat[, i] != "white", 1], "blue", "red")
  }
  mat[mat[, 1]!="white",1] <- apply(mat[mat[, 1]!="white",], 1, function(x) ifelse(any(x=="blue"), "blue", "red"))
  
  plot.new()
  rasterImage(t(mat), 0, 0, 1, 1)
  mtext(paste0(.SD[dir=="F", enhF], " x ", .SD[dir=="F", enhR]), line = 1)
  mtext(paste0(.SD[dir=="R", enhF], " x ", .SD[dir=="R", enhR]), side = 1, line = 1)
  abline(h= seq(0, 1, length.out= ncol(mat)+1))
  y <- seq(1, 0, length.out= ncol(mat)*2+1)
  y <- y[(seq(y) %% 2)==0]
  axis(side = 2, at = y, labels = colnames(mat), las= 1)
  
  # Add forward features
  feat[, {
    .c <- as.data.table(matchPattern(DNAString(sequence), 
                                     subject = paste0(align$refseq, collapse= ""))@ranges)
    if(nrow(.c)>0)
    {
      arrows(.c$start/nrow(mat), 1.1, 
             .c$end/nrow(mat), 1.1, 
             lwd= 5, 
             col = Cc, 
             xpd= T, 
             lend= 2, 
             length = 0.025)
    }
  }, (feat)]
  # Add reverse features
  feat[, {
    .c <- as.data.table(matchPattern(reverseComplement(DNAString(sequence)), 
                                     subject = paste0(align$refseq, collapse= ""))@ranges)
    if(nrow(.c)>0)
    {
      arrows(.c$end/nrow(mat), 1.1, 
             .c$start/nrow(mat), 1.1, 
             lwd= 5, 
             col = Cc, 
             xpd= T, 
             lend= 2, 
             length = 0.025)
    }
  }, (feat)]
  print("DONE")
}, cdition]
dev.off()
