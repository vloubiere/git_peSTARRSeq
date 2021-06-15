setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(BSgenome.Dmelanogaster.UCSC.dm3)

int <- fread("Rdata/selected_introns.txt")

# add 20nt on each side
int[,c("start", "end"):= .(start-10, end+10)]

# extract sequence
seq <- as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(int)))

# Fusion PCR overhang/spacer/MaubI site
L <- rep(NA, length(seq))
i <- 20
while(any(is.na(L)))
{
  pattern <- paste0("(^[A-Z]{", i, "}).*")
  seqs <- gsub(pattern, "\\1", seq)
  bool <- sapply(seqs, function(x) vl_oligo_Tm(x)[[1]]>58)
  bool[!is.na(L)] <- F
  L[bool] <- seqs[bool]
  i <- i+1
}

R <- rep(NA, length(seq))
i <- 2
while(any(is.na(R)))
{
  pattern <- paste0(".*([A-Z]{", i, "})$")
  seqs <- sapply(gsub(pattern, "\\1", seq), function(x) as.character(reverseComplement(DNAString(x))))
  bool <- sapply(seqs, function(x) vl_oligo_Tm(x)[[1]]>58)
  bool[!is.na(R)] <- F
  R[bool] <- seqs[bool]
  i <- i+1
}

int[, c("primer_L", "primer_R"):= .(L, R)]
int <- melt(int, id.vars = c("FBgn", "strand", "size"), measure.vars = patterns("primer"))
int[, name:= paste0(c(FBgn, size, strand), collapse= "_"), (int)]
int[variable=="primer_L", name:= paste0(name, "_L"), name]
int[variable=="primer_R", name:= paste0(name, "_R"), name]
