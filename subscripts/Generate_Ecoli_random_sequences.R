library(BSgenome.Ecoli.NCBI.20080805)
library(GenomicRanges)

#length
my.size <- 249

#initialise list to store chromosome sizes
my_chr_size <- list()
for (i in Ecoli@seqinfo@seqnames){
  my_chr_size[[i]] <- length(Ecoli[[i]])
}

#initialise vectors for storing random coordinates
my_random_start  <- vector()
my_random_end    <- vector()

#loop through number of regions
chr_random_region_list <- lapply(names(my_chr_size), function(chr){
  my_max <- my_chr_size[[chr]]-my.size
  Random <- data.frame(Chr=chr,
                       Start=seq(1, my_max, by = my.size),
                       End=seq(my.size, (my_max+my.size), by=my.size)
  )
  
  # remove overlapped regions with peaks (use 2kb window to check the overlap)
  gr_Random <- GRanges(seqnames = Random$Chr,
                       ranges = IRanges(start = Random$Start, end = Random$End))
  
  # confirm there are no sequences with "NN..."
  Sequence <- getSeq(Ecoli, gr_Random)
  if(length(grep("N", Sequence))>0) gr_Random <- gr_Random[-grep("N", Sequence),]
  
  return(gr_Random)
  
})

# join all regions in one GRanges object
all_random_region_list <- do.call("c", chr_random_region_list)

# select "n" sequences
n <- 40
set.seed(1234)
Ecoli_backbones <- all_random_region_list[sample(1:length(all_random_region_list), n)]
# chose strands
strand(Ecoli_backbones) <- "+"

# get sequences
getSeq(Ecoli, Ecoli_backbones)

