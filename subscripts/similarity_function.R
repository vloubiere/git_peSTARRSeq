require(GenomicRanges)
require(BSgenome)
require(stringdist)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(BSgenome.Ecoli.NCBI.20080805)

similar <- function(file1, file2= NULL, genome1= "dm3", genome2= "dm3")
{
  if(is.null(file2))
  {
    print("Only one file provided: reports idx of rows for which another row")
    print("has similar sequence (<2 diff)")
    file2 <- file1
  }
  
  ### Specify genome1
  if(!(genome1 %in% c("dm3", "Ecoli")) | !(genome2 %in% c("dm3", "Ecoli")))
  {
    stop("unknown genome specified")
  }
  
  if(genome1 == "dm3"){
    genome1 <- BSgenome.Dmelanogaster.UCSC.dm3
  }else if(genome1 == "Ecoli"){
    genome1 <- BSgenome.Ecoli.NCBI.20080805
  }
  
  ### Import sequences
  plus1 <- file1
  strand(plus1) <- "+"
  fw1 <- as.character(getSeq(genome1, plus1))
  start_fw1 <- substring(fw1, 1, 10)
  
  minus1 <- file1
  strand(minus1) <- "-"
  rev1 <- as.character(getSeq(genome1, minus1))
  start_rev1 <- substring(rev1, 1, 10)
  
  if(identical(file1, file2))
  {
    similar <- mclapply(seq(start_fw1), function(i) 
    {
      dist_1 <- stringdist(start_fw1[i], start_fw1)
      dist_2 <- stringdist(start_rev1[i], start_rev1)
      same <- any(dist_1[-i] < 2 | dist_2[-i] < 2)
      return(same)
    })
  }else
  {
    print("Two files provided: reports idx of rows of file 1 for which ")
    print("a row from the file 2 has a similar sequence (<2 diff)")
    
    ### Specify genome2
    if(genome2 == "dm3"){
      genome2 <- BSgenome.Dmelanogaster.UCSC.dm3
    }else if(genome2 == "Ecoli"){
      genome2 <- BSgenome.Ecoli.NCBI.20080805
    }
    
    ### Import sequences
    plus2 <- file2
    strand(plus2) <- "+"
    fw2 <- as.character(getSeq(genome2, plus2))
    start_fw2 <- substring(fw2, 1, 10)
    
    minus2 <- file2
    strand(minus2) <- "-"
    rev2 <- as.character(getSeq(genome2, minus2))
    start_rev2 <- substring(rev2, 1, 10)
    
    similar <- mclapply(seq(start_fw1), function(i) 
    {
      dist_1 <- stringdist(start_fw1[i], start_fw2)
      dist_2 <- stringdist(start_rev1[i], start_rev2)
      same <- any(dist_1 < 2 | dist_2 < 2)
      return(same)
    })
  }
  return(unlist(similar))
}