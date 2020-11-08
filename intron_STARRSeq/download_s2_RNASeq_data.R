# setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
setwd("/home/vloubiere/public/")
# sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
sapply(list.files("/home/vloubiere/functions/", ".R$", full.names = T), source)
require(R.utils)

#---------------------------------------#
# download_s2_PROSeq_Duarte_2016 
#---------------------------------------#
if(!dir.exists("/home/vloubiere/public/PROSeq_s2/"))
{
  dir.create("/home/vloubiere/public/PROSeq_s2/")
}

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2055nnn/GSM2055132/suppl/GSM2055132_S2_NHS_PRO_seq_norm_minus.bedgraph.gz", 
              destfile = "/home/vloubiere/public/PROSeq_s2/GSM2055132_S2_NHS_PRO_seq_norm_minus.bedgraph.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2055nnn/GSM2055132/suppl/GSM2055132_S2_NHS_PRO_seq_norm_plus.bedgraph.gz", 
              destfile = "/home/vloubiere/public/PROSeq_s2/GSM2055132_S2_NHS_PRO_seq_norm_plus.bedgraph.gz")
lapply(list.files("PROSeq_s2/", pattern = ".gz", full.names = T), function(x){
  .c <- import(gunzip(x, remove= F, skip= T))
  .c$score <- abs(.c$score)
  info <- Seqinfo(genome="dm3")
  seqlengths(.c) <- seqlengths(info)[match(names(seqlengths(.c)), seqnames(info))]
  export.bw(.c, gsub(".bedgraph.gz", ".bw", x))})


#---------------------------------------#
# download_s2_RNASeq_Gan_2019 
#---------------------------------------#
if(!dir.exists("/home/vloubiere/public/RNASeq_S2/"))
{
  dir.create("/home/vloubiere/public/RNASeq_S2/")
}

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM480nnn/GSM480160/suppl/GSM480160_GA0840_Drosophila_S2_RNAseq.bed.gz", 
              destfile = "/home/vloubiere/public/RNASeq_S2/GSM480160_GA0840_Drosophila_S2_RNAseq.bed.gz")
.c <- import(gunzip("/home/vloubiere/public/RNASeq_S2/GSM480160_GA0840_Drosophila_S2_RNAseq.bed.gz", remove= F, skip= T))
.cov <- coverage(.c)
export.bw(.cov, "RNASeq_S2/GSM480160_GA0840_Drosophila_S2_RNAseq.bw")




