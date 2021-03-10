# alignment ####
aln <- data.table(file= list.files("db/fastq/", "fq.gz", full.names = T))
aln[, bam:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq_2/db/bam/", gsub(".fq.gz", ".bam", basename(file)))]
aln[, lib:= gsub("^[^_]*_(.*)", "\\1", basename(bam))]
aln[, lib:= tstrsplit(lib, "_", keep= 1)]
aln[, idx:= grep(paste0(lib, ".*index$"), list.dirs("/groups/stark/vloubiere/genomes/Rsubread_genome/dm3"), value = T), lib]
aln[, idx:= unique(gsub("_idx.*$", "_idx", list.files(idx, full.names = T))), idx]
aln[, txt:= gsub(".bam$", ".txt", bam)]
aln[, {
  print(paste("START", file))
  stats <- capture.output(align(index = idx,
                                readfile1 = file, 
                                type= "dna", 
                                output_file = bam, 
                                maxMismatches= 1, 
                                unique= T, 
                                nTrim5 = 15, 
                                nthreads= getDTthreads()-1))
  writeLines(stats, gsub(".bam$", "_stats.txt", bam))
  cmd <- paste("module load build-env/2020; module load samtools/1.9-foss-2018b; samtools view -@", 
               getDTthreads()-1, bam, ">", txt)
  system(cmd)
  print(paste(bam, "DONE!"))
}, .(file, bam, txt, idx)]
