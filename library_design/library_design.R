setwd("~/Dropbox (VBC)/R_db/dm3/0002_library_design/")
pack <- c("knitr", "colorRamps",  "plyr", "dplyr", "rtracklayer", "GenomicRanges", "data.table", "Hmisc",
          "yarrr", "Biostrings", "BSgenome", "stringdist", "BSgenome.Dmelanogaster.UCSC.dm3")
for(p in pack){require(p, character.only = T, quietly = T)}
source("../../functions/my_plots.R")
source("scripts/function_compute_similarity.R")

################################################################################
#### Use Bernardo's twist library to design paire-element STARR-Seq library ####
################################################################################
# import twist library Bernardo
twist <-  GRanges(readRDS("../Bernardo_fanny/Rdata/twist_enhancers_strongest_strand.rds"))

# Low input cutoff
input <- mcols(twist)[, grep("input", colnames(mcols(twist)))]
par(mfrow= c(3,3))
my_empty_plot(ylim=c(0,1), xlim=c(0, 10000))
invisible(lapply(seq(ncol(input)), function(i){lines(ecdf(input[,i]), col= i)}))
legend("topleft", col = seq(ncol(input)), legend= colnames(input), bty= "n", lty= 1)
cutoff <- 1500
abline(v= cutoff)
twist_clean <- twist[apply(input, 1, function(x) all(x>cutoff))]

# Select only canonical chromosomes and put everything on positive strand
seqlevels(twist_clean, pruning.mode="coarse") <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")
strand(twist_clean) <- "+"

# Save clean twist
save_clean <- twist_clean
names <- paste0(twist_clean$enhancer_group_detail, "_", round(twist_clean$dev_log2FoldChange, 2))
mcols(save_clean) <- data.frame(name= names)
export(save_clean, "Rdata/clean_twist.bed")

#####################################################################
### Compute N dev enhancers for all genes in order to select loci of interests
# Import genes
genes <- get(load("../../../genomes/dm3/dmel-all-no-analysis-r5.57_simple.Rdata"))
genes <- genes[genes$type == "gene"]
# Extract dev enhancers from twist library
dev_enhancers <- twist_clean[twist_clean$enhancer_group_detail %in% c("dev", "shared")]
seqlevels(genes, pruning.mode="coarse") <- seqlevels(dev_enhancers) 
# Count enhancers ovelrapping each gene
N_overlapping_enhancers <- countOverlaps(genes, dev_enhancers, ignore.strand= T)
genes$N_overlapping_enhancers <- N_overlapping_enhancers

# Extract non-overlapping remaining enhancers and assign non-overlapping enhancers to closest TSS
unassigned_dev_enhancers <- subsetByOverlaps(dev_enhancers, genes, invert = TRUE, ignore.strand = TRUE)
gene_TSSs <- promoters(genes, upstream = 0, downstream = 1)
hits <- nearest(unassigned_dev_enhancers, gene_TSSs, ignore.strand = TRUE, select = "all")
counts <- as.data.table(subjectHits(hits))
counts <- counts[,.N, V1]
genes$N_closest_non_overlapping_enhancers <- 0
genes$N_closest_non_overlapping_enhancers[counts$V1] <- counts$N

# Promoters loci of interest
loci <- genes[genes$Name %in% c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo")]
POI <- split(promoters(loci, upstream = 0, downstream = 1), seq(loci))

#####################################################################
### Select 600 representative dev enhancers with a wide range of strenght + expected to collaborate endogenously
dev_enhancers <- twist_clean[twist_clean$enhancer_group_detail %in% c("dev", "shared")]

# Select ~274 dCP which are found around a set of loci of interest
dev_sel <- mclapply(POI, function(x)
{
  dist <- as.data.table(distanceToNearest(dev_enhancers, x, ignore.strand= T))
  setorder(dist, distance)
  sel <- dev_enhancers[dist[1:30, queryHits]]
  return(sel)
})
dev_sel <- unique(Reduce(c, dev_sel))
# remove similar sequences
dev_sel <- dev_sel[!similar(dev_sel)]

# Use the remaining 326 slots to equilibrate the library with a wide range of strenght 
# Lowly active TWIST enhancers slightly disadvanteged !!!
dev_breaks <- c(-Inf,4,6,8,10)
dev_quant <- cut(dev_enhancers$dev_log2FoldChange, dev_breaks, include.lowest = T, labels = 1:4)

counts <- cut(dev_sel$dev_log2FoldChange, dev_breaks, include.lowest = T, labels = c("inactive", "weak", "medium", "strong"))
freq_sel <- table(counts)
print(table(counts))

# Add n unique very strong dev between 8 and 10 (take them all)
set.seed(1)
strong <- dev_enhancers[which(dev_quant==4 & !(dev_enhancers$ID %in% dev_sel$ID))]
strong <- strong[!similar(strong) & !similar(strong, dev_sel)]
dev_sel <- Reduce(c, GRangesList(dev_sel, strong))

# Add n enhancers with strength between 4 and 6 (reach 197 total)
set.seed(1)
weak <- dev_enhancers[which(dev_quant==2 & !(dev_enhancers$ID %in% dev_sel$ID))]
weak <- weak[!similar(weak) & !similar(weak, dev_sel)]
weak <- weak[sample(seq(weak), 197-freq_sel[2])]
dev_sel <- Reduce(c, GRangesList(dev_sel, weak))

# Add n enhancers with strength between 4 and 6 (reach 197 total)
set.seed(1)
medium <- dev_enhancers[which(dev_quant==3 & !(dev_enhancers$ID %in% dev_sel$ID))]
medium <- medium[!similar(medium) & !similar(medium, dev_sel)]
medium <- medium[sample(seq(medium), 197-freq_sel[3])]
dev_sel <- Reduce(c, GRangesList(dev_sel, medium))

# Add group columns and add to library file
mcols(dev_sel) <- data.frame(BA_ID= dev_sel$ID, group= "dev", genome= "dm3",
                             detail= cut(dev_sel$dev_log2FoldChange, dev_breaks, include.lowest = T, labels = c("inactive", "weak", "medium", "strong")))
vl_library <- dev_sel

#####################################################################
### Select 100 representative hk enhancers with a wide range of strenght + expected to collaborate endogenously
hk_enhancers <- twist_clean[twist_clean$enhancer_group_detail == "hk"]
hk_breaks <- quantile(hk_enhancers$hk_log2FoldChange, seq(0, 1, length.out= 5))
hk_quant <- cut(hk_enhancers$hk_log2FoldChange, hk_breaks, include.lowest = T, labels = 1:4)

counts <- cut(hk_enhancers$hk_log2FoldChange, hk_breaks, include.lowest = T, labels = 1:4)
print(table(counts))

# First selection of 50 inactive hk enhancers (1st quant), than refine to 25
# 2 steps approach avoids removing a many sequences because they are duplicated in unused ones
set.seed(1)
hk_sel <- hk_enhancers[sample(which(hk_quant==1), 50)]
set.seed(1)
hk_sel <- hk_sel[!similar(hk_sel) & !similar(hk_sel, vl_library)][1:25]

# Select 25 weak hk enhancers (2nd quant)
weak <- hk_enhancers[hk_quant==2]
set.seed(1)
weak <- weak[sample(which(!similar(weak)  & !similar(weak, vl_library)), 25)]
hk_sel <- Reduce(c, GRangesList(hk_sel, weak))

# Select 25 medium hk enhancers (3rd quant)
medium <- hk_enhancers[hk_quant==3]
set.seed(1)
medium <- medium[sample(which(!similar(medium)  & !similar(medium, vl_library)), 25)]
hk_sel <- Reduce(c, GRangesList(hk_sel, medium))

# Select 25 strong hk enhancers (4th quant)
strong <- hk_enhancers[hk_quant==4]
set.seed(1)
strong <- strong[sample(which(!similar(strong)  & !similar(strong, vl_library)), 25)]
hk_sel <- Reduce(c, GRangesList(hk_sel, strong))

# Add group columns
mcols(hk_sel) <- data.frame(BA_ID= hk_sel$ID, group= "hk",  genome= "dm3",
                            detail= cut(hk_sel$hk_log2FoldChange, hk_breaks, include.lowest = T, labels = c("inactive", "weak", "medium", "strong")))
vl_library <- c(vl_library, hk_sel)

#####################################################################
### First selection of 150 flat genomic regions controls, than refine to 100
### 2 steps approach avoids removing a many sequences because they are duplicated in unused ones
ctls <- twist_clean[twist_clean$enhancer_group_detail == "Controls"]
set.seed(1)
ctl_sel <- ctls[sample(seq(ctls), 150)]
set.seed(1)
ctl_sel <- ctl_sel[!similar(ctl_sel) & !similar(ctl_sel, vl_library)][1:100]
mcols(ctl_sel) <- data.frame(BA_ID= ctl_sel$ID, genome= "dm3", group= "control", detail= "flat_genomic_region")
vl_library <- c(vl_library, ctl_sel)

#####################################################################
##### Define ecd and hs inducible enhancers #########################
#####################################################################
### Import enhancers bernardo
ecd_enhancers <- GRanges(fread("../Bernardo_fanny/Rdata/Ecdysone_inducible_final_table.txt"))
ecd_enhancers <- ecd_enhancers[ecd_enhancers$ID %in% twist_clean$ID & ecd_enhancers$enhancer_group_detail == "Ecdysone"]
strand(ecd_enhancers) <- "+"

ecd_sel <- ecd_enhancers[ecd_enhancers$dev_log2FoldChange < 2]
ecd_sel <- ecd_sel[!similar(ecd_sel) & !similar(ecd_sel, vl_library)]
ecd_sel <- ecd_sel[order(ecd_sel$dev_Ecdysone_log2FoldChange, decreasing = T)][1:50]

# clean table
mcols(ecd_sel) <- data.frame(BA_ID= ecd_sel$ID,  genome= "dm3", group= "ecdysone", detail= "low_dCP")
vl_library <- c(vl_library, ecd_sel)

#####################################################################
### Select 50 HS-inducible enhancers with low basal activity in dCP STARR-Seq
hs_enhancers <- GRanges(fread("../Bernardo_fanny/Rdata/Heatshock_inducible_enhancers_final_table.txt"))
hs_enhancers <- hs_enhancers[hs_enhancers$ID %in% twist_clean$ID & hs_enhancers$enhancer_group_detail == "HeatShock"]
strand(hs_enhancers) <- "+"

hs_sel <- hs_enhancers[hs_enhancers$dev_log2FoldChange < 2]
hs_sel <- hs_sel[!similar(hs_sel) & !similar(hs_sel, vl_library)]
hs_sel <- hs_sel[order(hs_sel$dev_HS_log2FoldChange, decreasing = T)][1:50]

# clean table
mcols(hs_sel) <- data.frame(BA_ID= hs_sel$ID,  genome= "dm3", group= "heatshock", detail= "low_dCP")
vl_library <- c(vl_library, hs_sel)

################################################################################
####          Compute other sets of controls                                ####
################################################################################
### First selection of 100 OSC-specific enhancer from the list generated by Bernardo, than refine to 50
### 2 steps approach avoids removing a many sequences because they are duplicated in unused ones
source("scripts/OSC_specific_enhancers_bernardo.R")
OSC_enhancers <- GRanges(get(load("../Bernardo_fanny/Rdata/OSC_specific_enhancers.Rdata")))
seqlevels(OSC_enhancers) <- seqlevels(twist_clean)
strand(OSC_enhancers) <- "+"

OSC_sel <- OSC_enhancers[order(OSC_enhancers$enrichment, decreasing= T)][1:100]
OSC_sel <- OSC_sel[!similar(OSC_sel) & !similar(OSC_sel, vl_library)][1:50]

# clean table
mcols(OSC_sel) <- data.frame(BA_ID= NA, genome= "dm3", group= "OSC", detail= "OSC_specific")
vl_library <- c(vl_library, OSC_sel)

### First selection of 100 exons, than refine to 25
### 2 steps approach avoids removing a many sequences because they are duplicated in unused ones
genes <- get(load("../../../genomes/dm3/dmel-all-no-analysis-r5.57_simple.Rdata"))
mRNA_names <- unlist(genes$Name[genes$type=="mRNA"])
exons <- genes[genes$type=="exon"]
mRNA_exons <- exons[which(unlist(tstrsplit(exons$Name, ":", keep=1)) %in% unlist(tstrsplit(mRNA_names, "-", keep=1)))]
strand(mRNA_exons) <- "+"

n <- 100
set.seed(1)
exon_sel <- exons[sample(seq(exons), n)]
exon_sel <- resize(exon_sel, fix= "center", width = 249)
exon_sel <- exon_sel[!similar(exon_sel) & !similar(exon_sel, vl_library)][1:25]

mcols(exon_sel) <- data.frame(BA_ID= NA, genome= "dm3", group= "control", detail= "exon")
vl_library <- c(vl_library, exon_sel)

################################################################################
####                        Ecoli sequences                                 ####
################################################################################
### First selection of 50 Ecoli sequences than refine to 50
### 2 steps approach avoids removing a many sequences because they are duplicated in unused ones
require(BSgenome.Ecoli.NCBI.20080805)
info <- as.data.frame(seqinfo(BSgenome.Ecoli.NCBI.20080805))
set.seed(1)
chr <- data.table(seqnames= sample(rownames(info), 50, replace= T))[,.N, seqnames]

Ecoli_sel <- mclapply(seq(nrow(chr)), function(i)
{
  c_chr <-chr$seqnames[i]
  c_N <- chr$N[i]
  set.seed(1)
  coor <- sample(info[rownames(info)==c_chr, "seqlengths"], c_N)
  sub <- GRanges(c_chr, IRanges(coor, coor), strand= "+")
  return(sub)
})
suppressWarnings(Ecoli_sel <- Reduce(c, Ecoli_sel))
Ecoli_sel <- resize(Ecoli_sel, width = 249, fix = "center")
Ecoli_sel <- Ecoli_sel[!similar(Ecoli_sel, genome1 = "Ecoli") & !similar(Ecoli_sel, vl_library, genome1= "Ecoli")]
sequence <- as.character(getSeq(BSgenome.Ecoli.NCBI.20080805, Ecoli_sel))
rm <- unique(grep("K|N", sequence))
if(length(rm)>0)
{
  Ecoli_sel <- Ecoli_sel[-rm]
}
Ecoli_sel <- Ecoli_sel[1:25]

mcols(Ecoli_sel) <- data.frame(BA_ID= NA,  genome= "Ecoli_20080805", group= "control", detail= "Ecoli")
genome(Ecoli_sel) <- "2008/08/05"

################################################################################
####         Extract sequences, merge and add linkers                       ####
################################################################################
vl_library$enh_sequence <- as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, vl_library))
Ecoli_sel$enh_sequence <- as.character(getSeq(BSgenome.Ecoli.NCBI.20080805, Ecoli_sel))
vl_library <- c(vl_library, Ecoli_sel)
if(any(grep("N", vl_library$enh_sequence))){stop("'N' within nucleotide sequences!!!")}
if(any(grep("K", vl_library$enh_sequence))){stop("'N' within nucleotide sequences!!!")}

### Add linker seauences
vl_library$linker_ID <- NA

vl_library$linker_ID[sample(which(vl_library$group == "dev" & vl_library$detail == "inactive"), 25)] <- "A"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "weak"), 25)] <- "A"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "medium"), 25)] <- "A"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "strong"), 5)] <- "A"
vl_library$linker_ID[sample(which(vl_library$group == "control" & is.na(vl_library$linker_ID) & vl_library$detail == "flat_genomic_region"), 20)] <- "A"

vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "inactive"), 70)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "weak"), 70)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "medium"), 70)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "dev" & is.na(vl_library$linker_ID) & vl_library$detail == "strong"), 20)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "hk" & is.na(vl_library$linker_ID) & vl_library$detail == "inactive"), 10)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "hk" & is.na(vl_library$linker_ID) & vl_library$detail == "weak"), 10)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "hk" & is.na(vl_library$linker_ID) & vl_library$detail == "medium"), 10)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "hk" & is.na(vl_library$linker_ID) & vl_library$detail == "strong"), 10)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "ecdysone" & is.na(vl_library$linker_ID)), 25)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "heatshock" & is.na(vl_library$linker_ID)), 25)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "OSC" & is.na(vl_library$linker_ID)), 25)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "control" & is.na(vl_library$linker_ID) & vl_library$detail == "flat_genomic_region"), 35)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "control" & is.na(vl_library$linker_ID) & vl_library$detail == "exon"), 10)] <- "B"
vl_library$linker_ID[sample(which(vl_library$group == "control" & is.na(vl_library$linker_ID) & vl_library$detail == "Ecoli"), 10)] <- "B"

vl_library$linker_ID[is.na(vl_library$linker_ID)] <- "C"
overall_design <- as.data.table(mcols(vl_library))
overall_design <- overall_design[,.N,.(group, detail, linker_ID)]
overall_design <- overall_design[order(group, detail, linker_ID)]
overall_design
sum(overall_design$N)

# fw linker sequences
vl_library$fw_linker <- "TTGACAGTGAGCGCGTCTCTCACCG" # Forward primer zuber lab
vl_library$fw_linker <- paste0("G", vl_library$fw_linker ) # Add random G to the 5'end of the oligo to reach 26 nt
# Reverse linker sequences
vl_library$rev_linker <- NA
#rev complement "Sensor_reverse mod"  PCR primer:             CCTAGGATCGACGCGGACAA
vl_library$rev_linker[vl_library$linker_ID == "A"] <- paste0("CCTAGGATCGACGCGGACAACAC", "CG") # Rev complement (seq after EcoRI) of "Sensor_reverse mod" 
                                                                                              # + random CG to the 3'end to reach 25 nt
#rev complement "PGKFmod" PCR primer:                         GCTTTTGAAGCGTGCAGAAT
vl_library$rev_linker[vl_library$linker_ID == "B"] <- paste0("GCTTTTGAAGCGTGCAGAATGAA", "TA") # Rev complement (seq after EcoRI) of "PGKFmod"
                                                                                              # + random TA to the 3'end to reach 25 nt
#rev complement "Angpt2-1" PCR primer:                        CCTGTGCCTTAGACAGCAGCTGAG
vl_library$rev_linker[vl_library$linker_ID == "C"] <- paste0("CCTGTGCCTTAGACAGCAGCTGA", "GT") # Rev complement (seq after EcoRI) of "Angpt2-1"
                                                                                              # + GT to the 3'end to match primer sequence while reaching 25 nt
vl_library$oligo_full_sequence <- paste0(vl_library$fw_linker, vl_library$enh_sequence, vl_library$rev_linker)
vl_library$ID_vl <- paste0(vl_library$group, "_", vl_library$detail, "_", vl_library$linker_ID, "_", sprintf("%05d", seq(vl_library)))

################################################################################
####         SAVE SAVE SAVE SAVE SAVE SAVE SAVE SAVE                        ####
################################################################################
save(vl_library, file="/Volumes/stark/vloubiere/library_design_112019/vl_library_112019.Rdata")
save(vl_library, file= "Rdata/vl_library_112019.Rdata")

saveRDS(vl_library, file="/Volumes/stark/vloubiere/library_design_112019/vl_library_112019.rds")
saveRDS(vl_library, file="Rdata/vl_library_112019.rds")

write.table(vl_library, file= "/Volumes/stark/vloubiere/library_design_112019/vl_library_112019.txt")
write.table(vl_library, file= "Rdata/vl_library_112019.txt")

bed <- vl_library
mcols(bed) <- data.frame(name= vl_library$ID_vl)
export(bed, con= "Rdata/vl_library_112019.bed")

groups <- unique(vl_library$group)
mclapply(groups, function(c_group)
{
  sub <- vl_library[vl_library$group == c_group]
  mcols(sub) <- data.frame(name= sub$ID_vl)
  export(sub, paste0("Rdata/bed_library_groups/", c_group, "_112019_library.bed"))
})

##################################
####        Last checks        ####
##################################
print(paste("N similar sequences=", length(which(similar(vl_library[vl_library$genome == "dm3"])))))

#### Check plots
twist <-  GRanges(readRDS("../Bernardo_fanny/Rdata/twist_enhancers_strongest_strand.rds"))
twist <- as.data.frame(mcols(twist))
twist <- twist[,-(4:6)]

ecd <- GRanges(fread("../Bernardo_fanny/Rdata/Ecdysone_inducible_final_table.txt"))
ecd <- as.data.table(mcols(ecd))
cols <- colnames(ecd)[c(1, grep("Ecd", colnames(ecd)))]
ecd <- ecd[, ..cols]

hs <- GRanges(fread("../Bernardo_fanny/Rdata/Heatshock_inducible_enhancers_final_table.txt"))
hs <- as.data.table(mcols(hs))
cols <- colnames(hs)[c(1, grep("HS", colnames(hs)))]
hs <- hs[, ..cols]

dat <- merge(twist, ecd, by.x= "ID", by.y= "ID")
dat <- merge(dat, hs, by.x= "ID", by.y= "ID")
lib <- as.data.frame(mcols(vl_library))
dat <- as.data.table(merge(lib, dat, by.x= "BA_ID", by.y= "ID"))
setkey(dat, group)

boxplot(dat$dev_log2FoldChange~dat$group, drop=T, main= "all", ylab= "dCP STARR-Seq")
boxplot(dat$hk_log2FoldChange~dat$group, drop=T, main= "all", ylab= "hkCP STARR-Seq")
boxplot(dat[c("control","dev"), dev_log2FoldChange]~dat[c("control","dev"), detail], drop=T, main= "dev + ctls", ylab= "dCP STARR-Seq")
boxplot(dat[c("control","hk"), hk_log2FoldChange]~dat[c("control","hk"), detail], drop=T, main= "hk + ctls", ylab= "hkCP STARR-Seq")
boxplot(dat["ecdysone", .(dev_log2FoldChange, dev_Ecdysone_log2FoldChange)], drop=T, main= "Ecdysone", ylab= "dCP STARR-Seq without/with Ecd")
boxplot(dat["heatshock", .(dev_log2FoldChange, dev_HS_log2FoldChange)], drop=T, main= "Heatshock", ylab= "dCP STARR-Seq without/with HS")
