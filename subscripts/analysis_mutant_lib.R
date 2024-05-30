setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(ggseqlogo)

# Import data ----
dat <- readRDS("Rdata/vl_library_twist015_112022.rds")[sublib=="A"]

# Make sequence ID unique (numbers are only unique within one group) ----
dat[, c("group", "mut", "enh"):= tstrsplit(ID, "_", keep= c(1,2,4))]

# Twist sequences ----
twist <- merge(dat[grepl("^twist", group) & mut=="WT", .(enh, enh_sequence)],
               dat[grepl("^twist", group) & mut=="mutTwist", .(enh, enh_sequence)],
               by= "enh",
               suffixes= c("_wt", "_mut"))
twist$pos <- vl_motif_pos(twist$enh_sequence_wt,
                          sel = "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
                          genome = "dm3",
                          collapse.overlapping = T,
                          p.cutoff = 1e-4)
# twist$pos <- vl_motif_pos(twist$enh_sequence_mut,
#                           sel = "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
#                           genome = "dm3",
#                           collapse.overlapping = T,
#                           p.cutoff = 1e-4)

pdf("pdf/draft/reviewer_motif_mutation.pdf", width = 8, height = 2)
plot.new()
x <- seq(0, 1, length.out= 250)
twist[1, {
  rect(x[pos[[1]]$start], 0, x[pos[[1]]$end], 1)
  sapply(1:249, function(i) vl_plotLetter(substr(enh_sequence_wt, i, i),
                                          xleft = x[i],
                                          ytop= 1,
                                          width = x[2]-x[1],
                                          height = .5))
  sapply(1:249, function(i) vl_plotLetter(substr(enh_sequence_mut, i, i),
                                          xleft = x[i],
                                          ytop= .5,
                                          width = x[2]-x[1],
                                          height = .5))
}]
dev.off()

twist <- twist[, {
  .c <- pos[[1]]
  .c <- .c[, {
    wt <- substr(enh_sequence_wt, start, end)
    # ifelse(strand=="-", vl_revComp(wt), wt)
    mut <- substr(enh_sequence_mut, start, end)
    # ifelse(strand=="-", vl_revComp(mut), mut)
    .(wt, mut)
  # }, .(start, end, strand)]
  }, .(start, end)]
  .c
}, .(enh, enh_sequence_wt)]


set.seed(1)
pwm <- consensusMatrix(DNAStringSet(twist$wt), as.prob = T)
p <- ggseqlogo::ggseqlogo(pwm[1:4,], method='p') + ggtitle("WT")
plot(p)
pwm <- consensusMatrix(DNAStringSet(twist$mut), as.prob = T)
p <- ggseqlogo::ggseqlogo(pwm[1:4,], method='p') + ggtitle("WT")
plot(p)
