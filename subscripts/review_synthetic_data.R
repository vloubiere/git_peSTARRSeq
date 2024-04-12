setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data
dat <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
synt <- CJ(indL= seq(-1.5, 7.5, .1),
           indR= seq(-1.5, 7.5, .1))
synt[, additive:= log2(2^indL+2^indR-1)]
synt[, multiplicative:= indL+indR]

# The number of transcripts goes as follows:


synt[indL==0 & indR==0, .(indL,
                          indR,
                          `Exp. add`= paste0(additive, " (n=", 2^additive, " transcripts)"), 
                          `Exp. mult`= paste0(multiplicative, " (n=", 2^multiplicative, " transcripts)"))]
synt[indL==1 & indR==1, .(indL,
                          indR,
                          `Exp. add`= paste0(additive, " (n=", 2^additive, " transcripts)"), 
                          `Exp. mult`= paste0(multiplicative, " (n=", 2^multiplicative, " transcripts)"))]


plot(dat$additive,
     dat$multiplicative)
abline(0, 1)
abline(h= 0)
abline(v= 0)
plot(dat$additive,
     dat$multiplicative,
     xlim= c(-0.25, .25),
     ylim= c(-0.25, .25))
abline(0, 1)
abline(h= 0)
abline(v= 0)