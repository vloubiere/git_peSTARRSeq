setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Metadata
dat <- fread("Rdata/")

# Large WT oligo pool activities
dat <- readRDS("db/FC_tables/vllib002_DESeq2.rds")
tmp <- tempfile("WT_oligo_pool_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Mutated oligo pool activities
dat <- readRDS("db/FC_tables/vllib029_DESeq2.rds")
tmp <- tempfile("Mutated_oligo_pool_activities_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Focused oligo pool hkCP activities
dat <- readRDS("db/FC_tables/vllib016_DESeq2.rds")
tmp <- tempfile("Focused_oligo_pool_activities_hkCP_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

# Focused oligo pool dCP activities
dat <- readRDS("db/FC_tables/vllib015_DESeq2.rds")
tmp <- tempfile("Focused_oligo_pool_activities_dCP_", fileext = ".txt")
fwrite(dat,
       tmp,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)
vl_dropbox_upload(tmp, ".")

