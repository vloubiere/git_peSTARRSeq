# Load metadata
tmp <- tempfile(fileext = ".R")
writeLines(readLines("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/pipeline_3.0.R")[1:48],
           tmp)
source(tmp)

# Compute idx
meta[, idx:= .I]

# Identify groups that should be run together
parallel_meta <- rbind(meta[!is.na(DESeq2_group), c("g_type", "grouping"):= .("DESeq2", DESeq2_group)],
                       meta[is.na(DESeq2_group), c("g_type", "grouping"):= .("single", output_prefix)])
setorderv(parallel_meta, "grouping")

# run groups in parallel
parallel_meta[, {
  if(length(.I)>0)
  {
    script <- readLines("/groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/subscripts/pipeline_3.0.R")
    cmd <- c(script[1:48], # First part pipeline
             "\n",
             # Restrict to current group
             rep("#---------------------------------------------------------------------------------------------------------------------------------------------------------------#"),
             rep("#----RESTRICT TO SUBGROUP-----#"),
             deparse(substitute(meta <- meta[idx])),
             rep("#---------------------------------------------------------------------------------------------------------------------------------------------------------------#"),
             rep("\n", 2),
             # End of the pipeline
             script[49:length(script)])
    # Save as a .R script
    tmp <- tempfile(fileext = ".R")
    writeLines(cmd, tmp)
    # Bsub
    Rcmd <- paste("module load build-env/2020; module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript", tmp)
    if(g_type=="DESeq2")
      bsub_cmd <- paste("/groups/stark/software-all/shell/bsub",
                        "-C 12", # N cpus
                        "-m 32", # memory
                        paste("-n", grouping), #name
                        "-T '08:00:00'", #name
                        "-o /groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
                        "-e /groups/stark/vloubiere/projects/pe_STARRSeq/logs/") #stde
    if(g_type=="single") # Less ressources
      bsub_cmd <- paste("/groups/stark/software-all/shell/bsub",
                        "-C 4", # N cpus
                        "-m 8", # memory
                        paste("-n", grouping), #name
                        "-T '02:00:00'", #name
                        "-o /groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
                        "-e /groups/stark/vloubiere/projects/pe_STARRSeq/logs/") #stde
    # Wrap and Submit
    bsub_cmd <- paste0(bsub_cmd, " \"", Rcmd, "\"")
    Sys.unsetenv("SBATCH_RESERVATION")
    Sys.unsetenv("SBATCH_WCKEY")
    job_ID <- system(bsub_cmd, intern = T)
    # Return bsub ID
    unlist(job_ID[2])
  }
}, .(g_type, grouping)]
