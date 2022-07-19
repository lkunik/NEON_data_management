#######################
#' 1_preprocess_NEON.r
#' Description: Define all arguments to the NEON "preprocess" function
#'   to be called sequentially or in parallel
#'   NOTE: This is the master script which calls the sub-functions for NEON
#'   eddy covariance download/extract steps
#######################

rm(list=ls()) #Clear environment to avoid issues

#########################################
# Begin Main
#########################################

# Important: run in parallel on CHPC using rslurm?
run_chpc <- FALSE

# source function to be called
source("1_preprocess_NEON_wrapper.r")

data_dir <- "/scratch/general/nfs1/u6017162/neon/"
working_dir <- "/uufs/chpc.utah.edu/common/home/lin-group5/ltk/eddy_covariance/NEON/data_proc/"
#sites <- c("TEAK", "SOAP", "SJER")
#years <- 2017:2022
sites <- c("TEAK", "TALL")
years <- c(2019, 2020)
overwrite_preQC <- TRUE
NEON_token <- Sys.getenv("LTK_NEON_TOKEN")

# define the site-years to use as inputs.
siteyears <- expand.grid(sites, years)
nsiteyears <- nrow(siteyears)

if(run_chpc){

  library(rslurm) #assumes rslurm is installed on CHPC
  message("submitting NEON preprocessing job to slurm")

  #include SLURM options
  slurm_options <- list(time = '5:00:00', #set limit of each task at 5 hours
                        account = 'lin-kp', #indicates that code should be run on kingspeak on john lin's nodes
                        partition = 'lin-kp')
  job <- paste0('NEON_preQC')

  # For the rest of the inputs, repeat
  # values to fill out the rest of the arguments list
  args.df <- data.frame(site = siteyears[,1],
                        year = siteyears[,2],
                        data_dir = rep(data_dir, nsiteyears),
                        working_dir = rep(working_dir, nsiteyears),
                        overwrite_preQC = rep(overwrite_preQC, nsiteyears),
                        NEON_token = rep(NEON_token, nsiteyears))

  # send the list of jobs off to SLURM with slurm_apply
  rslurm::slurm_apply(f = preprocess_NEON,
                      params = args.df,
                      jobname = job,
                      nodes = 4,
                      cpus_per_node = 10,
                      pkgs = 'base',
                      slurm_options = slurm_options)
} else{ #if running all site-years sequentially, not on CHPC

  #TEST!! Remove when done. sink output to file.
  filecon <- file("/uufs/chpc.utah.edu/common/home/lin-group5/ltk/eddy_covariance/NEON/data_proc/_rslurm_NEON_preQC/output.20220714.log", open = "wt")
  sink(filecon, type = "output")
  sink(filecon, type = "message")

  # loop through siteyears and call the wrapper function for each site-year pair
  for(ii in 1:nsiteyears){

    # Call the function with the provided arguments:
    preprocess_NEON(siteyears[ii,1], data_dir, working_dir, siteyears[ii,2],
                    overwrite_preQC, NEON_token)

  }

  sink() #TEST
  sink() #TEST
}
