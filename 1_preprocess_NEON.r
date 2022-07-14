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
run_chpc <- TRUE

# source function to be called 
source("1_preprocess_NEON_wrapper.r")

data_dir <- "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
working_dir <-"/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/"
sites <- c("TEAK", "SOAP", "SJER") 
years <- 2017:2022
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
  
  # loop through siteyears and call the wrapper function for each site-year pair
  for(ii in 1:nsiteyears){
    # Call the function with the provided arguments:
    preprocess_NEON(siteyears[ii,1], data_dir, working_dir, siteyears[ii,2], 
                    overwrite_preQC, NEON_token)
  }
}




