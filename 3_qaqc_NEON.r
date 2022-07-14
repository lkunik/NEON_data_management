#######################
#' 3_qaqc_NEON.r
#' Description: Define all arguments to the NEON "QAQC" function
#'   to be called sequentially or in parallel 
#'   NOTE: This is the master script which calls the sub-functions for NEON
#'   eddy covariance QAQC and plot steps
#######################

rm(list=ls()) #Clear environment to avoid issues

#########################################
# Begin Main
#########################################

# Important: run in parallel on CHPC using rslurm?
run_chpc <- TRUE

# source function to be called 
source("3_qaqc_NEON_wrapper.r")

data_dir <- "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
working_dir <-"/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/"
sites <- c("TEAK", "SOAP", "SJER")
nsites <- length(sites)
overwrite_postQC <- TRUE
plot_QA_diagnostics <- FALSE
plot_QA_extra_diagnostics <- FALSE
plot_flux_results <- TRUE

if(run_chpc){
  
  library(rslurm) #assumes rslurm is installed on CHPC
  message("submitting NEON QAQC job to slurm")
  
  #include SLURM options
  slurm_options <- list(time = '5:00:00', #set limit of each task at 5 hours
                        account = 'lin-kp', #indicates that code should be run on kingspeak on john lin's nodes
                        partition = 'lin-kp')
  job <- paste0('NEON_QC')
  
  # For the rest of the inputs, repeat 
  # values to fill out the rest of the arguments list
  args.df <- data.frame(site = sites, 
                        data_dir = rep(data_dir, nsites),
                        working_dir = rep(working_dir, nsites),
                        overwrite_postQC = rep(overwrite_postQC, nsites),
                        plot_QA_diagnostics = rep(plot_QA_diagnostics, nsites),
                        plot_QA_extra_diagnostics = rep(plot_QA_extra_diagnostics, nsites),
                        plot_flux_results = rep(plot_flux_results, nsites))
  
  # send the list of jobs off to SLURM with slurm_apply
  rslurm::slurm_apply(f = qaqc_NEON,
                      params = args.df,
                      jobname = job,
                      nodes = 4,
                      cpus_per_node = 10,
                      pkgs = 'base',
                      slurm_options = slurm_options)
} else{ #if running all site-years sequentially, not on CHPC
  
  # loop through siteyears and call the wrapper function for each site-year pair
  for(ii in 1:nsites){
    
    # Call the function with the provided arguments:
    qaqc_NEON(sites[ii], data_dir, working_dir, 
              overwrite_postQC, 
              plot_QA_diagnostics,
              plot_QA_extra_diagnostics,
              plot_flux_results)
  }
}
