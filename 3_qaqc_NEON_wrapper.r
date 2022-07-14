#######################
#' 3_qaqc_NEON_wrapper.r
#' Description: Wrapper function to execute Univ of Utah (UU) QAQC routine on a 
#'   given NEON site's eddy covariance data
#'   NOTE: This function is designed to download and extract data for one site
#'   at a time (works on multiple years) - this is called from a master script 
#'   which can be run in parallel on CHPC
#######################

#########################################
# Function defined below
#########################################

qaqc_NEON <- function(site, data_dir, working_dir, 
                            overwrite_postQC = TRUE,
                            plot_QA_diagnostics = FALSE,
                            plot_QA_extra_diagnostics = FALSE,
                            plot_flux_results = TRUE,
                            turb_axis_min = -10, turb_axis_max = 5,
                            stor_axis_min = -10, stor_axis_max = 5){
  
  # Run UU qaqc procedure
  source(paste0(working_dir, "3a_qaqc.r"))
  site_avail <- qaqc_NEON_data(site = site, data_dir = data_dir,  
                               plot_diagnostics = plot_QA_diagnostics,
                               plot_extra_diagnostics = plot_QA_extra_diagnostics,
                               overwrite_postQC = overwrite_postQC)
  
  # The function above returns the given site if QAQC was successful
  if(length(site_avail) == 0){
    message(paste0("no data files available for ", site))
    return() # if the site had no data then exit function now
  }
  
  # if desired, move to the last step and plot NEE, turbulent and storage fluxes
  if(plot_flux_results){
    source(paste0(working_dir, "3b_plot_fluxes.r"))
    plot_fluxes(site = site, data_dir = data_dir, turb_axis_min = turb_axis_min,
                turb_axis_max = turb_axis_max, stor_axis_min = stor_axis_min,
                stor_axis_max = stor_axis_max)
  }
  
}
