#######################
#' 1_preprocess_NEON_wrapper.r
#' Description: wrapper function to define NEON data download parameters and
#'   execute data extraction steps
#'   NOTE: This function is designed to download and extract data for one site and 
#'   one year at a time - this is called from a master script which can be run
#'   in parallel on CHPC
#######################

#########################################
# Function defined below
#########################################

preprocess_NEON <- function(site, data_dir, working_dir, year = NA, 
                            overwrite_preQC = TRUE,
                            NEON_token){
    
  # Download hdf5 data bundles for this site
  source(paste0(working_dir, "1a_download_hdf5.r"))
  year_avail <- manage_local_EC_archive(sites = site, data_dir = data_dir, 
                                        years = year)
  
  # The function above returns the given year if download was successful or files already exist
  if(length(year_avail) == 0){
    message(paste0("no data files available for ", site, " - ", year))
    return() # if the year had no data (e.g. you input 2050) then exit function now
  }
    
  
  source(paste0(working_dir, "1b_extract_data.r"))
  extract_NEON_data(site = site, data_dir = data_dir, year = year, 
                    NEON_token = NEON_token,
                    overwrite_preQC = overwrite_preQC)

}
