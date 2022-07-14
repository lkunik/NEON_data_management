#######################
#' 2_clean_NEON_files.r
#' Description: clean up duplicate files 
#'   NOTE: This function is designed to download and extract data for one site and 
#'   one year at a time - this is called from a master script which can be run
#'   in parallel on CHPC
#######################

rm(list=ls()) #Clear environment to avoid issues

#########################################
# Load packages
#########################################


#########################################
# Begin Main
#########################################

data_dir <- "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
working_dir <-"/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/"
remove_hdf5_duplicates <- TRUE


source(paste0(working_dir, "1a_download_hdf5.r"))

# Remove duplicate HDF5 files in local data dir
# note, "sites" and "years" arguments don't matter for this - files from all sites and 
# years will be searched for within the function and only the most recent files 
# for a given date will be kept
# ALSO! "dry_run" must be set to "FALSE" for this to work
year_avail <- manage_local_EC_archive(sites = "all", data_dir = data_dir, 
                                      get = FALSE, trim = remove_hdf5_duplicates,
                                      dry_run = FALSE)
