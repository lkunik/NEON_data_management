#######################
#' 1_download_hdf5.r
#' Description: Download NEON HDF5 data bundles for select sites
#' 
#'  Based on Dave Bowling modifications to Rich Fiorella's NEONiso package
#' function named "manage_local_EC_archive()"
#' Changes made:
#'   Accepts multiple site codes as function arguments
#'   Note DB modification to include "expanded" datasets (rather than "basic")
#'   NOTE! for the "trim" functionality, and for the functionality of the 
#'    steps that follow this, it is important that only "expanded"
#'    files exist in your local data directories. This code is not written with
#'    the flexibility to identify between the types of file, and may not
#'    reference the intended file for data extraction.
#######################


#########################################
# Define Function(s)
#########################################

#' manage_local_EC_archive
#'
#' Utility function to help retrieve new EC data and/or prune duplicates,
#' as NEON provisions new data or re-provisions data for an existing site
#' and month. NOTE: CURRENTLY THE TRIM FUNCTIONALITY IS TURNED OFF BY DEFAULT, 
#' AND MUST BE MANUALLY INVOKED WITH trim=TRUE.
#' Also, note that dry_run must be FALSE in order to actually delete files.
#' Defaults true to protect against unintended data loss.
#'
#' @author Rich Fiorella \email{rich.fiorella@@utah.edu} with modifications
#'   from Dave Bowling and Lewis Kunik
#'
#' @param sites atomic vector of every 4-letter site code you want to evaluate
#'              for eddy covariance (EC) data updates
#' @param data_dir Specify the root directory where the local EC store is kept.
#' @param years (Optional) specify which years to refine your search
#' @param get Pull down data from NEON API that does not exist locally?
#' @param trim Search through local holdings, and remove older file where
#'             there are duplicates?
#' @param dry_run List files identified as duplicates, but do not actually
#'             delete them? Default true to prevent unintended data loss.
#' @export
#'

manage_local_EC_archive <- function(sites, data_dir, years = NA, get = TRUE,
                                    trim = FALSE, dry_run = TRUE) {

  
  #DEBUG
  # sites = c("SOAP")
  # data_dir = "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
  # get = TRUE
  # trim = FALSE
  # dry_run = TRUE
  # years = NA
  
  # Make sure that the "1_raw" dir exists, and if not, create it now
  raw_data_dir <- paste0(data_dir, "1_raw/") #note, for this to work, data_dir MUST end in slash
  if(!dir.exists(raw_data_dir)){
    message(paste0(raw_data_dir, " directory not found - creating dir..."))
    dir.create(raw_data_dir)
  }
  
  # define this variable now - will be returned even if no years of data are downloaded
  return_years <- NULL 
  
  
  if (get == TRUE) {

    # Define data/api variables
    data_product <- "DP4.00200.001"
    neon_api_address <- "http://data.neonscience.org/api/v0/products/"

    # see what sites have data
    data_request <- httr::GET(paste0(neon_api_address, data_product))
    # parse JSON object into an R list
    available <- httr::content(data_request, as = "parsed")
    # get number of sites.
    nsites <- length(available$data$siteCodes)

    # loop through sites, and download data.
    for (i in 1:nsites) {
      
      # get site name
      site_name <- available$data$siteCodes[[i]]$siteCode

      # check to see if site [i] is a core/relocatable site
      if (!(site_name %in% sites)) {
        print(paste("Site name", site_name,
                    "is not a requested site...skipping..."))
        next
      } else {
        message(paste0("Downloading raw NEON data bundles for ", site_name, ", year = ", 
                       ifelse(all(is.na(years)), "All available", 
                              ifelse(length(years) == 1, years, 
                                     paste0(years[1], " - ", tail(years, 1))))))
      }

      # get a vector of site months available for site i
      site_months_all <- unlist(available$data$siteCodes[[i]]$availableMonths)
      
      # if distinct years are specified, filter for them here
      if(all(is.na(years))){
        site_months <- site_months_all #default is "NA" meaning get all available years
      } else{
        isite_months <- rep(F, length(site_months_all))
        
        # cumulatively 
        for(yy in years){
          isite_months <- isite_months | grepl(yy, site_months_all) 
        }
        site_months <- site_months_all[isite_months]
      }
      
      if(length(site_months) == 0)
        return(NULL)

      return_years <- list()
      
      # check to see if data folder exists for site, otherwise create.
      ifelse(!dir.exists(paste0(raw_data_dir, site_name)),
             dir.create(paste0(raw_data_dir, site_name)), FALSE)

      # okay - now loop through months and get the data files.
      if (!is.null(length(site_months))) {

        for (j in 1:length(site_months)) {

          # re-query api w/ given site code and month.
          sitemonth_urls_json <- httr::GET(
            unlist(available$data$siteCodes[[i]]$availableDataUrls[isite_months][j]))

          # list files returned.
          sitemonth_urls_parsed <- httr::content(sitemonth_urls_json,
                                                 as = "parsed")

          # extract just file names and URLs
          fnames <- sapply(sitemonth_urls_parsed$data$files, "[[", "name")
          furls  <- sapply(sitemonth_urls_parsed$data$files, "[[", "url")

          # get basic zipfile for now, but should kick out to a
          # function argument later on.
          #fnames_basic <- (grepl("basic", fnames) & grepl("h5.gz", fnames))

          # DB change
          # expanded zipfile
          fnames_expanded <- (grepl("expanded", fnames) & grepl("h5.gz", fnames))

          # check to see if files already exist, and download if missing.
          dl_names <- fnames[fnames_expanded]
          dl_urls  <- furls[fnames_expanded]

          for (k in 1:length(dl_names)) {
            print(dl_names[k])
            if (!length(dl_names[k])==0) {
              if (!is.na(dl_names[k])) {
                
                # At this point, there is either an existing file for this
                # date (in which case we skip and move on) or there is a file
                # to be downloaded from NEON's server. So, we can add this year
                # to the return_years list to indicate that data do exist for this
                # year
                thisyear <- as.numeric(substring(site_months[j], 1, 4))
                return_years <- c(return_years, thisyear)
                
                # check to see if file exists in folder
                if (file.exists(paste0(raw_data_dir, site_name, "/", dl_names[k]))) {
                  print(paste(dl_names[k], "exists...skipping..."))
                  next
                } else { #doesn't exist, so download it.
                  print(paste("Downloading", dl_names[k]))
                  httr::GET(url = dl_urls[k],
                            httr::write_disk(paste0(raw_data_dir,
                                                    site_name,
                                                    "/",
                                                    dl_names[k]),
                                             overwrite = TRUE))
                  
                } # end if(file exists)
              } # end if(filename is NA)
            } # end if(at least one name provided)
          } # end k loop (download files)
        } # end j loop (months available for this site)
      } # end if is.null

      # uncomment if API throttling becomes an issue
      #Sys.sleep(100) 

    } # end i loop (sites)
  } # end if(get new files).

  
  ##########################################
  ## remove duplicate files
  ##########################################
  
  if (trim == TRUE) {
    # list the files in raw_data_dir
  
    files <- list.files(path = raw_data_dir,
                        pattern = "*.h5",
                        recursive = TRUE,
                        full.names = TRUE)

    #need to extract sites, months from file names.
    file_pieces <- strsplit(files, split = ".", fixed = TRUE)

    #Unless there's a very unusual file path, the
    #site name should be element 3 of resulting list,
    #and site year/month/day should be element 8. We'll
    #also want element 10, which is either 'h5' in the old
    #file naming convention, or is the publication date/time.

    fsites <- sapply(file_pieces, "[[", 3) #4-digit site code
    yrmnd  <- sapply(file_pieces, "[[", 8) #year-month-day
    fdiff <- sapply(file_pieces, "[[", 10) #file publication timestamp

    # print(head(fsites))
    # print(head(yrmnd))
    # print(head(fdiff))

    site_list <- unique(fsites)

    for (i in seq_along(site_list)) {
      # get list of files where site == site[i]
      isite <- fsites == site_list[i]

      # check to see if there are duplicates of yrmnd
      yrmnd_isite <- yrmnd[isite] #4-digit site code
      fdiff_isite <- fdiff[isite] #year-month-day
      files_isite <- files[isite] #file publication timestamp
      
      # if there are any unzipped ".h5" files (*not* ".h5.gz") remove them here
      #note in the pattern, the "\\." specifies that only files with the extension
      # of .h5 will be searched, the $ specifies that ".h5" is the end of the string
      # (so ".h5.gz" files will be excluded from this list)
      h5_files <- list.files(path = paste0(raw_data_dir, site_list[i]),
                             pattern = "\\.h5$", full.names = TRUE)
      if(length(h5_files) > 0){
        message(paste0("removing uncompressed .h5 files from ", 
                       raw_data_dir, site_list[i]))
        lapply(h5_files, FUN = function(x) system(paste0("rm ", x)))
      }

      # test to see if there are any duplicates
      if (sum(duplicated(yrmnd_isite)) > 0) {

        # get list of duplicated months
        #print(paste(site_list[i],yrmnd_isite[duplicated(yrmnd_isite) | duplicated(yrmnd_isite, fromLast = TRUE)]))

        # print list of files that are duplicated?
        dups <- duplicated(yrmnd_isite) | duplicated(yrmnd_isite, fromLast = TRUE)

        # narrow list to just list of candidate duplicate files.
        dup_candidates <- files_isite[dups]
        dup_yrmnd       <- yrmnd_isite[dups]
        dup_fdiff      <- fdiff_isite[dups]

        # check to see if one site has an 'h5' fdiff and the duplicates
        # have a date.
        # [LTK] note, not sure what this actually does because h5 file
        # extensions (probably) won't show up in the fdiff slot,
        # those are publication timestamps
        if (!is.null(dup_candidates) & any(fdiff_isite[dups] == 'h5')) {
          
          h5files <- fdiff_isite[dups] == 'h5'
          print(paste('Removing:',dup_candidates[h5files]))
          
          if (!dry_run) {
            file.remove(dup_candidates[h5files]) # remove files.
          }
          
        } else { # none are simply h5, so need to determine which is the most recent file.
          
          # Get publication timestamps of all files
          h5_times <- as.POSIXct(dup_fdiff, format = "%Y%m%dT%H%M%SZ")
          
          for (j in 1:length(unique(dup_yrmnd))) {
            # get publication timestamps associated w/ particular duplicate.
            thisdate_h5_times <- as.POSIXct(dup_fdiff[dup_yrmnd == unique(dup_yrmnd)[j]], format = "%Y%m%dT%H%M%SZ")
            # determine which files are not the most recent.
            # get file names for only this yrmnd.
            dups_yrmnd <- dup_candidates[(dup_yrmnd == unique(dup_yrmnd)[j]) & (h5_times != max(thisdate_h5_times))]
            # print which files to remove
            print(paste('Removing:',dups_yrmnd))
            if (!dry_run) {
              file.remove(dups_yrmnd)
            }
          }
        }
      }
    }
  } #end if(trim) 

  
  #############################
  ## Final return call
  #############################
  return(unique(return_years))
  
} #end function definition
