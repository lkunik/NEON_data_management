#######################
#' 2a_qaqc.r
#' Description: Univ of Utah (UU) QAQC routine for NEON eddy covariance data
#' 
#'  Based on Rich Fiorella's NEON extraction functions
#######################


#########################################
# Define Function(s)
#########################################

#' qaqc_NEON_data()
#' 
#' Performs data filtering by removing outlier data based on 2D histograms of  
#' turbulent and storage fluxes. Creates diagnostic plots if desired
#' 
#' @param site string, 4-letter NEON site code
#' @param data_dir string, main directory where data folders are stored
#'         (e.g. parent directory to the "1_raw/" and "2_extracted/" dirs)
#' @param plot_diagnostics boolean, Should diagnostic plots be created? (plots go to 
#'  standard output and are not saved to files)
#' @param plot_extra_diagnostics boolean, Should the extra diagnostic plots be created?
#'  (plots go to standard output and are not saved to files)
#' @param overwrite_postQC boolean (T/F) should the postQC outfile be 
#'         overwritten if it already exists? (default = FALSE)
#' @export
#'

qaqc_NEON_data <- function(site, data_dir, plot_diagnostics = TRUE,
                           plot_extra_diagnostics = FALSE,
                           overwrite_postQC = FALSE){
  
  #for debug
  # data_dir <- "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
  # working_dir <-"/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/"
  # site <- "TEAK"
  # overwrite_postQC <- TRUE
  # #plot_data <- TRUE
  # plot_extra_diagnostics = FALSE
  
  library(ggplot2)
  library(plotly)
  library(patchwork)
  
 
  
  inpath <- paste0(data_dir, "2_extracted/") # path to "preQC" input data
  outpath <- paste0(data_dir, "3_qaqc") # path to "postQC" output data (note - must not end with tailing slash)
  if(!dir.exists(outpath)){
    message(paste0("creating dir ", outpath))
    dir.create(outpath)
  }
  
  # site-specific path has the "NEO" identifier in the sub-folder
  site_path <- sprintf("%sNEO_%s/", inpath, site)
  
  # get files in folder (these are saved by individual years)
  files <- list.files(site_path)
  
  # sanity check - are there any input files to use for this step?
  if(length(files) == 0){
    message(paste0("no input files found for ", site, " site"))
    return(NULL)
  }
  
  # read and concatenate annual files
  for (i in 1:length(files)) {
    df <- read.table(paste0(site_path, files[i]), header = TRUE)
    if (i==1) {
      df_all_years <- df
    }
    else {
      df_all_years <- rbind.data.frame(df_all_years, df)
    }
  }
  
  #re-assign to "df" object now that "df_all_years" contains all years of data
  df <- df_all_years 
  rm(df_all_years) 
  
  # find start and end year for the data
  years <- unique(df$Year)
  start_year <- min(years)
  end_year <- max(years)
  
  # Define the name of the output file
  outfile <- sprintf("%s/NEO_%s-%4d-%4d-postQC.txt", outpath, site, start_year, end_year)
  
  # Check if any output files already exist for this site
  outfiles <- list.files(outpath, pattern = paste0("NEO_", site), full.names = T)
  
  # if there any existing output files, decide what to do with them
  if(length(outfiles) > 0){
    
    # loop through all the existing output files
    for(check_file in outfiles){
      
      # check if this filename is a duplicate of the one we are about to write
      # i.e. the start and end year match the input data we are about to process
      if(check_file == outfile){
        
        # Delete this file *if* the overwrite_postQC input argument is set to TRUE
        if(overwrite_postQC){
          message(paste0(outfile, " found, removing so file can be overwritten"))
          system(paste0("rm ", check_file))
        } else{
          # otherwise, skip this file
          message(paste0(outfile, " exists and overwrite_postQC = FALSE, skipping..."))
          return(site) #return this site name to indicate that the old file exists and is ok to use in future steps
        }
      } else{ #if the file to check != the output file we were going to write (i.e. contains fewer years of data)
        message(paste0("found old existing file ", check_file, ", deleting..."))
        system(paste0("rm ", check_file))
      } #end if(checkfile == outfile) statement 
    } #end checkfile for-loop (existing files in output dir)
  } #end if(any existing files in output dir) stmt
  
  #---------------
  # parse data file
  fracyr <- df$Year + (df$DoY + df$Hour/24)/365 #fractional year
  dd <- df$DoY + df$Hour/24 #decimal day of year
  fracday <- df$Hour/24 #fractional day
  
  # 8 plots per page
  layout(matrix(c(1,2,3,4,5,6,7,8), nrow=2, ncol=4, byrow=TRUE))
  
  #TURB###############################################################
  #---------------
  # QC for turbulent flux using light response during day
  
  PAR_threshold <- 5     # this is used to distinguish night versus day
  
  # light response curve for turbulent flux (all data)
  if (plot_extra_diagnostics) {
    message("plotting light response of turbulent fluxes")
    title = paste0(site, ': light response turb flux (PAR > PAR_threshold)')
    plot(df$PAR, df$data.fluxCo2.turb.flux, type="p", ylim=c(-40,40),
         main=title, pch=16, col="black")
    test <- !is.na(df$data.fluxCo2.turb.flux)
    sum_total <- sum(test, na.rm=TRUE) # this is used below to calculate percentages
    text(500, 37, sprintf("n=%d", sum_total, col="red"))
  
    # split based on NEON flags
    test <- df$qfqm.fluxCo2.turb.qfFinl == 1
    plot(df$PAR[test], df$data.fluxCo2.turb.flux[test], type="p", ylim=c(-40,40),
         main='qfqm.fluxCo2.turb.qfFinl == 1 bad', pch=16, col="dodgerblue",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    text(500, 37, sprintf("n=%d", sum(test, na.rm=TRUE)), col="red")
  
    test <- df$qfqm.fluxCo2.turb.qfFinl == 0
    plot(df$PAR[test], df$data.fluxCo2.turb.flux[test], type="p", ylim=c(-40,40),
         main='qfqm.fluxCo2.turb.qfFinl == 0 good', pch=16, col="forestgreen",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    pct_removed <- 100-100*(sum_total-sum(test, na.rm=TRUE))/sum_total
    text(500, 37, sprintf("n=%d, %1.2f %% removed", sum(test, na.rm=TRUE), pct_removed), col="blue")
  
  }
  
  # now parse based on 2D histogram
  message(paste0("parsing data based on 2D turbulent flux histogram, site: ", site))
  
  # initial data set, x and y must match in dimension
  x <- df$PAR
  y <- df$data.fluxCo2.turb.flux
  
  # filter for only datatime data (force data below the PAR_threshold to NA)
  test <- df$PAR >= PAR_threshold
  x[!test] <- NA
  y[!test] <- NA
  
  # omit extremes in x
  test <- x > 2500
  x[test] <- NA
  
  # omit extremes in y
  flux_ylim <- 40
  test <- y < -1*flux_ylim | y > flux_ylim
  y[test] <- NA
  
  # bins for x and y
  nbins <- 75
  seq_x <- seq(0,2500, length=nbins)
    seq_x[1] <- PAR_threshold     # omit lowest PAR in first bin (night)
  seq_y <- seq(-1*flux_ylim, flux_ylim, length=nbins)
  
  # #-------top---------------------------------------
  # # this code shows the 2D histogram (not needed for removal of outliers)
  # # calculate frequency of occurrence (count) in each bin
  # bin_count <- matrix(data=NA, nbins, nbins)  # bin count, linked to each bin (freq of occur in 2D histogram)
  # bin_count_xy <- rep(NA, length(x))     # the count of the bin in which each xy pair resides
  # for (i in 1:(length(seq_x)-1)) {
  #   for (j in 1:(length(seq_y)-1)) {
  #     test <- (x>=seq_x[i]&x<seq_x[i+1]) & (y>=seq_y[j]&y<seq_y[j+1])
  #     bin_count[i,j] <- sum(test, na.rm=TRUE)
  #     bin_count_xy[test] <- bin_count[i,j]       # bin count, linked to xy pairs
  #   }
  # }
  #
  # tmp <- bin_count[!is.na(bin_count) & bin_count != 0]
  # quant_50 <- quantile(tmp, probs = c(0.50), na.rm=TRUE)
  #
  # test <- tmp > quant_50
  # tmp_50 <- tmp[test]
  #
  # plot(ecdf(tmp))
  # points(ecdf(tmp_50), col="red")
  #
  # hist(tmp)
  #
  # df_2D_hist <- cbind.data.frame(seq_x, seq_y, bin_count)
  # fig <- plot_ly(dat=df_2D_hist, x=seq_y, y=seq_x, z=bin_count, type="surface")
  # fig <- fig %>% layout(
  #   title = "",
  #   scene = list(
  #     xaxis = list(title = "turbulent flux"),
  #     yaxis = list(title = "PAR"),
  #     zaxis = list(title = "# of occurrences")
  #   ))
  
  #-------bot------------------------------------------
  
  # calculate frequency of occurrence (count) in each bin
  bin_count <- rep(NA, length(df$PAR)) # same size as x and y, all NA
  for (i in 1:(length(seq_x)-1)) {
    for (j in 1:(length(seq_y)-1)) {
      test <- (x>=seq_x[i]&x<seq_x[i+1]) & (y>=seq_y[j]&y<seq_y[j+1])
      bin_count[test] <- sum(test, na.rm=TRUE)       # bin count, linked to xy pairs
    }
  }
  
  # calculate quantiles for bin counts, first omitting bin_count = 0 or NA
  # any xy pairs that are outside the bins will retain a bin_count=NA
  tmp <- bin_count[!is.na(bin_count) & bin_count != 0]
  quants <- quantile(tmp, probs = c(0.01, 0.02, 0.03), na.rm=TRUE)
  
  # plot all data, and those in each of the 3 quantiles
  if (plot_extra_diagnostics) {
    message("plotting turbulent flux quantiles")
    plot(x, y, pch=16, col="grey90", main=sprintf("%s: data.fluxCo2.turb.flux", site))
      test_01 <- bin_count <= quants[1]
      points(x[test_01], y[test_01], pch=16, col="red")
  
      test_02 <- bin_count <= quants[2] & bin_count > quants[1]
      points(x[test_02], y[test_02], pch=16, col="green")
  
      test_03 <- bin_count <= quants[3] & bin_count > quants[2]
      points(x[test_03], y[test_03], pch=16, col="purple")
  
      # some statistics
      text(500, 38, sprintf("p=.03 %1.2f %% removed",
                            100*sum(test_03|test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
      text(500, 35, sprintf("p=.02 %1.2f %% removed",
                            100*sum(test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
      text(500, 32, sprintf("p=.01 %1.2f %% removed",
                            100*sum(test_01, na.rm=TRUE)/sum_total), col="blue")
  
      rm(test_01, test_02, test_03)
  
  }
  
  # this test indicates the good data
  test <- bin_count > quants[3]
  
  # save this test for ANDing below
  # this should have same number of rows as the original dataframe df
  test_after_LR_correction <- test
  
  #---------------
  # QC for turbulent flux using T response at night
  message(paste0("parsing turbulent flux data with nighttime Temp response, site: ", site))
  
  # T response curve for turbulent flux (all data with PAR < threshold)
  # filter for only at night (force day part of y and flag columns to NA)
  x <- df$Tair
  y <- df$data.fluxCo2.turb.flux
  NEON_flag <- df$qfqm.fluxCo2.turb.qfFinl
  
  # for daytime periods force x and y to NA (maintains the original dimension)
  test <- df$PAR > PAR_threshold
  x[test] <- NA
  y[test] <- NA
  NEON_flag[test] <- NA
  
  # omit extremes in y
  test <- y < -1*flux_ylim | y > flux_ylim
  y[test] <- NA
  
  # limits for x axis
  xlim_min <- min(x, na.rm=TRUE)
  xlim_max <- max(x, na.rm=TRUE)
  
  if (plot_extra_diagnostics) {
    message("plotting Temp response of nighttime turbulent flux")
    title = paste0(site, ': T response turb flux (PAR < 5)')
    plot(x, y, type="p", ylim=c(-40,40), xlim=c(xlim_min, xlim_max),
         main=title, pch=16, col="black")
    test <- !is.na(y)
    sum_total <- sum(test, na.rm=TRUE)
    text(xlim_min + (xlim_max-xlim_min)*0.25, 37, sprintf("n=%d", sum_total, col="red"))

    # split based on NEON flags
    test <-  NEON_flag == 1
    plot(x[test], y[test], type="p", ylim=c(-40,40), xlim=c(xlim_min, xlim_max),
         main='qfqm.fluxCo2.turb.qfFinl == 1 bad', pch=16, col="dodgerblue",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    text(xlim_min + (xlim_max-xlim_min)*0.25, 37, sprintf("n=%d", sum(test, na.rm=TRUE)), col="red")

    test <- NEON_flag == 0
    plot(x[test], y[test], type="p", ylim=c(-40,40), xlim=c(xlim_min, xlim_max),
         main='qfqm.fluxCo2.turb.qfFinl == 0 good', pch=16, col="forestgreen",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    pct_removed <- 100*(sum_total-sum(test, na.rm=TRUE))/sum_total
    text(xlim_min + (xlim_max-xlim_min)*0.25, 37, sprintf("n=%d, %1.2f %% removed", sum(test, na.rm=TRUE), pct_removed), col="blue")

    rm(NEON_flag)
  }

  # now parse based on 2D histogram
  message("reparsing 2D histogram after flitering for Temp response")
  
  # initial data set, x and y must match in dimension
  
  # bins for x and y (nbins is specified above)
  seq_x <- seq(-40, 50, length=nbins)
  seq_y <- seq(-1*flux_ylim, flux_ylim, length=nbins)
  
  # calculate frequency of occurrence (count) in each bin
  bin_count <- rep(NA, length(df$PAR)) # same size as x and y, all NA
  for (i in 1:(length(seq_x)-1)) {
    for (j in 1:(length(seq_y)-1)) {
      test <- (x>=seq_x[i]&x<seq_x[i+1]) & (y>=seq_y[j]&y<seq_y[j+1])
      bin_count[test] <- sum(test, na.rm=TRUE)       # bin count, linked to xy pairs
    }
  }
  
  # calculate quantiles for bin counts, first omitting bin_count = 0 or NA
  tmp <- bin_count[!is.na(bin_count) & bin_count != 0]
  quants <- quantile(tmp, probs = c(0.01, 0.02, 0.03), na.rm=TRUE)
  
  # plot all data, and those in each of the 3 quantiles
  if (plot_extra_diagnostics) {
    message("plotting turbulent flux quantiles after 2nd round of parsing")
    plot(x, y, pch=16, col="grey90", main=sprintf("%s: data.fluxCo2.turb.flux", site),
         xlim=c(xlim_min, xlim_max),
         ylim = c(-40, 40))
    test_01 <- bin_count <= quants[1]
    points(x[test_01], y[test_01], pch=16, col="red")
  
    test_02 <- bin_count <= quants[2] & bin_count > quants[1]
    points(x[test_02], y[test_02], pch=16, col="green")
  
    test_03 <- bin_count <= quants[3] & bin_count > quants[2]
    points(x[test_03], y[test_03], pch=16, col="purple")
  
    # some statistics
    xloc <- xlim_min + (xlim_max-xlim_min)*0.25
    text(xloc, 38, sprintf("p=.03 %1.2f %% removed",
                          100*sum(test_03|test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
    text(xloc, 35, sprintf("p=.02 %1.2f %% removed",
                          100*sum(test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
    text(xloc, 32, sprintf("p=.01 %1.2f %% removed",
                          100*sum(test_01, na.rm=TRUE)/sum_total), col="blue")
    rm(test_01, test_02, test_03)
  
  }
  
  message("determining UU turbulent flux flags")
  
  # this test indicates the good data
  test <- bin_count > quants[3]
  
  # save this test for ANDing below
  test_after_TR_correction <- test
  rm(test)
  
  # merge data that passed both light response (LR) and temperature response (TR) tests
  test_both <- test_after_TR_correction | test_after_LR_correction
  test_both[is.na(test_both)] <- FALSE
  
  # add a column for final turb flux quality flag
  df <- mutate(df, test_both)
  names(df)[names(df)=="test_both"] <- "UU_turb_flux_flag"
  
  # layout(matrix(c(1), nrow=1, ncol=1, byrow=TRUE))
  # test <- df$UU_turb_flux_flag
  # plot(df$PAR, df$data.fluxCo2.turb.flux, pch=16, col="black", ylim=c(-40,40))
  # points(df$PAR[test], df$data.fluxCo2.turb.flux[test], pch=16, col="red")
  
  rm(test_after_TR_correction, test_after_LR_correction, test_both)

  ###############################################################
  # END TURBULENT FLUX QC
  ###############################################################
  
  
  
  
  
  
  ###############################################################
  # BEGIN STORAGE FLUX QC
  ###############################################################
  
  # QC for storage flux using diel pattern (storage flux versus hour of day)
  
  # diel pattern of storage flux (all data)
  if (plot_extra_diagnostics) {
    message("plotting diel storage flux patterns")
    title = paste0(site, ': storage flux')
    plot(fracday, df$data.fluxCo2.stor.flux, type="p", ylim=c(-40,40),
         main=title, pch=16, col="black")
    test <- !is.na(df$data.fluxCo2.stor.flux)
    sum_total <- sum(test, na.rm=TRUE)
    text(0.1, 37, sprintf("n=%d", sum_total, col="red"))
  
    # split based on NEON flags
    test <- df$qfqm.fluxCo2.stor.qfFinl == 1
    plot(fracday[test], df$data.fluxCo2.stor.flux[test], type="p", ylim=c(-40,40),
         main='qfqm.fluxCo2.stor.qfFinl == 1 bad', pch=16, col="dodgerblue",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    text(0.1, 37, sprintf("n=%d", sum(test, na.rm=TRUE)), col="red")
  
    test <- df$qfqm.fluxCo2.stor.qfFinl == 0
    plot(fracday[test], df$data.fluxCo2.stor.flux[test], type="p", ylim=c(-40,40),
         main='qfqm.fluxCo2.stor.qfFinl == 0 good', pch=16, col="forestgreen",
         sub=sprintf('n=%d', sum(test, na.rm = TRUE)))
    pct_removed <- 100-100*(sum_total-sum(test, na.rm=TRUE))/sum_total
    text(0.2, 37, sprintf("n=%d, %1.2f %% removed", sum(test, na.rm=TRUE), pct_removed), col="blue")
  
  }
  
  # now parse based on 2D histogram
  message(paste0("parsing data based on 2D storage flux histogram, site: ", site))
  
  # initial data set, x and y must match in dimension
  x <- fracday
  y <- df$data.fluxCo2.stor.flux
  # omit extremes in y
  test <- y < -1*flux_ylim | y > flux_ylim
  y[test] <- NA
  
  # bins for x and y
  nxbins <- 48
  nybins <- 75
  seq_x <- seq(0,1, length=nxbins)
  seq_y <- seq(-1*flux_ylim, flux_ylim, length=nybins)
  
  # calculate frequency of occurrence (count) in each bin
  bin_count <- rep(NA, length(x)) # same size as x and y, all NA
  for (i in 1:(length(seq_x)-1)) {
    for (j in 1:(length(seq_y)-1)) {
      test <- (x>=seq_x[i]&x<seq_x[i+1]) & (y>=seq_y[j]&y<seq_y[j+1])
      bin_count[test] <- sum(test, na.rm=TRUE)       # bin count, linked to xy pairs
    }
  }
  
  # calculate quantiles for bin counts, first omitting bin_count = 0 or NA
  tmp <- bin_count[!is.na(bin_count) & bin_count != 0]
  quants <- quantile(tmp, probs = c(0.01, 0.02, 0.03), na.rm=TRUE)
  
  # plot all data, and those in each of the 3 quantiles
  if (plot_extra_diagnostics) {
    plot(x, y, pch=16, col="grey90", main=sprintf("%s: data.fluxCo2.stor.flux", site))
    test_01 <- bin_count <= quants[1]
    points(x[test_01], y[test_01], pch=16, col="red")
  
    test_02 <- bin_count <= quants[2] & bin_count > quants[1]
    points(x[test_02], y[test_02], pch=16, col="green")
  
    test_03 <- bin_count <= quants[3] & bin_count > quants[2]
    points(x[test_03], y[test_03], pch=16, col="purple")
  
    # some statistics
    text(0.2, 38, sprintf("p=.03 %1.2f %% removed",
                           100*sum(test_03|test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
    text(0.2, 35, sprintf("p=.02 %1.2f %% removed",
                           100*sum(test_02|test_01, na.rm=TRUE)/sum_total), col="blue")
    text(0.2, 32, sprintf("p=.01 %1.2f %% removed",
                           100*sum(test_01, na.rm=TRUE)/sum_total), col="blue")
    rm(test_01, test_02, test_03)
  
  }
  
  # flag to indicate good storage flux data
  message("determining UU storage flux flags")
  
  # this should have same number of rows of the original dataframe df
  test <- bin_count > quants[3]
  test_after_stor_correction <- test
  test_after_stor_correction[is.na(test_after_stor_correction)] <- FALSE
  
  # add a column for final stor flux quality flag
  df <- mutate(df, test_after_stor_correction)
  names(df)[names(df)=="test_after_stor_correction"] <- "UU_stor_flux_flag"
  
  ###############################################################
  # END STORAGE FLUX QC
  ###############################################################
  
  # add a column for NEE, based on the UU quality flags
  test <- df$UU_turb_flux_flag & df$UU_stor_flux_flag
  UU_NEE <- rep(NA, length(df$Year))
  UU_NEE[test] <- df$data.fluxCo2.turb.flux[test] + df$data.fluxCo2.stor.flux[test]
  df <- mutate(df, UU_NEE)
  
  # save output file containing new flag columns
  message(paste0("writing postQC file: ", outfile))
  write.table(df, outfile, append=FALSE)
  
  # Check - do we want to plot diagnostics of light/temperature response for
  # turbulent and storage fluxes?
  if(plot_diagnostics){
    
    message("plotting basic diagnostics of light and temperature responses of turbulent and storage fluxes")
    
    # initialize layout for plots
    layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE))
    
    # plot all data and those that pass the turb flux test
    # light response
    plot(df$PAR, df$data.fluxCo2.turb.flux, type="p", ylim=c(-40,40),
         main=sprintf("%s: turb flux light response", site), pch=16, col="grey80",
         ylab="turbulent CO2 flux", xlab="PAR")
    points(df$PAR[df$UU_turb_flux_flag], df$data.fluxCo2.turb.flux[df$UU_turb_flux_flag], type="p", ylim=c(-40,40),
           pch=16, col="green3")
    # T response (all data, not just low PAR)
    plot(df$Tair, df$data.fluxCo2.turb.flux, type="p", ylim=c(-40,40), pch=16, col="grey80",
         main=sprintf("%s: turb flux T response", site),
         ylab="turbulent CO2 flux", xlab="air temperature")
    points(df$Tair[df$UU_turb_flux_flag], df$data.fluxCo2.turb.flux[df$UU_turb_flux_flag], type="p", ylim=c(-40,40),
           pch=16, col="forestgreen")
    
    # plot all data and those that pass the storage flux test
    plot(fracday, df$data.fluxCo2.stor.flux, type="p", ylim=c(-40,40),
         main=sprintf("%s: storage CO2 flux diel pattern", site), pch=16, col="grey80",
         ylab="storage CO2 flux", xlab="hour of day")
    points(fracday[df$UU_stor_flux_flag], df$data.fluxCo2.stor.flux[df$UU_stor_flux_flag], type="p", ylim=c(-40,40),
           pch=16, col="darkgoldenrod2")
    
    #---------------
    # time series of net, turbulent and storage fluxes (all data)
    layout(matrix(1:3, nrow=3, ncol=1, byrow=TRUE))
    
    plot(fracyr, df$data.fluxCo2.nsae.flux, type="p", ylim=c(-40,40),
         main=sprintf('%s: net flux', site), pch=16, col="grey80")
    
    plot(fracyr, df$data.fluxCo2.turb.flux, type="p", ylim=c(-40,40),
         main='turbulent flux', pch=16, col="grey80")
    points(fracyr[df$UU_turb_flux_flag], df$data.fluxCo2.turb.flux[df$UU_turb_flux_flag], pch=16, col="blue")
    
    plot(fracyr, df$data.fluxCo2.stor.flux, type="p", ylim=c(-40,40),
         main='storage flux', pch=16, col="grey80")
    points(fracyr[df$UU_stor_flux_flag], df$data.fluxCo2.stor.flux[df$UU_stor_flux_flag], pch=16, col="blue")
    
    #---------------
    # summary of clean data
    
    # 12 plots per page
    layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE))
    
    # monthly (30-day) diel pattern of turbulent flux (boxplot)
    for (i in seq(0, 330, by=30)) {
      y <- df$data.fluxCo2.turb.flux[df$UU_turb_flux_flag]
      hour_tmp <- df$Hour[df$UU_turb_flux_flag]
      dd_tmp <- dd[df$UU_turb_flux_flag]
      startday <- i
      endday <- i + 30
      if (endday == 360) {endday <- 365}
      period <- dd_tmp > startday & dd_tmp <= endday # time period selected (day of year)
      tmp <- cbind(y[(hour_tmp == 0 | hour_tmp == 0.5) & period],
                   y[(hour_tmp == 1 | hour_tmp == 1.5) & period],
                   y[(hour_tmp == 2 | hour_tmp == 2.5) & period],
                   y[(hour_tmp == 3 | hour_tmp == 3.5) & period],
                   y[(hour_tmp == 4 | hour_tmp == 4.5) & period],
                   y[(hour_tmp == 5 | hour_tmp == 5.5) & period],
                   y[(hour_tmp == 6 | hour_tmp == 6.5) & period],
                   y[(hour_tmp == 7 | hour_tmp == 7.5) & period],
                   y[(hour_tmp == 8 | hour_tmp == 8.5) & period],
                   y[(hour_tmp == 9 | hour_tmp == 9.5) & period],
                   y[(hour_tmp == 10 | hour_tmp == 10.5) & period],
                   y[(hour_tmp == 11 | hour_tmp == 11.5) & period],
                   y[(hour_tmp == 12 | hour_tmp == 12.5) & period],
                   y[(hour_tmp == 13 | hour_tmp == 13.5) & period],
                   y[(hour_tmp == 14 | hour_tmp == 14.5) & period],
                   y[(hour_tmp == 15 | hour_tmp == 15.5) & period],
                   y[(hour_tmp == 16 | hour_tmp == 16.5) & period],
                   y[(hour_tmp == 17 | hour_tmp == 17.5) & period],
                   y[(hour_tmp == 18 | hour_tmp == 18.5) & period],
                   y[(hour_tmp == 19 | hour_tmp == 19.5) & period],
                   y[(hour_tmp == 20 | hour_tmp == 20.5) & period],
                   y[(hour_tmp == 21 | hour_tmp == 21.5) & period],
                   y[(hour_tmp == 22 | hour_tmp == 22.5) & period],
                   y[(hour_tmp == 23 | hour_tmp == 23.5) & period])
      #[LTK] Adding extra check for empty boxplot data:
      if(length(tmp) == 0){
        boxplot(NULL, xlab = "hour of day", ylab="turbulent CO2 flux",
                ylim=c(-10,10),
                main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
      } else{
        boxplot(tmp, xlab = "hour of day", ylab="turbulent CO2 flux",
                ylim=c(-10,10),
                main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
        tmp_median <- apply(tmp, 2, FUN = median, na.rm=TRUE) # apply median function to each column
        lines(1:24, tmp_median, type="b", col="magenta", pch=16)
        lines(c(0,25), c(0,0), lty=2, col="black")
        grid()
      }
      rm(tmp)
    }
    
    #---------------
    # monthly (30-day) diel pattern of storage flux (boxplot)
    for (i in seq(0, 330, by=30)) {
      y <- df$data.fluxCo2.stor.flux[df$UU_stor_flux_flag]
      hour_tmp <- df$Hour[df$UU_stor_flux_flag]
      dd_tmp <- dd[df$UU_stor_flux_flag]
      startday <- i
      endday <- i + 30
      if (endday == 360) {endday <- 365}
      period <- dd_tmp > startday & dd_tmp <= endday # time period selected (day of year)
      tmp <- cbind(y[(hour_tmp == 0 | hour_tmp == 0.5) & period],
                   y[(hour_tmp == 1 | hour_tmp == 1.5) & period],
                   y[(hour_tmp == 2 | hour_tmp == 2.5) & period],
                   y[(hour_tmp == 3 | hour_tmp == 3.5) & period],
                   y[(hour_tmp == 4 | hour_tmp == 4.5) & period],
                   y[(hour_tmp == 5 | hour_tmp == 5.5) & period],
                   y[(hour_tmp == 6 | hour_tmp == 6.5) & period],
                   y[(hour_tmp == 7 | hour_tmp == 7.5) & period],
                   y[(hour_tmp == 8 | hour_tmp == 8.5) & period],
                   y[(hour_tmp == 9 | hour_tmp == 9.5) & period],
                   y[(hour_tmp == 10 | hour_tmp == 10.5) & period],
                   y[(hour_tmp == 11 | hour_tmp == 11.5) & period],
                   y[(hour_tmp == 12 | hour_tmp == 12.5) & period],
                   y[(hour_tmp == 13 | hour_tmp == 13.5) & period],
                   y[(hour_tmp == 14 | hour_tmp == 14.5) & period],
                   y[(hour_tmp == 15 | hour_tmp == 15.5) & period],
                   y[(hour_tmp == 16 | hour_tmp == 16.5) & period],
                   y[(hour_tmp == 17 | hour_tmp == 17.5) & period],
                   y[(hour_tmp == 18 | hour_tmp == 18.5) & period],
                   y[(hour_tmp == 19 | hour_tmp == 19.5) & period],
                   y[(hour_tmp == 20 | hour_tmp == 20.5) & period],
                   y[(hour_tmp == 21 | hour_tmp == 21.5) & period],
                   y[(hour_tmp == 22 | hour_tmp == 22.5) & period],
                   y[(hour_tmp == 23 | hour_tmp == 23.5) & period])
      #[LTK] Adding extra check for empty boxplot data:
      if(length(tmp) == 0){
        boxplot(NULL, xlab = "hour of day", ylab="turbulent CO2 flux",
                ylim=c(-10,10),
                main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
      } else{
        boxplot(tmp,
              xlab = "hour of day", ylab="storage CO2 flux",
              ylim=c(-10,10),
              main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
        tmp_median <- apply(tmp, 2, FUN = median, na.rm=TRUE) # apply median function to each column
        lines(1:24, tmp_median, type="b", col="magenta", pch=16)
        lines(c(0,25), c(0,0), lty=2, col="black")
        grid()
      }
      rm(tmp)
    }
    
    #---------------
    
    # light response of turbulent flux (boxplot)
    for (i in seq(0, 330, by=30)) {
      y <- df$data.fluxCo2.turb.flux[df$UU_turb_flux_flag]
      dd_tmp <- dd[df$UU_turb_flux_flag]
      PAR_tmp <- df$PAR[df$UU_turb_flux_flag]
      startday <- i
      endday <- i + 30
      if (endday == 360) {endday <- 365}
      period <- dd_tmp > startday & dd_tmp <= endday # time period selected (day of year)
      tmp <- cbind(y[PAR_tmp >= 0 & PAR_tmp < 100 & period],
                   y[PAR_tmp >= 100 & PAR_tmp <200 & period],
                   y[PAR_tmp >= 200 & PAR_tmp <300 & period],
                   y[PAR_tmp >= 300 & PAR_tmp <400 & period],
                   y[PAR_tmp >= 400 & PAR_tmp <500 & period],
                   y[PAR_tmp >= 500 & PAR_tmp <600 & period],
                   y[PAR_tmp >= 600 & PAR_tmp <700 & period],
                   y[PAR_tmp >= 700 & PAR_tmp <800 & period],
                   y[PAR_tmp >= 800 & PAR_tmp <900 & period],
                   y[PAR_tmp >= 900 & PAR_tmp <1000 & period],
                   y[PAR_tmp >= 1000 & PAR_tmp <1100 & period],
                   y[PAR_tmp >= 1100 & PAR_tmp <1200 & period],
                   y[PAR_tmp >= 1200 & PAR_tmp <1300 & period],
                   y[PAR_tmp >= 1300 & PAR_tmp <1400 & period],
                   y[PAR_tmp >= 1400 & PAR_tmp <1500 & period],
                   y[PAR_tmp >= 1500 & PAR_tmp <1600 & period],
                   y[PAR_tmp >= 1600 & PAR_tmp <1700 & period],
                   y[PAR_tmp >= 1700 & PAR_tmp <1800 & period],
                   y[PAR_tmp >= 1800 & PAR_tmp <1900 & period],
                   y[PAR_tmp >= 1900 & PAR_tmp <2000 & period],
                   y[PAR_tmp >= 2000 & PAR_tmp <2100 & period],
                   y[PAR_tmp >= 2100 & PAR_tmp <2200 & period],
                   y[PAR_tmp >= 2200 & period])
      #[LTK] Adding extra check for empty boxplot data:
      if(length(tmp) == 0){
        boxplot(NULL, xlab = "hour of day", ylab="turbulent CO2 flux",
                ylim=c(-10,10),
                main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
      } else{
        boxplot(tmp,
              xlab = "PAR",
              ylab= "turbulent CO2 flux", ylim=c(-10,10),
              #names=c("200","400","600","800","1000","1200","1400","1600","1800","2000","2200","2400"),
              main=sprintf('%s: days %d - %d', site, startday, endday), col="grey85")
        grid()
        lines(c(0,2300), c(0,0), lty=2, col="black")
      }
      rm(tmp)
    }
  } #end if(plot_diagnostics)

  # return this site code if successful
  return(site)
  
} #end function definition





