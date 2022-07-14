#######################
#' 1b_extract_data.r 
#' Description: extract data needed from NEON data bundles using NEONutilities 
#' package functions
#######################


#########################################
# Define Function(s)
#########################################

#' extract_NEON_data()
#' 
#' uses functions from the NEONutilities package to extract data from the NEON
#' HDF5 data bundles, then supplements 
#' 
#' Inputs:
#' @param site string, 4-letter NEON site code 
#' @param data_dir string, main directory where data folders are stored
#'         (e.g. parent directory to the "1_raw/" and "2_extracted/" dirs)
#' @param year 4-digit year over which to extract fluxes
#' @param NEON_token string, unique personal NEON API token generated at neonscience.org
#' @param overwrite_preQC boolean (T/F) should the preQC outfile be 
#'         overwritten if it already exists? (default = FALSE)
#'
#' @export
#' 
#' Outputs:
#'  None, but an output "preQC" file is saved with desired site and year's worth of pre-QC/QC data


# NOTE!  This requires a NEON API TOKEN as an argument. This can be obtained from neonscience.org
#        This token is accessed when neonUtilities functions are called (for example, search in this file for stackEddy or loadByProduct)
#     details at https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-tokens-tutorial

extract_NEON_data <- function(site, data_dir, year, NEON_token,
                              overwrite_preQC = FALSE) {
  
  # For Debug
  # site <- "SOAP"
  # data_dir <- "/Users/lkunik/Documents/Eddy_Covariance/NEON_data_mgmt/data/"
  # year <- 2019
  # overwrite_preQC <- TRUE
  # NEON_token <- Sys.getenv("LTK_NEON_TOKEN")
  
  
  
  # list required packages
  require(rhdf5)
  require(neonUtilities)
  require(xts)
  require(lubridate)
  require(tidyverse)
  require(pracma)
  #require(UUfluxproc)
 

  # set working directory
  #setwd('/Users/lkunik/Documents/Eddy_Covariance/NEON/')
  
  pracma::tic() #start timer to see how long this takes to run
  
  # location of NEON data files
  raw_data_dir <- paste0(data_dir, "1_raw/")
  
  # file to write output to (check first if file exists)
  outpath <- paste0(data_dir, "2_extracted/NEO_", site, "/")
  outfile <- paste0(outpath, "NEO_", site, "-", year, "-preQC.txt")
  
  if(file.exists(outfile) & !overwrite_preQC){
    message(paste0(outfile, " exists, and overwrite option is disabled... skipping extraction."))
    return()
  }
  
  # added this- need to create a directory with this year of downloaded data
  # for the "one site year" case
  if(!dir.exists(paste0(raw_data_dir, site, "/", year))){
    dir.create(paste0(raw_data_dir, site, "/",year))
  }
  
  # temporarily copy the gzipped h5 files to this year/ directory
  system(paste0("cp ", raw_data_dir, site, "/*.nsae.", year, "-*.gz ", raw_data_dir, 
                site, "/", year))
   
  #----------------------------------
  # validate site
  sites <- c("BONA","CLBJ","CPER","GUAN","HARV","KONZ",
             "NIWO","ONAQ","ORNL","OSBS","PUUM","SCBI",
             "SJER","SRER","TALL","TOOL","UNDE","WOOD",
             "WREF","YELL","ABBY","BARR","BART","BLAN",
             "DELA","DSNY","GRSM","HEAL","JERC","JORN",
             "KONA","LAJA","LENO","MLBS","MOAB","NOGP",
             "OAES","RMNP","SERC","SOAP","STEI","STER",
             "TEAK","TREE","UKFS","DCFS","DEJU")
  
  if (!(site %in% sites)) {
    stop("Check site argument - not a current 4-letter NEON code for a terrestrial site.")
  }
  
  # stack flux data for one year
  message(paste0("Running stackEddy to extract data from hdf5 files - ", site, " - ", year))
  
  # This will unzip the hdf5 data bundles
  fluxes <- neonUtilities::stackEddy(paste0(raw_data_dir, site,"/",year),level="dp04")
  message(paste0("Finished running stackEddy - ", site, " - ", year))
  
  fluxes.flat <- fluxes[[site]] # flatten list structure.
  
  # extract required variables
  # detailed description can be found in NEON.DOC.004571vA.pdf (the Algorithm Theoretical Basis Document
  # for the eddy-covariance data products bundle)
  # note these are not piped to the rename function as in the older code (commented out below)
  fluxes.reduced <- fluxes.flat %>%
    dplyr::select(timeBgn,timeEnd, # time variables
                  # CO2 fluxes
                  data.fluxCo2.nsae.flux,
                  data.fluxCo2.stor.flux,
                  data.fluxCo2.turb.flux,
                  data.fluxCo2.turb.fluxCor,
                  data.fluxCo2.turb.fluxRaw,
                  # H2O fluxes
                  data.fluxH2o.nsae.flux,
                  data.fluxH2o.stor.flux,
                  data.fluxH2o.turb.flux,
                  data.fluxH2o.turb.fluxCor,
                  data.fluxH2o.turb.fluxRaw,
                  # sensible heat fluxes
                  data.fluxTemp.nsae.flux,
                  data.fluxTemp.stor.flux,
                  data.fluxTemp.turb.flux,
                  # friction velocity
                  data.fluxMome.turb.veloFric,
                  # quality flags and metrics (there are lots more)
                  # for CO2 fluxes
                  qfqm.fluxCo2.nsae.qfFinl,
                  qfqm.fluxCo2.nsae.qfFinlStor,
                  qfqm.fluxCo2.nsae.qfFinlTurb,
                  qfqm.fluxCo2.stor.qfFinl,
                  qfqm.fluxCo2.turb.qfFinl,
                  # for H2O fluxes
                  qfqm.fluxH2o.nsae.qfFinl,
                  qfqm.fluxH2o.nsae.qfFinlStor,
                  qfqm.fluxH2o.nsae.qfFinlTurb,
                  qfqm.fluxH2o.stor.qfFinl,
                  qfqm.fluxH2o.turb.qfFinl)
  
  # convert "NEON" time to POSIXct
  fluxes.reduced$timeBgn <- as.POSIXct(fluxes.reduced$timeBgn,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")
  fluxes.reduced$timeEnd <- as.POSIXct(fluxes.reduced$timeEnd,format="%Y-%m-%dT%H:%M:%S.%OSZ",tz="UTC")
  
  # convert fluxes to xts.
  fluxes.notime <- fluxes.reduced %>%
    dplyr::select(-timeBgn,timeEnd)
  
  flux.xts <- as.xts(fluxes.notime,order.by=fluxes.reduced$timeBgn)
  
  # get lat/lon and timezone of station.
  slist <- list.files(path=paste0(raw_data_dir, site,"/", year, "/"),pattern=".h5",
                      full.names=TRUE,recursive=TRUE)
  
  attrs <- h5readAttributes(slist[[1]],site)
  
  lat <- as.numeric(attrs$LatTow)
  lon <- as.numeric(attrs$LonTow)
  tzone <- as.character(attrs$ZoneTime)
  
  #------------------------------------------------------------
  # load met data.
  #------------------------------------------------------------
  
  # set start month.
  start.mon <- paste0(year,"-01")
  end.mon   <- paste0(year,"-12")
  
  # load met data
  # Rh.tmp <- loadByProduct("DP1.00098.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  # Rg.tmp <- loadByProduct("DP1.00023.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  # Ta.tmp <- loadByProduct("DP1.00003.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  # PAR.tmp <- loadByProduct("DP1.00024.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  # Tsoil.tmp <- loadByProduct("DP1.00041.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  # SWC.tmp <- loadByProduct("DP1.00094.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F)
  
  Rh.tmp <- loadByProduct("DP1.00098.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                          token=NEON_token)
  Rg.tmp <- loadByProduct("DP1.00023.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                          token=NEON_token)
  Ta.tmp <- loadByProduct("DP1.00003.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                          token=NEON_token)
  PAR.tmp <- loadByProduct("DP1.00024.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                           token=NEON_token)
  Tsoil.tmp <- loadByProduct("DP1.00041.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                             token=NEON_token)
  SWC.tmp <- loadByProduct("DP1.00094.001",site=site,startdate=start.mon,enddate=end.mon,avg=30,check.size=F,
                           token=NEON_token)
  
  # parse met data
  Rh <- Rh.tmp$RH_30min %>%
    dplyr::filter(as.numeric(horizontalPosition) == 0) %>%
    dplyr::filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    dplyr::select(RHMean,startDateTime)
  Ta <- Ta.tmp$TAAT_30min %>%
    dplyr::filter(as.numeric(horizontalPosition) == 0) %>%
    dplyr::filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    dplyr::select(tempTripleMean,startDateTime)
  Rg <- Rg.tmp$SLRNR_30min %>%
    dplyr::filter(as.numeric(horizontalPosition) == 0) %>%
    dplyr::filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    dplyr::select(inSWMean,outSWMean,inLWMean,outLWMean,startDateTime)
  PAR <- PAR.tmp$PARPAR_30min %>%
    dplyr::filter(as.numeric(horizontalPosition) == 0) %>%
    dplyr::filter(as.numeric(verticalPosition) == max(as.numeric(verticalPosition))) %>%
    dplyr::select(PARMean,startDateTime)
  Tsoil_501 <- Tsoil.tmp$ST_30_minute %>%
    dplyr::filter(as.numeric(horizontalPosition) == 1) %>%
    dplyr::filter(as.numeric(verticalPosition) == 501) %>%
    dplyr::select(soilTempMean,startDateTime)
  Tsoil_502 <- Tsoil.tmp$ST_30_minute %>%
    dplyr::filter(as.numeric(horizontalPosition) == 1) %>%
    dplyr::filter(as.numeric(verticalPosition) == 502) %>%
    dplyr::select(soilTempMean,startDateTime)
  #  ... more depths available - be very careful adding depths, there are lots of places where the 501 etc need to be changed ...
  SWC_501 <- SWC.tmp$SWS_30_minute %>%
    dplyr::filter(as.numeric(horizontalPosition) == 1) %>%
    dplyr::filter(as.numeric(verticalPosition) == 501) %>%
    dplyr::select(VSWCMean,startDateTime)
  SWC_502 <- SWC.tmp$SWS_30_minute %>%
    dplyr::filter(as.numeric(horizontalPosition) == 1) %>%
    dplyr::filter(as.numeric(verticalPosition) == 502) %>%
    dplyr::select(VSWCMean,startDateTime)
  #  ... more depths available ...
  
  names(Rh) <- c("Rh","time")
  names(Ta) <- c("Ta","time")
  names(Rg) <- c("SWdown","SWup","LWdown","LWup","time")
  names(PAR) <- c("PAR","time")
  names(Tsoil_501) <- c("Tsoil_501","time")
  names(Tsoil_502) <- c("Tsoil_502","time")
  names(SWC_501) <- c("SWC_501","time")
  names(SWC_502) <- c("SWC_502","time")
  
  # the SWC_502 seems to be missing at SJER - quick fix, replace with timestamp from 501 and NA
  if (nrow(SWC_502)==0) {
    SWC_502 <- SWC_501     # duplicate SWC_501
    names(SWC_502) <- c("SWC_502","time")     # rename correctly
    SWC_502$SWC_502 <- rep(NA, length(SWC_502$SWC_502))   # fill with NA
  }
  
  # change time vars to character, then to posix ct.
  Rg$time <- as.POSIXct(Rg$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Rh$time <- as.POSIXct(Rh$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Ta$time <- as.POSIXct(Ta$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  PAR$time <- as.POSIXct(PAR$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  
  Tsoil_501$time <- as.POSIXct(Tsoil_501$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  Tsoil_502$time <- as.POSIXct(Tsoil_502$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  SWC_501$time <- as.POSIXct(SWC_501$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  SWC_502$time <- as.POSIXct(SWC_502$time,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")
  
  # convert data to xts objects and then merge.
  Rg.xts <- xts(Rg[,1:4],order.by=Rg$time)
  Rh.xts <- xts(Rh$Rh,order.by=Rh$time)
  Ta.xts <- xts(Ta$Ta,order.by=Ta$time)
  PAR.xts <- xts(PAR$PAR,order.by=PAR$time)
  Tsoil_501.xts <- xts(Tsoil_501$Tsoil_501,order.by=Tsoil_501$time)
  Tsoil_502.xts <- xts(Tsoil_502$Tsoil_502,order.by=Tsoil_502$time)
  SWC_501.xts <- xts(SWC_501$SWC_501,order.by=SWC_501$time)
  SWC_502.xts <- xts(SWC_502$SWC_502,order.by=SWC_502$time)
  
  names(Rg.xts) <- c("SWdown","SWup","LWdown","LWup")
  names(PAR.xts) <- "PAR"
  names(Rh.xts) <- "Rh"
  names(Ta.xts) <- "Ta"
  names(Tsoil_501.xts) <- "Tsoil_501"
  names(Tsoil_502.xts) <- "Tsoil_502"
  names(SWC_501.xts) <- "SWC_501"
  names(SWC_502.xts) <- "SWC_502"
  
  minTime <- as.POSIXct(min(c(min(index(Rg.xts)),min(index(Rh.xts)),min(index(Ta.xts)),
                              min(index(Tsoil_501.xts)),min(index(Tsoil_502.xts)),
                              min(index(SWC_501.xts)),min(index(SWC_502.xts)),
                              min(index(flux.xts)))),origin="1970-01-01")
  maxTime <- as.POSIXct(max(c(max(index(Rg.xts)),max(index(Rh.xts)),max(index(Ta.xts)),
                              max(index(Tsoil_501.xts)),max(index(Tsoil_502.xts)),
                              max(index(SWC_501.xts)),max(index(SWC_502.xts)),
                              max(index(flux.xts)))),origin="1970-01-01")
  
  dummy.ts <- seq.POSIXt(minTime,maxTime,by=1800)
  dummy.data <- rep(NA,length(dummy.ts))
  
  dummy.xts <- xts(dummy.data,order.by=dummy.ts)
  
  # create bound xts
  all.data <- merge.xts(dummy.xts,Rg.xts,Rh.xts,Ta.xts,PAR.xts,Tsoil_501.xts,Tsoil_502.xts,SWC_501.xts,SWC_502.xts,flux.xts)
  
  # calculate vpd from Tair and Rh
  all.data$vpd <- (100-all.data$Rh)/100*ifelse(all.data$Ta < 0,
                                               exp(23.33086-6111.72784/(all.data$Ta+273.15)+0.15215*log(all.data$Ta+273.15)), # vapor over ice
                                               exp(53.67957-6743.769/(all.data$Ta+273.15)-4.8451*log(all.data$Ta+273.15))) # vapor over liquid
  
  all.data <- all.data[paste0(start.mon,"/",end.mon)]
  
  if (tzone == "PST") {
    tzone(all.data) <- "Etc/GMT+8"
  } else if (tzone == "MST") {
    tzone(all.data) <- "Etc/GMT+7"
  } else if (tzone == "CST") {
    tzone(all.data) <- "Etc/GMT+6"
  } else if (tzone == "EST") {
    tzone(all.data) <- "Etc/GMT+5"
  } else if (tzone == "AKST") {
    tzone(all.data) <- "Etc/GMT+9"
  }
  
  #------------------------------------------------------------
  # this is an old and likely obsolete comment
  # load more data if we're going for the expanded package.
  # convert to MPI required format.
  out.data <- coredata(all.data)
  time.tmp <- index(all.data)
  
  # pull out time variables as required.
  yr <- lubridate::year(time.tmp)
  doy <- lubridate::yday(time.tmp)
  hr <- as.numeric(lubridate::hour(time.tmp)+lubridate::minute(time.tmp)/60)
  
  head1 <- c("Year","DoY","Hour",
             # CO2 fluxes
             "data.fluxCo2.nsae.flux",    # the sum of the storage and turbulent fluxes
             "data.fluxCo2.stor.flux",    # storage flux
             "data.fluxCo2.turb.flux",    # turbulent flux
             "data.fluxCo2.turb.fluxCor", # drift-corrected (mean) turbulent flux
             "data.fluxCo2.turb.fluxRaw", # raw turbulent flux without drift correction
             # H2O fluxes
             "data.fluxH2o.nsae.flux",    # the sum of the storage and turbulent fluxes
             "data.fluxH2o.stor.flux",    # storage flux
             "data.fluxH2o.turb.flux",    # turbulent flux
             "data.fluxH2o.turb.fluxCor", # column present but no water drift correction at present
             "data.fluxH2o.turb.fluxRaw", # raw turbulent flux without drift correction
             # sensible heat fluxes
             "data.fluxTemp.nsae.flux",
             "data.fluxTemp.stor.flux",
             "data.fluxTemp.turb.flux",
             # friction velocity (moved with the met data)
             # quality flags and metrics (there are lots more)
             # for CO2 fluxes
             "qfqm.fluxCo2.nsae.qfFinl",
             "qfqm.fluxCo2.nsae.qfFinlStor",
             "qfqm.fluxCo2.nsae.qfFinlTurb",
             "qfqm.fluxCo2.stor.qfFinl",
             "qfqm.fluxCo2.turb.qfFinl",
             # for H2O fluxes
             "qfqm.fluxH2o.nsae.qfFinl",
             "qfqm.fluxH2o.nsae.qfFinlStor",
             "qfqm.fluxH2o.nsae.qfFinlTurb",
             "qfqm.fluxH2o.stor.qfFinl",
             "qfqm.fluxH2o.turb.qfFinl",
             # met
             "SWup","SWdown","LWup","LWdown","Tair","rH","VPD",
             "Ustar",
             "PAR",
             "Tsoil_501","Tsoil_502",
             "SWC_501","SWC_502")
  
  data.out <- data.frame(yr,doy,hr,
                         # CO2 fluxes
                         out.data[,"data.fluxCo2.nsae.flux"],
                         out.data[,"data.fluxCo2.stor.flux"],
                         out.data[,"data.fluxCo2.turb.flux"],
                         out.data[,"data.fluxCo2.turb.fluxCor"],
                         out.data[,"data.fluxCo2.turb.fluxRaw"],
                         # H2O fluxes
                         out.data[,"data.fluxH2o.nsae.flux"],
                         out.data[,"data.fluxH2o.stor.flux"],
                         out.data[,"data.fluxH2o.turb.flux"],
                         out.data[,"data.fluxH2o.turb.fluxCor"],
                         out.data[,"data.fluxH2o.turb.fluxRaw"],
                         # sensible heat fluxes
                         out.data[,"data.fluxTemp.nsae.flux"],
                         out.data[,"data.fluxTemp.stor.flux"],
                         out.data[,"data.fluxTemp.turb.flux"],
                         # friction velocity (moved in with the met data below)
                         # quality flags and metrics (there are lots more)
                         # for CO2 fluxes
                         out.data[,"qfqm.fluxCo2.nsae.qfFinl"],
                         out.data[,"qfqm.fluxCo2.nsae.qfFinlStor"],
                         out.data[,"qfqm.fluxCo2.nsae.qfFinlTurb"],
                         out.data[,"qfqm.fluxCo2.stor.qfFinl"],
                         out.data[,"qfqm.fluxCo2.turb.qfFinl"],
                         # for H2O fluxes
                         out.data[,"qfqm.fluxH2o.nsae.qfFinl"],
                         out.data[,"qfqm.fluxH2o.nsae.qfFinlStor"],
                         out.data[,"qfqm.fluxH2o.nsae.qfFinlTurb"],
                         out.data[,"qfqm.fluxH2o.stor.qfFinl"],
                         out.data[,"qfqm.fluxH2o.turb.qfFinl"],
                         # met
                         out.data[,"SWup"],out.data[,"SWdown"],out.data[,"LWup"],out.data[,"LWdown"],
                         out.data[,"Ta"],
                         out.data[,"Rh"],out.data[,"vpd"],
                         out.data[,"data.fluxMome.turb.veloFric"], # ustar
                         out.data[,"PAR"],
                         out.data[,"Tsoil_501"],out.data[,"Tsoil_502"],
                         out.data[,"SWC_501"],out.data[,"SWC_502"])
  
  # return a data frame.
  names(data.out) <- head1
  
  attr(data.out,"lat") <- lat
  attr(data.out,"lon") <- lon
  attr(data.out,"tzone") <- tzone
  
  # Write data to file
  dir.create(outpath)
  write.table(data.out, file = outfile ,sep = "\t", row.names = FALSE, append = FALSE)
  
  # clean up the copied data bundles for this year so we just have the *.gz files 
  system(paste0("rm -R ", raw_data_dir, site, "/", year, "/"))
  
  
  pracma::toc() #stop timer - print out total function elapsed time
  
}










