#######################
#' 2b_plot_fluxes.r
#' Description: create plots of turbulent and storage flux data based on QAQC step
#' requires "postQC.txt" file saved by the QAQC step for given site
#'  
#######################


#########################################
# Define Function(s)
#########################################

#' plot_fluxes()
#' 
#' Loads the "postQC" data and creates turbulent/storage flux plots.
#' Saves to the "data/4_plots/" directory
#' 
#' @param site string, 4-letter NEON site code
#' @param data_dir string, main directory where data folders are stored
#'         (e.g. parent directory to the "1_raw/" and "2_extracted/" dirs)
#' @param turb_axis_min numeric, what should be the minimum value on the turbulent
#'  flux axis (in units of umol CO2 m-2 s-1, default = -10)
#' @param turb_axis_max numeric, what should be the maximum value on the turbulent
#'  flux axis (in units of umol CO2 m-2 s-1, default = 5)
#' @param stor_axis_min numeric, what should be the minimum value on the storage
#'  flux axis (in units of umol CO2 m-2 s-1, default = -10)
#' @param stor_axis_max numeric, what should be the maximum value on the storage
#'  flux axis (in units of umol CO2 m-2 s-1, default = 5)
#' @export
#'


plot_fluxes <- function(site, data_dir, turb_axis_min = -10,
                        turb_axis_max = 5, stor_axis_min = -10,
                        stor_axis_max = 5){
  

  graphics.off()  # close any existing plots
  
  library(ggplot2)
  library(plotly)
  library(patchwork)
  
  # Define output directory for plots, print a message to console
  outpath <- paste0(data_dir, "4_plots/")
  message(paste0("creating flux plots for site: ", site))
  
  # Define input data path, get files within that folder which match this site
  inpath <- paste0(data_dir, "3_qaqc/")
  files <- list.files(inpath, pattern = site)
  
  # If previous steps are run as intended, there should be only one file 
  # matching this site, but check file handling here anyways
  if(length(files) > 1){
    warning(paste0("Found more than 1 file matching site: ", site, 
                    ", using first file in list"))
  }
  
  # read input "postQC" file
  df <- read.table(paste0(inpath, files[1]), header = TRUE)

  # parse data file
  fracyr <- df$Year + (df$DoY + df$Hour/24)/365 #fractional year
  dd <- df$DoY + df$Hour/24 #decimal day of year
  fracday <- df$Hour/24 #fractional day

  # create a factor for PAR bins to be used in plotting
  PARfactor <- rep(NA, length(df$PAR))
  PARinc <- 250 #define the bin increments
  
  # loop through data and assign a bin number to all datapoints
  for (i in seq(0, 2250, by=PARinc)) {
    test <- df$PAR >= i & df$PAR < i+PARinc
    PARfactor[test] = i
  }
  PARfactor <- as.factor(PARfactor)

  # create a factor for Tair bins to be used in plotting
  TAIRfactor <- cut_interval(df$Tair, 20)

  # create a factor for HOUR bins to be used in plotting (24 instead of 48 half hours as in the data files)
  HOURfactor <- rep(NA, length(df$Hour))
  for (i in seq(0, 23, by=1)) {
    test <- df$Hour >= i & df$Hour < i+1
    HOURfactor[test] = i
  }
  HOURfactor <- as.factor(HOURfactor)

  # create a factor for 30-day bins (roughly 12 months) to be used in plotting
  #month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  MONTHfactor <- rep(NA, length(df$DoY))
  DoYinc <- 30
  month_count <- 0
  for (i in seq(0, 330, by=30)) {
    test <- df$DoY >= i & df$DoY < i+30
    month_count <- month_count + 1
    # this bit lumps in DoY 330+ as the final "month"
    if (month_count == 12) {
      test <- df$DoY >= 330
    }
    MONTHfactor[test] = month_count
    #MONTHfactor[test] = month_names[month_count]
  }
  MONTHfactor <- as.factor(MONTHfactor)

  # create dataframe with additional columns useful for ggploting
  # note that some columns are manipulated (added, or ANDed) in the call to cbind.data.frame()
  dat <- cbind.data.frame(fracyr, df$PAR, df$Tair, df$data.fluxCo2.turb.flux, df$data.fluxCo2.stor.flux,
                          df$data.fluxCo2.nsae.flux, (df$data.fluxCo2.turb.flux + df$data.fluxCo2.stor.flux),
                          as.logical(df$qfqm.fluxCo2.nsae.qfFinl), (df$UU_turb_flux_flag & df$UU_stor_flux_flag),
                          as.logical(df$qfqm.fluxCo2.turb.qfFinl), as.logical(df$qfqm.fluxCo2.stor.qfFinl),
                          df$UU_turb_flux_flag, df$UU_stor_flux_flag,
                          PARfactor, HOURfactor, MONTHfactor, TAIRfactor)
  names(dat) <- c("fracyr", "PAR", "Tair", "turb.flux", "stor.flux",
                  "NEON_NEE", "UU_NEE",
                  "NEON_NEE_flag", "UU_NEE_flag",
                  "NEON_turb_flux_flag", "NEON_stor_flux_flag",
                  "UU_turb_flux_flag", "UU_stor_flux_flag",
                  "PARfactor", "HOURfactor", "MONTHfactor","TAIRfactor")

  # invert the NEON flags so TRUE = good data
  # (original NEON flags = 0 indicate good data, but the opposite convention is used for plotting)
  dat$NEON_NEE_flag <- !dat$NEON_NEE_flag
  dat$NEON_turb_flux_flag <- !dat$NEON_turb_flux_flag
  dat$NEON_stor_flux_flag <- !dat$NEON_stor_flux_flag

  #---------------------------------------------------------------------------
  # T response at night, with E0 lines
  # dat2 <- subset(dat, PAR<5 & UU_NEE_flag)
  #
  # T0 <- -46.02
  # Tref <- 0
  # Rref <- 1
  #
  # E0_1 <- 0
  # E0_2 <- 100
  # E0_3 <- 500
  # E0_4 <- -100
  #
  # R1 <- Rref * exp(E0_1*( (1/(Tref - T0)) -  (1/(dat2$Tair-T0)))  )
  # R2 <- Rref * exp(E0_2*( (1/(Tref - T0)) -  (1/(dat2$Tair-T0)))  )
  # R3 <- Rref * exp(E0_3*( (1/(Tref - T0)) -  (1/(dat2$Tair-T0)))  )
  # R4 <- Rref * exp(E0_4*( (1/(Tref - T0)) -  (1/(dat2$Tair-T0)))  )
  #
  # dat3 <- mutate(dat2, R1) # note the first one uses dat2
  # dat3 <- mutate(dat3, R2) # but the others dat3
  # dat3 <- mutate(dat3, R3)
  # dat3 <- mutate(dat3, R4)
  # p1 <- ggplot() +
  #   geom_point(data=dat3, aes(x=Tair, y=UU_NEE), color="grey30") +
  #   geom_line(data=dat3, aes(x=Tair, y=R1), color="black", size=2) +
  #   geom_line(data=dat3, aes(x=Tair, y=R2), color="blue", size=2) +
  #   geom_line(data=dat3, aes(x=Tair, y=R3), color="green", size=2) +
  #   geom_line(data=dat3, aes(x=Tair, y=R4), color="red", size=2) +
  #   ylim(-1,10)
  #
  # rm(dat2, dat3)


  #---------------------------------------------------------------------------
  # T response at night , split by month
  
  message(paste0("plotting T response to ", 
                 sprintf("NEO_%4s_5_monthly_T_response_turb_flux_night.png",site)))
  
  dat2 <- subset(dat, PAR<5 & UU_NEE_flag)
  #p1 <- ggplot() +
    #geom_point(data=dat2, aes(x=Tair, y=turb.flux)) +
    #geom_hline(yintercept = 0, color = "black", linetype="dashed") +
  p1 <- ggplot(dat, aes(x=TAIRfactor, y=UU_NEE)) +
    geom_boxplot() +
    ylab("turbulent flux") + ylim(-2,10) +
    xlab("air temperature") +
    ggtitle(paste0(site,": monthly temperature response of UU_NEE (NIGHT only)")) +
    facet_wrap(~MONTHfactor) +
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
  )

  ggsave(sprintf("NEO_%4s_5_monthly_T_response_turb_flux_night.png",site), path=outpath,
         device = "png", width=11, height=8.5, units="in")

  rm(dat2)

  #---------------------------------------------------------------------------
  # time series of NEE, turbulent flux, storage flux
  # using NEON flags
  
  message(paste0("plotting NEE, turbulent, storage fluxes with NEON flags to ", 
                 sprintf("NEO_%4s_1_0_NEON_flux_time_series.png",site)))
  
  count_all_data <- sum(!is.na(dat$NEON_NEE))
  count_good_data <- sum(dat$NEON_NEE_flag, na.rm=TRUE)
  p1 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=NEON_NEE, color=NEON_NEE_flag)) +
    ylab("net ecosystem flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "net ecosystem exchange", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    ggtitle(paste0(site,": flux time series, with NEON flags")) +
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
      legend.position = "right"
    )
  # turbulent flux
  count_all_data <- sum(!is.na(dat$turb.flux))
  count_good_data <- sum(dat$NEON_turb_flux_flag, na.rm=TRUE)
  p2 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=turb.flux, color=NEON_turb_flux_flag)) +
    ylab("turbulent flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "turbulent flux", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    theme(
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )
  # storage flux
  count_all_data <- sum(!is.na(dat$stor.flux))
  count_good_data <- sum(dat$NEON_stor_flux_flag, na.rm=TRUE)
  p3 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=stor.flux, color=NEON_stor_flux_flag)) +
    ylab("storage flux") + ylim(stor_axis_min,stor_axis_max) +
    xlab("year") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "storage flux", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    theme(
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )

  p <- p1 + p2 + p3 + plot_layout(ncol = 1, nrow=3, widths = c(1, 2))
  #print(p)

  ggsave(sprintf("NEO_%4s_1_0_NEON_flux_time_series.png",site), path=outpath, device = "png", width=11, height=8.5, units="in")


  #---------------------------------------------------------------------------
  # time series of NEE, turbulent flux, storage flux
  # using UU flags
  
  message(paste0("plotting NEE, turbulent, storage fluxes with UU flags to ", 
                 sprintf("NEO_%4s_1_1_UU_flux_time_series.png",site)))
  
  count_all_data <- sum(!is.na(dat$UU_NEE))
  count_good_data <- sum(dat$UU_NEE_flag, na.rm=TRUE)
  p1 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=UU_NEE, color=UU_NEE_flag)) +
    ylab("net ecosystem flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "net ecosystem exchange", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    ggtitle(paste0(site,": flux time series, with UU flags")) +
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
      legend.position = "right"
    )
  # turbulent flux
  count_all_data <- sum(!is.na(dat$turb.flux))
  count_good_data <- sum(dat$UU_turb_flux_flag, na.rm=TRUE)
  p2 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=turb.flux, color=UU_turb_flux_flag)) +
    ylab("turbulent flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "turbulent flux", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    theme(
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )
  # storage flux
  count_all_data <- sum(!is.na(dat$stor.flux))
  count_good_data <- sum(dat$UU_stor_flux_flag, na.rm=TRUE)
  p3 <- ggplot() +
    geom_point(data=dat, aes(x=fracyr, y=stor.flux, color=UU_stor_flux_flag)) +
    ylab("storage flux") + ylim(stor_axis_min,stor_axis_max) +
    xlab("year") + xlim(2017,2022) +
    annotate("text", x=2017, y=10, label= "storage flux", size=6, hjust=0) +
    annotate("text", x=2017, y=5, label=sprintf("all data, n=%d",count_all_data), size=6, hjust=0) +
    annotate("text", x=2017, y=0, label=sprintf("good data, n=%d", count_good_data), size=6, hjust=0) +
    annotate("text", x=2017, y=-5, label=sprintf("%3.1f %%", 100*count_good_data/count_all_data), size=6, hjust=0) +
    theme(
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )

  p <- p1 + p2 + p3 + plot_layout(ncol = 1, nrow=3, widths = c(1, 2))
  #print(p)

  ggsave(sprintf("NEO_%4s_1_1_UU_flux_time_series.png",site), path=outpath, device = "png", width=11, height=8.5, units="in")


  #---------------------------------------------------------------------------
  # light response boxplot, split by month
  
  message(paste0("plotting light response by month to ", 
                 sprintf("NEO_%4s_4_monthly_light_response_turb_flux.png",site)))
  
  p1 <- ggplot(dat, aes(x=PARfactor, y=turb.flux)) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_boxplot() +
    ylab("turbulent flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("PAR") +
    ggtitle(paste0(site,": monthly light response of turbulent flux")) +
    facet_wrap(~MONTHfactor) +
    scale_x_discrete(labels = c("0","","500","","1000","","1500","","2000","","NA"))
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )

  ggsave(sprintf("NEO_%4s_4_monthly_light_response_turb_flux.png",site), path=outpath,
         device = "png", width=11, height=8.5, units="in")


    
    #---------------------------------------------------------------------------
    # light response boxplot, NEE filtered by NEON flags, split by month
  
  message(paste0("plotting light response by month (NEE filtered with NEON flags) to ", 
                 sprintf("NEO_%4s_4_monthly_light_response_NEON_NEE.png",site)))
  
    p1 <- dat %>%
      filter(NEON_NEE_flag == 1) %>%
      ggplot(aes(x=PARfactor, y=NEON_NEE)) +
      geom_hline(yintercept = 0, color = "black", linetype="dashed") +
      geom_boxplot() +
      ylab("NEON_NEE") + ylim(turb_axis_min,turb_axis_max) +
      xlab("PAR") +
      ggtitle(paste0(site,": monthly light response of NEON NEE")) +
      facet_wrap(~MONTHfactor) +
      scale_x_discrete(labels = c("0","","500","","1000","","1500","","2000","","NA"))
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )
    #print(p1)

    ggsave(sprintf("NEO_%4s_4_monthly_light_response_NEON_NEE.png",site), path=outpath,
           device = "png", width=11, height=8.5, units="in")


    #---------------------------------------------------------------------------
    # light response boxplot, NEE filtered by UU flags, split by month
    
    message(paste0("plotting light response by month (NEE filtered with UU flags) to ", 
                   sprintf("NEO_%4s_4_monthly_light_response_UU_NEE.png",site)))
    
    p1 <- dat %>%
      filter(UU_NEE_flag == 1) %>%
      ggplot(aes(x=PARfactor, y=UU_NEE)) +
      geom_hline(yintercept = 0, color = "black", linetype="dashed") +
      geom_boxplot() +
      ylab("UU_NEE") + ylim(turb_axis_min,turb_axis_max) +
      xlab("PAR") +
      ggtitle(paste0(site,": monthly light response of UU NEE")) +
      facet_wrap(~MONTHfactor) +
      scale_x_discrete(labels = c("0","","500","","1000","","1500","","2000","","NA"))
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )

    ggsave(sprintf("NEO_%4s_4_monthly_light_response_UU_NEE.png",site), path=outpath,
           device = "png", width=11, height=8.5, units="in")

    
    
        
  #---------------------------------------------------------------------------
  # diel pattern of NEE, boxplot, split by month
    
  message(paste0("plotting diurnal NEE by month to ", 
                 sprintf("NEO_%4s_1_3_monthly_UU_NEE.png",site)))  
    
  p1 <- ggplot(dat, aes(x=HOURfactor, y=UU_NEE)) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_boxplot() +
    ylab("net ecosystem exchange") + ylim(turb_axis_min,turb_axis_max) +
    xlab("hour of day") +
    ggtitle(paste0(site,": monthly diel pattern of UU_NEE")) +
    facet_wrap(~MONTHfactor) +
    scale_x_discrete(labels = c("0","","2","","4","","6","","8","","10","","12",
                                "","14","","16","","18","","20","","22",""))
  theme(
    plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
    axis.line=element_line(size=0.75),
    axis.title.x=element_text(size=12,color="black"),
    axis.text.x=element_text(size=12,color="black"),
    axis.title.y=element_text(size=12,color="black"),
    axis.text.y=element_text(size=12,color="black"),
    panel.background =element_blank(),
    panel.grid.major=element_line(linetype = "dotted", color="gray60"),
  )

  ggsave(sprintf("NEO_%4s_1_3_monthly_UU_NEE.png",site), path=outpath, device = "png", width=11, height=8.5, units="in")


  #---------------------------------------------------------------------------
  # diel pattern of turbulent flux, boxplot, split by month
  
  message(paste0("plotting diurnal turbulent fluxes by month to ", 
                 sprintf("NEO_%4s_2_monthly_turb_flux.png",site)))  
  
  p1 <- ggplot(dat, aes(x=HOURfactor, y=turb.flux)) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_boxplot() +
    ylab("turbulent flux") + ylim(turb_axis_min,turb_axis_max) +
    xlab("hour of day") +
    ggtitle(paste0(site,": monthly diel pattern of turbulent flux")) +
    facet_wrap(~MONTHfactor) +
    scale_x_discrete(labels = c("0","","2","","4","","6","","8","","10","","12",
                                "","14","","16","","18","","20","","22",""))
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
    )

  ggsave(sprintf("NEO_%4s_2_monthly_turb_flux.png",site), path=outpath, device = "png", width=11, height=8.5, units="in")


  #---------------------------------------------------------------------------
  # diel pattern of storage flux, boxplot, split by month
  
  message(paste0("plotting diurnal storage fluxes by month to ", 
                 sprintf("NEO_%4s_3_monthly_stor_flux.png",site)))  
  
  p1 <- ggplot(dat, aes(x=HOURfactor, y=stor.flux)) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_boxplot() +
    ylab("storage flux") + ylim(stor_axis_min,stor_axis_max) +
    xlab("hour of day") +
    ggtitle(paste0(site,": monthly diel pattern of storage flux")) +
    facet_wrap(~MONTHfactor) +
    scale_x_discrete(labels = c("0","","2","","4","","6","","8","","10","","12",
                                "","14","","16","","18","","20","","22",""))
    theme(
      plot.title = element_text(color="purple", hjust = 0.5, size=18, face="bold"),
      axis.line=element_line(size=0.75),
      axis.title.x=element_text(size=12,color="black"),
      axis.text.x=element_text(size=12,color="black"),
      axis.title.y=element_text(size=12,color="black"),
      axis.text.y=element_text(size=12,color="black"),
      panel.background =element_blank(),
      panel.grid.major=element_line(linetype = "dotted", color="gray60"),
      panel.grid.minor=element_line(linetype = "dotted", color="gray80"),
    )

  ggsave(sprintf("NEO_%4s_3_monthly_stor_flux.png",site), path=outpath, device = "png", width=11, height=8.5, units="in")
 
} #end function definition
