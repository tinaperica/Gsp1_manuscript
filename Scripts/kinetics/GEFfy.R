## this script fits the GEF-mediated nucleotide exchange FRET data from the HiTek Synergy H1 plate reader
### it does either an exponential or a linear fit, using the provided starting parameters
## it outputs the initial velocities and finally the Michaelis-Menten fits
#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)
library(minpack.lm)
library(ggrepel)

## READ and FORMAT the raw ASCII data from Synergy H1
####### a function to read in all the biotek files, gather the data and reformat the time columns
read_and_gather <- function(file) {
  print(file)
  raw_in <-  read_lines(file)
  first_data_row <- grep(raw_in, pattern = "Time\tT")
  last_data_row <- grep(raw_in, pattern = "Results")
  if (length(last_data_row) > 0) {
    n_rows_to_read <- last_data_row - first_data_row - 2
    data_in <-  read_tsv(file, col_names = T, skip = (first_data_row - 1),
                           n_max = n_rows_to_read, locale = locale(encoding = 'windows-1250'))
  } else {
    data_in <-  read_tsv(file, col_names = T, skip = (first_data_row - 1),
                           locale = locale(encoding = 'UTF-8'))
  }
  data_gathered <- data_in %>% 
    select(., -2) %>% 
    gather(., key = well, value = fluorescence, -Time) %>% 
    mutate("Time" = as.numeric(hms(Time)))  #### hms gives a warning when parsing 00:00:00
  dataset <- index %>% 
    filter(., data_file == file) %>% 
    inner_join(., data_gathered, by = "well") %>% 
    mutate("condition" = str_c(date, sample, well, conc, sep = "-"),
           "row" = str_sub(well, 1, 1), "column" = str_sub(well, 2)) %>% 
    mutate("cutoff_time" = ifelse(is.na(cutoff_time), max(Time, na.rm = T), cutoff_time)) %>% 
    filter(Time < cutoff_time)
  return(dataset)
}

# combine all the files into one tibble
(files <- dir("Data/RanGEF_assay/data", pattern = "GEF_FRET_assay", full.names = T))
outfile <- "Data/RanGEF_assay/data_parsed.txt"
### load the index file (has conditions per well)
(index <- read_tsv("Data/RanGEF_assay/data_index.txt", col_names = T))
files <- index %>% 
  pull(data_file) %>% unique()

### read in the data files, join them with the info from the index file and make them tidy
#data_points_to_discard <- read_tsv("GEF_assay/2018_data/data_to_discard.txt")
( dataset <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
write_tsv(dataset, path = outfile)

run_nls <- function(data, debug = FALSE) {
  
  c0 <- data$conc[1]
  GEF_conc <- data$GEF_conc[1]
  k_est <- data$k_est[1]
  span1_est <- data$span1[1]
  span2_est <- data$span2[1]
  deadtime <- data$deadtime[1]
  f_plateau_est <- data$f_plateau[1]
  back_signal <- data$back_signal[1]
  if (back_signal == "exp") {
    min_k_back <- 1e-5
    #min_k_back <- 5e-5
  } else {
    min_k_back <- 0
  }
  if (deadtime > 0) {
    data <- data %>% 
      mutate("Time" = Time + deadtime)
  }
  #### debug mode: if debug is TRUE, print out condition being fit 
  if (debug) { print(unique(data$condition)) }
  start <- list(
                f_plateau = f_plateau_est,
                span1 = span1_est,
                span2 = span2_est,
                k = k_est,
                k_background = 1e-4)
  lower <- c(
                f_plateau_est - 0.2 * f_plateau_est,
                span1_est - 0.05*span1_est,
               0,
              5e-4,
              min_k_back)
  upper <- c(
              f_plateau_est + 0.2 * f_plateau_est,
             span1_est + 0.05*span1_est,
              max(data$observed, na.rm = T),
              0.1,
              3e-4)

  out <- nlsLM(observed ~ span1 * exp(-k * Time) + span2 * (exp(-k_background * Time)) + f_plateau,
               data = data,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 500))
  
  #### save optimal parameters
  f_plateau <- coef(out)[1]
  span1 <- coef(out)[2]
  span2 <- coef(out)[3]
  k <- coef(out)[4]
  print(str_c(data$condition[1], " ", k, " ", f_plateau, " ", span1))
  k_background <- coef(out)[5]

  data$predicted <- span1 * exp(-k * data$Time) + span2 * (exp(-k_background * data$Time)) + f_plateau
  
  data$exchange <- span1 * exp(-k * data$Time)  + f_plateau + span2
  data$background <-  span2 * (exp(-k_background * data$Time)) + f_plateau #+ span1

  #### save optimal parameters in the data table
  data$max_flo <- max(data$observed, na.rm = T)
  data$k <- k
  data$f_plateau <- f_plateau
  data$k_background <- k_background
  data$span1 <- span1
  data$span2 <- span2
  f0 <- span1 + span2 + f_plateau
  data$f_mid <- f0 - span1
  data$f0 <- f0
  
  ## get the fit statistics
  # Chi Square in kaleidagraph is the total sum of the squared errors (sum((y-f(x))/sigma)^2)
  # R in kaleidagraph is the Pearson's R (root 1- chi^2/sum(sig*(y-mean(y))^2)
  sigma <- summary(out)$sigma
  resid <- summary(out)$residuals
  chisq <- out$m$deviance()
  data$chisq <- chisq
  data$pearsonr <- sqrt(1-chisq/sum(sigma*resid^2))
  data$span1_pval <- log10(summary(out)$coefficients[,4][2])
  data$k_pval <- log10(summary(out)$coefficients[,4][4])
    ### calculate the initial rate
  data$vf0 <- (span1 * k * exp(k * 0)) / (GEF_conc*0.001)   ### initial rate in fluorescence units
  return(data)
}

run_linear <- function(data) {
  max_fluo <- max(data$observed, na.rm = T)
  est_span <- 0.6 * max_fluo
  initial_fraction <- data$linear_fit_fraction[1]
  cutoff_fluo <- max_fluo - (initial_fraction * est_span)
  #initial_change <- initial_fraction * (max_fluo - min(data$observed, na.rm = T))
  #cutoff_fluo <- max_fluo - initial_change 
  initial_data <- data %>% filter(observed > cutoff_fluo)
  fit <- lm(observed ~ Time, data = initial_data)
  data %>% ggplot(aes(x = Time, y = observed)) + geom_point() + 
   geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1])
  slope <- fit$coefficients[2]
  data$slope <- slope
  data$intercept <- fit$coefficients[1]
  data$vf0 <- - slope / (data$GEF_conc[1] * 0.001)
  data$f0 <- max_fluo
  data <- data %>% 
    mutate("initial_linear" = ifelse(observed > cutoff_fluo, "initial_linear", "all"))
  return(data)
}
run_linear_with_fixed_intercept <- function(data) {
  print(unique(data$condition))
  fixed_intercept <- data$intercept[1]
  est_span <- 0.6 * fixed_intercept
  initial_fraction <- data$linear_fit_fraction[1]
  cutoff_fluo <- fixed_intercept - (initial_fraction * est_span)
  cutoff_times <- data %>% filter(observed > cutoff_fluo) %>% 
    arrange(observed) %>% pull(Time)
  cut_time <- cutoff_times[1]
  #initial_data <- data %>% filter(observed > cutoff_fluo)
  initial_data <- data %>% filter(Time < cut_time)
  #fit <- lm(observed ~ Time, data = initial_data)
  fit <- lm(I(observed - fixed_intercept) ~ Time - 1, data = initial_data)
  slope <- fit$coefficients
  data %>% ggplot(aes(x = Time, y = observed)) + geom_point() + 
    geom_abline(slope = slope, intercept = fixed_intercept)
  
  data$slope <- slope
  data$vf0 <- - slope / (data$GEF_conc[1] * 0.001)
  data$f0 <- fixed_intercept
  data <- data %>% 
    mutate("initial_linear" = ifelse(Time < cut_time, "initial_linear", "all"))
  return(data)
}
fit_conversion_factor <- function(data) {
  sample <- data$sample[1]
  date <- data$date[1]
  c_plateau_values <- list()
  f_plateau_values <- list()
  conditions <- data %>% 
    #filter(fit == "exp") %>% 
    pull(condition) %>% unique()
  for (i in seq_along(conditions)) {
    data_temp <- data %>% filter(condition == conditions[i])
    c_plateau_values[[i]] = 0.995*data_temp$conc[1]  # account 1/200 molecules being GDP instead of mGTP at steady-state
    f_plateau_values[[i]] = data_temp$f0[1]
  }
  c_plateau_values <- unlist(c_plateau_values)
  f_plateau_values <- unlist(f_plateau_values)
  tib <- tibble("conc" = c_plateau_values, "fluor" = f_plateau_values)
  out <- lm(c_plateau_values ~ f_plateau_values)
  tib %>% ggplot(aes(fluor, conc)) + geom_point() + 
    geom_abline(slope = out$coeff[2], intercept = out$coeff[1]) +
    ggtitle(str_c(sample, " ", date))
  data$conversion_ratio = out$coef[2] # need to divide by the factor that accounts for trp signal decrease for mant bound
  #data$v0 = data$vf0 * data$conversion_ratio
  return(data)
}
fit_poly_conversion_factor <- function(data) {
  sample <- data$sample[1]
  date <- data$date[1]
  c_plateau_values <- list()
  span1_values <- list()
  conditions <- data %>% 
    filter(fit == "exp") %>% 
    pull(condition) %>% unique()
  for (i in seq_along(conditions)) {
    data_temp <- data %>% filter(condition == conditions[i])
    c_plateau_values[[i]] = 0.995*data_temp$conc[1]  # account 1/200 molecules being GDP instead of mGTP at steady-state
    span1_values[[i]] = data_temp$k_background[1]
  }
  c_plateau_values <- unlist(c_plateau_values)
  span1_values <- unlist(span1_values)
  fluor_conc_tib <- tibble("f_mid" = span1_values, "conc" = c_plateau_values)
  linear.fit <- lm(f_mid ~ conc, data = fluor_conc_tib[fluor_conc_tib$conc < 5,])
  fluor_conc_tib %>% ggplot(aes(x = conc, y = f_mid)) + 
    geom_point() +
    geom_abline(slope = linear.fit$coefficients[2], intercept = linear.fit$coefficients[1])
  tib <- tibble("conc" = c_plateau_values, "span1" = span1_values, "span1_sq" = span1_values^2)
  quadratic_fit <- lm(conc ~ span1 + span1_sq, data = tib)
  predicted_conc <- predict(quadratic_fit, list(span1 = span1_values, span1_sq = span1_values^2))
  predicted_line <- tibble("span1" = span1_values, "conc" = predicted_conc)
  predicted_line %>% 
    ggplot(aes(x = span1, y = conc, color = as.character(date))) + 
    geom_line() + 
    geom_point(data = data, aes(x = span1, y = conc))
  vf0_values <- data %>% pull(vf0) %>% unique()
  predicted_v0 <- predict(quadratic_fit, list(span1 = vf0_values, span1_sq = vf0_values^2))
  predicted_v0_tib <- tibble("vf0" = vf0_values, "v0_poly" = predicted_v0)
  predicted_v0_tib %>% 
    ggplot(aes(x = vf0, y = v0_poly)) + 
    geom_point() 
  plot_values <- data %>% select(conc, vf0) %>% 
    unique() %>%
    inner_join(., predicted_v0_tib)
  data <- data %>% 
    inner_join(., predicted_v0_tib, by = "vf0")
  fit <- lm(v0_poly ~ vf0, data = plot_values[plot_values$conc < 5,])
  plot_values %>% ggplot(aes(x = vf0, y = v0_poly)) + geom_point() +
    geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1])
  plot_values %>% ggplot(aes(x = conc, y = vf0)) + geom_point()
  return(data)
  
  
}
plot_raw_data <- function(data, by_sample = FALSE, output = getwd()) {
  
  today = gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
  
  if (by_sample) {
    samples = unique(data$sample)
    for (j in seq_along(samples)) {
        data_mutant <- data %>% filter(sample == samples[j])
        plots <- list()
        conditions = unique(data_mutant$condition)
        for (i in seq_along(conditions)) {
          data_to_plot <- data_mutant %>% arrange(conc) %>% filter(condition == conditions[i]) 
          plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = fluorescence)) +
            geom_point(color = "black") +
            ggtitle(conditions[i])
        }
        pdf(paste0(output, paste(samples[j], 'raw_data.pdf', sep = '_')))
        print(plots)
        dev.off()
    }
  if (!by_sample) {
    plots <- list()
    conditions = unique(data$condition)
    for (i in seq_along(conditions)) {
      data_to_plot <- data %>% filter(condition == conditions[i]) 
      plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed)) +
        geom_point(color = "black") +
        ggtitle(conditions[i])
    }
    pdf(paste0(output, 'raw_data.pdf'))
    print(plots)
    dev.off()
    }
  }
}
plot_fits_show_bkgrnd <- function(data, output = getwd()) {
  samples <- data %>% pull(sample) %>% unique()
  for (j in seq_along(samples)) {
    data_mutant <- data %>% filter(sample == samples[j])
    plots <- list()
    conditions <- data_mutant %>% 
      arrange(conc, condition) %>% pull(condition) %>% unique()
    for (i in seq_along(conditions)) {
      condition_data <- data_mutant %>% 
        filter(condition == conditions[i])
      fit <- condition_data %>% select(fit) %>% unique()
      
      if (fit == "exp") {
        title <- condition_data %>% 
          select(condition, pearsonr, k, span1) %>% 
          unique() %>%
          mutate("title" = str_c("exponential fit\n", condition, "R =", round(pearsonr, 3), "k =", round(k, 4), "span =", round(span1, 0), sep = " ")) %>% pull(title)
        data_to_plot <- condition_data %>% 
          select(Time, observed, predicted, exchange, background) %>% 
          gather(`curve type`, value, -Time)
        plots[[i]] <- ggplot(data_to_plot[data_to_plot[["curve type"]] != "observed", ], mapping = aes(x = Time, y = value, color = `curve type`)) +
          geom_point(data = data_to_plot[data_to_plot[["curve type"]] == "observed", ], color = "black", alpha = 0.5) + 
          geom_line() + ylab("fluorescence") +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          ggtitle(title)
      } else if (fit == "lin") {
        title <- condition_data %>% 
          select(condition, slope, intercept) %>% 
          mutate("title" = str_c("linear fit\n", condition,  "slope =", round(slope, 1), sep = " ")) %>% 
          pull(title) %>% unique()
        slope <- condition_data %>% pull(slope) %>% unique()
        intercept <- condition_data %>% pull(intercept) %>% unique()
        data_to_plot <- condition_data %>% 
          filter(Time < 3000) %>% 
          select(Time, initial_linear, observed) %>% 
          mutate("data" = initial_linear)
        plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = Time, y = observed, color = data)) +
          geom_point(alpha = 0.5) + 
          geom_abline(slope = slope, intercept = intercept) + ylab("fluorescence") +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          ggtitle(title)
      }
    }
    pdf(paste0(output, str_c(samples[j], 'fits_show_bkgrnd.pdf', sep = '_')))
    print(plots)
    dev.off()
  }
}
fit_MM <- function(data) { ### fits Michaelis-Menten
  sample <- data %>% pull(sample) %>% unique()
  root_n <- data %>% pull(condition) %>% unique() %>% length() %>% sqrt()
  print(sample)
  sample_conc <- data %>% pull(conc) %>% unique()
  data_to_fit <- data %>% 
    #filter(Time == 0)
    select(condition, sample, conc, v0, date) %>% 
    unique()
  #start <- list(Vmax = 5, Km = 2)
  start <- list(Vmax = max(data_to_fit$v0, na.rm = T), Km = max(sample_conc)/5)
  lower <- c(0.8 * max(data_to_fit$v0, na.rm = T), 0.5)
  upper <- c(10 * max(data_to_fit$v0), max(sample_conc))
  data_to_fit %>% ggplot(., aes(conc, v0)) + geom_point()
  out <- nlsLM(v0 ~ (conc * Vmax) / (Km + conc),
               data = data_to_fit,
               start = start,
               lower = lower,
               upper = upper,
               control = nls.lm.control(maxiter = 1000))
  Vmax <- coef(out)[1]
  Km <- coef(out)[2]
  kcat <- Vmax
  data$kcat <- kcat
  data$Km <- Km
  kcat_se <- summary(out)$parameters[1,2]
  data$kcat_se <- kcat_se
  kcat_sd <- kcat_se * root_n
  data$kcat_sd <- kcat_sd
  Km_se <- summary(out)$parameters[2,2]
  data$Km_se <- Km_se
  Km_sd <- Km_se * root_n
  data$Km_sd <- Km_sd
  data$kcat_Km <- kcat / Km
  #kcat_Km_se <- sqrt(kcat_se^2) + sqrt(Km_se^2)
  kcat_Km_se <- sqrt( (kcat_se/kcat)^2 + (Km_se/Km)^2 ) * (kcat/Km)
  kcat_Km_sd <- sqrt( (kcat_sd/kcat)^2 + (Km_sd/Km)^2 ) * (kcat/Km)
  data$kcat_Km_se <- kcat_Km_se
  data$kcat_Km_sd <- kcat_Km_sd
  data$predicted_v0 <- (data$conc * Vmax) / (Km + data$conc)
  return(data)
}   
plot_MM <- function(data, output = getwd()) {
  plots <- list()
  samples = unique(data$sample)
  for (i in seq_along(samples)) {
    data_to_plot <- data %>% filter(sample == samples[i]) #%>% filter(Time == 0)
    plots[[i]] <- ggplot(data_to_plot, mapping = aes(x = conc, y = v0, color = as.character(date))) + 
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0, color = sample)) +
      ggtitle(paste(samples[i], paste("Km:", round(data_to_plot$Km[1],2), sep = " "), paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "), sep = "\n"))
  }
  pdf(paste0(output, 'MM.pdf'), width = 15)
  print(plots)
  dev.off()
}
plot_MM_bins <- function(data, output = getwd()) {
  plots <- list()
  samples = unique(data$sample)
  merged_data <- data.frame()
  for (i in seq_along(samples)) {
    
    data_to_plot <- data %>%
      filter(sample == samples[i]) %>%
      #filter(Time == 0) %>% 
      mutate("floor_conc" = ifelse(conc > 1, floor(conc), conc))
    
   data_to_plot <- data_to_plot %>%
      group_by(floor_conc) %>%
      summarise("mean_v0" = mean(v0)) %>%
      inner_join(data_to_plot, by = "floor_conc")

    data_to_plot <- data_to_plot %>%
      group_by(floor_conc) %>%
      summarise("sd_v0" = sd(v0)) %>%
      inner_join(data_to_plot, by = "floor_conc")
    merged_data <- rbind(merged_data, data_to_plot)

    plots[[i]] <- ggplot(data_to_plot, aes(x = floor_conc, y = mean_v0)) +
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0)) +
      geom_errorbar(aes(ymin = mean_v0 - sd_v0, ymax = mean_v0 + sd_v0)) +
      ggtitle(paste(samples[i],
                paste("Km:", round(data_to_plot$Km[1],2), sep = " "),
                paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "),
              sep = "\n"))
  }
  
  pdf(paste0(output, 'MM_bins.pdf'), width = 15)
  print(plots)
  dev.off()
  
  for (i in seq_along(samples)) {
    
    data_to_plot <- data %>%
      filter(sample == samples[i]) 

    plots[[i]] <- ggplot(data_to_plot, aes(x = conc, y = v0, color = as.character(date))) +
      geom_point() + 
      geom_line(aes(x = conc, y = predicted_v0)) +
      ggtitle(paste(samples[i],
                paste("Km:", round(data_to_plot$Km[1],2), sep = " "),
                paste("kcat:", round(data_to_plot$kcat[1],2), sep = " "),
              sep = "\n"))
  }
  pdf(paste0(output, 'MM.pdf'), height = 4, width = 7.5)
  print(plots)
  dev.off()
  
}


#### set output directory
today <- gsub('-', '', today(tzone="US/Pacific")) # set date, for filenaming purposes
output <- str_c('Data/RanGEF_assay/', today, '_GEF_output', '/')
dir.create(output, showWarnings = FALSE)

### first plot raw data for everything
#plot_raw_data(dataset, by_sample = T, output = output)  ### do this only once

#### Fit data assuming photobleaching decay is exponential, with observations in fluorescence units
exp_fits <- dataset %>% 
  filter(fit == "exp") %>% pull(sample) %>% unique()
if (length(exp_fits) > 0) {
  processed.data_exp <- dataset %>%
    filter(fit == "exp") %>% 
    mutate("observed" = fluorescence) %>%
    group_by(sample, condition) %>%
    do(run_nls(., debug = T)) %>%  # fit curve
    ungroup()
}

lin_fits <- dataset %>% 
  filter(fit == "lin") %>% pull(sample) %>% unique()
if (length(lin_fits) > 0) {
processed.data_lin <- dataset %>% 
  filter(fit == "lin") %>% 
  mutate("observed" = fluorescence) %>% 
  group_by(sample, condition) %>%
  do(run_linear_with_fixed_intercept(.)) %>%  # linear fit
  ungroup() 
}

if (length(lin_fits) > 0 & length(exp_fits) > 0) {
  processed.data <- full_join(processed.data_exp, processed.data_lin) 
} else if (length(exp_fits) > 0) {
  processed.data <- processed.data_exp
} else if (length(lin_fits) > 0) {
  processed.data <- processed.data_lin
}

processed.data <- processed.data %>% 
  group_by(sample, date) %>%
  do(fit_conversion_factor(.)) %>% 
  mutate("v0" = vf0 * conversion_ratio) %>% 
  ungroup()

# plot the fits - don't do this every time - too slow
#plot_fits_show_bkgrnd(processed.data, output = output)

if (length(exp_fits) > 0) {
  fit.parameters <- processed.data %>% 
    select(- data_file, -fluorescence, -observed, -predicted, -exchange) %>% 
    unique() %>% 
    arrange(sample, conc) 
} else {
    fit.parameters <- processed.data %>% 
      select(-observed, -data_file) %>% 
      unique() %>% 
      arrange(sample, conc)
}

# fit.parameters %>%
#    ggplot(., aes(x = conc, y = v0 , color = row, shape = as.character(date))) + 
#    geom_point() + facet_grid(~ sample)

initial_rates <- processed.data %>%
  select(sample, v0, conc, condition, date) %>% 
  unique() 
initial_rates %>% write_tsv('Data/RanGEF_assay/GEF_assay_initial_rate_data.txt')
MM.data <- processed.data %>%
  select(sample, v0, conc, condition, date) %>% 
  unique() %>% 
  separate(sample, c("num", "mutant"), remove = F) %>% 
  group_by(mutant) %>%
  do(fit_MM(.)) %>% # fit michaelis-menten
  ungroup() 

plot_MM_bins(MM.data, output = output)

MM.data_to_save <- MM.data %>% 
  select(sample, mutant, kcat, kcat_se, kcat_sd, Km, Km_se, Km_sd, kcat_Km, kcat_Km_se, kcat_Km_sd) %>% unique()
write_tsv(MM.data_to_save, file.path(output, "20200506_GEF_MM_data.txt"))
