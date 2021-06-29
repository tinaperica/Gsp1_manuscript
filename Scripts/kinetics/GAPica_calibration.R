#### import libraries
library(tidyverse)
library(viridis)
library(lubridate)

## READ and FORMAT the raw ASCII data from Synergy H1
####### a function to read in all the biotek files, gather the data and reformat the time columns
read_and_gather <- function(file) {
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
  cal_dataset <- index %>% 
    filter(., data_file == file) %>% 
    inner_join(., data_gathered, by = "well") %>% 
    mutate("row" = str_sub(well, 1, 1), "column" = str_sub(well, 2)) %>% 
    mutate("condition" = str_c(date, row, column, Pi_conc, sensor_conc, sep = "-"))
  return(cal_dataset)
}

# combine all the files into one tibble
outfile <- "GTPase_assay/calibration_data_parsed.txt"
index <- read_tsv("GTPase_assay/2018_data/calibration_index.txt", col_names = T)
files <- index %>% 
  pull(data_file) %>% unique()

### read in the data files, join them with the info from the index file and make them tidy
( calibration_data <- do.call("bind_rows", lapply(files, FUN = read_and_gather)) )
calibration_data %>% pull(Pi_conc) %>% unique()
write_tsv(calibration_data, path = outfile)


blank <- calibration_data %>% 
  filter(Pi_conc == 0) %>% 
  group_by(sensor_conc) %>% 
  summarize("mean_blank" = mean(fluorescence)) %>% 
  ungroup()

calibration_data <- calibration_data %>%
  inner_join(., blank, by = "sensor_conc") %>% 
  mutate("norm_fluor" = fluorescence - mean_blank)

calibration_data %>% ggplot(aes(x = norm_fluor, y = Pi_conc, color = as.character(sensor_conc))) + 
  geom_point() +
  scale_color_viridis(discrete = TRUE) +
  theme_bw()

calibration_data %>% ggplot(aes(x = fluorescence, y = Pi_conc, color = as.character(sensor_conc))) + 
  geom_point() +
  scale_color_viridis(discrete = TRUE) +
  theme_bw()


calibration_data %>% 
  filter(Time == 100) %>% 
  ggplot(aes(x = fluorescence, y = Pi_conc, color = as.character(sensor_conc))) + 
  geom_point() +
  scale_color_viridis(discrete = TRUE) +
  theme_bw()


low_sensor_fit <- lm(Pi_conc ~ norm_fluor, data = calibration_data[calibration_data$sensor_conc == 10 & calibration_data$Pi_conc < 3,])
high_sensor_fit <- lm(Pi_conc ~ norm_fluor, data = calibration_data[calibration_data$sensor_conc == 20 & calibration_data$Pi_conc < 13,])
higher_sensor_fit <- lm(Pi_conc ~ norm_fluor, data = calibration_data[calibration_data$sensor_conc == 50 & calibration_data$Pi_conc < 40,])
highest_sensor_fit <- lm(Pi_conc ~ norm_fluor, data = calibration_data[calibration_data$sensor_conc == 70 & calibration_data$Pi_conc < 50,])

calibration_data <- tibble("sensor_conc" = c(20, 70), 
                   "fit_slope" = c(higher_sensor_fit$coefficients[2], highest_sensor_fit$coefficients[2]),
                   "fit_intercept" = c(higher_sensor_fit$coefficients[1], highest_sensor_fit$coefficients[1])) %>% 
  inner_join(., calibration_data, by = "sensor_conc")
calibration_data %>% 
  filter(sensor_conc == 20) %>% 
  group_by(Pi_conc) %>% 
  mutate("mean_fluorescence" = mean(norm_fluor)) %>% 
  ungroup() %>% 
  ggplot(aes(x = fluorescence, y = Pi_conc, color = as.character(sensor_conc))) + 
  geom_point() + 
  #geom_abline(aes(slope = fit_slope, intercept = fit_intercept)) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw()

