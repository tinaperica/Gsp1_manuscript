library(tidyverse)
library(scales)
library(minpack.lm)
source('ucsf_colors.R')

# Parsing function, which returns a list of two dataframes per experiment: data and metadata
parse_CD_file = function(filepath) {
  
  connection = file(filepath, "r")
  metadata <- data.frame(param=character(), value=character(), stringsAsFactors=F) 
  last_param = ''
  while (T) {
    line <- strsplit(readLines(connection, n = 1), split = '\t')
    last_param <- line[[1]][1]
    if ( length(line) == 0 | last_param == 'XYDATA') {
      break
    }
    metadata <- rbind(metadata, data.frame('param'=line[[1]][1],
                                           'value'=line[[1]][2],
                                           stringsAsFactors=F))
  }
  
  X <- as.character(subset(metadata, param == 'XUNITS')$value)
  Y <- as.character(subset(metadata, param == 'YUNITS')$value)
  Y2 <- as.character(subset(metadata, param == 'Y2UNITS')$value)
  
  data <- data.frame(X=character(), Y=character(), Y2=character(), stringsAsFactors=F) 

  while (T) {
    line <- strsplit(readLines(connection, n = 1), split = '\t')
    if ( length(line) == 0 ) {
      break
    }
    data <- rbind(data, data.frame('X'=as.numeric(line[[1]][1]),
                                   'Y'=as.numeric(line[[1]][2]),
                                   'Y2'=as.numeric(line[[1]][3])))
  }
  colnames(data) <- c(X, Y, Y2)
  close(connection)
  
  # Scale CD data from ellipticity (millidegrees) to molar ellipticity (deg*cm*dmol-1)
  # pathlength = 0.2 cm, concentration = 2e-6 M
  data <- mutate(data, 'Molar ellipticity' = 100 * (`CD[mdeg]`/1000) / (0.2 *2e-6))
  
  return(list(data = data,
              metadata = metadata,
              name = substr(filepath, 25, nchar(filepath) - 9)))
}

# Function for fitting a sigmoid to a CD melt curve
fit_sigmoid_to_melt <- function(data) {
  x <- data$`Temperature [C]`
  y <- data$`Molar ellipticity`
  approx_Tm <- x[which.min(abs(y - mean(y)))]
  nls(y ~ ((yf + mf*x) + (yu + mu*x)*exp(m*(1/Tm - 1/x))) / (1 + exp(m*(1/Tm - 1/x))),
        start = list(m = 2000, Tm = approx_Tm, mf = 0, mu = 0, yf = min(y), yu = max(y)))
}

# Plotting functions, for melts and scans
plot_melt <- function(expt) {

  sigmoid_model <- fit_sigmoid_to_melt(expt$data)
  data_to_plot <- mutate(expt$data, 'fit' = predict(sigmoid_model, `Temperature [C]`)) 
  Tm <- summary(sigmoid_model)$coefficients[2,1]
  x <- data_to_plot$`Temperature [C]`
  y <- data_to_plot$fit

  plot <-
    ggplot(data_to_plot) +
    geom_point(aes_string(x = '`Temperature [C]`', y = '`Molar ellipticity`'), size = 0.1) +
    geom_line(aes_string(x = '`Temperature [C]`', y = 'fit'), size = 0.2, color = ucsf_colors$pink1) +
    scale_y_continuous(breaks = pretty_breaks(), label=scientific_format()) +
    geom_vline(xintercept = Tm, color = ucsf_colors$blue1, linetype = 'dashed', size = 0.2) +
    ggtitle(paste(expt$name, paste0('approx Tm = ', round(Tm, digits = 1)), sep='\n')) +
    xlab('Temperature (°C)') +
    ylab('Molar Ellipticity (deg * cm / dmol)') +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.ticks = element_line(size = 0.05),
      axis.ticks.length = unit(0.05, 'cm'),
      axis.line = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5, size = 6)
    )
  
  return(list(plot = plot, Tm_data = c(expt$name, round(Tm, digits = 1))))
}

plot_scan <- function(expt) {
  ggplot(expt$data, aes_string(x = '`NANOMETERS`', y = '`Molar ellipticity`')) +
    geom_line(size = 0.2) +
    scale_y_continuous(breaks= pretty_breaks()) +
    geom_hline(yintercept = 0, color = ucsf_colors$pink1, linetype = 'dashed', size = 0.2) + 
    ggtitle(expt$name) +
    xlab('Wavelength (nm)') +
    ylab('Molar Ellipticity (deg * cm / dmol)') +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.ticks = element_line(size = 0.05),
      axis.ticks.length = unit(0.05, 'cm'),
      axis.line = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5, size = 6)
    )
}

# set data directory
directory_path = 'Data/Circular_Dichroism'

# import melts
melts_filenames = dir(path = directory_path, pattern = '*Melt.txt', full.names = T)
melts <- lapply(melts_filenames, parse_CD_file)
names(melts) <- lapply(melts_filenames, function(x) substr(x, 25, nchar(x) - 9))

# plot melts and make a table of apparent Tms
melt_output <- lapply(melts, plot_melt)

plot_list <- lapply(melt_output, function(x) x[[1]])
pdf('Supplemental_Figures/CD_Melt_Plots.pdf', height = 1.5, width = 2, onefile = TRUE)
invisible(print(plot_list))
dev.off()

Tm_list <- lapply(melt_output, function(x) x[[2]])
Tm_df <- do.call(rbind.data.frame, Tm_list)
colnames(Tm_df) <- c('Gsp1 variant','Apparent Tm')
write.csv(Tm_df, 'Supplementary_Data_Tables/Supp_Table9_CD_apparent_Tm.csv')

# import scans
scans_filenames = dir(path = directory_path, pattern = '*Scan.txt', full.names = T)
scans <- lapply(scans_filenames, parse_CD_file)
names(scans) <- lapply(scans_filenames, function(x) substr(x, 25, nchar(x) - 9))

# plot scans
pdf('Supplemental_Figures/CD_Scan_Plots.pdf',
    height = 1.5, width = 2, onefile = TRUE)
invisible(print(lapply(scans, plot_scan)))
dev.off()







