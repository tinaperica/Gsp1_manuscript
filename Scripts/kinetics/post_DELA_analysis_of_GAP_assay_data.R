library(tidyverse)
library(ggrepel)
library(pracma)
outfiles <- "Data/RanGAP_assay/"

#### analyse and parse parameters from DELA
parameters <- read_tsv("Data/RanGAP_assay/20200720_parameters_summary.txt") %>% 
  select(data, kcat, Km, "conc" = `[S]`, basS) %>% 
  separate(col = data, into = c("date", "sample", "well", "input_conc", "GAP_conc", "sensor_conc"), sep = "-", convert = T, remove = T) %>% 
  mutate("aprox_loading_eff" = conc/input_conc) %>% 
  separate(col = sample, into = c("PE", "mutant")) %>% 
  arrange(mutant, date) %>% 
  mutate('individual_kcat_Km' = kcat/Km)
interval_table <- read_tsv('Data/RanGAP_assay/interval_table.txt') 
read_tsv("Data/RanGAP_assay/20200720_parameters_summary.txt")%>% 
  inner_join(., interval_table, by = c('data' = 'condition')) %>% 
  write_tsv('Data/RanGAP_assay/parameters_with_interval.txt')
avg_parameters <- parameters %>% 
  group_by(mutant) %>% 
  summarize("mean_kcat" = mean(kcat), "mean_Km" = mean(Km), 
            "kcat_sd" = sd(kcat), "Km_sd" = sd(Km), 
            "kcat_se" = std_err(kcat), "Km_se" = std_err(Km), 
            'mean_kcat_Km' = mean(individual_kcat_Km),    #### take the mean of individual kcat/Km
            'sd' = sd(individual_kcat_Km),  ### this sd reports the sd of the kcat/Km ratio between individual measurements
            'se' = std_err(individual_kcat_Km),
            "combined_sd" =  (sqrt( (kcat_sd/mean_kcat)^2 + (Km_sd/mean_Km)^2 )) * mean_kcat_Km,  ### but use the more strict way for getting the std.dev
            "combined_se" =  (sqrt( (kcat_se/mean_kcat)^2 + (Km_se/mean_Km)^2 )) * mean_kcat_Km,  ### but use the more strict way for getting the std.dev
            # 'log_kcat_over_Km' = log(mean_kcat_Km),
            # 'log_sd' = sd/mean_kcat_Km
            ) %>% 
  ungroup() %>% 
  select(mutant, 'mean_kcat' = mean_kcat, 'mean_Km' = mean_Km, mean_kcat_Km, everything())

GAP_parameters <- parameters %>% 
  inner_join(., avg_parameters, by = 'mutant') %>% 
  select(-date, -PE, -well, -input_conc, -aprox_loading_eff, -basS)
write_tsv(GAP_parameters, path = file.path(outfiles, "GAP_kinetics_MichaelisMenten_parameters.txt"))
# avg_parameters %>% 
#   select(mutant,'kcat' = mean_kcat, 'Km' = mean_Km, everything()) %>% 
#   write_tsv(., path = file.path(outfiles, "GAP_MM_parameters_clean.txt"))
#### kcat vs Km scatterplot
avg_parameters %>% 
  ggplot(aes(x = mean_kcat, y = mean_Km, label = mutant)) + 
  geom_point(aes(size = sd), alpha = 0.5) +
  geom_point(data = avg_parameters[avg_parameters$mutant == "WT",], 
             aes(size = sd, color = "red"), alpha = 0.5) +
  geom_text_repel() +
  xlab(label = "kcat / s-1") +
  ylab(label = "Km / uM") +
  labs(size = "standard deviation\nof kcat/Km ") +
  ggtitle("Michaelis Menten parameters of SpGAP facilitated GTP hydrolysis of Gsp1 mutants")
#ggsave(filename = file.path(outfiles, "kcat_vs_Km_scatterplot.pdf"), height = 7, width = 12)
gap_table <- avg_parameters %>% 
  select(mutant, mean_kcat, kcat_se, mean_Km, Km_se, mean_kcat_Km, se) 
write_tsv(gap_table, path = 'Revisions/Supplementary_Files/Supplementary_File_1/GAP_table_raw.txt')
write_tsv(gap_table, path = '~/Desktop/test.txt')

parameters <- parameters %>%
  mutate('row' = substr(well, 1,1)) %>%
  mutate('exp' = str_c(date, row, mutant))
intrinsic_hydrolysis_data <- read_tsv('Data/RanGAP_assay/intrinsic_data/revisions_intrinsic_hydrolysis_data.txt') %>% 
  mutate('exp' = str_c(date, row, mutant))

intrinsic_replicates_to_remove <- intrinsic_hydrolysis_data %>%   #### remove the intrinsic measurements for many of the test wt runs, where final concncetration was probably not well estimated
  filter(mutant == 'WT' & !(exp %in% unique(parameters$exp))) %>% pull(exp) %>% unique()

intrinsic_hydrolysis_data <- intrinsic_hydrolysis_data %>% 
  filter(! exp %in% intrinsic_replicates_to_remove)

intrinsic_hydrolysis_data %>% 
  filter(mutant %in% c('T34Q', 'T34E')) %>% 
  ggplot(aes(x = Time, y = fluorescence, color = condition)) + 
  geom_point() +
  facet_grid(~mutant)
### remove the intrinsic measurements from 20190606 and the T34Q one from 20190606 because that was a faulty batch of sensor
### also so far remove PE19_R78K from 20190315 because it looks way to high but need to add more repeats to confirm
intrinsic_hydrolysis <- intrinsic_hydrolysis_data %>% 
  select(condition, sample, mutant, final_product_conc, date, slope, intercept, rel_rate) %>% 
  filter(! (date == '20190606' | 
              (date == '20190605' & sample == 'PE13_T34Q') |
              (date == '20190315' & (sample == 'PE19_R78K')) |
              (date == '20200531' & sample == 'PE1_WT') |
              (date == '20190715' & (sample == 'PE2_T34A' | sample == 'PE1_WT')))) %>% 
  unique() %>% 
  arrange(rel_rate)
intrinsic_hydrolysis_summary <- intrinsic_hydrolysis %>% 
  mutate('rel_rate_log' = log(rel_rate)) %>% 
  select(sample, mutant, condition, rel_rate, rel_rate_log) %>% 
  unique() %>% 
  group_by(sample, mutant) %>% 
  summarize('mean_rel_rate' = mean(rel_rate), 'sd_rel_rate' = sd(rel_rate), 'se_rel_rate' = std_err(rel_rate),  
            'mean_log_rel_rate' = mean(rel_rate_log),
            'sd_log_rel_rate' = sd_rel_rate/mean_rel_rate,
            'se_log_rel_rate' = se_rel_rate/mean_rel_rate
            )

intrinsic_hydrolysis %>% 
  inner_join(., intrinsic_hydrolysis_summary, by = 'mutant') %>% 
  select(mutant, final_product_conc, rel_rate, mean_rel_rate, sd_rel_rate, se_rel_rate, mean_log_rel_rate, sd_log_rel_rate, se_log_rel_rate) %>% 
  arrange(mean_rel_rate) %>% write_tsv(path = str_c(outfiles, '/intrinsic_hydrolysis.txt'))

write_tsv(intrinsic_hydrolysis_summary, file.path(outfiles,'intrinsic_hydrolysis_parameters.txt')) 


