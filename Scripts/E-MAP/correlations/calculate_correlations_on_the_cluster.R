remove(list = ls())
options(stringsAsFactors = F)

weighted_pearson <- function (x, y, w) {
  df <- data.frame(x, y, w)
  df <- df[df$w != 0, ]
  x <- df$x
  y <- df$y
  w <- df$w
  ### Compute the weighted means
  mean_x <- sum(x * w) / sum(w)
  mean_y <- sum(y * w) / sum(w)
  ### Compute the weighted variance
  vx <- sum( w * (x - mean_x)^2 ) / sum(w)
  vy <- sum( w * (y - mean_y)^2 ) / sum(w)
  ### Compute the covariance
  vxy <- sum( w * (x - mean_x) * (y - mean_y) ) / sum(w)
  # Compute the correlation
  w_correlation <- (vxy / sqrt(vx * vy))
  return(w_correlation)
}
weighted_pearson_wrapper <- function( score_x, score_y, weight_x, weight_y) {
  ## pick pairwise observable
  score_df <- data.frame(score_x, score_y)
  weight_df <- data.frame(weight_x, weight_y)
  weight_df[, "max"] <- apply(weight_df[,], 1, FUN = function(x) {max(x, na.rm = T)})
  weight_df[, "min"] <- apply(weight_df[,], 1, FUN = function(x) {min(x, na.rm = T)})
  weight_df[, "mean"] <- apply(weight_df[,], 1, FUN = function(x) {mean(x, na.rm = T)})
  df <- cbind(score_df, weight_df)
  df <- df[complete.cases(df),]
  x <- df$score_x
  y <- df$score_y
  normal_pearson <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  corr_test <- cor.test(x, y, alternative = 'two.sided', method = 'pearson')
  two_sided_p_value <- corr_test$p.value
  corr_test <- cor.test(x, y, alternative = 'g', method = 'pearson')
  greater_p_value <- corr_test$p.value
  weighted_min_pearson <- weighted_pearson(x = x, y = y, w = df$min)
  #weighted_max_pearson <- weighted_pearson(x = x, y = y, w = df$max)
  #weighted_mean_pearson <- weighted_pearson(x = x, y = y, w = df$mean)
  values <- c(normal_pearson, two_sided_p_value, greater_p_value, weighted_min_pearson)
  return(values)
}


task_n <- as.numeric( Sys.getenv( "SGE_TASK_ID" ) ) #
taskfilename <- paste0("all_correlations_task_info/", task_n, "_task_info.RData")
outputfilename <- paste0("20190403_", task_n, "_correlations.txt")
load(taskfilename)
load("spitzemapko_for_corr.rda")

output_file_path <- file.path("Output/correlations", outputfilename)

correlations_df <- data.frame("query_uniq1" = task.info[["pairs"]][1, ],
                              "query_uniq2" = task.info[["pairs"]][2, ], 
                              "pearson" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double"),
                              "two_sided_p_value" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double"),
                              "greater_p_value" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double"),
                              "min_weighted_pearson" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double")
                              #"max_weighted_pearson" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double"),
                              #"mean_weighted_pearson" = vector(length = length(task.info[["pairs"]][1, ]), mode = "double")
)

for ( p in seq_along(task.info[["pairs"]][1,]) ) {
  query1 <- task.info[["pairs"]][1, p]
  query2 <- task.info[["pairs"]][2, p]
  spitzemap.query1 <- spitzemapko_for_corr[spitzemapko_for_corr[["query_allele_name"]] == query1, ]
  spitzemap.query2 <- spitzemapko_for_corr[spitzemapko_for_corr[["query_allele_name"]] == query2, ]
  merged.data <- merge(spitzemap.query1, spitzemap.query2, by = "array_ORF")
  merged.data <- merged.data[complete.cases(merged.data),]
  merged.data <- merged.data[complete.cases(merged.data),]
  values <- round(weighted_pearson_wrapper
                  (score_x = merged.data$score.x, score_y = merged.data$score.y,
                    weight_x = merged.data$weight.x, weight_y = merged.data$weight.y), 4)
  correlations_df[p, 3:6] <- values
}

write.table(correlations_df, file = output_file_path, sep = "\t", quote = F, row.names = F)

