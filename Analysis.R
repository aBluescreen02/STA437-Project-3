library(tidyverse)

df <- read.csv('sdss_data.csv')

cat("Objects:", nrow(df), ", Variables:", ncol(df), "\n")
table(df$class)

# Photometric features for multivariate analysis
photo_cols <- c('u', 'g', 'r', 'i', 'z')
X <- as.matrix(df[, photo_cols])
y <- df$class
