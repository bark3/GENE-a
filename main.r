# GENE-e 
# 02/13/2021
# Courtney Barkley and Zach Barkley
# Description: Analyzing 12K data

# Set environment
library(readxl)
library(tidyverse)

# Read in lung patient samples
failed_donor <- read_xlsx("/Users/zacharybarkley/Documents/Courtney R Project/Raw Data/LC Disparities Data Analysis-Zach edit.xlsx")
luad <- read_xlsx("/Users/zacharybarkley/Documents/Courtney R Project/Raw Data/LC Disparities Data Analysis-Zach edit.xlsx", sheet = 2)
lusc <- read_xlsx("/Users/zacharybarkley/Documents/Courtney R Project/Raw Data/LC Disparities Data Analysis-Zach edit.xlsx", sheet = 3)

# Clean
failed_donor <- select(failed_donor, -c('Delta C(t)', 'RQ', 'Avg RQ', 'STDEV', 'SEM'))
luad <- select(luad, -c('Delta C(t)', 'RQ', 'Avg RQ', 'STDEV', 'SEM'))
lusc <- select(lusc, -c('Delta C(t)', 'RQ', 'Avg RQ', 'STDEV', 'SEM'))

# Calculate Delta C(t)
failed_donor["delta_ct"] <- failed_donor['C(t)']-failed_donor['Ctrl C(t)']
failed_donor["RQ"] <- 2^-(failed_donor['delta_ct'])

boogy woogy woogy