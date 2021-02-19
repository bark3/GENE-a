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

# Calculate Delta C(t) and RQ for the each patient status
add_dctrq <- function(dt){
  dt["delta_ct"] <- dt['C(t)']-dt['Ctrl C(t)']
  dt["RQ"] <- 2^-(dt['delta_ct'])
  
  return(dt)
}

failed_donor <- add_dctrq(failed_donor)
luad <- add_dctrq(luad)
lusc <- add_dctrq(lusc)

# Calculate failed donor vs agg lung cancer
names(failed_donor) <- str_replace_all(names(failed_donor), c(" " = "." , "," = "" ))
fdVsAgControl <- failed_donor %>% group_by(Target.Name) %>% summarize(avgRQ = mean(RQ))

