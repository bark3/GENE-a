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
names(luad) <- str_replace_all(names(luad), c(" " = "." , "," = "" ))
names(lusc) <- str_replace_all(names(lusc), c(" " = "." , "," = "" ))

agCan <- bind_rows(luad,lusc)

miRNAs <- unique(agCan["Target.Name"])


# Loop to calculate t.test for each micro-RNA
size = dim(miRNAs)

i <- 1
test <- tibble(pValue = 1.1:size[1], muX = 1.1:size[1], muY =1.1:size[1])
while (i <= size[1]){
  xGroup <- filter(failed_donor, Target.Name == miRNAs[[i,1]])
  x <- xGroup["RQ"]
  yGroup <- filter(agCan, Target.Name == miRNAs[[i,1]])
  y <- yGroup["RQ"]
  iM <- t.test(x=x,y=y)
  test[i,1] <- unname(iM[[3]])
  test[i,2] <- unname(iM[[5]][[1]])
  test[i,3] <- unname(iM[[5]][[2]])
  i = i+1
}

testResults <- bind_cols(miRNAs, test)

testResults["foldChange"] <- testResults["muY"]/testResults['muX']

accept <- filter(testResults, pValue <.05, foldChange > 2 | foldChange < .5)


  
  
  
  
  
  