
#### Aim of the script: Analyse STR-MPS stutter artefacts

# Parameters: Minimum read count is 10 (T = 10) reads and Hb (Heterozygote balance) > 0.3

## 1. Load packages, source helpfunctions and export and create the initial data frame:

# 1.1. Load packages:
library(ggplot2)
library(tidyr)
library(knitr)
library(kableExtra)
library(stringi)
library(ggpubr)
library(readxl)
library(stringr)
library(lattice)
library(betareg)
library(viridis)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(gtools)
theme_set(theme_pubr())

# 1.2. Load data and create initial dataframe:

# Set working directory and path for reading input file:
setwd("C:/your_path/")
path <- "C:/your_path/"

# Read input file:
initial_df <- read.table(paste0(path,"datafile.txt"), sep="\t", header = TRUE)

# Create new variable with only digits from SampleID, create Run variable with values 3 and 4:
initial_df <- initial_df %>%
  mutate(SampleNumber = as.numeric(gsub("[^[:digit:]]", "", SampleID))) %>%
  select(-SampleNumber)

# File for creating variable STR_complexity (step number 8):
STR_complexity_df <- read_excel("STR_Complexity.xlsx", sheet = 2)

# 1.3. Load data and create initial dataframe:

source("getMotifStutterInfo.R")
source("getMotifReps.R")
source("1_Helpfun_InitialDataframe.R")
source("2_Helpfun_qualityCheckAllele2Version2.R")
source("3_Helpfun_stutterTabAnalysis.R")
source("4_Helpfun_StutterFilterAllele2StutterPos.R")
source("5_Helpfun_AddVariables.R")
source("6_Helpfun_stutterCriteria.R")
source("7_Helpfun_LUSCalculation.R")
source("8_Helpfun_ComplexATContVar.R")

## 2. Initial working data frame with minimum reads of 10 (T = 10) and Hb = 30. Definition of true alleles.

# Parameters
MinimumReads <- 10
Hb_Threshold <- 30
# Initial working data frame
dfList_initialdf <- setNames(helpfun_initialDataframe(initial_df,MinimumReads,Hb_Threshold), c("initial_df2","noise_underMinimumReads")) # Define the names of the data frames inside of the list
list2env(dfList_initialdf,.GlobalEnv) # list2env: From a named list x, create an environment containing all list components as objects, or "multi-assign" from x into a pre-existing environment
# Allele 2 drop outs quality check
dfList_allele2DropOuts <- setNames(helpfun_qualityCheckAllele2(initial_df2, Hb_Threshold), c("possibleAllele2_DropOuts","initial_df2"))
list2env(dfList_allele2DropOuts,.GlobalEnv)
# Calculate noise sequences frequencies by reads:
table(noise_underMinimumReads$Reads)
nrow(noise_underMinimumReads)

## 3. Generate the stutter table:

# Stutter table
dfListStutterTab <- setNames(helpfun_stutterTabAnalysis(initial_df2), c("stutterInfoTable","included_SequencesDF","excluded_SequencesDF")) # Define the names of the data frames inside of the list
list2env(dfListStutterTab,.GlobalEnv)
## 4.Remove Allele 2 in n-1 position and final filter stutter table:

# Parameters
Hb_min <- Hb_Threshold
Hb_max <- 50
# Final stutter table
dfListFinalStutterTab <- setNames(helpfun_StutterFilterAllele2StutterPos(stutterInfoTable,initial_df2,Hb_min, Hb_max), c("allele1StutterPosMin_MaxHb","stutterInfoTable_Final", "initial_df2"))
list2env(dfListFinalStutterTab,.GlobalEnv)

## 5. Add variables to the final data frame for summarising the data:

# Final stutter table with new variable:
dfList_Variables <- setNames(helpfun_AddVariables(initial_df2,included_SequencesDF,excluded_SequencesDF), c("TotalExcluded","TotalIncluded","ExludedPriori","initial_df2"))
list2env(dfList_Variables,.GlobalEnv)

## 6. Additional filtering of stutters.

# Obtains unambiguous stutters, frequent and infrequent, where:
# a. Unambiguous stutters are those stutters that can only be originated from one parental allele.
# b. Frequent stutters are the most observed stutter types: n-1, n+1, n-2, n+2 and *n0?*
# c. Infrequent stutters are complicated stutters eg. n-3 or n-4

# Stutter filtering. Obtains new data frames:
dfList_StutterFiltering <- setNames(helpfun_StutterFilterAdd(stutterInfoTable_Final), c("ambiguous_stutter","unambiguous_stutter", "infrequent_stutters", "frequent_stutters")) # Define the names of the data frames inside of the list
list2env(dfList_StutterFiltering,.GlobalEnv) # list2env: From a named list x, create an environment containing all list components as objects, or "multi-assign" from x into a pre-existing environment

## 7. Add variables:

# LUS_Rep: contains LUS calculation
# parentMotif_length1: first parental uninterrupted stretch (PTUS1)
# parentMotif_length2: second parental uninterrupted stretch (PTUS2)

final_StutterTable <- getLUS_lengthReps(frequent_stutters)

## 8. Add variables Complexity and AT-content:

final_StutterTable <- helpfun_ComplexATContVar(final_StutterTable)


