
### 1_Helpfun_InitialDataframe.R
### Aim: Function for step 2 in script 1_StutterTableMR10_Hb30.Rmd. Build initial data frame with filters: 
### set minimum reads for removing noise (non-stutters) / Set Hb for definition of true alleles.

### Additional Comments:
### With this strategy, filtering of the minimum reads is applied after the sequence classification (eg. Allele 1, etc...). ProHb's are build before application of this filtering, 
### which will lead to extreme ProHb values for those cases where we can only obtain one of the parental alleles (either due to HMZ or drop out?). 
### Check table: Version2AT_extremeHbs_example.csv

helpfun_initialDataframe <- function(df, MinimumReads,Hb_Threshold) {
  # a. Calculate per locus (depth of coverage?) the number of mapped reads and append a new column to original df:
  df2 <- transform(df, perlocus_reads = ave(df$Reads, list(df$SampleID,df$Locus), FUN = sum))
  
  # b. Calculate parameters for analysis (Per Sample and per Marker):
  
  # First Max Value and second max value. Percent reads, Heterozygote balance calculation and AT:
  df2 <- df2 %>% 
    mutate(SeqForw_length = nchar(Forward_Strand_Sequence)) %>%
    group_by(SampleID, Locus) %>%
    arrange(desc(Reads), .by_group= TRUE) %>% # within each group (group is define by a particular combination of Sample with Locus)
    mutate(
      Index = row_number(), # Creates a variable with an index per group, where the group is every possible combination of sample with marker (SampleID with Locus)
      Max_reads = Reads[[1]], # First max read counts
      Sec_max_reads = ifelse(length(Reads) == 1, NA, Reads[[2]]), # Second max read counts (It takes always the second largest value of the Reads column, if nrow reads column = 1, then we do not have a sec max reads, add NA)
      Percent_reads = Reads*100/Max_reads, # Percentage of sequence reads compared to Max_reads variable (100%)
      ProHb = (Sec_max_reads/Max_reads)*100 # Percent proportion of Sec_max_reads and Max_reads (Sec_max_reads / Max_reads = X/100)
    )
  
  # Comment HMZ alleles: If we have one single allele (eg. Allele 1), we are uncertain about the second alleles. Cases:
    # 1. Homozygous allele, Max_reads = 1/2 Max_reads and Sec_max_reads = 1/2 Max_reads. Then, ProHb would be 100% (Hb = 1). But we then, can have possibility 2.
    # 2. If we assign the second allele by using the ProHb value, we can fall in the possibility that the ProHb (high unbalance) is beyond the limit and therefore, the second allele is excluded: Drop-out    
    # 3. Mutation in primer site, no amplification of the allele.
  
  # c. Add unique identification label for each instance by creating a new variable with the row number (I want the identification of the samples after 
  # I have ordered the sequences within each group)
  df2$SequenceID <- as.numeric(rownames(df2)) # We need the ID variable to be numeric
  # Move the column (last position) to the start:
  df2 <- df2 %>% dplyr::select(SequenceID, everything())

  # d. Add new column with sequence/allele type. Classification of Allele 1, Allele2, stutter and noise:
  
  # Comments:
  # Intralocus balance (As a rule of thumb and validation results, 60% minimum intralocus balance, under this percentage allele observation can 
  # be due to stochasticity). In this case we will use the following criteria to avoid losing highly imbalanced parental sequences during the analysis:
  # Reduction from 60% to 30% the minimum intralocus balance percentage (assessment of 4 Positive Controls gave a minimum intralocus balance of 31% for capturing all parental sequences)
  # Noise definition: PCR errors (not leading to stutter) / Sequencing errors 
  
  # We define true alleles with Hb > 30 as obtained in previous analysis. We have to consider all the possible cases and apply rules for the cases. Can we build a general rule? 
  df2 <- df2 %>%
    group_by(SampleID, Locus) %>%
    mutate(
      sequence_type = case_when(
        Index == 1 ~ "Allele1", # Accounts for regular HTZ cases and also when Max_reads and Sec_max_reads 
        Index == 2 & ProHb > Hb_Threshold ~ "Allele2", # Accounts for regular HTZ cases and also when Max_reads and Sec_max_reads 
        TRUE ~ "Stutter"
      ))
  
  # e. Write a function for creating the variable with the AT/GC content:
  calculate_ATContent <- function(df2) {
    count_A <- str_count(df2$Forward_Strand_Sequence, "A")
    count_T <- str_count(df2$Forward_Strand_Sequence, "T")
    AT_content <-round(((count_A+count_T) / (df2$SeqForw_length))*100, 2)
    return(AT_content)
  }
  df2$AT_contentVar <- calculate_ATContent(df2) 
  
  # f. Filter for the minimum number of reads that are accepted for analysis This filtering will primarily
  # remove noise non-stutter sequences eg. base substitution sequencing errors.
  
  df3 <- df2[df2$Reads >= MinimumReads,] # Including the reads imposed by MinimumReads
  
  # Build a data frame with the noise that has been removed (Complement)
  noise_underMinimumReads <- df2[df2$Reads < MinimumReads,]
  
  # Convert sequence_type variable into 
  df3$sequence_type <- factor(df3$sequence_type,
                                      levels <- c("Allele1", "Allele2", "Stutter"),
                                      labels <- c(1, 2, 3))

  return(list(df3, noise_underMinimumReads))
}


 
