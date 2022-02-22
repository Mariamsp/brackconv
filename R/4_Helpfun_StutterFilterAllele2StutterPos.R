
### 4_Helpfun_StutterFilterAllele2StutterPos
### Aim: Function for step 4: Remove n-1 Allele 2 and final filter stutter table ###
### Remove parental Allele2 that falls in stutter position n-1  when comparing Allele1 (1) and Allele2 (2)
### Keep only those rows where Seq1 == Allele1 or Allele2 (1 or 2) and Seq2 == Stutter (3)
### With this script we will indicate in the initial / general data frame those Allele 2 that are part of the stutter Table and are in n-1 position. 
### They will be labelled as Allele2_bwd1Position

helpfun_StutterFilterAllele2StutterPos <- function(df1,df_ext,Hb_min, Hb_max) {

  ## 1. Add rownames as ID column for additional filtering:
  df1 <- cbind(ID = rownames(df1), df1)
  
  ## 2. Remove parental Allele2 that falls in stutter position n-1 / n+1 when comparing Allele1 (1) and Allele2 (2):
  
  ### Step1. Add column with Hb proportion variable with the percentage proportion of reads of Seq2 compared to Seq1:
  df1$Hb_proportion <- df1$`Seq2 counts`*100 / df1$`Seq1 counts`  
  
  ### Step2. Construct a data frame that contains the instances we want to filter out later from our main data frame: 
  
  ### We want to filter out those sequences that fulfills the following criteria for removing sequences:
  
  # A) # * Sequence Seq1 is parental Allele 1 and sequence Seq2 is parental Allele 2
  # * The difference found between both sequences is of the class: stutter type n-1. In which case the parental allele 2 is in a stutter position,
  # and this is a confounder factor for determining number of reads from the true parental and number of reads from a potential stutter in that position
  # * We also get a more realistic data set when framing those parental alleles 2 with a low Heterozygous balance (high intralocus imbalance), and we determine that the 
  # optimal Heterozygous balance lower an upper threshold  for filtering out sequences is 30% (given that in the 4 Positive control assay there
  # were the two true allele differences in read number of down to 31%) and we go up to 50% (discussion with ?yvind). Previously, we define a minimum
  # heterozygous balance of 30%, in which case, parental alleles will be only defined when Hb_proportion > 30% and no observations will be found under this
  # value. This is the reason why we filter >30% and <50%.
  
  # B) ### Code for filtering (dplyr): # Check point ### !!!! Check if there is any Hb < 30% for these sequences /// The key stutter to look
  ### at is the n-1 stutter (the most frequent type). And we see that the stutter proportions are high approx 0.3 (30%)
  df2 <- df1 %>% 
    group_by(Samples, Locus) %>%
    filter(sequence_TypeSeq1 == 1 & sequence_TypeSeq2 == 2 & Hb_proportion > Hb_min & Hb_proportion  < Hb_max & stutter_type == "bwd-1")
  
  # C) Generate a new variable in initial_df2 and classify Allele 2 sequences in n-1 position as Allele2_stutter / real Allele 2 into Allele2_real:
  df_ext$Allele2_type[df_ext$SequenceID %in% df2$SequenceIDSeq2 & df_ext$sequence_type == 2 & df_ext$DropOut_type == "NoDropOutAllele2"] = "Allele2_bwd1Position"
  df_ext$Allele2_type[!df_ext$SequenceID %in% df2$SequenceIDSeq2 & df_ext$sequence_type == 2 & df_ext$DropOut_type == "NoDropOutAllele2"] = "Allele2"
  df_ext$Allele2_type[is.na(df_ext$Allele2_type)] <- "Not_Allele2" # Replace NAs with the label "Not_Allele2" for other sequences that are not Allele2

  # D) ### Create a temporal data frame with rows where Seq1 (origin of stutter) is parental Allele 2 in stutter position n-1 (Observed in allele1StutterPos)
  ### Important: We only consider parental Allele 2 with Hb between 30% and 50% parental Allele1 and falling in stutter position, because for these instances
  ### sequences can be a mixed of true parental allele 2 and stutter of parental allele 1.
  stutterInfoTable_SubsetNOT <- df1[unlist(sapply(1:NROW(df2), function(i)
    which(df1$Samples == df2$Samples[i] & df1$Seq1 == df2$Seq2[i]))),]
  
  # E)  Final Data frame. Subset and NOT to include rows from stutterInfoTable_SubsetNOT:   
  df_final <- df1[!df1$ID %in% stutterInfoTable_SubsetNOT$ID,]
  
  ## 3. Keep only those rows where Seq1 == Allele1 or Allele2 (1 or 2) and Seq2 == Stutter (3):
  df_final <- df_final[df_final[,9] == 1 & df_final[,10] == 3 | df_final[,9] == 2 & df_final[,10] == 3  ,, drop = FALSE]
  # View(stutterInfoTable_Final)
  
  return(list(df2,df_final,df_ext))
  
}

