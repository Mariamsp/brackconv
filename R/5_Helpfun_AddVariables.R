
### 5_Helpfun_AddVariables
### Aim: Function for step 4: Remove n-1 Allele 2 and final filter stutter table ###
### Remove parental Allele2 that falls in stutter position n-1  when comparing Allele1 (1) and Allele2 (2)
### Keep only those rows where Seq1 == Allele1 or Allele2 (1 or 2) and Seq2 == Stutter (3)
### With this script we will indicate in the initial / general data frame those Allele 2 that are part of the stutter Table and are in n-1 position. 
### They will be labelled as Allele2_bwd1Position

helpfun_AddVariables <- function(inputDF_stutterTable, in_df, ex_df1) {
 
  # PART 1. Check point of excluded and/or included sequences in the stutter table (stutter analysis):
  # Sequences excluded from the stutter table (a posteriori script) + sequences excluded from stutter analysis (a priori script)
  TotalExcluded <- inputDF_stutterTable[!inputDF_stutterTable$SequenceID %in% in_df$SequenceID,] 
  
  # Equals to sequences included in the stutter analysis
  TotalIncluded <- inputDF_stutterTable[!inputDF_stutterTable$SequenceID %in% TotalExcluded$SequenceID,]
   
  # Sequences excluded from stutter analysis (a priori script)
  ExludedPriori <- TotalExcluded[!TotalExcluded$SequenceID %in% ex_df1$SequenceID,]
  
  # PART 2. Create variable inclusion_stutterTab that will tell us whether a sequence has been included or not in the stutter analysis:
  
  inputDF_stutterTable$inclusion_stutterTab <- ""
  inputDF_stutterTable$inclusion_stutterTab[inputDF_stutterTable$SequenceID %in% in_df$SequenceID & inputDF_stutterTable$sequence_type == 1] = "ParentalAllele1"
  inputDF_stutterTable$inclusion_stutterTab[inputDF_stutterTable$SequenceID %in% in_df$SequenceID & inputDF_stutterTable$sequence_type == 2] = "ParentalAllele2"
  inputDF_stutterTable$inclusion_stutterTab[inputDF_stutterTable$SequenceID %in% in_df$SequenceID & inputDF_stutterTable$sequence_type == 3] = "Stutter"
  inputDF_stutterTable$inclusion_stutterTab[inputDF_stutterTable$SequenceID %in% ex_df1$SequenceID] = "Excluded_stutterTable"
  inputDF_stutterTable$inclusion_stutterTab[inputDF_stutterTable$SequenceID %in% ExludedPriori$SequenceID] = "Excluded_stutterTable"

  return(list(TotalExcluded,TotalIncluded,ExludedPriori,inputDF_stutterTable))
}
