
### 2_Helpfun_qualityCheckAllele2
### Aim: Function for step 2 (Quality Check; Allele 2 drop-outs) in script  1_StutterTableMR10_Hb30.Rmd
### In the first part we obtain the summary classification of the sequences into Allele 1 and Allele 2 (coded as 1 for the first and 2 for the second)
### In the second part we do a second analysis only of the allele 2 sequences. We check drop-outs of allele 2 due to removal of sequences because the
### Hb values falls out of the set threshold (inputHb)

helpfun_qualityCheckAllele2 <- function(initial_df2, InputHb) {
  InputHb = Hb_Threshold
  
  # Part 1. Quality check. Check for the correctness of the classification of sequences into Allele 1 /  Allele 2. Check ADO / Locus drop-outs / Drop-ins
  
  # 1.1. Classification into Allele 1
  # True Allele 1----------------------------------------------------------------------------
  parental_allele1_TestThreshold <- aggregate(sequence_type ~ SampleID + Locus,initial_df2, function(x) sum(x==1)) # Sums up number of Allele 1 observations per SampleID and Locus
  isTRUE(parental_allele1_TestThreshold$sequence_type == "1") # Check if for all the rows the value = 1, if false means there are samples/locus rows with more than one allele 1 designation or 0 allele 1 designation
  
  # Summarise the counts in parental_allele1_100:
  allele1_summary <- table(parental_allele1_TestThreshold$sequence_type)
  
  # 1.2. Classification into Allele 2
  # True Allele 2---------------------------------------------------------------------------
  parental_allele2_TestThreshold <- aggregate(sequence_type ~ SampleID + Locus,initial_df2, function(x) sum(x==2)) # Sums up number of Allele 1 observations per SampleID and Locus
  isTRUE(parental_allele2_TestThreshold$sequence_type == "2") # Check if for all the rows the value = 2, if false means there are samples/locus rows with more than one allele 2 designation or 0 allele 2 designation
  
  # Summarise the sequence types in parental_allele2_TestThreshold:
  allele2_summary <- table(parental_allele2_TestThreshold$sequence_type)
  
  
  # Part 2. Further analysis of quality check of Allele 2
  
  outliersParental2_TestThresholdtype0 <- parental_allele2_TestThreshold[which(parental_allele2_TestThreshold$sequence_type == 0),]# Indices of those rows with Allele 2 = 0 
  #outliersParental2_TestThresholddf1type0 <- parental_allele2_TestThreshold[outliersParental2_TestThresholdtype0,] # dim(outliersParental2_TestThresholddf1type0) # Selected rows 
  outliersParental2_TestThresholddf2type0 <- outliersParental2_TestThresholdtype0[,1:2] # dim(outliersParental2_TestThresholddf2type0) # Sample and Locus columns
  
  # We have a subset data frame with all the sequences for the sample/marker where allele 2 has drop out or is non-existent:
  outliersParental2_TestThresholddfFinaltype0 <- initial_df2[unlist(sapply(1:NROW(outliersParental2_TestThresholddf2type0), function(i)
    which(initial_df2$SampleID == outliersParental2_TestThresholddf2type0$SampleID[i] & initial_df2$Locus == outliersParental2_TestThresholddf2type0$Locus[i]))),]

  # Data frame with possible allele 2 drop-outs:
  possibleAllele2_DropOuts <- subset(outliersParental2_TestThresholddfFinaltype0,outliersParental2_TestThresholddfFinaltype0$Reads == outliersParental2_TestThresholddfFinaltype0$Sec_max_reads)
  
  # Build new variable for assigning the multiple types of drop-outs:
  # DropOut_LowHb = Possible allele 2 (index == 2) but with Hb < 30 so it will be categorised as a stutter
  # NoDropOutAllele2 = confirmed Allele 1  / confirmed Allele 2 / confirmed stutter
  possibleAllele2_DropOuts$DropOut_type <- ifelse(possibleAllele2_DropOuts$ProHb <= InputHb, "DropOut_LowHb","NoDropOutAllele2")
                                           
  # Create simple data frame with DropOut_type variable and SequenceID:
  DropOuts_df <- possibleAllele2_DropOuts[,c("SequenceID","DropOut_type")] 
  
  # Add DropOuts to initial_df2 and assign NA the value "NoDropOutAllele2":
  initial_df2 <- merge(initial_df2, DropOuts_df, all.x = TRUE) # LEFT JOIN the data frame with the DropOut variable and the initial df
  initial_df2$DropOut_type[is.na(initial_df2$DropOut_type)] <- "NoDropOutAllele2"
 
  return(list(possibleAllele2_DropOuts,initial_df2))
}


