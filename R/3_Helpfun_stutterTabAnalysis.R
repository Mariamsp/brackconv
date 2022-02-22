
### 4_Helpfun_stutterTabAnalysis
### Aim: Function for step 1.6. Creation of stutter table in script Version1AnalyticalThresholdMinFilter10.R

helpfun_stutterTabAnalysis <- function(inputDF_stutterTable) {
  
  ## Samples and markers vectors:
  #inputDF_stutterTable = initial_df2
  samples_test <- unique(inputDF_stutterTable$SampleID)
  markers <- unique(inputDF_stutterTable$Locus)
  
  ## Empty container. Initialisation of stutterInfo:
  stutterInfoTable <- c()
  excluded_SequencesDF <- c()
  included_SequencesDF <- c()
  ## Iteration through all markers first and all the samples afterwards. It will 
  ## group observations by marker across all the samples:
  for (marker in markers) {
    
    #print(marker) ## If we print here we will be able to check how the R loop goes 
    # through the markers and we can see where it stops, if it stops at some point
    for (sample in samples_test) {
      ## Subset the data into the sample and the marker:
      newdata <- subset(inputDF_stutterTable,inputDF_stutterTable$SampleID == sample & inputDF_stutterTable$Locus == marker)
      #newdata <- subset(initial_df2,initial_df2$SampleID == "Sample43" & initial_df2$Locus == "D19S433")
      
      if (nrow(newdata)==1) next # Eg. marker PENTA E for sample A4PC, it has one sequence only, so there are no other 
      # sequences to compare with for this marker
      
      ## 3.1. Input vector for helpfunction
      seqs2 <- rbind(newdata$Reads, newdata$sequence_type, newdata$SequenceID, newdata$AT_contentVar)
      colnames(seqs2) <- newdata$Forward_Strand_Bracketed_form
      
      
      ## 3.2. Create stutter table info 
      stutterInfo <- getMotifStutterInfo(seqs2)
      
      
      ## 3.3. Dat frame with excluded sequences:
      excluded_Sequences <- subset(newdata, !(newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq1 | newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq2))
      
      included_Sequences <- subset(newdata, newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq1 | newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq2)
      
      ####################### End of Test area Check sequences that are not included in stutter table ########################################################################################
      
      if(length(stutterInfo)==0) next
      
      if(nrow(stutterInfo)==0) next ## Avoids stopping the loop when it encounters an empty 
      ## stutterInfo data frame. "Error in data.frame(..., check.names = FALSE) : 
      ## arguments imply differing number of rows: 1, 0. stutterInfo is empty and
      ## it does not want to append the empty data frame to stutterInfoTable
      
      ## 3.3. Merge all the matrices (as many as markers should be)
      stutterInfoTable <- rbind(stutterInfoTable, cbind(sample, marker, stutterInfo))
      #View(stutterInfoTable)
      
      # Merge data frames with included sequences for checking later:
      included_SequencesDF <- rbind(included_SequencesDF,included_Sequences)
      
      # Merge data frames with excluded sequences for checking later:
      excluded_SequencesDF <- rbind(excluded_SequencesDF, excluded_Sequences)
    }
  }  
  
  ## 3.4. Conversion into data frame:
  stutterInfoTable <- data.frame(stutterInfoTable, stringsAsFactors = FALSE)
  colnames(stutterInfoTable) <- c("Samples", "Locus", "Seq1", "Seq2","MotifDifference", "Seq1 counts", "Seq2 counts", 
                                  "sequence_TypeSeq1", "sequence_TypeSeq2", "SequenceIDSeq1","SequenceIDSeq2","ATcontentSeq1", "ATcontentSeq2", "Stutter proportion")
  
  ## 3.5. Create stutter_type variable for re-coding the stutter type (MotifDifference - stutter_type):
  
  stutterInfoTable <- stutterInfoTable %>%
    mutate(
      stutter_type = case_when(
        grepl("^(?!.*/.*).*\\-1$", stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "bwd-1",
        grepl("^(?!.*/.*).*\\-2$", stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "bwd-2",
        grepl("^(?!.*/.*).*\\+1$", stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "fwd+1",
        grepl("^(?!.*/.*).*\\+2$", stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "fwd+2",
        grepl("^[[:alpha:]]+:\\+1/[[:alpha:]]+:\\-1$",  stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "fwd+1_bwd-1",
        grepl("^[[:alpha:]]+:\\-1/[[:alpha:]]+:\\+1$",  stutterInfoTable$MotifDifference, perl = TRUE) == TRUE ~ "bwd-1_fwd+1",
        TRUE ~ "Infrequent" ## Includes all those stutters that are not in the above categories
      ))
  
  return(list(stutterInfoTable,included_SequencesDF,excluded_SequencesDF))
}

