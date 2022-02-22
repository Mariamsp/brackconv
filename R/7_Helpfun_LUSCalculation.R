
### 7_Helpfun_LUSCalculation.R ###
### Aim: Function for step 7: generate new variables: LUS, parental motif length 1 and parental motif length 2 ###
### Extract information from the parental sequences for building a beta regression model ###
### Pitfalls of this approach:
### 1. LUS_Rep column: gives the LUS repetition values only when:
    # * max element for Seq1 is equal to Seq2
    # * The motif difference between the stutter and the parental are the same
### 2. parentMotif_length column: gives the repeat length of the parental element that stutters:
    # * find differences between parental and stutter 
    # * Select only those vectors with 1 element (two will be n0)

## getLUSRepeats function (extended version) for creating LUS_reps, parentMotif_length1 and parentMotif_length2 
## variables ------------------------------------------------------------------------------------------------------

getLUS_lengthReps <- function (df) {
  LUS_Rep <- c() # empty vector for LUS reps
  parentMotif_length1 <- c() # empty vector for parental length motifs 1
  parentMotif_length2 <- c() # empty vector for parental length motifs 2
  isLUS <- c()

  # Iterate through every row in stutter table:
  for (i in 1:nrow(df)) { 
    # Preliminary requirements:
    splittedSeqs1 <- strsplit(df$Seq1," ")[[i]] # split up Seq1 into motifs and reps
    splittedSeqs2 <- strsplit(df$Seq2," ")[[i]] # split up Seq2 into motifs and reps
    MotifRepsSeqs1 <- list() # Create empty list for motif reps from Seq1
    MotifRepsSeqs2 <- list() # Create empty list for motif reps from Seq2
    
    # Iterate through every split element Seq1
    for( j in 1:length(splittedSeqs1)) {
      # Seq1
      nMotifReps1 <- sapply(splittedSeqs1,getMotifReps) # apply getMotifReps to all elements of splittedSeqs1 list
      nreps1 <- as.integer(nMotifReps1[2,]) # Number of reps for each motif of the list
      names(nreps1) <- nMotifReps1[1,] # Name of each of the motifs
      MotifRepsSeqs1 [[j]] <- nreps1 # List with motif name + reps
      # Obtain max rep number and associated motif for determining LUS sequence
      MotifRepsSeqs1_vector <- unlist(MotifRepsSeqs1) # Vector with motif name + reps
      max_element_MotifRepsSeqs1 <- which.max(MotifRepsSeqs1_vector) # Index of the max element of the vector MotifRepsSeqs1_vector
      motifName_compareSeq1 <- names(max_element_MotifRepsSeqs1) # Name of the max element of the vector MotifRepsSeqs1_vector
      motifRep_forlusSeq1 <- as.integer(MotifRepsSeqs1_vector[max_element_MotifRepsSeqs1]) # Number of reps of the max element of the vector MotifRepsSeqs1_vector
      
      for (k in 1:length(splittedSeqs2)) {
        # Seq2
        nMotifReps2 <- sapply(splittedSeqs2,getMotifReps)
        nreps2 <- as.integer(nMotifReps2[2,])
        names(nreps2) <- nMotifReps2[1,]
        MotifRepsSeqs2 [[j]] <- nreps2 
        # Obtain max rep number and associated motif for determining LUS sequence
        MotifRepsSeqs2_vector <- unlist(MotifRepsSeqs2)
        max_element_MotifRepsSeqs2 <- which.max(MotifRepsSeqs2_vector) # Index of the max element of the vector MotifRepsSeqs1
        motifName_compareSeq2 <- names(max_element_MotifRepsSeqs2)
        motifRep_forlusSeq2 <- as.integer(MotifRepsSeqs2_vector[max_element_MotifRepsSeqs2])
        
      } # k element loop closure (inner 2)
      
    } # j element loop closure (inner 1)
    
    # Add columns to table with parental motif length info (1 AND 2):
    diff_vector = which(splittedSeqs1 != splittedSeqs2)# Vector containing the indexes of the elements that are different
    df$parentMotif_length1[i] <- ifelse(length(diff_vector) == 1 | length(diff_vector) == 2, as.integer(MotifRepsSeqs1_vector[diff_vector[1]]), 0) # If vector has only 
    # one element (we exclude N0 with this approach) OR vector of differences has two elements (we include N0 stutters) then select by the index (motif that diff should have the same index in parental and in 
    # stutter) the 1st element of the vector and take the integer (number of reps)
    df$parentMotif_length2[i] <- ifelse(length(diff_vector) == 2, as.integer(MotifRepsSeqs1_vector[diff_vector[2]]), 0) # If vector has 2 elements
    # (we include N0 with this approach), then select by the index (motif that diff should have the same index in parental and in 
    # stutter) the second element of the vector and take the integer (number of reps)
    
    # Add LUS column for Seq1:
    df$LUS_Rep[i] <- motifRep_forlusSeq1
    
    # if (max_element_MotifRepsSeqs1 == max_element_MotifRepsSeqs2) {
    #   # Find motifDifference parental - stutter and compare them:
    #   motifDiff_split <- strsplit(df$MotifDifference, ":")[[i]]
    #   motifDiff_stutter <- unlist(motifDiff_split)[1]
    #   df$LUS_Rep [i] <- ifelse (motifName_compareSeq1 == motifDiff_stutter,motifRep_forlusSeq1, 0)  # LUS column can have value 0 if 
    #   # MotifName_compareSeq1 (max reps motif parental) != motifDiff_stutter ( motif dif with stutter taken from var MotifDifference stutter)
    # } else {df$LUS_Rep[i] <- 0} # LUS column can have value 0 if max_element_MotifRepsSeqs1 != max_element_MotifRepsSeqs2 
    
    # Add isLUS variable:
    # This is a binary variable, where:
      # 0 if LUS != parental motif 
      # 1 if LUS == parental motif 1 
    # NB! Only valid for non n0 stutters (comparison is only between LUS and parental motif 1). We compare the length and the motif of the LUS with the parental:
    df$isLUS[i] <- as.numeric(ifelse(df$LUS_Rep[i] == df$parentMotif_length1[i] & motifName_compareSeq1 == names(MotifRepsSeqs1_vector[diff_vector[1]]), 1, 0))
    
  } # Outer i element loop closure
  return(df)
} # Function closure

