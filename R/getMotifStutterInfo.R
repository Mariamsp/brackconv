
################################################################################
# HELPFUNCTION

##### Function for creating the final table with stutter information. It will
##### traverse through all the sequences (short bracket format) and compare each
##### sequence to the remaining sequences for stutter determination.

#' @title getMotifStutterInfo
#'
#' @author Oyvind Bleka <oyvble.at.hotmail.com>
#'
#' @description Creates final table with stutter information
#'
#' @details Help function. Takes as input the short bracket format sequences.
#' Compares one sequence with the remaining sequences by extracting and comparing
#' their repetitive motifs. If equal, then number of repetitions of the motifs
#' will be compared and variations will be considered as stutters.
#'
#' NOTE: STUTTER PRODUCTS ARE ASSUMED TO NOT EXCEED THE COMPARING SEQUENCE
#'
#' @param seqBrack matrix with STRs sequence converted into bracket format, reads, allele_type,
#' Sequence ID and AT-content
#' @param addSeq1Coverage TRUE condition for adding Seq1 read counts ("coverage") to the
#' final stutter table
#' @param addSeq2Coverage TRUE condition for adding Seq2 read counts to the final
#' stutter table
#' @param addStutterProportion TRUE condition for adding the stutter proportion to the
#' stutter table. Proportion of read counts Seq2 compare to the total (Seq1 + Seq2)
#'
#' @return data frame with sequences Seq1 and Seq2, number of different motifs (motifDifference),
#' read counts (coverage) for Seq1 and Seq2, and stutter proportion.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bracketformat <- c("[ATCT]12","[ATCT]10", "[ATCT]11", "[ATCT]9")
#' reads <- c(869, 60, 28, 23)
#' allele_type <- c(1,2,3,3)
#' Sequence_ID <- c(227760, 227761, 227762,227763)
#' AT_contentVar <- c(75.00, 33.33, 64.00, 77.33)
#' seqBrack <- rbind(reads, allele_type, Sequence_ID, AT_contentVar)
#' colnames(seqBrack) <- bracketformat
#' getMotifStutterInfo(seqBrack)}

#-------------------------------------------------------------------------------
getMotifStutterInfo = function(seqBrack,addSeq1Coverage=TRUE,addSeq2Coverage=TRUE,addSeq1allele_type= TRUE,
                               addSeq2allele_type=TRUE,addSeq1SequenceID=TRUE,addSeq2SequenceID=TRUE,
                               addSeq1ATcontent=TRUE,addSeq2ATcontent=TRUE,addStutterProportion=TRUE) {

  #Comparing two sequences and do following steps:
  #1) For each sequence we extract motif info from bracket format
  #2) For each aligned block we compare the motif string between the sequence, and if identical, we compare the repeats
  #Input:
  #seqBrack: sequence vector with bracket info

  colname = colnames(seqBrack)
  coverage = NULL
  if(!is.null(colname)) {
    coverage = seqBrack[1,]
    allele_type = seqBrack[2,]
    SequenceID = seqBrack[3,]
    ATcontent = seqBrack[4,]
    seqBrack = colnames(seqBrack)
  } #seqBrack is already bracket

  stuttTab = numeric() #init table
  for(xind in seq_len(length(seqBrack))) { #traverse each sequence
    # xind = 1

    #Extracting info from 'point-of-view' sequence (parent allele)
    xVector = strsplit(seqBrack[xind]," ")[[1]] #split seq to list
    nBlocks = length(xVector) #number of blocks

    for(yind in seq_len(length(seqBrack))) { #traverse each other sequence
           # yind = 2

      #SKIPPING COMPARING SEQUENCE:
      if(xind==yind) next #skip since same sequence
      if(!is.null(coverage) && coverage[yind]>=coverage[xind]) next #skip if comparing sequence is STRICTLY LARGER

      #extract info from comparing sequence
      yVector = strsplit(seqBrack[yind]," ")[[1]] #split seq to list
      if(length(xVector) != length(yVector)) next #don't consider comparison if different number of blocks

      #Need to infer difference between two sequences
      diffVector = numeric() #init difference vector
      isComparable = TRUE #indicate whether the sequences can be compared
      for(block in seq_len(nBlocks)) {#loop through each block
        #      block=1
        xmotif = getMotifReps(xVector[block])
        ymotif = getMotifReps(yVector[block])
        if(xmotif[1]!=ymotif[1]) {
          isComparable = FALSE #The two sequences should not be compared
          break #skip if not same motif
        }

        motifdiff = as.integer(ymotif[2]) - as.integer(xmotif[2])  #number of motifs of y compared to x
        if(motifdiff==0) next #don't record if identical number of motifs

        if(motifdiff>0) motifdiff = paste0("+",motifdiff)
        diffVector = c(diffVector, paste0(xmotif[1],":",motifdiff))
      }
      if(!isComparable || length(diffVector)==0) next
      diffVector = paste0(diffVector,collapse="/") #collapse if several differences

      newrow = c(xind,yind,diffVector)
      stuttTab = rbind(stuttTab , newrow) #add table
    }
  }
  if(length(stuttTab)==0) return(NULL)
  colnames = c("Seq1","Seq2","MotifDifference")

  #whether to insert coverage of seq1:
  if(addSeq1Coverage) {
    colnames = c(colnames,"Seq1Counts")
    stuttTab = cbind(stuttTab,NA) #add column
    ncol(stuttTab)
    if(!is.null(coverage)) stuttTab[,ncol(stuttTab)] = coverage[as.integer(stuttTab[,1])]
  }

  #whether to insert coverage of seq2:
  if(addSeq2Coverage) {
    colnames = c(colnames,"Seq2Counts")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(coverage)) stuttTab[,ncol(stuttTab)] = coverage[as.integer(stuttTab[,2])]
  }

  #whether to insert allele_type of seq1:
  if(addSeq1allele_type) {
    colnames = c(colnames,"Seq1Allele_Type")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(allele_type)) stuttTab[,ncol(stuttTab)] = allele_type[as.integer(stuttTab[,1])]
  }

  #whether to insert allele_type of seq2:
  if(addSeq2allele_type) {
    colnames = c(colnames,"Seq2Allele_Type")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(allele_type)) stuttTab[,ncol(stuttTab)] = allele_type[as.integer(stuttTab[,2])]
  }

  #whether to insert SequenceID of seq1:
  if(addSeq1SequenceID) {
    colnames = c(colnames,"Seq1SequenceID")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(SequenceID)) stuttTab[,ncol(stuttTab)] = SequenceID[as.integer(stuttTab[,1])]
  }

  #whether to insert SequenceID of seq2:
  if(addSeq2SequenceID) {
    colnames = c(colnames,"Seq2SequenceID")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(SequenceID)) stuttTab[,ncol(stuttTab)] = SequenceID[as.integer(stuttTab[,2])]
  }

  #whether to insert ATcontent of seq1:
  if(addSeq1ATcontent) {
    colnames = c(colnames,"Seq1ATcontent")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(ATcontent)) stuttTab[,ncol(stuttTab)] = ATcontent[as.numeric(stuttTab[,1])]
  }

  #whether to insert ATcontent of seq2:
  if(addSeq2ATcontent) {
    colnames = c(colnames,"Seq2ATcontent")
    stuttTab = cbind(stuttTab,NA) #add column
    if(!is.null(ATcontent)) stuttTab[,ncol(stuttTab)] = ATcontent[as.numeric(stuttTab[,2])]
  }

  #whether to insert stutter proportions:
  if(addStutterProportion) {
    colnames = c(colnames,"StutterProportion")
    stuttTab = cbind(stuttTab,NA) #add column

    if(!is.null(coverage)) {
      for(row in seq_len(nrow(stuttTab)) ) {
        seq1cov = coverage[as.integer(stuttTab[row,1])]
        seq2cov = coverage[as.integer(stuttTab[row,2])]
        stuttTab[row,ncol(stuttTab)] = seq2cov/(seq1cov + seq2cov) #stutter proportions
        stuttTab
      }
    }
  }

  #Obtain Motif names (not necessary unique)
  stuttTab[,1] = seqBrack[as.integer(stuttTab[,1])]
  stuttTab[,2] = seqBrack[as.integer(stuttTab[,2])]
  colnames(stuttTab) = colnames
  rownames(stuttTab)=NULL

  #Final data frame
  df = data.frame(stuttTab)#,stringsAsFactors = FALSE)
  if(ncol(df)>3) { #if coverage is given
    for(i in 4:ncol(df)) df[,i] = as.numeric(df[,i]) #convert datatype to numeric
    #df$StutterProportion = as.numeric(df$StutterProportion)
  }
  return(df)
}

