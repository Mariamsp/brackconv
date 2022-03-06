
#' @title initialStutterTable
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Help function for obtaining the initial table with the stutter sequences
#'
#' @details Obtains the initial stutter table where pairs of parental and stutter sequences are
#' stored for further analysis. In parallel, two tables with the included and excluded sequences are returned.
#'
#' @param inputDF_stutterTable Input data frame with all parental and stutter sequences
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#'
#' @return Three data frames: 1) Pairs of parental and stutter 2) All sequences included in the analysis
#' 3) Sequences excluded from the analysis
#'
#' @export
#'

initialStutterTable <- function(inputDF_stutterTable) {

  ## 1. Samples and markers vectors:
  samples_test <- unique(inputDF_stutterTable$SampleID)
  markers <- unique(inputDF_stutterTable$Locus)

  ## 2. Initialisation of empty containers:
  stutterInfoTable <- c()
  excluded_SequencesDF <- c()
  included_SequencesDF <- c()

  ## 3. Iteration per marker and per sample:
  for (marker in markers) {
    for (sample in samples_test) {
      ## Subset the data into the sample and the marker:
      newdata <- subset(inputDF_stutterTable,inputDF_stutterTable$SampleID == sample & inputDF_stutterTable$Locus == marker)

      if (nrow(newdata)==1) next # Eg. marker PENTA E for sample A4PC, it has one sequence only, so there are no other
      # sequences to compare with for this marker

      ## a. Input vector for help function getMotifStutterInfo.R
      seqs2 <- rbind(newdata$Reads, newdata$sequence_type, newdata$SequenceID, newdata$AT_contentVar)
      colnames(seqs2) <- newdata$Forward_Strand_Bracketed_form

      ## b. Create stutter table info:
      stutterInfo <- getMotifStutterInfo(seqs2)

      ## c. Data frame with excluded sequences:
      excluded_Sequences <- subset(newdata, !(newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq1 | newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq2))

      included_Sequences <- subset(newdata, newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq1 | newdata$Forward_Strand_Bracketed_form %in% stutterInfo$Seq2)

      ## Control flow. Avoids stopping the loop when it encounters an empty data frame:
      if(nrow(stutterInfo)==0) next

      ## d. Merge all the matrices (as many as markers should be)
      stutterInfoTable <- rbind(stutterInfoTable, cbind(sample, marker, stutterInfo))
      #View(stutterInfoTable)

      # Merge data frames with included sequences for checking later:
      included_SequencesDF <- rbind(included_SequencesDF,included_Sequences)

      # Merge data frames with excluded sequences for checking later:
      excluded_SequencesDF <- rbind(excluded_SequencesDF, excluded_Sequences)
    }
  }

  ## 4. Conversion into data frame:
  stutterInfoTable <- data.frame(stutterInfoTable, stringsAsFactors = FALSE)
  colnames(stutterInfoTable) <- c("Samples", "Locus", "Seq1", "Seq2","MotifDifference", "Seq1 counts", "Seq2 counts",
                                  "sequence_TypeSeq1", "sequence_TypeSeq2", "SequenceIDSeq1","SequenceIDSeq2","ATcontentSeq1", "ATcontentSeq2", "Stutter proportion")

  ## 5. Create stutter_type variable for re-coding the stutter type (MotifDifference - stutter_type):

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

