
#' @title AddVariables
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Help function for adding variable "inclusion_stutterTab" to input data frame
#'
#' @details Variable "inclusion_stutterTab" indicates if a sequence has been included in the stutter analysis or not
#'
#' @param inputDF_stutterTable Input data frame with pairs of parental and stutter sequences
#' @param in_df Data frame with sequences included in the analysis
#' @param ex_df1 Data frame with sequences excluded in the analysis
#'
#' @return Four data frames: 1) Data frame with all excluded sequences 2) Data frame with all included sequences
#' 3) Data frame sequences excluded before the stutter analysis (first analytical step) 4) Input data frame with pairs of parental
#' and stutter sequences (inlcuidng new variable)
#'
#' @export
#'

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
