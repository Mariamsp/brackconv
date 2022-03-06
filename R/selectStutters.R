
#' @title selectStutters
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Help function for selecting unambiguous stutters
#'
#' @details Obtains stutter sequences that can only be originated from one parental allele.
#' ### A second filtering step is applied and only most common stutters plus n0
#' ### (n-1, n-2, n+1, n+2 and n0) are further analysed. Infrequent stutters are also stored in a separate data frame.
#'
#' @param stutterInit_df Input data frame with pairs of parental and stutter sequences
#'
#' @return Four data frames: 1) Data frame with ambiguous stutters 2) Data frame with unambiguous stutters
#' 3) Data frame with infrequent stutters 4) Data frame with frequent stutters
#'
#' @export
#'

selectStutters <- function(stutterInit_df) {
  # 1. Change names PENTA E and PENTA D for avoiding errors:
  stutterInit_df$Locus <- gsub("PENTA D", "PentaD", stutterInit_df$Locus)
  stutterInit_df$Locus <- gsub("PENTA E", "PentaE", stutterInit_df$Locus)

  # 2. Remove duplicated stutters:
  # 2.1. Function for removing the duplicated stutters including the parental. Apply the function per sample and per locus:
  remove_duplicates <- function(stutterInit_df) {
    ## It excludes duplicated rows searching from the first value and after from the last
    new_df <- stutterInit_df[!(duplicated(stutterInit_df$Seq2) | duplicated(stutterInit_df$Seq2, fromLast = TRUE)), ]
    return(new_df)
  }

  # 2.2. Create data frame with unambiguous stutters (unambiguous definition = stutters which possible origin is one single parental sequence):
  markers <- unique(stutterInit_df$Locus)
  samples <- unique(stutterInit_df$Samples)

  unambiguous_stutter <- c()
  for (marker in markers) {
    for (sample in samples) {
      sliced_df <- subset(stutterInit_df, stutterInit_df$Samples == sample & stutterInit_df$Locus == marker )
      new_stutter <- remove_duplicates(sliced_df)
      if(nrow(new_stutter)==0) next
      unambiguous_stutter <- rbind(unambiguous_stutter,new_stutter)
    }
  }

 # 3. Create data frame for the ambiguous stutters:
  ambiguous_stutter <- stutterInit_df[!stutterInit_df$ID %in% unambiguous_stutter$ID,]

 # 4. Filtering of the unambiguous stutters (unambiguous_stutter data frame (df)):

  # Split the initial data frame containing all the unambiguous stutters into two data
  # frames, one containing the infrequent (rare stutters) and the other data frame will
  # contain the most observed/frequent stutters (n+1, n-1, n+2, n-2 and n0):
  # 4.1. Create "infrequent stutters" df:

  # 4.1.1. Selects one motif changes and two motif changes with a repetition number larger than 2 (selects > n + or - 2 stutters):
  infrequent_stutters <- unambiguous_stutter[grepl("[3-9]|\\d{2,}|^.*:[[:punct:]]2/.*:[[:punct:]]1|^.*:[[:punct:]]1/.*:[[:punct:]]2", unambiguous_stutter[["MotifDifference"]]),]
  # 4.1.2. Selects two motif changes with a repetition number equal to one but same sign for both changes (+/+ or -/-) (selects non n0!):
  infrequent_stutters2 <- unambiguous_stutter[grepl("^.*:\\+1/.*:\\+1|^.*:\\-1/.*:\\-1", unambiguous_stutter[["MotifDifference"]]),]
  # 4.1.3. Merge by row both data frames:
  infrequent_stutters <- rbind(infrequent_stutters,infrequent_stutters2)

  # 4.2. Create "frequent stutters" df by filtering out infrequent stutters from unambiguous_stutter df:
  frequent_stutters <- subset(unambiguous_stutter, !(unambiguous_stutter$ID %in% infrequent_stutters$ID))

  # 5. Convert categorical variable stutter_type into factors:
  levelsType <- unique(frequent_stutters$stutter_type)
  frequent_stutters$stutter_type <- factor(frequent_stutters$stutter_type)

  return(list(ambiguous_stutter,unambiguous_stutter, infrequent_stutters, frequent_stutters))
}
