
################################################################################
# HELPFUNCTION

##### Compares short bracket format to original
##### full length sequence

#' @title getMotifReps
#'
#' @author Oyvind Bleka <oyvble.at.hotmail.com>
#'
#' @description Get motifs and repetitions from short bracket format
#'
#' @details Help function. Obtains motif and repeats from short bracket format
#'
#' @param tmp A string with bracket format e.g. "[ATCG]6"
#' @param brksign Bracket signs to use
#'
#' @return Vector with 2 elements: Extracted motif and corresponding repetitions
#'
#' @export
#'
#' @examples
#' \dontrun{
#' getMotifReps("[ATCG]6")
#' }

getMotifReps = function( tmp, brksign=c("[","]") ) {
  #input: tmp = ""
  if( grepl("\\[",tmp) ) {
    tmp2 = strsplit(tmp,"\\]")[[1]]
    motif = gsub("\\[","",tmp2[1]) #motif string (without brackets)
    motifReps = tmp2[2] #convert number of repititions to number
  } else {
    motif = tmp
    motifReps = 1 #only 1 motif repetition
  }
  return(c(motif,motifReps)) #return motif string and number of motif
}
