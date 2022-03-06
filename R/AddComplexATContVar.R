
#' @title AddComplexATContVar
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Help function for adding complexity and AT-content variables to data frame
#'
#' @details # 1. Motif type as a factor variable and establish groups (do it per marker)
#' # 2. STR complexity as factor variable:
#' * Simple
#' * Compound
#' * Complex
#'
#' @param df Data frame output from \code{\link{getLUS}}
#'
#' @importFrom dplyr left_join
#'
#' @return Two data frames: 1) final data frame with stutter sequences and 2) data frame with noise sequences
#' @export
#'

helpfun_ComplexATContVar <- function(df, STR_complexity_df) {
  # 1. Add Motif type variable:
  df$Motif_Type <- factor(df$MotifDifference)

  # 2. STR complexity variable:
  # Create object with STR compelxity classes:
  STR_complexity_df <- system.file("extdata", "STR_Complexity.xlsx", package = "brackconv")

  # Merge df with table containing STR complexity classes:
  df <- left_join(df, STR_complexity_df)

  return(df)
  }
