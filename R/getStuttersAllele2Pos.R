
#' @title getStuttersAllele2Pos
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Help function for obtaining Allele 2 sequences in n-1 position
#'
#' @details Identifies Allele 2 sequences in n-1 position when Hb is < 0.3 and removes them from
#' stutter table
#'
#' @param df1 Input data frame with pairs of parental and stutter sequences
#' @param df_ext Input data frame with all parental and stutter sequences
#' @param Hb_min Minimum Hb value
#' @param Hb_max Maximum Hb value
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#'
#' @return Three data frames: 1) Data frames with stutter sequences that need to be removed 2) New stutter table with removed sequences
#' 3) Input data frame with all parental and stutter sequences and new column with Allele2_type variable
#'
#' @export
#'

getStuttersAllele2Pos <- function(df1,df_ext,Hb_min, Hb_max) {

  ## 1. Add rownames as ID column for additional filtering:
  df1 <- cbind(ID = rownames(df1), df1)

  ## 2. Remove parental allele2 that falls in stutter position n-1 respective to allele 1:

  ### Step1. Add column with Hb proportion variable with the percentage proportion of reads of Seq2 compared to Seq1:
  df1$Hb_proportion <- df1$`Seq2 counts`*100 / df1$`Seq1 counts`

  ### Step2. Construct a data frame that contains the instances we want to filter out later from our main data frame:

  ### We want to filter out those sequences that fulfills the following criteria for removing sequences:

  # A) ### * Sequence Seq1 is parental Allele 1 and sequence Seq2 is parental allele 2.
  # * The difference found between both sequences is of the class: stutter type n-1. In which case the parental allele 2 is in a stutter position,
  # and this is a confounder for determining number of reads from the true parental and number of reads from a potential stutter in that position.
  # * We also get a more realistic data set when framing those parental alleles 2 with a low Heterozygous balance (high intralocus imbalance),
  # and we determine that the optimal Heterozygous balance lower an upper threshold  for filtering out sequences is 0.3 and 0.5
  # (initial value based on a positive control assay).
  # * Parental alleles will be only defined when Hb_proportion > 0.3 and stutters will be found under this
  # value.

  # B) ### Code for filtering: # Check point ### !!!! Check if there is any Hb < 0.3 for these sequences /// The main stutters we want
  # to look at is the n-1 stutter (the most frequent type). And we see that the stutter proportions are high, approx 0.3.
  df2 <- df1 %>%
    group_by(Samples, Locus) %>%
    filter(sequence_TypeSeq1 == 1 & sequence_TypeSeq2 == 2 & Hb_proportion > Hb_min & Hb_proportion  < Hb_max & stutter_type == "bwd-1")

  # C) Generate a new variable in initial_df2 and classify Allele 2 sequences in n-1 position:
  df_ext$Allele2_type[df_ext$SequenceID %in% df2$SequenceIDSeq2 & df_ext$sequence_type == 2 & df_ext$DropOut_type == "NoDropOutAllele2"] = "Allele2_bwd1Position"
  df_ext$Allele2_type[!df_ext$SequenceID %in% df2$SequenceIDSeq2 & df_ext$sequence_type == 2 & df_ext$DropOut_type == "NoDropOutAllele2"] = "Allele2"
  df_ext$Allele2_type[is.na(df_ext$Allele2_type)] <- "Not_Allele2" # Replace NAs with the label "Not_Allele2" for other sequences that are not Allele2

  # D) ### Create a temp data frame with rows where Seq1 (origin of stutter) is parental Allele 2 in stutter position n-1 (Observed in allele1StutterPos)
  ### Important: We only consider parental Allele 2 with Hb between 0.3 and 0.5 parental Allele1 and falling in stutter position, because for these instances
  ### sequences can be a mixed of true parental allele 2 and stutter of parental allele 1.
  stutterInfoTable_SubsetNOT <- df1[unlist(sapply(1:NROW(df2), function(i)
    which(df1$Samples == df2$Samples[i] & df1$Seq1 == df2$Seq2[i]))),]

  # E)  Final Data frame. Subset and NOT to include rows from stutterInfoTable_SubsetNOT:
  df_final <- df1[!df1$ID %in% stutterInfoTable_SubsetNOT$ID,]

  ## 3. Keep only those rows where Seq1 == Allele1 or Allele2 (1 or 2) and Seq2 == Stutter (3):
  df_final <- df_final[df_final[,9] == 1 & df_final[,10] == 3 | df_final[,9] == 2 & df_final[,10] == 3  ,, drop = FALSE]
  # View(stutterInfoTable_Final)

  return(list(df2,df_final,df_ext))

}

