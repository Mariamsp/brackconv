
#' @title getInitialDataframe
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Helpfunction for building working data frame
#'
#' @details set minimum reads for removing noise sequences and
#' set Hb for definition of true alleles.
#'
#' @param MinimumReads A numerical value for minimum read counts
#' @param Hb_Threshold Minimum vaue of heterozygote balance
#' @param df Imported data frame from lusSTR
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom stringr str_count
#'
#' @return Two data frames: 1) final data frame with stutter sequences and 2) data frame with noise sequences
#' @export
#'

getInitialDataframe <- function(df, MinimumReads,Hb_Threshold) {
  # a. Calculate per locus (depth of coverage?) the number of mapped reads and append a new column to original df:
  df2 <- transform(df, perlocus_reads = ave(df$Reads, list(df$SampleID,df$Locus), FUN = sum))

  # b. Calculate parameters for analysis (Per Sample and per Marker):

  # First max value and second max value. Percent reads, Heterozygote balance calculation and T:
  df2 <- df2 %>%
    mutate(SeqForw_length = nchar(Forward_Strand_Sequence)) %>%
    group_by(SampleID, Locus) %>%
    arrange(desc(Reads), .by_group= TRUE) %>%
    mutate(
      Index = row_number(), # Creates a variable with an index per group
      Max_reads = Reads[[1]], # First max read counts
      Sec_max_reads = ifelse(length(Reads) == 1, NA, Reads[[2]]), # Second max read counts
      Percent_reads = Reads*100/Max_reads, # Percentage of sequence reads compared to Max_reads variable (100%)
      ProHb = (Sec_max_reads/Max_reads)*100 # Percent proportion of Sec_max_reads and Max_reads (Sec_max_reads / Max_reads = X/100)
    )

  # c. Add unique identification label for each instance
  df2$SequenceID <- as.numeric(rownames(df2)) # We need the ID variable to be numeric
  # Move the column (last position) to the start:
  df2 <- df2 %>% dplyr::select(SequenceID, everything())

  # d. Add new column with sequence/allele type. Classification of Allele 1, Allele2, stutter and noise:

  # We define true alleles with Hb > 0.3 as obtained in previous analysis.
  df2 <- df2 %>%
    group_by(SampleID, Locus) %>%
    mutate(
      sequence_type = case_when(
        Index == 1 ~ "Allele1", # Accounts for regular HTZ cases and also when Max_reads and Sec_max_reads
        Index == 2 & ProHb > Hb_Threshold ~ "Allele2", # Accounts for regular HTZ cases and also when Max_reads and Sec_max_reads
        TRUE ~ "Stutter"
      ))

  # e. Write a function for creating the variable with the AT/GC content:
  calculate_ATContent <- function(df2) {
    count_A <- str_count(df2$Forward_Strand_Sequence, "A")
    count_T <- str_count(df2$Forward_Strand_Sequence, "T")
    AT_content <-round(((count_A+count_T) / (df2$SeqForw_length))*100, 2)
    return(AT_content)
  }
  df2$AT_contentVar <- calculate_ATContent(df2)

  # f. Filter for the minimum number of reads that are accepted for analysis This filtering will primarily
  # remove noise non-stutter sequences eg. base substitution sequencing errors.

  df3 <- df2[df2$Reads >= MinimumReads,] # Including the reads imposed by MinimumReads

  # Build a data frame with the noise that has been removed (Complement)
  noise_underMinimumReads <- df2[df2$Reads < MinimumReads,]

  # Convert sequence_type variable into
  df3$sequence_type <- factor(df3$sequence_type,
                                      levels <- c("Allele1", "Allele2", "Stutter"),
                                      labels <- c(1, 2, 3))

  return(list(df3, noise_underMinimumReads))
}



