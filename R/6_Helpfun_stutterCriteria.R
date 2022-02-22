
### 6_Helpfun_stutterCriteria.R
### Aim: Function for step 6: additional filtering of stutters ###

### Obtain unambiguous stutters:
### We want to frame those stutters that can only be originated from one parental allele. ###
### This stutters will be unique per sample. ###
### We will also keep the stutters that potentially can be originated from both parental alleles in a
### separate data frame. ###
### A second filtering step will be applied and we will only keep the most common stutters plus n0 
### (bwd-1, bwd-2, fwd+1, fwd+2 and n0), infrequent stutters will be set aside in a different data 
### frame, and we can take them for other studies later ###

helpfun_StutterFilterAdd <- function(stutterInit_df) {
  #stutterInit_df = stutterInfoTable_Final
  # 1. Change names PENTA E and PENTA D for avoiding errors:
  stutterInit_df$Locus <- gsub("PENTA D", "PentaD", stutterInit_df$Locus)
  stutterInit_df$Locus <- gsub("PENTA E", "PentaE", stutterInit_df$Locus)
  
  # 2. Remove duplicated stutters:
  # 2.1. Function for removing the duplicated stutters including the parental. Apply the function per sample and per locus:
  remove_duplicates <- function(stutterInit_df) {
    new_df <- stutterInit_df[!(duplicated(stutterInit_df$Seq2) | duplicated(stutterInit_df$Seq2, fromLast = TRUE)), ] ## It excludes duplicated rows searching from the first value and after from the last
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
  # frames, one containing the infrequent (complicated stutters) and the other data frame will 
  # contain the most observed stutters (+1, -1 or +1/-1 and n0):
  
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