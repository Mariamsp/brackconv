
#' @title readlusSTRFile
#'
#' @author Maria Martin Agudo <maagud.at.ous-hf.no>
#'
#' @description Read output file from luSTR and prepare it for later analysis
#'
#' @details Import output file obtained from lusSTR
#'
#' @param path Add path to lusSTRFile
#' @param lusSTRFile name with extension of txt input file obtained from lusSTR
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
#' @return Data frame object prepared analysis in \code{\link{getInitialDataframe}}
#' @export
#'
#' @example
#' \dontrun{
#' path = "C:/Users/username/lusSTR/Project_Test"
#' lusSTRFile = "STRaitRazor_forenseq_final.txt"
#' initial_df <- readlusSTRFile(path, lusSTRFile)}
#'
readlusSTRFile <- function(path,lusSTRFile) {
  initial_df <- read.table(paste0(path,lusSTRFile), sep="\t", header = TRUE)

  # Create variable SampleNumber with only digits from SampleID:
  initial_df <- initial_df %>%
    mutate(SampleNumber = as.numeric(gsub("[^[:digit:]]", "", SampleID))) %>%
    select(-SampleNumber)
}
