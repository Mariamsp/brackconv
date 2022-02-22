
### 8_Helpfun_ComplexATContVar.R ###
### Aim: Function for step 8: generate new variables: Motif type and STR complexity ###

# 1. Motif type as a factor variable and establish groups (do it per marker)
# 2. STR complexity as factor variable
# * Simple
# * Compound
# * Complex

helpfun_ComplexATContVar <- function(df) {
  # 1. Add Motif type variable:
  df$Motif_Type <- factor(df$MotifDifference)
  
  # 2. STR complexity variable:
  # Merge the table finalStutterTable with the file containing the table with the STR complexity classes :
  df <- left_join(df, STR_complexity_df)

  return(df)
  }