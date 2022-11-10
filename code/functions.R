# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# FUNCTIONS
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/OneDrive/Documents/Homework/OSU/R Homework")

# SIMPLE-RATIO INDEX ------------------------------------------------------
SRI<- function(df){
  z<-ncol(df) # Name number of columns
  var.names = names(df)  # Make a placeholder for variable names
  sum.table = c() # Make a place for summary table
    x<- 
    yA<-
    yB<-
    yAB<- 
  SRIAB<- x/x+yA+yB+yAB
  for(i in 1:z) {
    sum.table<-
        rbind(sum.table,c(mean(x[,i]),sd(x[,i]),min(x[,i]),max(x[,i])))
  }
}

