# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# MCMC GLMM TESTS
###########################################################################

# Set working directory here
setwd("../data")

## load all necessary packages
library(ade4) # Look at Dai Shizuka/Jordi Bascompte
library(asnipe) # mrqap.dsp
library(assortnet) # associative indices
library(kinship2) # genetic relatedness
library(ggplot2) # Visualization
library(doParallel) # For faster coding
library(abind) # array
library(statip) # dbern
library(MCMCglmm) # MCMC models
library(brms) # Baysian
library(nimble) # For MCMC
library(mcmcplots) # For MCMC plots
library(MCMCvis)
library(tidyverse)
source("../code/attach.nimble_v2.R")
source("../code/functions.R") # nxn

# Read in social association matrix and listed data
dist_BG <- readRDS("dist_BG.RData") # BG Sim Matrix
dist_FG <- readRDS("dist_FG.RData") # FG Sim Matrix
dist_SD <- readRDS("dist_SD.RData") # SD Sim Matrix
kov <- readRDS("kov.RDS")  # Home range overlap
nxn <- readRDS("nxn_ovrlap.RData") # Association Matrix
list_years <- readRDS("list_years_ovrlap.RData") # Data listed into periods
SE_array <- readRDS("SE_array.RData")

# Now make a sex and age data frame
ILV_df <- list_years[[2]][!duplicated(list_years[[2]][, "Code"]), c("Code", "Sex", "Age")]
ILV_df$Sex <- ifelse(ILV_df$Sex == "Female", 0, 
       ifelse(ILV_df$Sex == "Male", 1, NA))

# Estimate unknowns
ILV_df$Sex <- ifelse(is.na(ILV_df$Sex), rbinom(n = nrow(ILV_df), size = 1, prob = 0.5), ILV_df$Sex)
ILV_df$Age <- ifelse(is.na(ILV_df$Age), floor(runif(n = nrow(ILV_df), min = 1, max = 56)), ILV_df$Age) # uniform probability for all ages

# Make sex sim and age diff matrix
sex_sim <- matrix(0, nrow = nrow(nxn[[1]]), ncol = ncol(nxn[[1]]))
age_diff <- matrix(0, nrow = nrow(nxn[[1]]), ncol = ncol(nxn[[1]]))

for (i in 1:(nrow(sex_sim) - 1)) {
  for (j in (i + 1):ncol(sex_sim)) {
    sex_sim[i, j] <- as.numeric(ILV_df$Sex[i] == ILV_df$Sex[j])
    age_diff[i, j] <- abs(ILV_df$Age[i] - ILV_df$Age[j])
    
    # Since sex_sim and age_diff are symmetric matrices, update the corresponding values
    sex_sim[j, i] <- sex_sim[i, j]
    age_diff[j, i] <- age_diff[i, j]
  }
}

# Smaller sample data
dist_BG <- lapply(dist_BG, function(ls) ls[c(1:20),c(1:20)])
dist_FG <- lapply(dist_FG, function(ls) ls[c(1:20),c(1:20)])
dist_SD <- lapply(dist_SD, function(ls) ls[c(1:20),c(1:20)])
kov <- kov[c(1:20),c(1:20)]
nxn <- lapply(nxn, function(ls) ls[c(1:20),c(1:20)])
ILV_df <- ILV_df[c(1:20),]
nxn <- abind(nxn, along=3)
SE_array <- SE_array[1:20, 1:20, 1:2]

################################################################################
# PART 1: Create HI Disimilarity Matrix  ------------------------------------------------

# Genetic relatedness matrix
GR <- lapply(list_years, function(df) {
  
  pedigree_matrix <- list_years_sexage[, c("Code", "Mom", "Dad")]
  kin_matrix <- kinship(pedigree_matrix)
  
  return(kin_matrix)
  
})

# Extract specific columns from each data frame in list_years
aux_data <- function(list_years) {
aux <- lapply(list_years, function(df) {
  data.frame(Code = df$Code,
    Behaviors = df$Behaviors,
    HumanInteraction = df$HumanInteraction,
    ConfHI = df$ConfHI)})

# Add the 'Foraging' variable to each data frame in the 'aux' list
aux <- lapply(aux, function(df) {
  df$Foraging <- "Other"
  df$Foraging[grepl(pattern = 'Feed', x = df$Behaviors, ignore.case = FALSE)] <- "Feed"
  df
})
return(aux)
}

aux <- aux_data(list_years)

# Categorize ID to Foraging
ID_forg <- function(aux_data) {
IDbehav <- lapply(aux_data, function(df) {
  df <- table(df$Code, df$Foraging)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- df[, c(1, 3)]
  colnames(df) <- c("Code", "Forg_Freq")
  df <- aggregate(. ~ Code, data = df, sum)
  df
})
return(IDbehav)
}

IDbehav <- ID_forg(aux)

# HI behaviors should be partitioned into 3 different types---------------------
#' BG = Beg: F, G
#' SD = Scavenge and Depredation: B, C, D, E
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("B", "C", "D", "E"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", "None")))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(aux)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
rawHI_diff <- lapply(aux_data, function(df) {
  table_df <- as.data.frame(table(df$Code, df$DiffHI))
  colnames(table_df) <- c("Code", "DiffHI", "Freq")
  return(table_df)
})}

rawHI_diff <- diff_raw(aux)

# Create a frequency count for each HI behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    HI_freq <- rawHI_diff_data[[i]]$Freq[rawHI_diff_data[[i]]$DiffHI == HI]
    df$HI <- HI_freq[match(df$Code, rawHI_diff_data[[i]]$Code)]
    colnames(df) <- c("Code", "Foraging", "HI")
    df
  })
}

# Including zeros
IDbehav_BG <- get_IDHI("BG", IDbehav, rawHI_diff)
IDbehav_FG <- get_IDHI("FG", IDbehav, rawHI_diff)
IDbehav_SD <- get_IDHI("SD", IDbehav, rawHI_diff)

# Clump all the HI behaviors together
clump_behav <- function(aux_data) {
for (i in seq_along(aux_data)) {
  aux_data[[i]]$ConfHI <- ifelse(aux_data[[i]]$ConfHI != "0", 1, 0)}

# Categorize ConfHI to IDs
rawHI <- lapply(aux_data, function(df) {
  # Sum up the frequencies of HI by code
  aggregated_df <- aggregate(ConfHI ~ Code, data = df, sum)
  unique_codes_df <- data.frame(Code = unique(df$Code))
  # Merge the unique codes data frame with the aggregated data frame
  merged_df <- merge(unique_codes_df, aggregated_df, by = "Code", all.x = TRUE)
  # Fill missing Freq values (if any) with 0
  merged_df$ConfHI[is.na(merged_df$ConfHI)] <- 0
  return(merged_df)
})
return(rawHI)
}

rawHI <- clump_behav(aux)

# Get HI Freq
create_IDbehav_HI <- function(IDbehav_data, rawHI_data){
IDbehav_HI <- lapply(seq_along(IDbehav_data), function(i) {
      df <- IDbehav_data[[i]]
      df$HI <- rawHI_data[[i]]$ConfHI
      colnames(df) <- c("Code", "Foraging", "HI")
      df
    })
return(IDbehav_HI)
}

IDbehav_HI <- create_IDbehav_HI(IDbehav, rawHI)

# Proportion of time Foraging spent in HI
Prop_HI <- function(IDbehav) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HIprop <- as.numeric(df$HI) / as.numeric(df$Foraging)
    df$HIprop[is.na(df$HIprop)] <- 0
    # Keep only 'Code' and 'HIprop' columns
    df <- df[, c('Code', 'HIprop')]
    df
  })
}

prob_HI <- Prop_HI(IDbehav_HI)
prob_BG <- Prop_HI(IDbehav_BG)
prob_SD <- Prop_HI(IDbehav_SD)
prob_FG <- Prop_HI(IDbehav_FG)

# Dissimilarity of HI proportion among individual dolphins, using Euclidean distance
dis_matr <- function(Prop_HI) {
  dissimilarity_HI <- list()
  for (i in seq_along(Prop_HI)) {
    fake_HIprop <- Prop_HI[[i]]$HIprop
    dissimilarity_HI[[i]] <- as.matrix(dist(matrix(fake_HIprop), method = "euclidean"))
  }
  return(dissimilarity_HI)
}

dist_BG <- dis_matr(prob_BG)
dist_SD <- dis_matr(prob_SD)
dist_FG <- dis_matr(prob_FG)

saveRDS(dist_BG, "dist_BG.RData")
saveRDS(dist_SD, "dist_SD.RData")
saveRDS(dist_FG, "dist_FG.RData")

################################################################################
# PART 2: Diagnostics ----------------------------------------------------------

# Check for collinearity 
# Check if it is based off zeros
## Create a list of predictor matrices
year <- 1
predictor_matrices <- list(dist_BG[[year]], dist_FG[[year]], dist_SD[[year]])

## Calculate correlation matrix
num_predictors <- length(predictor_matrices)
correlation_matrix <- matrix(NA, nrow = num_predictors, ncol = num_predictors)

for (i in 1:num_predictors) {
  for (j in (1+i):num_predictors) {
    mtest <- mantel.rtest(as.dist(predictor_matrices[[i]]), as.dist(predictor_matrices[[j]]), nrepet=999)
    correlation_matrix[i, j] <- mtest$obs
  }
}

## Print the correlation matrix
print(correlation_matrix) # It seems that BEG and HI are highly correlated


###########################################################################
# PART 3: Create MCMC GLMMs  ------------------------------------------------

# Write a Nimble model: SRI ~ HRO + SEX + AGE + GR + HAB(BP + FG + SD)
model1 <- nimbleCode({
  
  #Priors
  
    ## ILV Effects
    HRO_Effect ~ dt(mu=0, sigma=1, df=1)
    SEX_Effect ~ dt(mu=0, sigma=1, df=1)
    AGE_Effect ~ dt(mu=0, sigma=1, df=1)
    #GR_Effect ~ dt(mu=0, sigma=1, df=1)
    Rand.Err ~ T(dt(mu=0, sigma=1, df=1), 0, )
    
    ## HI Effects
    for (p in 1:n.per) {
      BP_Effect[p] ~ dt(mu=0, sigma=1, df=1)
      FG_Effect[p] ~ dt(mu=0, sigma=1, df=1)
      SD_Effect[p] ~ dt(mu=0, sigma=1, df=1)
    } #p
    
    # Run through matrix
    for(i in 1:n.ind){
      
    # Estimate unknowns
    SEX[i] ~ dbern(prob = 0.5) # Same probability for female or male
    AGE[i] ~ dunif(0, 56) # uniform probability for all ages
    AGE.Est[i] <- floor(AGE[i]) # Get a whole number
    
    # Random effect term
    u[i] ~ dnorm(mean = 0, sd = Rand.Err)
    
    for (j in (1+i):n.ind) {
      
      # Intercept Prior for each ID
      IJ_Random[i, j] <- (u[i]+ u[j])
      
      #Impute missing sexes & ages:
      SEX.SIM[i, j] <- (SEX[i] == SEX[j])
      AGE_Diff[i, j] <- abs(AGE.Est[i] - AGE.Est[j])
      #AGE_Diff_Scaled[i, j] <- (AGE_Diff[i, j] - mean(AGE_DIFF[1:n.ind,1:n.ind],na.rm=T)/ sd(AGE_DIFF[1:n.ind,1:n.ind],na.rm=T))
      
      ## HI Effects
      for (p in 1:n.per) {
        for (y in 1:n.yr) {
          # SRI Observation Error
          NXN[i, j, y] <-  a[i, j, y] / (a[i, j, y] + b[i, j, y] + c[i, j, y])
          squared_diff[i, j, y, p] <- (NXN[i, j, y] - SRI[i, j, p])^2
        } #y
      
      Obs.Err[i, j, p] <- sqrt(sum(squared_diff[i, j, p]) / (n.yr - 1))/ sqrt(n.yr)
      
      # Process Model
      logit(SRI.Exp[i, j, p]) <- IJ_Random[i, j] + 
        HRO[i, j]*HRO_Effect + SEX.SIM[i, j]*SEX_Effect + AGE_Diff[i, j]*AGE_Effect + #GR[i, j]*GR_Effect +
        BP[i, j, p]*BP_Effect[p] +  FG[i, j, p]*FG_Effect[p] +  SD[i, j, p]*SD_Effect[p]
      
      # Observation Model (Likelihood) # a = mean^2*(1-mean)/sd^2-mean, b = mean*(1-mean)^2/sd^2+mean-1
      SRI[i, j, p] ~ dbeta(shape1 = (SRI.Exp[i, j, p]*(SRI.Exp[i, j, p]*(1-SRI.Exp[i, j, p])/(Obs.Err[i, j, p]^2)-1)), 
                           shape2 = ((1-SRI.Exp[i, j, p])*(SRI.Exp[i, j, p]*(1-SRI.Exp[i, j, p])/(Obs.Err[i, j, p]^2)-1)))
      }#p
    }#j
    
  }#i
  
  
})#model1




# Parameters monitored (are there any new parameters to include?)
parameters <- c("IJ_Random", 
                "HRO_Effect", "SEX_Effect", "AGE_Effect", #"GR_Effect", 
                "BP_Effect", "FG_Effect", "SD_Effect")

# MCMC Settings
ni <- 40000
nt <- 40
nb <- 20000
nc <- 3

# Array for the matrix

# Data
nimble.data = list(SRI = nxn,
                   HRO = kov,
                   SEX = ILV_df$Sex,
                   AGE = ILV_df$Age,
                   #GR = gr_list,
                   BP = abind(dist_BG, along=3),
                   FG = abind(dist_FG, along=3),
                   #SD = abind(dist_SD, along=3),
                   Obs.Err = SE_array
                   )

nimble.constants = list(n.ind = length(unique(ILV_df$Code)),
                        n.per = length(nxn))

saveRDS(nimble.constants, "../data/nimble.constants.RData")
saveRDS(nimble.data, "../data/nimble.data.RData")

n.cores <- detectCores() 
registerDoParallel(n.cores)
mcmc.output <- nimbleMCMC(code = model1,
                          data = nimble.data,
                          constants=nimble.constants,
                          monitors = parameters,
                          niter = ni,
                          nburnin = nb,
                          nchains = nc,
                          thin=nt,
                          summary=TRUE,
                          samplesAsCodaMCMC = TRUE)
### End parallel processing
stopImplicitCluster()



attach.nimble(mcmc.output$samples)

MCMCtrace(object = mcmc.output$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "AGE_Effect")

# Gelman-Rubin diagnostic (AKA RHat or PSRF)
gelman.diag(mcmc.output$samples)

# Visualize all of the relevant plots at the same time:
mcmcplot(mcmc.output$samples)

# Summarize the posterior distributions:
hist(AGE_Effect)


################################################################################

## Package for Bayes
# Prepare random effect for MCMC
num_nodes <- lapply(nxn, function(df) dim(df)[1])
node_names <- lapply(nxn, function(df) colnames(df))

# Separate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) matrix(rep(1:df, each = df), nrow = df, ncol = df))
node_ids_j <- lapply(node_ids_i, function(df) t(df))

df_list <- vector("list", length = length(nxn))
for (i in seq_along(nxn)) {
  upper_tri <- upper.tri(nxn[[i]], diag = TRUE)
  
  df_dolp <- data.frame(edge_weight = nxn[[i]][upper.tri(nxn[[i]])], 
                        age_difference = age_list[[i]][upper.tri(age_list[[i]])],
                        sex_difference = sex_list[[i]][upper.tri(sex_list[[i]])],
                        HI_differences = dist_HI_sexage[[i]][upper.tri(dist_HI_sexage[[i]])],
                        HRO = kov_sexage[[i]][upper.tri(kov_sexage[[i]])],
                        node_id_1 = factor(as.vector(node_names[[i]][node_ids_i[[i]][upper_tri]]), levels = node_names[[i]]),
                        node_id_2 = factor(as.vector(node_names[[i]][node_ids_j[[i]][upper_tri]]), levels = node_names[[i]]))
  
  df_list[[i]] <- df_dolp
}

# Calculate the mean value for each individual
combined_array <- lapply(seq_along(df_list), function(i) {
  upper_tri <- upper.tri(SE_array[,,i], diag = TRUE)
  df_list[[i]]$SE <- SE_array[,,i][upper_tri]
  return(df_list[[i]])
})

# Multimembership models in MCMCglmm
fit_mcmc <- MCMCglmm(edge_weight ~ HI_differences + HRO + age_difference + sex_difference, 
                     random=~mm(node_id_1 + node_id_2), data=df_list[[year]])
summary(fit_mcmc)


###########################################################################
# PART 5: Assortivity Index Based on HI Over Time  ------------------------------------------------

# Match Code with matrix and vector
get_HI_vector <- function(prop_HI) {
  HI_vector <- lapply(seq_along(nxn), function(i) {
    matrix_index <- match(rownames(nxn[[i]]), prop_HI[[i]]$Code)
    reordered_prob_HI <- prop_HI[[i]][matrix_index, ]
    return(reordered_prob_HI)
  })
  return(HI_vector)
}

# Get each combined and seperate HI
HI_vector <- get_HI_vector(prob_HI)
BG_vector <- get_HI_vector(prob_BG)
SD_vector <- get_HI_vector(prob_SD)
FG_vector <- get_HI_vector(prob_FG)

# Look at HI assortivity coefficient over periods
calculate_assortment <- function(HI_vector) {
  n.cores <- detectCores()
  registerDoParallel(n.cores)
  
  assort_HI <- NULL
  # se <- NULL
  for (i in seq_along(nxn)) {
    coeff <- assortment.continuous(nxn[[i]], HI_vector[[i]][, "HIprop"], SE = FALSE)
    assort_HI[i] <- coeff$r
    # se[i] <- coeff$se
  }
  
  # End parallel processing
  stopImplicitCluster()
  
  assort_HI_df <- data.frame(HI_assort = unlist(assort_HI), Year = c(1:2))
  return(assort_HI_df)
}

# Look at HI combined and separate
assort_HI <- calculate_assortment(HI_vector)
assort_BG <- calculate_assortment(BG_vector)
assort_SD <- calculate_assortment(SD_vector)
assort_FG <- calculate_assortment(FG_vector)

# Combine the assort dataframes and add a behavior column
assort_BG$Behavior <- "BG"
assort_SD$Behavior <- "SD"
assort_FG$Behavior <- "FG"

combined_assort <- rbind(assort_BG, assort_SD, assort_FG)
combined_assort$HI_assort <- ifelse(is.na(combined_assort$HI_assort), 0.75, 
                                    combined_assort$HI_assort)

# Calculate sample size
BG_1 <- length(unique(BG_vector[[1]]$Code[BG_vector[[1]]$HIprop > 0]))
BG_2 <- length(unique(BG_vector[[2]]$Code[BG_vector[[2]]$HIprop > 0]))
FG_1 <- length(unique(FG_vector[[1]]$Code[FG_vector[[1]]$HIprop > 0]))
FG_2 <- length(unique(FG_vector[[2]]$Code[FG_vector[[2]]$HIprop > 0]))
SD_1 <- length(unique(SD_vector[[1]]$Code[SD_vector[[1]]$HIprop > 0]))
SD_2 <- length(unique(SD_vector[[2]]$Code[SD_vector[[2]]$HIprop > 0]))

sample_sizes <- c(BG_1, BG_2, FG_1, FG_2, SD_1, SD_2)

# Create the combined plot with facets
ggplot(combined_assort, aes(x = Year, y = HI_assort, label = paste("n =", sample_sizes))) +
  geom_point() +
  geom_line() +
  labs(x = "Period", y = "HI assortment") +
  ggtitle("Whisker Plot of HI assortment") +
  theme_minimal() +
  facet_grid(Behavior ~ ., scales = "free_y", space = "free_y") +
  geom_text(size = 2, vjust = 0.5)

#' The assortivity of each behavior decreased. Begging and foraging around
#' fixed gear only slightly decreased while scavenging and depredation assortivity
#' greatly decreased.