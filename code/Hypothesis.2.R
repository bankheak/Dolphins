# Network Analysis of Anthropogenic Influences in Bottlenose dolphins

# Local Network Analysis Hypothesis 2 #

# Set working directory here
setwd("../data")

# Load all necessary packages
## Mapping
if(!require(GGally)){install.packages('GGally'); library(GGally)} # For mapping networks in ggplot
if(!require(ggraph)){install.packages('ggraph'); library(ggraph)}
if(!require(ggpattern)){install.packages('ggpattern'); library(ggpattern)}
if(!require(ggpattern)){install.packages('ggpattern'); library(ggpattern)} # geom_tile_pattern
if(!require(grid)){install.packages('grid'); library(grid)}
if(!require(viridis)){install.packages('viridis'); library(viridis)} # plot themes
if(!require(ggpattern)){install.packages('ggpattern'); library(ggpattern)} # heatmap hatches
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)}
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)}
## Network
if(!require(network)){install.packages('network'); library(network)} # For assigning coordinates to nodes %v%
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)} # For strength gradient network colors
if(!require(tnet)){install.packages('tnet'); library(tnet)} # For weights
if(!require(igraph)){install.packages('igraph'); library(igraph)} # Measure centrality here
if(!require(assortnet)){install.packages('assortnet'); library(assortnet)} # associative indices
source("../code/functions.R") # Matrix_to_edge_list
## LMM Test
if(!require(lme4)){install.packages('lme4'); library(lme4)} 
if(!require(nlme)){install.packages('nlme'); library(nlme)} 
if(!require(lmerTest)){install.packages('lmerTest'); library(lmerTest)} 
## Bayesian
if(!require(abind)){install.packages('abind'); library(abind)} # array
if(!require(brms)){install.packages('brms'); library(brms)} # For brm modellibrary(coda)
if(!require(bayesplot)){install.packages('bayesplot'); library(bayesplot)} # plot parameters
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Run parallel processing
if(!require(rstan)){install.packages('rstan'); library(rstan)} # To make STAN run faster
if(!require(tidybayes)){install.packages('tidybayes'); library(tidybayes)} # get_variables

# Read in full datasheet and list (after wrangling steps)
list_years <- readRDS("list_years.RData") # (1995-2000)/(2001-2006)/(2007-2012)
nxn <- readRDS("nxn.RData") # association matrix of list_years

###########################################################################
# PART 1: Wrangle Data ---------------------------------------------

# Read in aux data
aux <- readRDS("aux.RData")

# Look at how many individuals have HI
HI_1 <- unique(aux[[1]]$Code[aux[[1]]$DiffHI != "None"])
HI_2 <- unique(aux[[2]]$Code[aux[[2]]$DiffHI != "None"])
HI_3 <- unique(aux[[3]]$Code[aux[[3]]$DiffHI != "None"])
length(unique(c(HI_1, HI_2, HI_3)))

# Read in IDbehav data
IDbehav_BG <- readRDS("IDbehav_BG.RData")
IDbehav_FG <- readRDS("IDbehav_FG.RData")
IDbehav_SD <- readRDS("IDbehav_SD.RData")

# Get total number of HI individuals
BG_IDs <- lapply(IDbehav_BG, function (df) df$Code[df$Behav != 0])
FG_IDs <- lapply(IDbehav_FG, function (df) df$Code[df$Behav != 0])
SD_IDs <- lapply(IDbehav_SD, function (df) df$Code[df$Behav != 0])
ovrlap_IDs <- intersect(intersect(unlist(BG_IDs), unlist(FG_IDs)), unlist(SD_IDs))


###########################################################################
# PART 2: Calculate Local Metrics ---------------------------------------------

## Create social network
ig_func <- function(nxn) {
  ig <- lapply(nxn, function (df) {
    graph_from_adjacency_matrix(
      df,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE)})
  return(ig)}

ig <- ig_func(nxn)

# Set the node names based on row names
row_name_assign <- function(nxn, ig) {
  row_names <- lapply(nxn, function (df) {rownames(df)})
  for (i in seq_along(ig)) {
    V(ig[[i]])$name <- row_names[[i]]
  }
}

row_name_assign(nxn, ig)

# Save ig object
saveRDS(ig, "ig.RData")
ig <- readRDS("ig.RData")

# Edgelist: Nodes (i & j) and edge (or link) weight
n.cores <- detectCores()
registerDoParallel(n.cores)
el_years <- lapply(nxn, function (list) matrix_to_edgelist(list, rawdata = FALSE, idnodes = FALSE))
### End parallel processing
stopImplicitCluster()

# Save the el_list
saveRDS(el_years, "el_years.RData")
el_years <- readRDS("el_years.RData")

# step by step: connectance = length(which(as.dist(orca_hwi)!=0))/(N*(N-1)/2)
# Number of nodes (number of rows in the association matrix)
N = nrow(nxn[[1]])
# Number of possible links: Nodes*(Nodes-1)/2: (-1 removes the node itself; /2 removes repetitions)
total = N*(N-1)/2
# Number of realized links: all non-zero cells in the association matrix
real = length(which(as.dist(nxn[[1]])!=0))
# Connectance: realized/total
real/total # Connectance is low

# Set the node names based on row names
get_names <- function (matrix, metric) {
  row_names <- lapply(matrix, function (df) {rownames(df)})
  for (i in seq_along(metric)) {
    metric[[i]][,1] <- row_names[[i]]
  }
  return(metric)
}

# Degree and strength centrality 
strength <- lapply(el_years, function (df) {degree_w(df, measure=c("degree","output"), type="out", alpha=0.5)})
strength_diffs <- get_names(nxn, strength)

compare_strength <- merge(
  merge(strength_diffs[[1]], strength_diffs[[2]], by = "node"),
  strength_diffs[[3]], by = "node"
)

colnames(compare_strength) <- c("ID", "Before_HAB_degree", "Before_HAB_strength", 
                                "During_HAB_degree", "During_HAB_strength", 
                                "After_HAB_degree", "After_HAB_strength")

compare_strength[, c(2:7)] <- sapply(compare_strength[, c(2:7)], as.numeric)

# Look at all of the local metrics together
HI_data <- readRDS("aux.RData")

## Add a column containing HI type
names_BG <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "BG"]))})
names_SD <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "SD"]))})
names_FG <- lapply(HI_data, function (df) {
  as.vector(unique(df$Code[df$DiffHI == "FG"]))})

### Find all of the IDs that match to make sure below works
matching_unique_ids <- list()
for (i in 1:3) {
  matching_unique_ids[[i]] <- unique(c(names_BG[[i]], names_FG[[i]], names_SD[[i]]))
}

names_NF <- list()
for (i in 1:3) {
  unique_codes <- unique(HI_data[[i]]$Code)
  
# Check for codes that are not in any of the names_BG, names_FG, names_SD
names_NF[[i]] <- unique_codes[!(unique_codes %in% names_BG[[i]] | 
                                    unique_codes %in% names_FG[[i]] | 
                                    unique_codes %in% names_SD[[i]])]}


HI_list <- list(BG = names_BG, FG = names_FG, SD = names_SD, NF = names_NF)

saveRDS(HI_list, "HI_list.RData")
HI_list <- readRDS("HI_list.RData")

# Read in different HI data
prob_BG <- readRDS("prob_BG.RData")
prob_FG <- readRDS("prob_FG.RData")
prob_SD <- readRDS("prob_SD.RData")
prob_NF <- readRDS("prob_NF.RData")

# Order data
order_HI_data <- function(prob_HI) {
  
  order_rows <- compare_strength$ID
  
  # Now reorder the dataframe
  prob_HI <- lapply(prob_HI, function (df) {
    
    df <- data.frame(Code = order_rows,
                     PropHI = df$HIprop[match(order_rows, df$Code)])
    return(df)
    
  })
  return(prob_HI)
}

prob_BG <- order_HI_data(prob_BG)
prob_FG <- order_HI_data(prob_FG)
prob_SD <- order_HI_data(prob_SD)
prob_NF <- order_HI_data(prob_NF)

# Combine the data
local_metrics_HI <- data.frame(ID = compare_strength$ID,
                               Period = c(rep("1-Before_HAB", nrow(compare_strength)), 
                                          rep("2-During_HAB", nrow(compare_strength)),
                                          rep("3-After_HAB", nrow(compare_strength))),
                               Degree = c(compare_strength$Before_HAB_degree, 
                                          compare_strength$During_HAB_degree,
                                          compare_strength$After_HAB_degree),
                               Strength = c(compare_strength$Before_HAB_strength, 
                                            compare_strength$During_HAB_strength,
                                            compare_strength$After_HAB_strength),
                               Prop_BG = c(prob_BG[[1]]$PropHI, prob_BG[[2]]$PropHI,
                                           prob_BG[[3]]$PropHI),
                               Prop_FG = c(prob_FG[[1]]$PropHI, prob_FG[[2]]$PropHI,
                                           prob_FG[[3]]$PropHI),
                               Prop_SD = c(prob_SD[[1]]$PropHI, prob_SD[[2]]$PropHI,
                                           prob_SD[[3]]$PropHI),
                               Prop_NF = c(prob_NF[[1]]$PropHI, prob_NF[[2]]$PropHI,
                                           prob_NF[[3]]$PropHI))

# Add HI_type column
local_metrics_HI_1 <- local_metrics_HI[local_metrics_HI$Period == "1-Before_HAB", ]
local_metrics_HI_2 <- local_metrics_HI[local_metrics_HI$Period == "2-During_HAB", ]
local_metrics_HI_3 <- local_metrics_HI[local_metrics_HI$Period == "3-After_HAB", ]
local_list <- list(local_metrics_HI_1, local_metrics_HI_2, local_metrics_HI_3)

## Initialize a new dataframe to store the results
result_df <- data.frame()
## Initialize a counter
for (p in 1:3) {
  counter <- 0
  result_df_new <- data.frame()
  for (i in HI_list) {
    # Increment the counter
    counter <- counter + 1
    index <- local_list[[p]]$ID %in% i[[p]]
    # Create a new row for each ID that falls into four categories
    new_rows <- local_list[[p]][index, ]
    new_rows$HI <- names(HI_list[counter])
    # Append the new rows to the result dataframe
    result_df_new <- rbind(result_df_new, new_rows)
  }
  result_df <- rbind(result_df, result_df_new)
}

# Make period a binary variable and HI a factorial
result_df$During <- ifelse(result_df$Period == "2-During_HAB", 1, 0)
result_df$After <- ifelse(result_df$Period == "3-After_HAB", 1, 0)
result_df$HI <- as.factor(result_df$HI)

# Save dataset
saveRDS(result_df, "result_df.RData")

# Look at demographics of HI dolphins
result_df <- readRDS("result_df.RData")
ILV_dem <- read.csv("ILV_dem.csv") # Read in demographics
ILV_dem <- ILV_dem[, c("Alias", "Sex", "BirthYear")]
colnames(ILV_dem) <- c("ID", "Sex", "Age")

result_df <- readRDS("result_df.RData")
result_df <- merge(result_df, ILV_dem[, c("ID", "Sex", "Age")], by = "ID", all.x = TRUE)


## How many did two HC
result_df <- result_df[result_df$HI != "NF",]
# Step 1: Group by ID and count the unique HIs
HI_counts <- result_df %>%
  group_by(ID) %>%
  summarise(unique_HIs = n_distinct(HI))

# Step 2: Filter individuals with exactly two unique HIs
two_HI_individuals <- HI_counts %>%
  filter(unique_HIs == 2)

# Step 3: Count these individuals
num_two_HI_individuals <- nrow(two_HI_individuals)

# Print the result
num_two_HI_individuals

## Sex
length(result_df$ID[result_df$HI == "BG" & result_df$Sex == "Female"])
length(result_df$ID[result_df$HI == "BG" & result_df$Sex == "Male"])
length(result_df$ID[result_df$HI == "FG" & result_df$Sex == "Female"])
length(result_df$ID[result_df$HI == "FG" & result_df$Sex == "Male"])
length(result_df$ID[result_df$HI == "SD" & result_df$Sex == "Female"])
length(result_df$ID[result_df$HI == "SD" & result_df$Sex == "Male"])
length(result_df$ID[result_df$HI == "NF" & result_df$Sex == "Female"])
length(result_df$ID[result_df$HI == "NF" & result_df$Sex == "Male"])


# Make table of demographics
# Define the vectors for metrics
HI <- c("BG", "FG", "SD")
Per <- c("1-Before_HAB", "2-During_HAB", "3-After_HAB")
sx <- c("Male", "Female")

# Initialize the empty data frame
dem.table <- data.frame(HI = character(),
                        Period = character(),
                        Male = integer(),
                        Female = integer(),
                        stringsAsFactors = FALSE)

# Loop through each combination of HI and Per
for (hi in HI) {
  for (per in Per) {
    male_count <- length(result_df$ID[result_df$HI == hi & result_df$Period == per & result_df$Sex == "Male"])
    female_count <- length(result_df$ID[result_df$HI == hi & result_df$Period == per & result_df$Sex == "Female"])
    
    # Add a new row to the demographic table
    dem.table <- rbind(dem.table, data.frame(HI = hi, Period = per, Male = male_count, Female = female_count))
  }
}

# Melt the demographic table for ggplot
library(reshape2)  # For melting the data frame
dem.table.melted <- melt(dem.table, id.vars = c("HI", "Period"), variable.name = "Sex", value.name = "Count")

# Create the plot
ggplot(dem.table.melted, aes(x = Period, y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ HI) +
  labs(x = "Study Period",
       y = "Number of Individual Dolphins",
       fill = "Sex") +
  scale_x_discrete(labels = c("1-Before_HAB" = "Before", 
                              "2-During_HAB" = "During", 
                              "3-After_HAB" = "After")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines

###########################################################################
# PART 3: Look at HI on Group Sizes ---------------------------------------------

# Read in GBI
gbi <- readRDS("gbi.RData")

# Get the average group size for each ID
group_list <- lapply(gbi, function(group_matrix) {
  
  # Calculate group size for each group
  individual_group_size <- rowSums(group_matrix)
  
  # Create empty vectors to store results
  ids <- character()
  avg_group_sizes <- numeric()
  
  # Iterate through each individual in the group
  for (i in 1:ncol(group_matrix)) {
    
    # Get the individual ID
    individual_id <- colnames(group_matrix)[i]
    
    # Calculate the group size for the individual
    group_size <- ifelse(group_matrix[, individual_id] == 1, 
                         individual_group_size, 0)
    
    # Calculate the average group size for the individual
    group_size_non_zero <- group_size[group_size != 0]
    avg_group_size <- mean(group_size_non_zero)
    
    # Append the results to vectors
    ids <- c(ids, individual_id)
    avg_group_sizes <- c(avg_group_sizes, avg_group_size)
  }
  
  # Create a data frame for the current group
  group_data <- data.frame(ID = ids,
                           Average_Group_Size = avg_group_sizes)
  
  return(group_data)
})

# Add HI list
result_df <- readRDS("result_df.RData")
result_df$Group_size <- ifelse(result_df$Period == "1-Before_HAB", 
                               group_list[[1]]$Average_Group_Size[match(result_df$ID, group_list[[1]]$ID)], 
                               ifelse(result_df$Period == "2-During_HAB",
                                      group_list[[2]]$Average_Group_Size[match(result_df$ID, group_list[[2]]$ID)], 
                                      group_list[[3]]$Average_Group_Size[match(result_df$ID, group_list[[3]]$ID)]))

# Change the factor levels and add factor for Period
result_df$HI <- factor(result_df$HI, levels = c("NF", "BG", "FG", "SD"))
result_df$Period <- as.factor(result_df$Period)
write.csv(result_df, "result_df.csv")

# Plot the HI behaviors and group sizes for every year
ggplot(result_df, aes(x = HI, y = Group_size, fill = HI)) +
  geom_boxplot(outlier.shape = NA) + # Remove outliers
  geom_jitter(aes(color = HI), width = 0.2, alpha = 0.5) + # Set jitter points color and transparency
  facet_wrap(~ Period, labeller = labeller(Period = c("1-Before_HAB" = "Before", 
                                                      "2-During_HAB" = "During", 
                                                      "3-After_HAB" = "After"))) +
  labs(x = "Human-centric Behavior", y = "Individuals' Average Group Size") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank())
# Unequal variance

# Look at the difference in HI groups
## Check assumptions of anova
lm_model <- lm(Group_size ~ HI, data = result_df)
### Extract Residuals
residuals <- resid(lm_model)
### Scale Residuals
scaled_residuals <- scale(residuals)
hist(scaled_residuals) # normal

## Run LMM test for Group size
# Fit the model again to use lmerTest functions
f.model <- lmer(Group_size ~ HI * Period + (1 | ID), data = result_df)

# Check residuals and random effects
plot(f.model)
qqnorm(resid(f.model))
qqline(resid(f.model))

# Get the summary with p-values
summary(f.model)

###########################################################################
# PART 4: Run Model ---------------------------------------------

# Read in data
result_df <- read.csv("result_df.csv")

# Make ID numeric
result_df$numeric_ID <- as.numeric(factor(result_df$ID))

# Make sure there is only one ID in each period
result_df <- result_df[!duplicated(result_df[c("Period", "ID")]), ]

# # Standardize the variables
# result_df$Strength <- scale(result_df$Strength)
# result_df$Group_size <- scale(result_df$Group_size)
# 
# # Create a composite score (simple sum of standardized variables)
# result_df$Social_centrality <- result_df$Strength + result_df$Group_size

# Check var and hists
var(result_df$Prop_BG)
hist(result_df$Prop_BG[result_df$Period == "1-Before_HAB"])
hist(result_df$Prop_BG[result_df$Period == "2-During_HAB"])
hist(result_df$Prop_BG[result_df$Period == "3-After_HAB"])

var(result_df$Prop_FG)
hist(result_df$Prop_FG[result_df$Period == "1-Before_HAB"])
hist(result_df$Prop_FG[result_df$Period == "2-During_HAB"])
hist(result_df$Prop_FG[result_df$Period == "3-After_HAB"])

var(result_df$Prop_SD)
hist(result_df$Prop_SD[result_df$Period == "1-Before_HAB"])
hist(result_df$Prop_SD[result_df$Period == "2-During_HAB"])
hist(result_df$Prop_SD[result_df$Period == "3-After_HAB"])

# How many HI dolphins in each period?
length(result_df$ID[result_df$HI != "NF" & result_df$Period == "1-Before_HAB"])
length(result_df$ID[result_df$HI != "NF" & result_df$Period == "2-During_HAB"])
length(result_df$ID[result_df$HI != "NF" & result_df$Period == "3-After_HAB"])

# Check assumptions of model
## Visualize relationship
plot(result_df$Strength ~ result_df$Prop_BG)
plot(result_df$Strength ~ result_df$Prop_FG)
plot(result_df$Strength ~ result_df$Prop_SD)

## Check distributions
hist(result_df$Strength) # normal
## Fit the linear mixed-effects model with the specified variance structure
test_model <- lm(Strength ~ Prop_BG * Period + Prop_FG * Period + Prop_SD * Period,
                  data = result_df)
summary(test_model)
plot(fitted(test_model, resid(test_model)))
plot(resid(test_model))
## Check for variance among groups
bartlett.test(Strength ~ HI, data = result_df) # equal
### If I have unequal variance
#' fit_brm.0 <- brm(bf(Strength ~ 1 + (1 | numeric_ID), sigma ~ 0 + BG + FG + SD), 
#' chains = 3, family = gaussian, data = result_df)
## Independent
durbinWatsonTest(test_model) # not independent

# Help STAN run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set priors
full_priors <- c(
  # Prior for Prop_BG
  set_prior("normal(0, 1)", class = "b", coef = "Prop_BG"),
  # Prior for Prop_FG
  set_prior("normal(0, 1)", class = "b", coef = "Prop_FG"),
  # Prior for Prop_SD
  set_prior("normal(0, 1)", class = "b", coef = "Prop_SD")
  )

# Models in brms
fit_sc.0 <- brm(Strength ~ 1 + (1 | numeric_ID),
                 chains = 4, iter = 4000, warmup = 3000, family = gaussian, data = result_df)
fit_sc.1 <- brm(Strength ~ Prop_BG + Prop_FG + Prop_SD + Period + (1 | numeric_ID), 
                 chains = 4, iter = 4000, warmup = 3000, family = gaussian, data = result_df)
fit_sc.2 <- brm(Strength ~
                   Prop_BG * Period + 
                   Prop_FG * Period +
                   Prop_SD * Period + 
                   (1 | numeric_ID),
                 chains = 4, iter = 4000, warmup = 3000, 
                 family = gaussian, data = result_df, prior = full_priors)

loo(fit_sc.0, fit_sc.1, fit_sc.2, compare = T)
saveRDS(fit_sc.0, "fit_sc.0.RData")
saveRDS(fit_sc.1, "fit_sc.1.RData")
saveRDS(fit_sc.2, "fit_sc.2.RData")
fit_sc.2 <- readRDS("fit_sc.2.RData")
summary(fit_sc.2)

# Check for model convergence
model <- fit_sc.2
plot(model)
pp_check(model) # check to make sure they line up

# Find the significance
posterior_samples <- as.data.frame(as.matrix(posterior_samples(model)))
coefficients <- colnames(posterior_samples)
mean(posterior_samples$`b_Prop_BG` < 0)
mean(posterior_samples$`b_Prop_FG` > 0)
mean(posterior_samples$`b_Prop_SD` > 0)
mean(posterior_samples$`b_Period3MAfter_HAB` > 0)
mean(posterior_samples$`b_Prop_BG:Period2MDuring_HAB` < 0)
mean(posterior_samples$`b_Prop_BG:Period3MAfter_HAB` < 0)
mean(posterior_samples$`b_Period2MDuring_HAB:Prop_FG` < 0)
mean(posterior_samples$`b_Period3MAfter_HAB:Prop_FG` < 0)
mean(posterior_samples$`b_Period2MDuring_HAB:Prop_SD` < 0)
mean(posterior_samples$`b_Period3MAfter_HAB:Prop_SD` < 0)

# Plot the posterior distribution
get_variables(model) # Get the names of the parameters

## Period Centrality
theme_update(text = element_text(family = "sans"))

# Create mcmc_areas plot
mcmc_plot <- mcmc_intervals(
  as.array(model), 
  pars = c("b_Period3MAfter_HAB", "b_Period2MDuring_HAB"),
  prob = 0.95, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  labs(
    title = "Posterior parameter distributions",
    subtitle = "with medians and 95% intervals"
  ) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_Period2MDuring_HAB" = "During HAB", 
    "b_Period3MAfter_HAB" = "After HAB"
  )
)

## BG
# Create mcmc_areas plot
mcmc_plot <- mcmc_intervals(
  as.array(model), 
  pars = c("b_Prop_BG:Period2MDuring_HAB", "b_Prop_BG:Period3MAfter_HAB",
           "b_Prop_BG"),
  prob = 0.90, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  labs(
    title = "Posterior parameter distributions",
    subtitle = "with medians and 95% intervals"
  ) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_Prop_BG" = "Begging/Provisioning",
    "b_Prop_BG:Period2MDuring_HAB" = "BG: During", 
    "b_Prop_BG:Period3MAfter_HAB" = "BG: After"
  )
)

## FG
mcmc_plot <- mcmc_intervals(
  as.array(model), 
  pars = c("b_Period2MDuring_HAB:Prop_FG", "b_Period3MAfter_HAB:Prop_FG",
           "b_Prop_FG"),
  prob = 0.95, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  labs(
    title = "Posterior parameter distributions",
    subtitle = "with medians and 95% intervals"
  ) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_Prop_FG" = "Fixed Gear Foraging",
    "b_Period2MDuring_HAB:Prop_FG" = "FG: During", 
    "b_Period3MAfter_HAB:Prop_FG" = "FG: After"
  )
)

## SD
mcmc_plot <- mcmc_intervals(
  as.array(model), 
  pars = c("b_Period2MDuring_HAB:Prop_SD", "b_Period3MAfter_HAB:Prop_SD",
           "b_Prop_SD"),
  prob = 0.95, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  labs(
    title = "Posterior parameter distributions",
    subtitle = "with medians and 95% intervals"
  ) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )

mcmc_plot + scale_y_discrete(
  labels = c(
    "b_Prop_SD" = "Scavenging/Depredating",
    "b_Period2MDuring_HAB:Prop_SD" = "SD: During", 
    "b_Period3MAfter_HAB:Prop_SD" = "SD: After"
  )
)


###########################################################################
# PART 5: Circular heat map ---------------------------------------------

# Read in rank data
result_df <- readRDS("result_df.RData")

# Make sure there is only one ID in each period
unique_ids <- result_df[!duplicated(result_df[c("Period", "ID")]), ]

# Use the aggregate function to sum the Eigen values by ID
sum_by_id <- aggregate(Strength ~ ID, data = result_df, FUN = sum)

# Make HI groups
HI_by_id <- aggregate(HI ~ ID, data = result_df, FUN = function(x) paste(unique(x), collapse = ","))

# Merge the two data frames
rank_data <- merge(sum_by_id, HI_by_id, all = TRUE)

# Generate all combinations of BG, SD, and FG
all_combinations <- c("BG", "SD", "FG", "NF")

# Function to find the combination for each HI value
find_combination <- function(hi_value) {
  hi_elements <- unlist(strsplit(hi_value, ","))
  combination <- all_combinations %in% hi_elements
  paste(all_combinations[combination], collapse = ",")
}

# Apply the function to each HI value
rank_data$HI_combination <- sapply(rank_data$HI, find_combination)

# Make HI as a factor 
rank_data$HI_combination <- as.factor(rank_data$HI_combination)

# Split rows based on the number of HI categories each ID has
split_data <- rank_data[rep(1:nrow(rank_data), lengths(strsplit(rank_data$HI, ","))), ]

# Splitting the values in the HI column
hi_values <- strsplit(as.character(rank_data$HI), ",")
split_data$HI <- unlist(hi_values)

# Put this back into result_df
result_data <- data.frame(ID = split_data$ID,
                        Centrality = split_data$Strength,
                        HI = split_data$HI)

# Normalize the centrality values to fit in the circle
unique_ids$Strength <- (unique_ids$Strength - 
                                      min(unique_ids$Strength)) / 
                                      (max(unique_ids$Strength) - 
                                         min(unique_ids$Strength))

# Convert "Period" column to a factor to ensure correct ordering
unique_ids$Period <- factor(unique_ids$Period, levels = c("1-Before_HAB", "2-During_HAB", "3-After_HAB"))

# Reshape the data from wide to long format for plotting
df_long <- reshape2::melt(unique_ids, id.vars = c("ID", "Period"), measure.vars = c("Strength"))

# Make a Rank sum
rank_sum <- data.frame(ID = sum_by_id$ID, Period = "Rank-Sum",
                       variable = "Strength", 
                       value = sum_by_id$Strength)
 
# Normalize strength
rank_sum$value <- (rank_sum$value - min(rank_sum$value)) / (max(rank_sum$value) - min(rank_sum$value))

# Create different data frames for each HI behavior
## BG
df_long_BG <- df_long[df_long$ID %in% result_data$ID[result_data$HI == "BG"], ]
# Subset IDs that are not in before and/or during periods
change_behav_BG <- result_df[result_df$HI == "BG", c("ID", "Period", "Strength")]
# Add scaled rank sum values
rank_sum_BG <- rank_sum[rank_sum$ID %in% unique(change_behav_BG$ID), c("ID", "Period", "value")]
change_behav_BG <- merge(change_behav_BG, rank_sum_BG, all = T)
# Order data by Period
df_long_BG <- df_long_BG[order(df_long_BG$Period), ]

# Filter data for IDs where hatches should be added
## Initialize hatch column with FALSE
hatch_vect <- NULL
## Loop through each unique period in change_behav_BG
for (period in unique(change_behav_BG$Period)) {
  # Determine if ID is found in change_behav_BG for this period
  ids_in_period <- change_behav_BG$ID[change_behav_BG$Period == period]
  hatch <- df_long_BG$ID[df_long_BG$Period == period] %in% ids_in_period
  hatch_vect <- c(hatch_vect, hatch)
}

df_long_BG$hatch <- hatch_vect

## Filter data to include only rows where hatch is TRUE
df_hatched <- df_long_BG[df_long_BG$hatch,]

# Graph ranked BG individuals
ggplot(rank_sum_BG, aes(x = reorder(ID, value), y = value, fill = value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis(option = "plasma") +  # Using the plasma color palette from the viridis package
  labs(
    x = "ID",
    y = "Scaled Total Centrality"
  ) +
  theme(panel.background = element_blank())

# Plot with hatches for IDs found in each period
ggplot(df_hatched, aes(x = ID, y = Period, fill = value)) +
  geom_tile() +
  geom_tile_pattern(data = df_long_BG[!df_long_BG$hatch, ], 
                    aes(pattern = ID), 
                    pattern = "stripe", 
                    pattern_fill = "black", 
                    pattern_angle = 45, 
                    pattern_density = 0.1, 
                    pattern_spacing = 0.025,
                    color = "white") +
  scale_fill_viridis(option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ID", y = "Period", fill = "Strength") +
  coord_polar()

# Map out centrality over time
## Add proportion data
df_long_BG <- merge(df_long_BG, result_df[,c("ID", "Period", "Prop_BG")], 
                    by = c("ID", "Period"), all.x = TRUE)
## Calculate rank of Prop_BG
df_long_BG$Prop_BG_Rank <- rank(df_long_BG$Prop_BG)
## Plot
ggplot(df_long_BG, aes(x = Period, y = value, group = ID, color = Prop_BG_Rank)) +
  geom_line(color = "black") +  
  geom_point(size = 3) +  
  labs(x = "Period",
       y = "Social Centrality",
       color = "BG Engagement") +
  scale_color_viridis_c(name = "Prop_BG Rank") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

df_long_BG$Rank <- ifelse(df_long_BG$Prop_BG_Rank == 1, "None", 
                          ifelse(df_long_BG$Prop_BG_Rank < 6, "Low", "High"))
ggplot(df_long_BG,
       aes(x = Period, stratum = Rank, alluvium = ID,
           fill = Rank, label = Rank)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")

## FG
df_long_FG <- df_long[df_long$ID %in% result_data$ID[result_data$HI == "FG"], ]
# Subset IDs that are not in before and/or during periods
change_behav_FG <- result_df[result_df$HI == "FG", c("ID", "Period", "Strength")]
# Add scaled rank sum values
rank_sum_FG <- rank_sum[rank_sum$ID %in% unique(change_behav_FG$ID), c("ID", "Period", "value")]
change_behav_FG <- merge(change_behav_FG, rank_sum_FG, all = T)
# Order data by Period
df_long_FG <- df_long_FG[order(df_long_FG$Period), ]

# Filter data for IDs where hatches should be added
## Initialize hatch column with FALSE
hatch_vect <- NULL
## Loop through each unique period in change_behav_FG
for (period in unique(change_behav_FG$Period)) {
  # Determine if ID is found in change_behav_FG for this period
  ids_in_period <- change_behav_FG$ID[change_behav_FG$Period == period]
  hatch <- df_long_FG$ID[df_long_FG$Period == period] %in% ids_in_period
  hatch_vect <- c(hatch_vect, hatch)
}
df_long_FG$hatch <- hatch_vect

# Filter data to include only rows where hatch is TRUE
df_hatched <- df_long_FG[df_long_FG$hatch,]

## Graph ranked FG individuals
ggplot(rank_sum_FG, aes(x = reorder(ID, value), y = value, fill = value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis(option = "plasma") +  # Using the plasma color palette from the viridis package
  labs(
    x = "ID",
    y = "Scaled Total Centrality"
  ) +
  theme(panel.background = element_blank())

# Plot with hatches for IDs found in each period
ggplot(df_hatched, aes(x = ID, y = Period, fill = value)) +
  geom_tile() +
  geom_tile_pattern(data = df_long_FG[!df_long_FG$hatch, ], 
                    aes(pattern = ID), 
                    pattern = "stripe", 
                    pattern_fill = "black", 
                    pattern_angle = 45, 
                    pattern_density = 0.1, 
                    pattern_spacing = 0.025,
                    color = "white") +
  scale_fill_viridis(option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ID", y = "Period", fill = "Strength") +
  coord_polar()

# Map out centrality over time
## Add proportion data
df_long_FG <- merge(df_long_FG, result_df[,c("ID", "Period", "Prop_FG")], 
                    by = c("ID", "Period"), all.x = TRUE)
## Calculate rank of Prop_BG
df_long_FG$Prop_FG_Rank <- rank(df_long_FG$Prop_FG)
## Plot
ggplot(df_long_FG, aes(x = Period, y = value, group = ID, color = Prop_FG_Rank)) +
  geom_line(color = "black") +  
  geom_point(size = 3) +  
  labs(x = "Period",
       y = "Social Centrality",
       color = "FG Engagement") +
  scale_color_viridis_c(name = "Prop_FG Rank") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

## SD
df_long_SD <- df_long[df_long$ID %in% result_data$ID[result_data$HI == "SD"], ]
# Subset IDs that are not in before and/or during periods
change_behav_SD <- result_df[result_df$HI == "SD", c("ID", "Period", "Strength")]
# Add scaled rank sum values
rank_sum_SD <- rank_sum[rank_sum$ID %in% unique(change_behav_SD$ID), c("ID", "Period", "value")]
change_behav_SD <- merge(change_behav_SD, rank_sum_SD, all = T)
# Order data by Period
df_long_SD <- df_long_SD[order(df_long_SD$Period), ]

# Filter data for IDs where hatches should be added
## Initialize hatch column with FALSE
hatch_vect <- NULL
## Loop through each unique period in change_behav_SD
for (period in unique(change_behav_SD$Period)) {
  # Determine if ID is found in change_behav_BG for this period
  ids_in_period <- change_behav_SD$ID[change_behav_SD$Period == period]
  hatch <- df_long_SD$ID[df_long_SD$Period == period] %in% ids_in_period
  hatch_vect <- c(hatch_vect, hatch)
}
df_long_SD$hatch <- hatch_vect

# Filter data to include only rows where hatch is TRUE
df_hatched <- df_long_SD[df_long_SD$hatch,]

## Graph ranked SD individuals
ggplot(rank_sum_SD, aes(x = reorder(ID, value), y = value, fill = value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis(option = "plasma") +  # Using the plasma color palette from the viridis package
  labs(
    x = "ID",
    y = "Scaled Total Centrality"
  ) +
  theme(panel.background = element_blank())

# Plot with hatches for IDs found in each period
ggplot(df_hatched, aes(x = ID, y = Period, fill = value)) +
  geom_tile() +
  geom_tile_pattern(data = df_long_SD[!df_long_SD$hatch, ], 
                    aes(pattern = ID), 
                    pattern = "stripe", 
                    pattern_fill = "black", 
                    pattern_angle = 45, 
                    pattern_density = 0.1, 
                    pattern_spacing = 0.025,
                    color = "white") +
  scale_fill_viridis(option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ID", y = "Period", fill = "Strength") +
  coord_polar()

# Map out centrality over time
## Add proportion data
df_long_SD <- merge(df_long_SD, result_df[,c("ID", "Period", "Prop_SD")], 
                    by = c("ID", "Period"), all.x = TRUE)
## Calculate rank of Prop_BG
df_long_SD$Prop_SD_Rank <- rank(df_long_SD$Prop_SD)
## Plot
ggplot(df_long_SD, aes(x = Period, y = value, group = ID, color = Prop_SD_Rank)) +
  geom_line(color = "black") +  
  geom_point(size = 3) +  
  labs(x = "Period",
       y = "Social Centrality",
       color = "SD Engagement") +
  scale_color_viridis_c(name = "Prop_SD Rank") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

###########################################################################
# PART 6: Multinetwork Plots ---------------------------------------------

# Only show IDs of HI dolphins
HI_list <- readRDS("HI_list.RData")
HI_list <- HI_list[-4] # Get rid of natural foragers

# Read in network object and strength values
net <- readRDS("net.RData")
result_df <- readRDS("result_df.RData")

# Make sure there is only one ID in each period
result_df <- result_df[!duplicated(result_df[c("Period", "ID")]), ]

# Normalize strength
result_df$Strength <- (result_df$Strength - min(result_df$Strength)) / (max(result_df$Strength) - min(result_df$Strength))

# ---Plot network---
# Set up the plotting area with 1 row and 2 columns for side-by-side plots
# Initialize a list to store layout information for each graph
layout_list <- vector("list", length(net))

# Function to generate random layout
generate_random_layout <- function(net) {
  num_nodes <- network.size(net)
  layout <- matrix(runif(2 * num_nodes), ncol = 2)
  return(layout)
}

# Generate random layouts for each network in the list
for (i in 1:length(net)) {
  layout_list[[i]] <- generate_random_layout(net[[i]])
}

# Optionally, you can visualize the layouts to check the distribution
par(mfrow = c(1, length(net)))
for (i in 1:length(net)) {
  plot.network(net[[i]], coord = layout_list[[i]], main = paste("Network", i))
}

# Generate a color palette
result_df$colors <- viridis(length(result_df$Strength), 
                  option = "plasma")[as.numeric(cut(result_df$Strength, 
                                                    breaks = length(result_df$Strength)))]

# Create an empty list with dimensions num_i x num_j
plot_list <- vector("list", 3)
for (i in seq_along(plot_list)) {
  plot_list[[i]] <- vector("list", 3)
}

# Loop through the list of graphs and plot them side by side
for (j in 1:length(HI_list)) {  # Loop through columns first
  
  # Extract layout for this graph
  combined_layout <- layout_list[[3]]
  counter <- 0
  
  for (i in 1:length(net)) {  # Loop through rows
    
    counter <- counter + 1
    
    # Get nodes for each behavior
    labeled_nodes <- HI_list[[j]][[i]]  # Fixed index here
    
    # Filter the dataframe for the period
    period_val <- unique(result_df$Period)[i]
    filtered_df <- subset(result_df, Period == period_val)
    
    # Match the filtered dataframe to the vertex names in the graph object
    matched_indices <- match(net[[i]] %v% "vertex.names", filtered_df$ID)
    filtered_df <- filtered_df[matched_indices, ]
    
    # Match node colors
    node_colors <- filtered_df$colors
    
    # Filter the dataframe for the period and HI
    filtered_df <- subset(filtered_df, ID %in% labeled_nodes)
    matched_indices <- match(net[[i]] %v% "vertex.names", filtered_df$ID)
    
    # Handle NA values in matched_colors by assigning "black"
    node.colors <- ifelse(is.na(matched_indices), "black", node_colors)
    
    # Now make node sizes vector
    node.sizes <- ifelse(node.colors == "black", 0.1, 1)
    
    # Create the plot
    net_plot <- ggnet2(net[[i]],
           mode = combined_layout,
           edge.size = "weight", # edge thickness
           edge.color = "grey",
           size = node.sizes,
           node.label = labeled_nodes,
           label.color = "white", 
           label.size = 2,
           node.color = node.colors,
           edge.alpha = 0.5
           ) 
    
    plot_list[[j]][[i]] <- net_plot
    
  }
}

plot_list[[1]][[1]] # Before BG
plot_list[[1]][[2]] # During BG
plot_list[[1]][[3]] # After BG
plot_list[[2]][[1]] # Before FG 
plot_list[[2]][[2]] # During FG
plot_list[[2]][[3]] # After FG
plot_list[[3]][[1]] # Before SD
plot_list[[3]][[2]] # During SD
plot_list[[3]][[3]] # After SD

# Plot the density plots for each period
my_colors <- c("#FC4E07", "#009E73", "#00AFBB")

plots_list <- list()

for (i in 1:length(unique(result_df$Period))) {
  
  period_to_plot <- unique(result_df$Period)[i] # each period
  filtered_df <- subset(result_df, Period == period_to_plot) # Separate data by period
  
  mean_nf <- mean(filtered_df$Strength[filtered_df$HI == "NF"], na.rm = TRUE) # Calculate mean for HI=="NF"
  
  plot <- ggplot(filtered_df[filtered_df$HI != "NF", ], aes(x = HI, y = Strength, fill = HI)) +
    geom_violin(trim = FALSE, alpha = 0.4) + # Create violin plot
    geom_boxplot(width=0.1, color="black", alpha=0.2) +
    geom_jitter(width = 0.1, alpha = 0.6) + # Add jittered points for visibility
    geom_hline(yintercept = mean_nf, linetype = "dashed", color = "black", linewidth = 1.5) + # Add horizontal line
    scale_fill_manual(values = my_colors) + # Use custom color palette
    theme_ipsum() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) 
  
  plots_list[[i]] <- plot
}


# Output plots
plots_list[[1]] # Before
plots_list[[2]] # During
plots_list[[3]] # After


# Plot the density plots for each HI
## Define colors
grey_shades <- c("#DDDDDD", "#BBBBBB", "#999999")

plots_list_HI <- list()

# Loop through each HI category
for (i in 1:(length(unique(result_df$HI))-1)) {
  
  HI_to_plot <- unique(result_df$HI)[i] # each HI category
  filtered_df <- subset(result_df, HI == HI_to_plot) # Separate data
  
  # Reorder the data
  filtered_df$Period <- factor(filtered_df$Period,
                         levels = c('3-After_HAB', '2-During_HAB', '1-Before_HAB'),
                         ordered = TRUE)
  
  # Plot the graphs
  plot <- ggplot(filtered_df, aes(x = Period, y = Strength, fill = as.factor(Period))) +
    geom_violin(trim = FALSE, alpha = 0.4) + # Create violin plot
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    geom_jitter(width = 0.1, alpha = 0.6) + # Add jittered points for visibility
    scale_fill_grey(start = 0, end = .9) +
    theme_ipsum() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    guides(fill = guide_legend(title = "Period")) + # Add legend for Period
    coord_flip() 
  
  plots_list_HI[[i]] <- plot
}

# Output plots
plots_list_HI[[1]] # BG
plots_list_HI[[2]] # FG
plots_list_HI[[3]] # SD

