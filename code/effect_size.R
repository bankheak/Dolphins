# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# EFFECT SIZES
###########################################################################

# Set working directory here
setwd("C:/Users/bankh/My_Repos/Dolphins/data")

# Load all necessary packages
require(asnipe) # get_group_by_individual--Damien Farine
# Could do permutations
require(assocInd)
require(vegan)
# Run multiple cores for faster computing
require(doParallel)
require(parallel)
require(foreach)

source("OR3_NBDA code 1.2.15.r") # sources functions
source("OR2_sensitivity functions.r") # sources functions 

###########################################################################
# PART 1: Organize Raw Data ---------------------------------------------

# Read in & combine files
firstgen_data <- read.csv("firstgen_data.csv")
secondgen_data <- read.csv("secondgen_data.csv")
orig_data <- rbind(firstgen_data, secondgen_data)

# Make date into a date class
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))
write.csv(orig_data, "orig_data.csv")

# Group each individual by date and sighting
group_data <- cbind(orig_data[,c("Date","Sighting","Code","Year")]) 
group_data <- subset(group_data, subset=c(group_data$Code != "None"))
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- cbind(group_data[,3:5]) # Subset ID and group #

# Make a list of only one year per dataframe
years <- unique(group_data$Year)
list_years <- list()
for (i in 1:length(years)) {
  list_years[[i]] <- subset(group_data, subset=c(group_data$Year == years[i]))
}    

# Gambit of the group index
gxi <- list()
for (y in 1:length(years)) {
  gxi[[y]] <- get_group_by_individual(list_years[[y]][,c("Code", "Group")], data_format = "individuals")
}
saveRDS(gxi, file="gxi.RData")
gxi <- readRDS("gxi.RData")

# Create association matrix
source("../code/functions.R") # SRI & null permutation

n.cores <- detectCores()
system.time({
  registerDoParallel(n.cores)
  nxn <- list()
  for (i in 1:length(years)) {
    nxn[[i]] <- as.matrix(SRI.func(gxi[[i]]))
  }
})

# End parallel processing
stopImplicitCluster()

saveRDS(nxn, file="nxn_full.RData")

###########################################################################
# PART 2: Calculate Error in Association Sample Size ---------------------------------------------

source("../code/functions.R") # together_apart and NBDA_error

# Extract how many times each dyad has been seen together and how many times they have been seen apart
## These two matrices will be used to create a social network with observational error using a Bayesian approach
year <- 21
together_apart_output <- together_apart(sightings=gxi[[year]])
saveRDS(together_apart_output, "together_apart_output.RData")
write.csv(together_apart_output[,,1], file="../data/together_matrix.csv") # save 'together' matrix
write.csv(together_apart_output[,,2], file="../data/apart_matrix.csv") # save 'apart' matrix

together_apart_output <- readRDS("together_apart_output.RData")
together_matrix <- read.csv("together_matrix.csv")
apart_matrix <- read.csv("apart_matrix.csv")

# Determine how many times individuals have been seen and choose cut-off points
cutoff <- sort(unique(colSums(gxi[[year]]))) # extract colSums (=number of sightings)
cutoff <- cutoff[-length(cutoff)] # remove the maximum number of sightings  
cutoff <- cutoff[-length(cutoff)] # remove second last as well


# All done in the HPC ------------------------------------------------------

# Run a first simulation with dropping learners and s=8 (social learning)
  sim1 <- sensitivity_NBDA_ind_error(x=together_apart_output, # with s=8 and dropping learners
                                     sightings=gxi[[year]],
                                     cutoff=cutoff, 
                                     association_index="SRI", 
                                     iterations=10000, 
                                     s=8, 
                                     num_ind_learn=20, 
                                     cores=2,
                                     keep_learners = FALSE,
                                     delta_AICc = 2)

save(sim1, file="../data/sensitivity_NBDA_with_error_false neg_drop learners.RData") # save the sim1 object 
write.csv(sim1$raw,"../data/sensitivity_NBDA_with_error_false neg_drop learners_RAW.csv" )
write.csv(sim1$summary,"../data/sensitivity_NBDA_with_error_false neg_drop learners_SUMMARY.csv" )

sim1 <- readRDS("../data/sensitivity_NBDA_with_error.RData")

# Run a second simulation with keeping learners and s=8 (social learning)
sim2 <- sensitivity_NBDA_ind_error(x=together_apart_output, # with s=8 and keeping learners
                                   sightings=gxi[[year]],
                                   cutoff=cutoff, 
                                   association_index="SRI", 
                                   iterations=10000, 
                                   s=8, 
                                   num_ind_learn=20, 
                                   cores=2,
                                   keep_learners = TRUE,
                                   delta_AICc=2)

save(sim2, file="sensitivity_NBDA_with_error_false neg_keep learners.RData") # save the sim2 object 
write.csv(sim2$raw,"sensitivity_NBDA_with_error_false neg_keep learners_RAW.csv" )
write.csv(sim2$summary,"sensitivity_NBDA_with_error_false neg_keep learners_SUMMARY.csv" )

# Run a third simulation with dropping learners and s=0 (asocial learning)
sim3 <- sensitivity_NBDA_ind_error(x=together_apart_output, # with s=0 and dropping learners
                                   sightings=gxi[[year]],
                                   cutoff=cutoff, 
                                   association_index="SRI", 
                                   iterations=10000, 
                                   s=0, 
                                   num_ind_learn=20, 
                                   cores=2,
                                   keep_learners = FALSE,
                                   delta_AICc = 2)
save(sim3, file="sensitivity_NBDA_with_error_false pos_drop learners.RData") # save the sim3 object 
write.csv(sim3$raw,"sensitivity_NBDA_with_error_false pos_drop learners_RAW.csv" )
write.csv(sim3$summary,"sensitivity_NBDA_with_error_false pos_drop learners_SUMMARY.csv" )

# Run a fourth simulation with keeping learners and s=0 (asocial learning)
sim4 <- sensitivity_NBDA_ind_error(x=together_apart_output, # with s=0 and keeping learners
                                   sightings=gxi[[year]],
                                   cutoff=cutoff, 
                                   association_index="SRI", 
                                   iterations=10000, 
                                   s=0, 
                                   num_ind_learn=20, 
                                   cores=2,
                                   keep_learners = TRUE,
                                   delta_AICc = 2)
save(sim4, file="sensitivity_NBDA_with_error_false pos_keep learners.RData") # save the sim4 object 
write.csv(sim4$raw,"sensitivity_NBDA_with_error_false pos_keep learners_RAW.csv" )
write.csv(sim4$summary,"sensitivity_NBDA_with_error_false pos_keep learners_SUMMARY.csv" )

# Next take results from the HPC ------------------------------------------------
sim1 <- readRDS("neg_drop.RData")
sim2 <- readRDS("neg_keep.RData")
sim3 <- readRDS("pos_drop.RData")
sim4 <- readRDS("pos_keep.RData")

###########################################################################
# PART 3: Plot Simulated Results ---------------------------------------------

#----------------------------------------------------------------------------------------------------
# plot results of sim1 

windowsFonts(A = windowsFont("Times New Roman")) 

par(mfrow=c(1,2))

par(mar=c(5.1,4,3,1)) # set figure margins
plot(sim1$summary[,"cutoff"], sim1$summary[,"perc_social_learning_supported"], # plot the proportion of social models supported against the cut-off points. 
     type="b",
     ylim=c(0,100),
     xlab="Cut-off points", 
     ylab="Percent",
     family="A",
     main="a) dropping learners",
     cex.main=1)

par(new=TRUE) # prevent R from opening a new figure

plot(sim1$summary[,"cutoff"], sim1$summary[,"within_95%_CI"], # plot estimates for s falling within 95% CI into the same plot
     type="b",
     xaxt="n", 
     yaxt="n",
     xlab="",
     ylab="",
     col="red",
     pch=2,
     ylim=c(0,100),
     family="A") # set limits of y-axis to match your minimum and maximum estimates

par(family="serif")
legend(x=1.5,y=20, pch=c(1,2), legend=c("Social models supported","s withing 95% C.I."), # add a legend
       lty=c(1,4), col=c("black","red"),
       bty="n", 
       cex=0.9)
#----------------------------------------------------------------------------------------------------
# plot results of sim2
par(mar=c(5.1,2,3,3)) # set figure margins
plot(sim2$summary[,"cutoff"], sim2$summary[,"perc_social_learning_supported"], # plot the proportion of social models supported against the cut-off points. 
     type="b",
     ylim=c(0,100),
     xlab="Cut-off points", 
     ylab="",
     main="b) retaining learners",
     cex.main=1)

par(new=TRUE) # prevent R from opening a new figure

plot(sim2$summary[,"cutoff"], sim2$summary[,"within_95%_CI"], # plot estimates of s falling within 95% CI into the same plot 
     type="b",
     xaxt="n", 
     yaxt="n",
     xlab="",
     ylab="",
     col="red",
     pch=2,
     ylim=c(0,100)) # set limits of y-axis to match your minimum and maximum estimates 


par(family="serif")
legend(x=1.5,y=20, pch=c(1,2), legend=c("Social models supported", "s within 95% C.I."), # add a legend
       lty=c(1,4), col=c("black","red"),
       bty="n", 
       cex=0.9)
#----------------------------------------------------------------------------------------------------
# plot results of sim3
windowsFonts(A = windowsFont("Times New Roman")) 

par(mfrow=c(1,2))

par(mar=c(5.1,4,3,1)) # set figure margins
plot(sim3$summary[,"cutoff"], sim3$summary[,"perc_social_learning_supported"], # plot the proportion of social models supported against the cut-off points. 
     type="b",
     ylim=c(0,100),
     xlab="Cut-off points", 
     ylab="Percent",
     family="A",
     main="a) dropping learners",
     cex.main=1)

par(new=TRUE) # prevent R from opening a new figure

plot(sim3$summary[,"cutoff"], sim3$summary[,"within_95%_CI"], # plot estimates of s falling within 95% CI into the same plot 
     type="b",
     xaxt="n", 
     yaxt="n",
     xlab="",
     ylab="",
     col="red",
     pch=2,
     ylim=c(0,100),
     family="A") # set limits of y-axis to match your minimum and maximum estimates 

par(family="serif")
legend(x=1.5,y=50, pch=c(1,2), legend=c("Social models supported", "s within 95% C.I."), # add a legend
       lty=c(1,4), col=c("black","red"),
       bty="n", 
       cex=0.9)
#----------------------------------------------------------------------------------------------------
# plot results of sim4
par(mar=c(5.1,2,3,3)) # set figure margins
plot(sim4$summary[,"cutoff"], sim4$summary[,"perc_social_learning_supported"], # plot the proportion of social models supported against the cut-off points. 
     type="b",
     ylim=c(0,100),
     xlab="Cut-off points", 
     ylab="Percent",
     family="A",
     main="b) retaining learners",
     cex.main=1)

par(new=TRUE) # prevent R from opening a new figure

plot(sim4$summary[,"cutoff"], sim4$summary[,"within_95%_CI"], # plot estimates of s falling within 95% CI into the same plot 
     type="b",
     xaxt="n", 
     yaxt="n",
     xlab="",
     ylab="",
     col="red",
     pch=2,
     ylim=c(0,100),
     family="A") # set limits of y-axis to match your minimum and maximum estimates 

par(family="serif")
legend(x=1.5,y=50, pch=c(1,2), legend=c("Social models supported", "s within 95% C.I."), # add a legend
       lty=c(1,4), col=c("black","red"),
       bty="n", 
       cex=0.9)
#----------------------------------------------------------------------------------------------------


