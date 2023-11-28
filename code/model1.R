
# Load the parallel package
library(doParallel)
library(nimble)
source('../code/attach.nimble_v2.R')

# Read in social association matrix and listed data
nimble.data <- readRDS("../data/nimble.data.RData") 
nimble.constants <- readRDS("../data/nimble.constants.RData") 

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
    Obs.Err[p] ~ T(dt(mu=0, sigma=1, df=1), 0, ) 
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
        
        # Process Model
        logit(SRI.Exp[i, j, p]) <- IJ_Random[i, j] + 
          HRO[i, j]*HRO_Effect + SEX.SIM[i, j]*SEX_Effect + AGE_Diff[i, j]*AGE_Effect + #GR[i, j]*GR_Effect +
          BP[i, j, p]*BP_Effect[p] +  FG[i, j, p]*FG_Effect[p] +  SD[i, j, p]*SD_Effect[p]
        
        # Observation Model (Likelihood) # a = mean^2*(1-mean)/sd^2-mean, b = mean*(1-mean)^2/sd^2+mean-1
        SRI[i, j, p] ~ dbeta(shape1 = (SRI.Exp[i, j, p]*(SRI.Exp[i, j, p]*(1-SRI.Exp[i, j, p])/(Obs.Err[p]^2)-1)), 
                             shape2 = ((1-SRI.Exp[i, j, p])*(SRI.Exp[i, j, p]*(1-SRI.Exp[i, j, p])/(Obs.Err[p]^2)-1)))
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




attach.nimble(mcmc.output$samples)
