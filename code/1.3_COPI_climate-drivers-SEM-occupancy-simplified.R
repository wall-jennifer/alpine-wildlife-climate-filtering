# COPI_climate-drivers-SEM-occupancy-simplified


### Hierarchical occupancy model nested within a structural equation model with latent and observed variables


# Metadata:
# Visibility in km, >2 km visibility listed as 2
# Cloud: in % cloud cover of visible sky
# Wind: scaled from ________
# Talus.Score: talus scored on stability from 1 - 5
# Avg. Diameter: in cm, avg. scores, for small and large rocks, we chose predominant rock size (whichever was >50%)
# Vegetation: scored based on Daubenmire's scale-- 0 = 0, trace = 1, 1-5 = 2, 5-25 = 3, 25-50 = 4, 50-75 = 5, 75-95 = 6, 95-99 = 7, 100 = 8
# Species: X for not present, P for present, O for old, F for fresh


## JAGS occupancy model with different sign types
## By marginalizing the discrete occupancy state, Z, and including the likelihood functions,
## you can potentially improve processing time, effective sample sizes, 
## and make a model that’s transferable between JAGS and Stan

# Import libraries
library(tidyverse)
library(jagsUI)     # interfaces with JAGS
library(mcmcOutput) # interprets results

# Import and organize data
COPI <- read.csv("3.Analysis/Talussurveys_20212022_COPI_2025FinalVersion.csv") %>% 
  select(-X, -ID) %>% pivot_wider(names_from = "Type", values_from = "Obs") %>% 
  mutate(Scat = as.integer(ifelse(Scat >= 1, 1, 0)), row_num = row_number())

COPIdet_data <- as.data.frame(scale(COPI[,c("cloud", "avgdiam", "Wind", "Talus.Score", "gram", 
                                            "forb", "gramforb", "snow", "slope", "perm", "ros", 
                                            "sum.days.over.21", "wint.days.below.neg5", "max.sum",
                                            "min.wint", "wint.days.above.0", "evi")])) %>% 
  mutate(row_num = row_number()) %>% 
  left_join(COPI[,c("row_num", "Sp", "y", "x", "doy", "Year", "Month", "Day", "Start", 
                    "Indiv", "Scat", 'Vocal', "Haypile")], ., by = "row_num") %>% 
  filter(!is.na(avgdiam) & !is.na(gramforb))
# filtered out the rows with NA values (3 in total) left with 83 sites

data.COPI <- read.csv("3.Analysis/Talussurveys_20212022_COPI_2025_withCoOccurrences.csv")
scaled.COPI <- as.data.frame(scale(data.COPI[,c("gramforb", "avgdiam", "cloud")]))

# Fit the model
J <- 83
K <- 3

# Simplified SEM Model with the posterior predictive checks coded in
climate.COPI <- as.data.frame(scale(data.COPI[,c("snow", "slope", "perm", "heat")]))

COPI.climate.jags_data = list(J = J, K = K, y = data.COPI[, c("Scat", "Indiv", "Vocal")], 
                              Z = ifelse(rowSums(data.COPI[, c("Scat", "Indiv", "Vocal")], na.rm = T) > 0, 1, NA),
                              snow = climate.COPI$snow, slope = climate.COPI$slope,
                              perm = climate.COPI$perm, heat = climate.COPI$heat)


params.simp = c('a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6','g0', 'g1', 'g2',
                'Z', 'p', 'psi', 'D.obs', 'D.rep', 'N')

# Specify initial values
inits.simp <- function() {
  a4 <- rnorm(1, 0, .5)
  g0 <- rnorm(1, 0, .5)
  g1 <- rnorm(1, 0, .5)
  g2 <- rnorm(1, 0, .5)
  g3 <- rnorm(1, 0, .5)
  b0 <- rnorm(3, 0, .5)
  
  return(list(a4 = a4, g0 = g0, g1 = g1, g2 = g2, g3 = g3,
              bO1 = b0[1], bO2 = b0[2], bO3 = b0[3]))
}


sink("SEMOcc_SimplifiedSimulation_COPI2025.jags")
cat("

model{
  # Priors
  a4 ~ dnorm(-0.25, 0.5)  # perm ~ a4*heat

  b0[1] ~ dlogis(0,1)     # scat
  b0[2] ~ dlogis(0,1)     # individual sightings
  b0[3] ~ dlogis(0,1)     # vocalizations

  g0 ~ dlogis(0,1)        # psi ~ intercept
  g1 ~ dnorm(0.3, 0.5)    # psi ~ heat
  g2 ~ dnorm(0.25,0.15)   # psi ~ snow
  g3 ~ dnorm(0.25,0.15)   # psi ~ perm

  for (k in 1:1){
    sigma[k] ~ dgamma(1,1)
    tau[k] = 1/(sigma[k] * sigma[k])
  }
  
  # Likelihood
  for(i in 1:J){
    
    # model for permafrost given slope and heat
    perm[i] ~ dnorm(a4*heat[i], tau[1])
  
    # model for occupancy given heat and freezing
    Z[i] ~ dbern(psi[i])        # Z = 1 if true state occupied
    logit(psi[i]) <- g0 + g1*heat[i] + g2*snow[i] + g3*perm[i]
    
    # Generate Replicated Data
    Z.rep[i] ~ dbern(psi[i])

    # Calculate Residuals
    res[i] <- Z[i] - psi[i]
    res.new[i] <- Z.rep[i] - psi[i]
    
    # model for detection
      for(j in 1:K){
      logit(p[i,j]) = b0[j]
      y[i,j] ~ dbern(Z[i]*p[i,j])
    }
  }
  
  # Posterior predictive check for residuals
  D.obs <- sum(res[])
  D.rep <- sum(res.new[])
}
",fill = TRUE)
sink()

COPI_testocc_simplified <- jagsUI::jags(data = COPI.climate.jags_data, 
                                        inits=inits.simp,
                                        model.file = "SEMOcc_SimplifiedSimulation_COPI2025.jags",
                                        parameters.to.save = params.simp,
                                        n.chains = 25,              # chains to run
                                        n.iter = 15000,             # iterations including burn in
                                        n.thin = 10,                # thinning rate
                                        n.burnin = 5000,            # burn in
                                        parallel = T)

mcmcOutput::diagPlot(COPI_testocc_simplified) # diagnostic plots
print(COPI_testocc_simplified)

# check model fit
pp.check(COPI_testocc_simplified, observed = 'D.obs', simulated = 'D.rep')

# save the model output as R data
#save(COPI_testocc_simplified, file = "path-here.Rdata")

# also save the dataframe of summary stats
#write.csv(COPI_testocc_clim$summary, "path-here.csv")

