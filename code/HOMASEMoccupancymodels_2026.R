# HOMA_climate-filtering-occupancy


### Hierarchical single-species, single-season occupancy model of local, microclimate factors


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
HOMA <- read.csv("3.Analysis/Talussurveys_20212022_HOMA_2025FinalVersion.csv") %>% 
  select(-X, -ID) %>% pivot_wider(names_from = "Type", values_from = "Obs") %>% 
  mutate(Scat = as.integer(ifelse(Scat >= 1, 1, 0)), row_num = row_number())

HOMAdet_data <- as.data.frame(scale(HOMA[,c("cloud", "avgdiam", "Wind", "Talus.Score", "gram", 
                                          "forb", "gramforb", "snow.depth", "slope", "perm", "ros", 
                                          "sum.days.over.21", "wint.days.below.neg5", "max.sum",
                                          "min.wint", "wint.days.above.0", "evi")])) %>% 
  mutate(row_num = row_number()) %>% 
  left_join(HOMA[,c("row_num", "Sp", "y", "x", "doy", "Year", "Month", "Day", "Start", 
                       "Indiv", "Scat", 'Vocal')], ., by = "row_num") %>% 
  filter(!is.na(avgdiam) & !is.na(gramforb))
# filtered out the rows with NA values (3 in total) left with 83 sites

# Fit the model
J <- 83
K <- 3

# Organize the data
jags_data = list(J = J, K = K, det_data = HOMAdet_data[,c("Indiv", "Scat", 'Vocal')], 
                 avgdiam = HOMAdet_data$avgdiam, gramforb = vdet_data$gramforb,
                 cloud = HOMAdet_data$cloud)

HOMA.jags_data = list(J = J, K = K, y = HOMAdet_data[, c("Scat", "Indiv", "Vocal")], 
                      Z = ifelse(rowSums(HOMAdet_data[, c("Scat", "Indiv", "Vocal")], na.rm = T) > 0, 1, NA),
                      cloud = HOMAdet_data$cloud, avgdiam = HOMAdet_data$avgdiam, gramforb = HOMAdet_data$gramforb)

params = c('a0', 'a1', 'a2', 'a3', 'Z', 'p', 'psi', 'D.obs', 'D.rep')

# Specify initial values
inits <- function() {
  a0 <- rnorm(1, 0, .5)
  a1 <- rnorm(1, 0, .5)
  a2 <- rnorm(1, 0, .5)
  a3 <- rnorm(1, 0, .5)
  b0 <- rnorm(3, 0, .5)
  
  return(list(a0 = a0, a1 = a1, a2 = a2, a3 = a3, bO1 = b0[1], bO2 = b0[2], bO3 = b0[3]))
}

## JAGS code ______________________________________________________________________
# Now we fit the model
sink("HOMA_simpleoccupancy.jags")
cat("

model{

  # Priors
  a0 ~ dlogis(0,1)    # uniform non-informative
  a1 ~ dnorm(0.15, 0.5)
  a2 ~ dnorm(0.15, 0.5)
  a3 ~ dnorm(-0.25,0.5)
  
  b0[1] ~ dlogis(0,1)  # scat
  b0[2] ~ dlogis(0,1)  # individual sightings
  b0[3] ~ dlogis(0,1)  # vocalizations


  # Likelihood
  for(i in 1:J){
  
    # model for occupancy given heat and freezing
    Z[i] ~ dbern(psi[i])        # Z = 1 if true state occupied
    logit(psi[i]) <- a0 + a1*gramforb[i] + a2*avgdiam[i]

    # Generate Replicated Data
    Z.rep[i] ~ dbern(psi[i])

    # Calculate Residuals
    res[i] <- Z[i] - psi[i]
    res.new[i] <- Z.rep[i] - psi[i]
    
    # model for detection
    for(j in 1:K){
      logit(p[i,j]) = b0[j] + a3*cloud[i]
      y[i,j] ~ dbern(Z[i]*p[i,j])
    }
  }
  
  # Posterior predictive check for residuals
  D.obs <- sum(res[])
  D.rep <- sum(res.new[])
}
", fill = TRUE)
sink()



HOMA_testocc_sim <- jagsUI::jags(data = HOMA.jags_data, 
                                 inits=inits,
                                 model.file = "HOMA_simpleoccupancy.jags",
                                 parameters.to.save = params,
                                 n.chains = 25,              # chains to run
                                 n.iter = 15000,             # iterations including burn in
                                 n.thin = 10,                # thinning rate
                                 n.burnin = 5000,            # burn in
                                 parallel = T)

mcmcOutput::diagPlot(HOMA_testocc_sim) # diagnostic plots
print(HOMA_testocc_sim)

# check model fit
pp.check(HOMA_testocc_sim, observed = 'D.obs', simulated = 'D.rep')

# save the model output as R data
save(HOMA_testocc_sim, file = "~/UMT_2019-2021/PhD_Denali_Alpine_Wildlife/Chapter1_SpatiotemporalWildlifeOccupancy/ECY25-1309_Re-resubmission Materials/HOMA_occ_sim.Rdata")

# also save the dataframe of summary stats
write.csv(HOMA_testocc_sim$summary, "~/UMT_2019-2021/PhD_Denali_Alpine_Wildlife/Chapter1_SpatiotemporalWildlifeOccupancy/ECY25-1309_Re-resubmission Materials/HOMA_occ_sim_summary.csv")



