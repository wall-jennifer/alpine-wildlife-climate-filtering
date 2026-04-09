# HOMA_climate-drivers-SEM-occupancy


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


# Full SEM Model with the posterior predictive checks coded in

HOMA.climate.jags_data = list(J = J, K = K, y = HOMAdet_data[, c("Scat", "Indiv", "Vocal")], 
                              Z = ifelse(rowSums(HOMAdet_data[, c("Scat", "Indiv", "Vocal")], na.rm = T) > 0, 1, NA),
                              ros = HOMAdet_data$ros, slope = HOMAdet_data$slope,
                              evi = HOMAdet_data$evi, heat = HOMAdet_data$sum.days.over.21,
                              cold = HOMAdet_data$wint.days.below.neg5)

params.cl = c('a0', 'a1', 'a2', 'a3', 'a4', 'a5','g0', 'g1', 'g2', 'g3',
              'Z', 'p', 'psi', 'D.obs', 'D.rep', 'N')

# Specify initial values
inits.cl <- function() {
  a0 <- rnorm(1, 0, .5)
  a1 <- rnorm(1, 0, .5)
  a2 <- rnorm(1, 0, .5)
  a3 <- rnorm(1, 0, .5)
  a4 <- rnorm(1, 0, .5)
  a5 <- rnorm(1, 0, .5)
  g0 <- rnorm(1, 0, .5)
  g1 <- rnorm(1, 0, .5)
  g2 <- rnorm(1, 0, .5)
  g3 <- rnorm(1, 0, .5)
  b0 <- rnorm(3, 0, .5)
  
  return(list(a0 = a0, a1 = a1, a2 = a2, a3 = a3, a4 = a4, a5 = a5,
              g0 = g0, g1 = g1, g2 = g2, g3 = g3, bO1 = b0[1], bO2 = b0[2], bO3 = b0[3]))
}


sink("SEMOcc_ClimateSimulation_HOMA2025.jags")
cat("

model{
  # Priors
  a1 ~ dnorm(-0.15,0.25)  # heat ~ a1*slope
  a2 ~ dnorm(-0.15,0.25)  # evi ~ a3*slope + a4*heat
  a3 ~ dnorm(0.15, 0.5)   # evi ~ a3*slope + a4*heat
  a4 ~ dlogis(0,1)        # Freezing ~ a4*ros + a5*cold
  a5 ~ dlogis(0,1)        # Freezing ~ a4*ros + a5*cold
  
  b0[1] ~ dlogis(0,1)     # scat
  b0[2] ~ dlogis(0,1)     # individual sightings
  b0[3] ~ dlogis(0,1)     # vocalizations

  g0 ~ dlogis(0,1)        # psi ~ intercept
  g1 ~ dnorm(-0.25, 0.5)  # psi ~ evi
  g2 ~ dnorm(-0.15, 0.5)  # psi ~ heat
  g3 ~ dnorm(-0.25, 0.5)  # psi ~ freezing / winter temps 

  for (k in 1:3){
    sigma[k] ~ dgamma(1,1)
    tau[k] = 1/(sigma[k] * sigma[k])
  }
  
  # Likelihood
  for(i in 1:J){
  
    # model for heat given slope (no intercept bc covariates 0 centered)
    heat[i] ~ dnorm(a1*slope[i], tau[1])
    
    # model 
    evi[i] ~ dnorm(a2*slope[i] + a3*heat[i], tau[2])
    
    # model 
    F[i] <- a4*ros[i] + a5*cold[i] + tau[3]
  
    # model for occupancy given heat and freezing
    Z[i] ~ dbern(psi[i])        # Z = 1 if true state occupied
    logit(psi[i]) <- g0 + g1*evi[i] + g2*heat[i] + g3*F[i]
    
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

HOMA_testocc_clim <- jagsUI::jags(data = HOMA.climate.jags_data, 
                                  inits=inits.cl,
                                  model.file = "SEMOcc_ClimateSimulation_HOMA2025.jags",
                                  parameters.to.save = params.cl,
                                  n.chains = 25,              # chains to run
                                  n.iter = 15000,             # iterations including burn in
                                  n.thin = 10,                # thinning rate
                                  n.burnin = 5000,            # burn in
                                  parallel = T)

mcmcOutput::diagPlot(HOMA_testocc_clim) # diagnostic plots
print(HOMA_testocc_clim)

# check model fit
pp.check(HOMA_testocc_clim, observed = 'D.obs', simulated = 'D.rep')

# save the model output as R data
save(HOMA_testocc_clim, file = "~/UMT_2019-2021/PhD_Denali_Alpine_Wildlife/Chapter1_SpatiotemporalWildlifeOccupancy/ECY25-1309_Re-resubmission Materials/HOMA_SEMocc_clim.Rdata")

# also save the dataframe of summary stats
write.csv(HOMA_testocc_clim$summary, "~/UMT_2019-2021/PhD_Denali_Alpine_Wildlife/Chapter1_SpatiotemporalWildlifeOccupancy/ECY25-1309_Re-resubmission Materials/HOMA_SEMocc_clim_summary.csv")
