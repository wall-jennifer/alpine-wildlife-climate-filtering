# Pika and marmot occupancy paper title (Previously: Can marmots and pika take the heat?)
### Jennifer L. Wall and Jedediah F. Brodie

## Abstract
TBD. 

## Code
  1. [0.1_prior-estimation.r](/code/0.1_prior-estimation.r) Prior development and estimation.
  2. [Link to data]: Code for the simulation I ran of the various GLMs.
  3. [Link to data]: Simulated structural equation / occupancy model
  4. Collared pika *(Ochotona collaris)* models:
     
      - 1.1_COPI_climate-drivers-SEM-occupancy: Hierarchical occupancy model nested within a structural equation model with latent and observed variables.
      - 1.2_COPI_climate-filtering-occupancy: Hierarchical single-species, single-season occupancy model of local, microclimate factors, including: average rock diameter and vegetation scores, based on ...
      - 1.3_COPI_climate-drivers-SEM-occupancy-simplified:
      - 1.4_COPI_climate-drivers-SEM-Occupancy-cooccurrence"
    
  5. Hoary marmot *(Marmota caligata)* models: 
     
      - 2.1_HOMA_climate-drivers-SEM-occupancy: Hierarchical occupancy model nested within a structural equation model with latent and observed variables.
      - 2.2_HOMA_climate-filtering-occupancy: Hierarchical single-species, single-season occupancy model of local, microclimate factors, including: average rock diameter and vegetation scores, based on ...
      - 2.3_HOMA_climate-drivers-SEM-occupancy-simplified:
    
## Variables

- **slope**: Calculated from a digital elevation model from ArcticDEM (U.S. Geological Survey 2020) with 8 neighboring cells (Horn 1981).
- **max.sum.temp**: Maximum summer temperature (&deg;C), as calculated from MOD11A1 Version 6.1 daily Land Surface Temperature between May to August.
- **days.above.21**: Total summer days above 21&deg;C, as calculated from MOD11A1 Version 6.1 daily Land Surface Temperature between May to August.
- **min.wint**: Minimum winter temperature (&deg;C), as calculated from MOD11A1 Version 6.1 daily Land Surface Temperature between October of the previous year to March.
- **wint.days.below.neg5**: Total winter days below -5&deg;C, as calculated from MOD11A1 Version 6.1 daily Land Surface Temperature between October of the previous year to March.
- **wint.days.above.0**: Total winter days above 0&deg;C, as calculated from MOD11A1 Version 6.1 daily Land Surface Temperature between October of the previous year to March.
- **ros**: Total number of days with rain-on-snow events, as calculated from SnowModel where the 10 year average daily rainfall is ≥ 3mm on snow depths ≥ 1.5 cm.
- **snow.dens**: Average snow density (kg/m<sup>3</sup>), as calculated from the 10 year average SnowModel data between September of the previous year to April.
- **snow.depth**: Maximum snow depth (m), as calculated from the 10 year average SnowModel data between September of the previous year to April.
- **perm**: Continuous probability (%) of near-surface permafrost, based on models by Pastick 4326.
- **evi**: Average summer greenness from the previous summer, as measured via the enhanced vegetation index and calculated from MOD13Q1 v061 16-day EVI between June – August.

## Data

  1. Variables --
  2. 
