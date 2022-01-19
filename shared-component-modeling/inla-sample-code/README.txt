SampleData.csv - simulated data for 156 weeks and 4 strata with the following columns:
 - wk: week
 - strat: stratum id (1-4)
 - Nj: stratum population
 - yE: pathogen E disease counts in week wk and stratum j  (not observed)
 - yC: pathogen C disease counts in week wk and stratum j  (not observed)
 - yO: pathogen O disease counts in week wk and stratum j  (not observed)
 - yES: pathogen E severe disease counts in week wk and stratum j  (not observed)
 - yCS: pathogen C severe disease counts in week wk and stratum j  (not observed)
 - yOS: pathogen O severe disease counts in week wk and stratum j  (not observed)
 - yS: total number of severe cases in week wk and stratum j (observed)
 - yM: total number of mild cases in week wk and stratum j (observed)
 - kS: total number of severe cases subsampled for virology in week wk and stratum j (observed)
 - kM: total number of mild cases subsampled for virology in week wk and stratum j (observed)
 - zES: number of pathogen E severe cases in subsample in week wk and stratum j  (observed)
 - zCS: number of pathogen C severe cases in subsample in week wk and stratum j  (observed)
 - zOS: number of pathogen O severe cases in subsample in week wk and stratum j  (observed)
 - zEM: number of pathogen E mild cases in subsample in week wk and stratum j  (observed)
 - zCM: number of pathogen C mild cases in subsample in week wk and stratum j  (observed)
 - zOM: number of pathogen O mild cases in subsample in week wk and stratum j  (observed)
 - temp: average weekly temperature in week wk (observed)
 - logtheta.E: true value of log theta_E used to simulate this data
 - logtheta.C: true value of log theta_C used to simulate this data
 - logtheta.O: true value of log theta_O used to simulate this data

SampleCode.R contains the R code to estimate unobserved pathogen-specific disease counts and fit the joint-pathogen model using INLA with a single covariate, as described in "Time Series Modeling of Pathogen-Specific Disease Probabilities with Subsampled Data" by Fisher, L. et al. 