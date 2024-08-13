# zibbcp
Title: A Bayesian zero-inflated beta-binomial model for longitudinal data with group-specific changepoints

Journal: Statistics in Medicine

Author: Chun-Che Wen, Nathaniel Baker, Rajib Paul, Elizabeth Hill, Kelly Hunt, Hong Li, Kevin Gray, and Brian Neelon

This repository includes the simulation for the manuscript submitted in Statistics in Medicine. The script is the simulation for the zero-inflated beta-binomial (ZIBB) mixed model with the group-specific changepoints. The ZIBB mixed model can accommodate bounded and discrete data with overdispersion, repeated measures, and zero inflation. The group-specific changepoint can capture the dynamic change in treatment efficacy. 

This simulation includes the data generation, MCMC algorithm (Gibbs+MH steps), figures, traceplots.

## Files
zibbcp.R   =  R code simulation

zibbcp.Rda =  MCMC samples for Alpha/Beta/Rho/Kappa1/Kappa2/Random-Effect Variance

zibbcp-simulation.pdf =  R code simulation interpretation 



## Parameters
Alpha  = Fixed effects in binary component 

Beta   = Fixed effects in BB component

Rho    = Correlation parameter in BB distribution (0<rho<1)

Kappa1 = Placebo changepoint

Kappa2 = Treatment changepoint

Sigmab = Random effect variance (random intercept/slope/slope after cp)
