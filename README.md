# zibbcp
Title: A Bayesian zero-inflated beta-binomial model for longitudinal data with group-specific changepoints

Journal: Statistics in Medicine

Author: Chun-Che Wen, Nathaniel Baker, Rajib Paul, Elizabeth Hill, Kelly Hunt, Hong Li, Kevin Gray, and Brian Neelon

This repository includes the simulation for the manuscript submitted in Statistics in Medicine. The R script (ZIBBCP-Simulation-N500.R) is the simulation for the zero-inflated beta-binomial (ZIBB) mixed model with the group-specific changepoints. The ZIBB mixed model can accommodate bounded and discrete data with overdispersion, repeated measures, and zero inflation. The group-specific changepoint can capture the dynamic change in treatment efficacy. Two group-specific changepoints are included in both the binary and BB components of the zero-inflated model.  

This simulation includes the data generation, MCMC algorithm (Gibbs+MH steps), figures, traceplots.

## Files
ZIBBCP-Simulation-N500.R   =  R code simulation with sample size N=500

ZIBBCP-Simulation-N500.Rda =  MCMC samples for Alpha/Beta/Rho/Kappa1/Kappa2/Random-Effect Variance (generated from ZIBBCP-Simulation-N500.R script)

ZIBBCP-Simulation.pdf =  R code simulation markdown  

Figures folder = All figures from the manuscript and supplement 

## Parameters
Alpha  = Fixed effects in binary component 

Beta   = Fixed effects in BB component

Rho    = Correlation parameter in BB distribution (0<rho<1)

Kappa1 = Placebo changepoint

Kappa2 = Treatment changepoint

Sigmab = 6 x 6 matrix of Random effect variance (random intercept/slope/slope after cp)
