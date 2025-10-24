# Latent variable estimation with composite Hilbert space Gaussian processes
This is a repository for the codes used in the paper: Latent variable estimation with composite Hilbert space Gaussian processes. The codes are arranged as follows. Check 'indcompgpfns.R' to call necessary functions that can be used in all of the R scripts (for e.g: to customise any simulation conditions, fit models, prepare summary plots).

## Models
Our developed models: partial composite and partial derivative Hilbert space Gaussian processes (pcHSGPs and pdHSGPs) for Squared Exponential (SE) covariance functions are given as 'pcHSGP.stan' and 'pdHSGP.stan' to be found in 'Stan models' folder. Their exact GP counterparts, joint derivative GPs among others that are used as comparitive models in our experiments are also provided in the same folder.

## Simulation Study
See 'Simulation scripts' folder for the various simulation study scenarios. The script provides an output dataframe with all the results for post-processing. Our simulation experiment results are provided in the 'Simulation results' folders for the respective models that can be post-analysed using the appropriate R scripts.

## Model convergence
Check model convergence for all our simulation experiments using the 'ModelConvergence.R' script.

## Model Calibration
We use SBC (https://projecteuclid.org/journals/bayesian-analysis/volume-20/issue-2/Simulation-Based-Calibration-Checking-for-Bayesian-Computation--The-Choice/10.1214/23-BA1404.full) to check for our model calibration. See 'SBC_log_gamma.R' script for model calibration using the simulation results.

## Latent variable estimation
Check latent variable estimation accuracy for various simulation studies provided (or create your own) and compare between the different stan models available using the 'LatentxRecovery.R' script.

## GP hyperparameter estimates
Check the GP hyperparameter estimates and recovery of true simulated values using 'HyperparameterRecovery.R'.

## Case study
We apply our models on a single-cell biology data involving the maturation process of erythroid cells from progenitor cells (https://www.nature.com/articles/s41586-019-0933-9) to esimate latent cellular ordering. We compare the results between pcHSGP and pdHSGP as well as experimental time. See 'CaseStudy_pcHSGP.R' and 'CaseStudy_pdHSGP.R' to re-analyse the data provided in 'casestudydata' folder. The result figures can be generated using the script 'CaseStudyFigures.R'. For more details, check the 'Real-World Case Study' section of our paper.
