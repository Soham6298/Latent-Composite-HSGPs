#library
library(rstan)
library(posterior)
library(ggplot2)
library(bayesplot)

# Source fns 
source('indcompgpfns.R')
comphsgpmodel <- stan_model('indpcomphsgp_maternclass.stan')
compgpmodel <- stan_model('indpcompexactGP_maternclass.stan')

# Read data
dat_unsp <- read.csv('case study data/X_unspliced.csv', header = FALSE)
dat_sp <- read.csv('case study data/X_spliced.csv', header = FALSE)
dat_t <- read.csv('case study data/obs.csv')

