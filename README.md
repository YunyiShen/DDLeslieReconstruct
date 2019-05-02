# Leslie_Reconstruct_ChicagoDeer
## Introduction
This is a working repo for Yunyi's project, reconstruction of Chicago [white tailed deer](https://en.wikipedia.org/wiki/White-tailed_deer) (*Odocoileus virginianus*) population using culling data and Bayesian framework.
Basic ecology is Leslie growth with and without density dependency.

The goal for this project is reconstruction the population dynamics in the past decade, and try to figure out what kind of density dependency this population has. It is highly related to management of white tailed deer population, since if survival rate itself is negative related to density, control a lower density population can be much less efficient than a higher density one. Careful sensitivity analysis of the growth matrix is needed to figure out how to management the population when there exist density dependency.

## Method and the Model
Basic model follows a general Bayesian framework of population reconstruction. Follow the 4 level model of Mark C. Wheldon et al. (2013). Log survival rate, logit birth rate and logit culling proportion were assumed to follow Gaussian distribution, while the standard deviations follow Gamma distribution. Baseline population also follows Gaussian distribution. 

Every sample period, culling data was project according to the (non-)density dependent Leslie Matrix Model. Likelihood was calculated according to a log Gaussian noise distribution of population and all population parameters (e.g. survival rate) and their standard deviations follows gamma distribution. Priors are according to our best knowledge of the rates and population counts. 

## Simulation Study
### Simulation Population Setting
Coming soon
### Detection of Density Dependency
Coming soon
