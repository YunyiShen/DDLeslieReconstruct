

Introduction
============

White tailed deer (*Odocoileus virginianus*) is one of the most
important game species in north America. However, over abundance of this
species also cause problems, including over grazing one landscape,
conflict with human and spread diseases like Crown Waste Disease. In
suburban area, rapid growth of deer population increased the probability
of human-deer conflict. In Chicago suburban area, increasing deer
density promote concern of deer-vehicle collision (Etter et al. 2000, Jones and Witham 1995 ). In
1992, 2 human were killed and 145 injured due to deer-vehicle collision
(Etter et al. 2002, Etter 2001). Over growth of deer population also raise public health
concern of zoonoses in this area (Miller 2001, Etter 2001).

Researchers and Officials from Illinois Department of Natural Resources
started to try intensive culling as population control method in Cook
and DuPage Counties, Illinois (Complex 1) from 1992 to current.
Reconstruction of population dynamic is needed to evaluate this method
and its outcome.

Various method for reconstruction populations using multiple data
sources were proposed recently (Iijima et al. 2013, Wheldon et al. 2013). Bayesian method is one of
the most promising approach. Stochastic matrix projection model is one
of the most famous projection models for population with age structure
(Leslie 1945). Bayesian framework of reconstruction allow us to fuse multiple
data sources with different quality. Here we proposed and tested a
Bayesian framework reconstruction modified from Wheldon et al. 2013 in order to
understand the population’s dynamic under culling.

Method
======

Study Area and Environment
--------------------------

The study area encompassed 334,934 ha in the western suburbs of the
greater Chicago metropolitan area, in Cook and DuPage Counties,
Illinois. Land cover was dominated by urban/developed land (57.5%) and
associated urban grassland (14%). Other cover types included
forested/woodland (14.2%), cropland (4.9%), rural grassland (4.1%),
wetland (3.3%), open waters (1.8%), and barren/exposed land (0.2%;
Illinois Department of Natural Resources 1996).

Population of study occupied complex 1 included Argonne National
Laboratory (ANL), Waterfall Glen (WFG) and Woodridge (WDG) forest
preserves located in southeast DuPage County, and Black Partridge forest
preserve (directly adjacent to Woodridge) located in southwest Cook
County, Illinois. Habitat area of this complex remained constant through
culling years.

Intensive Culling
-----------------

Prior to 1992, no deer were legally harvested in Complex 1. Culling in
Complex 1 started in WFG in 1992. In 1993 and 1994, we culled deer from
WFG and WDG, and beginning in 1995 we culled deer from ANL and all other
forest preserve in Complex 1. Deers were culled during October–April by
sharpshooting, and capture and euthanasia as described by
(DeNicola and Williams 2008 Etter 2001). Antlerless deer were prioritized, but deer were
culled on a first opportunity basis.

Data Collection and Initial Estimation
--------------------------------------

Culled deers were aged into 8 age classes for female, 3 for male and
recorded from 1992 to 2006. Initial estimation and their uncertainty
which also serves as prior mean of vital rates and used to determine
hyperparameters came from previous study in the same area (Etter 2001).
Initial estimation of harvest rate was the ratio between previous
reconstruction study in the same population and total culling count,
while we assume harvest rate of fawns was half of harvest rate of adult
since the protocol of culling was strongly biased toward non-fawns.

Bayesian Reconstruction of Population Dynamics under Culling
------------------------------------------------------------

### Notations and Parameterization

Parameters of interest are time and age specific fecundity and survival,
as well as harvest rate and population size of female deers in study
area. We later on use <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{C}_{a,t}" title="cat" />   for culling counts at age a and time
t,  <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{s}_{a,t}" title="sat" /> 
for survival, <img src="https://latex.codecogs.com/svg.latex?\Large&space;f_{a,t}" title="fat" />
for fecundity, <img src="https://latex.codecogs.com/svg.latex?\Large&space;H_{a,t}" title="Hat" />
harvest rate, <img src="https://latex.codecogs.com/svg.latex?\Large&space;X_{a,t}" title="Xat" />
for latent living population size after culling,
<img src="https://latex.codecogs.com/svg.latex?\Large&space;Ae_{t}" title="Aet" /> 
for the post harvest aerial count estimation, 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;A_{t}" title="At" />
for the
aerial count detection, $SRB$ for sex ratio at birth and M for the
matrix projection model. Bold form of value is the the corresponding age
vector (e.g.<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{C}_{t}=(C_{1,t},C_{2,t}...)" title="Cteg" />
 ). Underline of
certain parameters means the best estimation and data we have currently
that will be used in the model. 

### Projection Model for Culling Dynamics

We assume a time inhomogeneous stochastic proportional harvest. Further,
harvest rate of fawns (age 0.5) is assumed to be different from yearling
and adults (age >0.5). We used a diagonal harvest matrix
<img src="https://latex.codecogs.com/svg.latex?\Large&space;H_{a,t}" title="Hat" />
to model the harvest in the projection to seperate
harvest and other mortality. Growth of living individuals are projected
using standard stochastic Leslie matrix model (Leslie 1945). Leslie
matrix contains <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{s}_{a,t}" title="sat" /> 
 and <img src="https://latex.codecogs.com/svg.latex?\Large&space;f_{a,t}" title="fat" />
 was noted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{H}_{a,t}" title="Hat" />.
Since harvest happened after reproduction, we left multiply the harvest
matrix. Vital rates’ distribution as described in the prior part of the
Bayesian framework.

Formally the projection model of harvest count vector<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{C}_{t}" title="\Large \mathbf{C}_{t}" /> 
is given by: 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{C}_{t+1}=\mathbf{H}_{t+1}\mathbf{L}_{t+1}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}" title="\Large \mathbf{C}_{t+1}=\mathbf{H}_{t+1}\mathbf{L}_{t+1}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t})" /> 

We also have aerial count data as estimation of post-harvest population
with imperfect detection. Model for aerial count is given below.

<img src="https://latex.codecogs.com/svg.latex?\Large&space;Ae_{t}=sum(A_{t}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t})" title="\Large Ae_{t}=sum(A_{t}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t})" />


In which I is identity matrix. Note that
<img src="https://latex.codecogs.com/svg.latex?\Large&space;(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}" title="\Large (\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}" /> 
 solves the living
individual after culling that undergone reproduction at time t+1. Also
note that baseline year should be trait differently since there is no
L0 needed. A graphical illustration of the dynamics is given in
Fig below.

![](https://github.com/YunyiShen/DDLeslieReconstruct/blob/uniform-aK0-prior/_figs_/LHD.jpg "Projection Model used in Bayesian Framework for Culling Dynamics")

### Bayesian Reconstruction

We generally followed Bayesian reconstruction framework of (Wheldon et al. 2013),
except we had a new set of parameters regarding the harvest rate and
aerial count. Further, likelihood used Poisson distribution rather than
log normal since we have 0 harvests. Reconstruction is equivalent of
estimating vital rates s, f, H and population population counts
X. We used the same 4 level setting to count for uncertainty of
initial estimation.

![](https://github.com/YunyiShen/DDLeslieReconstruct/blob/uniform-aK0-prior/_figs_/4level.jpg "Prior and Likelihood Part of the Bayesian Reconstruction Framework")

Relationship between data and parameters were shown in
Fig below. Note that alpha and beta are hyperparameters which encode the prior
knowledge we have for the uncertainty of certain vital rate or culling
count. Initial estimations of survival (logit transformed), fecundity
(log transformed), sex ratio at birth (SRB, logit transformed), harvest
rate (logit transformed), as well as aerial count probability (logit
transformed) were used as corresponding transformed prior mean of the
normal distribution whose sigma was invGamma distributed with
parameter alpha and beta. Culling
counts as well as aerial counts served in the likelihood part of the
model, prediction of the projection model served as expected value of
the Poisson distribution to evaluate the likelihood.

![](https://github.com/YunyiShen/DDLeslieReconstruct/blob/uniform-aK0-prior/_figs_/bayesian.jpg "Relationship between Various Data and Parameters in the Bayesian Reconstruction Model for Culling Dynamics")

### Determining the Hyperparameters

Determination of hyperparameters alpha and beta for vital rates
except for harvest rate were based on previous study’s error estimation,
use the same method in (Wheldon et al. 2013), but more conservative. Harvest rate’s
hyperparameter were set to be enough conservative that has .95 quantile
greater than 2. Detail hyperparameter setting is shown in


  parameters   |Survival  |Fecundity  |SRB  |Harvest  |Aerial detection
  -------------|----------|-----------|-----|---------|-----------------
  alpha        |   1      |    1      | 1   |  1      |       1
  beta         | .05      |  .01      |.05  | .05     |      .05


### Estimation

We draw samples from posterior distribution using Markov chain Monte
Carlo (MCMC) method and diagnose using R package `coda`
. The algorithm generally followed Wheldon et al. 2013. We
updated variance (sigmas) from conjugated full conditional
distributions as proposed density in Metropolis-Hastings algorithm.
Other vital rates including culling rate were updated using
Metropolis-Hastings steps and a symmetric normal proposal. Algorithm was
tuned by hand to achieve a reasonable acceptance rate. Follow Wheldon et
al., “iteration” was defined as one complete sweep through all age- and
time-specific parameters and variance parameters.

### Living Individuals

After we reconstruct the culling dynamics <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{C}_{t}" title="\Large \mathbf{C}_{t}" /> and harvest
rate <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{H}_{t}" title="\Large \mathbf{H}_{t}" />, we solve the living individual after culling
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{X}_{t}" title="\Large \mathbf{X}_{t}" /> using eqn:

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{X}_{t}=(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}" title="\Large \mathbf{X}_{t}=(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}" />


The model is implemented in R 3.6.0 modified from package
`popReconstruct` (Weldon et al. 2013) Source code is available on
[`GitHub`](https://github.com/YunyiShen/DDLeslieReconstruct) under MIT
license.

### Model Selection
In this frame work, we use [DIC](https://en.wikipedia.org/wiki/Deviance_information_criterion) as model selection criterion (Gelman 2002). Lower DIC means higher support of the model by data. 

Results
=======

Data Collection
---------------

Total 15 years of culling count was used in this study (1992-2006, 1992
as baseline). In 15 years total 3827 individuals were culled. Fecundity
and survival estimation was homogeneous for all age $>1.5$, however we
keep this. Average fecundity for adults is 1.86, yearling 1.53 and fawn
0.178. Mean survival rates are 0.85, 0.82 and 0.83 for fawn, yearling
and adults respectively (Etter 2008). Mean harvest rate for the population
is 0.5 (Etter et al. 2019 unpublished data).

  ***Culling counts***          |mean    |   Standard Error
  ------------------------------|--------|------------------
  Absolute Difference           |  3.38  |     1.65
  Posterior Standard Deviation  |  3.25  |     0.151
  ***Aerial counts***           |        | 
  Absolute Difference           |  2.25  |     0.33
  Posterior Standard Deviation  |  14.20 |     0.240


Model Selection
---------------
Four models were tested to chose the best supported:  
1) Four Harvest full model: Survival and Fecundity rates are time sex and age dependent, with 4 harvest rates: female fawn, femal adult male fawn and male adult which are time dependent.  
2) Three Harvest full model: Survival and Fecundity rates are time sex and age dependent, with 3 harvest rates: fawn, femal adult and male adult which are time dependent.  
3) Solely Density Dependent model: Survival and Fecundity rates are sex age and density but not time dependent, with 4 harvest rates, which are time dependent.  
4) Three age classes: All similar to 1 but fecundity and survival has only 3 age classes for both female and male.  

  Model | effective n parameters |  DIC |Best Model
  ------|------------------------|------|----------
  1     |187.11	                 |649.38|
  2     |196.3                   |631.48|*
  3     |108.81                  |835.74|
  4     |148.36                  |732.35|




Living Individuals After Culling
--------------------------------

After reconstruct the culling dynamic, we solved the total living
individuals’ posterior distribution as described by eqn.1.
Results were shown in Fig below. 

![](https://github.com/YunyiShen/DDLeslieReconstruct/blob/uniform-aK0-prior/_figs_/popsize.jpg "Reconstructed Population Size After Culling")

Density Dependency of Vital Rates
---------------------------------

Density dependency is one of the key features controlling population
growth (ueno 2010). To address whether there exist density dependency we
did linear regression between vital rates and reconstruct living
population size before reproduction (i.e. living individual counts after
culling of the previous year). Results were summarized in
Table below.

  Age              |     0.5    |     1.5    |     2.5    |     3.5    |     4.5    |     5.5    |      6.5     |      7.5
  -----------------| -----------| -----------| -----------| -----------| -----------| -----------| -------------| -------------
  ***Fecundity***  |     \*     |     \*     |     \*     |     \*     |     \*     |     \*     |      \*      | 
  p-value          |   0.00203  |   2.79e-5  |   2.27e-5  |   2.82e-4  |   7.71e-4  |   1.84e-3  |    0.0147    |    0.0343
  adj R^2          |  0.525     |  0.762     |  0.770     |  0.654     |  0.593     |  0.532     |   0.353      |   0.265
  beta             |-6.32e-05   | -1.01e-03  | -6.37e-04  | -4.96e-04  | -4.09e-04  | -3.28e-04  |  -2.70e-04   |  -2.49e-04
  SE beta          |1.61e-05    | 1.55e-04   | 9.540e-05  | 9.81e-05   | 9.16e-05   | 8.26e-05   |  9.49e-05    |  1.04e-04
  ***Survival***   |            |            |            |            |            |     \*     |      \*      | 
  p-value          |    0.365   |    0.150   |    0.886   |    0.130   |    0.514   |   0.00705  |    2.63e-4   |    0.0546
  adj R^2          |    0       | 0.0949     |    0       |  0.112     |    0       |  0.423     |0.657629475   |0.213787859
  beta             |-8.10e-05   | 2.64e-04   | -9.31e-06  | -1.69e-04  | 5.60e-05   | -4.42e-04  |  -2.12e-04   |  -9.33e-05
  SE beta          |8.61e-05    | 1.72e-04   | 6.39e-05   | 1.04e-04   | 8.32e-05   | 1.36e-04   |  4.16e-05    |  4.38e-05
  \*:  p:0.01      |            |            |            |            |            |            |              | 

For males:

  Age             |    Fawn    |  Yearling  |   Adult
  ----------------| -----------| -----------| ----------
  ***Survival***  |     \*     |            | 
  p-value         |  0.000441  |    0.402   |   0.0805
  adj R^2         |    0.628   |      0     |   0.169
  beta            |  -6.53e-04 |  7.731e-05 |  1.55e-04
  SE beta         |  1.36e-04  |  8.89e-05  |  8.17e-05
  \*:p:0.01       |            |            | 

  

Discussion
==========

By conducting intensive culling in Chicago’s suburban area, we
successfully controlled the overabundant deer population in this area.
Density dependency may reduce our ability to control population when
size is small. Our reconstruction shows that survival has no density
dependency probability since the source of mortality stays similar
before and after intensive culling. Fecundity shows a weak negative
density dependency for older individuals. This suggest in lower
population density, further culling of older individual is needed to
keep the same amount of control outcomes.

We assumed the population is closed to female which is generally
appropriate since deer is male dispersion. But if we are willing to
include sex ratio in this reconstruction we may need to consider
migration rates from the population for males. Sex ratio at birth is
also not considered here for a female reconstruction, however, it is
interesting to evaluate the sex ratio change under intensive culling.

The statistical model we modified can also be generalized to reconstruct
proportional harvesting dynamics using uncertain estimation of vital
rates and age-at-harvest data. It can also include relationships between
vital rates with environmental variables.

Supplementary Figures
=====================

All the supplementary figures and result summaries, as well as the
source code of this report are available
[here](https://github.com/YunyiShen/DDLeslieReconstruct/tree/uniform-aK0-prior/Main_analysis/figs/Poisson_harvest_longer_chain/4%20harvest).
