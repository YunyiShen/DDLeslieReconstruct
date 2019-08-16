

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
area. We later on use $C_{a,t}$ for culling counts at age $a$ and time
$t$, $s_{a,t}$ for survival, $f_{a,t}$ for fecundity, $H_{a,t}$ harvest
rate, $X_{a,t}$ for latent living population size after culling,
$Ae_{t}$ for the post harvest aerial count estimation, $A_{t}$ for the
aerial count detection, $SRB$ for sex ratio at birth and $M$ for the
matrix projection model. Bold form of value is the the corresponding age
vector (e.g. $\mathbf{C}_{t}=(C_{1,t},C_{2,t}...)$). Underline of
certain parameters means the best estimation and data we have currently
that will be used in the model. (see Fig.\[Fig.Bayes\])

### Projection Model for Culling Dynamics

We assume a time inhomogeneous stochastic proportional harvest. Further,
harvest rate of fawns (age 0.5) is assumed to be different from yearling
and adults (age $>0.5$). We used a diagonal harvest matrix
$\mathbf{H}_{t}$ to model the harvest in the projection to seperate
harvest and other mortality. Growth of living individuals are projected
using standard stochastic Leslie matrix model (Leslie 1945). Leslie
matrix contains $s_{a,t}$ and $f_{a,t}$ was noted by $\mathbf{L}_{t}$.
Since harvest happened after reproduction, we left multiply the harvest
matrix. Vital rates’ distribution as described in the prior part of the
Bayesian framework.

Formally the projection model of harvest count vector $\mathbf{C}_{t}$
is given by: $$\label{proj}
\mathbf{C}_{t+1}=\mathbf{H}_{t+1}\mathbf{L}_{t+1}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}$$

We also have aerial count data as estimation of post-harvest population
with imperfect detection \cite{}. Model for aerial count is given by
\[aeri\]. $$\label{aeri}
Ae_{t}=sum(A_{t}(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t})$$

In which $I$ is identity matrix. Note that
$(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}$ solves the living
individual after culling that undergone reproduction at time $t+1$. Also
note that baseline year should be trait differently since there is no
$L_{0}$ needed. A graphical illustration of the dynamics is given in
Fig.\[Fig.LHD\].

![Projection Model used in Bayesian Framework for Culling Dynamics](https://github.com/YunyiShen/DDLeslieReconstruct/tree/uniform-aK0-prior/_figs_/LHD.jpg)

### Bayesian Reconstruction

We generally followed Bayesian reconstruction framework of (Wheldon et al. 2013),
except we had a new set of parameters regarding the harvest rate and
aerial count. Further, likelihood used Poisson distribution rather than
log normal since we have 0 harvests. Projection model used in our study
was described by eqn.\[proj\]. Reconstruction is equivalent of
estimating vital rates $s$, $f$, $H$ and population population counts
$X$. We used the same 4 level setting to count for uncertainty of
initial estimation (Fig.\[Fig.Bayes4level\]).

![Prior and Likelihood Part of the Bayesian Reconstruction Framework](https://github.com/YunyiShen/DDLeslieReconstruct/tree/uniform-aK0-prior/_figs_/4level.jpg){width="0.8\linewidth"}

Relationship between data and parameters were shown in
Fig.\[Fig.Bayes\]. Note that $\alpha_{v}$ and $\beta_{v}$
$v\in\{s,f,H,A,SR\}$ are hyperparameters which encode the prior
knowledge we have for the uncertainty of certain vital rate or culling
count. Initial estimations of survival (logit transformed), fecundity
(log transformed), sex ratio at birth (SRB, logit transformed), harvest
rate (logit transformed), as well as aerial count probability (logit
transformed) were used as corresponding transformed prior mean of the
normal distribution whose $\sigma$ was invGamma distributed with
parameter $\alpha_{v}$ and $\beta_{v}$ $v\in\{s,f,H,A,SR\}$. Culling
counts as well as aerial counts served in the likelihood part of the
model, prediction of the projection model served as expected value of
the Poisson distribution to evaluate the likelihood.

![Relationship between Various Data and Parameters in the Bayesian Reconstruction Model for Culling Dynamics](https://github.com/YunyiShen/DDLeslieReconstruct/tree/uniform-aK0-prior/_figs_/bayesian.jpg){width="0.8\linewidth"}

### Determining the Hyperparameters

Determination of hyperparameters $\alpha$ and $\beta$ for vital rates
except for harvest rate were based on previous study’s error estimation,
use the same method in (Wheldon et al. 2013), but more conservative. Harvest rate’s
hyperparameter were set to be enough conservative that has .95 quantile
$>2$ ($\alpha=1$,$\beta=.1$). Detail hyperparameter setting is shown in
Table.\[tab:hyper\]

              Survival   Fecundity   SRB   Harvest   Aerial detection
  ---------- ---------- ----------- ----- --------- ------------------
   alpha$      1           1        1       1             1
   beta      .05         .01      .05     .05           .05

  : \[tab:hyper\]Hyperparameter Setting in This Study

### Estimation

We draw samples from posterior distribution using Markov chain Monte
Carlo (MCMC) method ([@RN14; @RN15]) and diagnose using R package `coda`
([@RN16]). The algorithm generally followed Wheldon et al. 2013. We
updated variance ($\sigma_{v}$s) from conjugated full conditional
distributions as proposed density in Metropolis-Hastings algorithm.
Other vital rates including culling rate were updated using
Metropolis-Hastings steps and a symmetric normal proposal. Algorithm was
tuned by hand to achieve a reasonable acceptance rate. Follow Wheldon et
al., “iteration” was defined as one complete sweep through all age- and
time-specific parameters and variance parameters.

### Living Individuals

After we reconstruct the culling dynamics $\mathbf{C}_{t}$ and harvest
rate $\mathbf{H}_{t}$, we solve the living individual after culling
$\mathbf{X}_{t}$ using eqn. \[eqn.living\]

$$\label{eqn.living}
\mathbf{X}_{t}=(\mathbf{H}_{t}^{-1}-\mathbf{I})\mathbf{C}_{t}$$

The model is implemented in R 3.6.0 modified from package
`popReconstruct` (Weldon et al. 2013) Source code is available on
[`GitHub`](https://github.com/YunyiShen/DDLeslieReconstruct) under MIT
license.

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

Model Checking
--------------

We use 2 distinct model checking indexes, Absolute Difference (AD)
defined as absolute value of difference between model predicted culling
counts and real counts which evaluates precision of the model. Posterior
Standard Deviation (PSD) defined as standard deviation of posterior
distribution of model predicted culling which evaluated uncertainty of
prediction.

                                  Mean    Standard Error
  ------------------------------ ------- ----------------
  ***Harvest***                          
  Absolute Difference             3.38         1.65
  Posterior Standard Deviation    3.25        0.151
  ***Aerial counts***                    
  Absolute Difference             2.25         0.33
  Posterior Standard Deviation    14.20       0.240

  : \[tab:check\]Model Checking Indexes for Reconstruction of Culling
  Data

Living Individuals After Culling
--------------------------------

After reconstruct the culling dynamic, we solved the total living
individuals’ posterior distribution as described by eqn.\[eqn.living\].
Results were shown in Fig.\[Fig.popsize\]. After 1996, posterior mean
showed the population was stabled around average 161 (sd=17) female
individuals which achieved our goal of controlling the population.
Harvest rate for fawns has mean of 0.26 and sd of 0.05, for other than
fawns has mean of 0.53 and sd of 0.13.

![Reconstructed Population Size After Culling](https://github.com/YunyiShen/DDLeslieReconstruct/tree/uniform-aK0-prior/_figs_/popsize.jpg)

Density Dependency of Vital Rates
---------------------------------

Density dependency is one of the key features controlling population
growth (ueno 2010). To address whether there exist density dependency we
did linear regression between vital rates and reconstruct living
population size before reproduction (i.e. living individual counts after
culling of the previous year). Results were summarized in
Table.\[tab:DDvital\_female\].

  Age                   0.5         1.5         2.5         3.5         4.5         5.5          6.5           7.5
  ----------------- ----------- ----------- ----------- ----------- ----------- ----------- ------------- -------------
  ***Fecundity***       \*          \*          \*          \*          \*          \*           \*       
  p-value             0.00203     2.79e-5     2.27e-5     2.82e-4     7.71e-4     1.84e-3      0.0147        0.0343
  adj R^2            0.525       0.762       0.770       0.654       0.593       0.532        0.353         0.265
  beta            -6.32e-05   -1.01e-03   -6.37e-04   -4.96e-04   -4.09e-04   -3.28e-04    -2.70e-04     -2.49e-04
  SE beta         1.61e-05    1.55e-04    9.540e-05   9.81e-05    9.16e-05    8.26e-05     9.49e-05      1.04e-04
  ***Survival***                                                                    \*           \*       
  p-value              0.365       0.150       0.886       0.130       0.514      0.00705      2.63e-4       0.0546
  adj R^2              0        0.0949         0         0.112         0         0.423     0.657629475   0.213787859
  \beta            -8.10e-05   2.64e-04    -9.31e-06   -1.69e-04   5.60e-05    -4.42e-04    -2.12e-04     -9.33e-05
  SE \beta         8.61e-05    1.72e-04    6.39e-05    1.04e-04    8.32e-05    1.36e-04     4.16e-05      4.38e-05
  \*:$p<0.01$                                                                                             

  : \[tab:DDvital\_female\]Linear Regressions of Female Vital Rates and
  Population Size after Culling

  Age                 Fawn      Yearling     Adult
  ---------------- ----------- ----------- ----------
  ***Survival***       \*                  
  p-value           0.000441      0.402      0.0805
  adj $R^2$           0.628         0        0.169
  $\beta$           -6.53e-04   7.731e-05   1.55e-04
  SE $\beta$        1.36e-04    8.89e-05    8.17e-05
  \*:$p<0.01$                              

  : \[tab:DDvital\_male\]Linear Regressions of Male Survival Rates and
  Population Size after Culling

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