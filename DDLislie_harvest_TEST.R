source('DDLislie_harvest.R')

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.1,0.05,0.05))
non0ferc = c(2,3)
nage = 3
K0 = 800
period = 7
Survival_inhomo = matrix( rep((Survival),period),nrow = nage,ncol = period )
Ferc_inhomo = matrix( rep((Ferc),period),nrow = nage,ncol = period )

### Test of helper functions
Lislie = getLislie(Survival,Ferc,T)
H = GetHarvest(Harvpar,nage)
H_prime = GetHarvest(0.05,nage) # age imspecific harvest rate
D = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),K0,T)
I = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),K0,T)
ProjectHarvest_helper(c(1,1,1), Lislie, H, global=T, E0=NULL, K0=K0, null = F)
data_homo = ProjectHarvest_homo(Survival=Survival,Harvpar=Harvpar,Ferc=Ferc,E0=NULL,K0=K0,global = T,null=T,bl=c(1,1,1),period = period,nage)
data_inhomo = ProjectHarvest_inhomo(Survival=Survival_inhomo,Harvpar=Harvpar,Ferc=Ferc_inhomo,E0=NULL,K0=K0,global = T,null=T,bl=as.matrix(c(1,1,1)),period = period,nage)
living_homo = getLivingIdividuals(H,data_homo)
log.lhood(log(data_homo),log(data_inhomo+matrix(rnorm(24,0,0.05),3,8)),1)

# skip the simplest functions

### Test of sampler function

mean.f = Ferc_inhomo
mean.s = Survival_inhomo
mean.b = as.matrix( c(1,1,1))
mean.K0 = as.matrix(K0)
mean.H = Harvpar
Harv.data = floor(data_homo+0.5) # this harv have density dependency
prop.vars = list(fert.rate = matrix(1,nrow = length(non0ferc),ncol = period),
                 surv.prop = matrix(1,nrow = nage, ncol = period),
                 H = matrix(1,nrow = nage,ncol = 1),
                 K0=1,
                 baseline.pop.count = matrix(1,nrow = nage,ncol = 1)) # should be a list, fert.rate for all age classes can give birth to new 
proj.periods = (ncol(Harv.data)-1)
ptm = proc.time()

### SIMULATION START HERE

SIMULATION_RES = HDDLislie.sampler(n.iter = 1000, burn.in = 10, mean.f = mean.f
                                   ,al.f = 1, be.f = 0.0001, al.s = 1, be.s = 0.0001
                                   , al.K0 = 1, be.K0 = 0.0001, al.n = 0.0001
                                   , be.n = 0.0001, al.H = 1, be.H = 0.0001
                                  , mean.s = mean.s, mean.b=mean.b,mean.K0 = mean.K0
                                  , mean.H = mean.H, Harv.data = (Harv.data+1e-4 * (Harv.data==0))
                                  , prop.vars = prop.vars)
