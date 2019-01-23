source('DDLislie_harvest_1overK0.r')

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.2,0.15,0.15))
non0ferc = c(2,3)
nage = 3
aK0 = 0
period = 10
#Survival_inhomo = matrix( rep((Survival),period),nrow = nage,ncol = period )
#Ferc_inhomo = matrix( rep((Ferc),period),nrow = nage,ncol = period )

### Test of helper functions, DONE 0122
Lislie = getLislie(Survival,Ferc,T)
H = GetHarvest(Harvpar,nage)
H_prime = GetHarvest(0.05,nage) # age imspecific harvest rate
D = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),aK0,F)
I = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),aK0,T)
ProjectHarvest_helper(c(1,1,1), Lislie, Harvpar, global=T, E0=NULL, aK0=aK0, null = F)
data_homo = 
  ProjectHarvest_homo(Survival=Survival,Harvpar=Harvpar,Ferc=Ferc,E0=NULL,aK0=aK0,global = T,null=F,bl=c(10,10,10),period = period,nage)
data_homo = data_homo + matrix(rnorm(nrow(data_homo)*ncol(data_homo)),nrow = nrow(data_homo),ncol(data_homo))
data_homo = round(data_homo)
#data_inhomo = ProjectHarvest_inhomo(Survival=Survival_inhomo,Harvpar=Harvpar,Ferc=Ferc_inhomo,E0=NULL,K0=K0,global = T,null=F,bl=as.matrix(c(10,10,10)),period = period,nage)
living_homo = getLivingIdividuals(Harvpar,data_homo)
log.lhood(log(data_homo),log(data_homo+matrix(rnorm(24,0,0.05),3,11)),1)

# skip the simplest functions

### Test of sampler function

mean.f = Ferc
mean.s = Survival
mean.b = as.matrix( c(10,10,10)/Harvpar)
mean.aK0 = 0
mean.H = Harvpar
Harv.data = floor(data_homo+0.5) # this harv have density dependency
prop.vars = list(fert.rate = matrix(.01,nrow = length(non0ferc),ncol = period),
                 surv.prop = matrix(.01,nrow = nage, ncol = period),
                 H = matrix(.01,nrow = nage,ncol = 1),
                 aK0=.0001,
                 baseline.pop.count = matrix(.01,nrow = nage,ncol = 1)) # should be a list, fert.rate for all age classes can give birth to new 
proj.periods = (ncol(Harv.data)-1)
ptm = proc.time()

### SIMULATION START HERE ###

 ## statistical inference of 1/K0
SIMULATION_RESDD = HDDLislie.sampler(n.iter = 50000, burn.in = 300, mean.f = mean.f
                                   ,al.f = 1, be.f = .1, al.s = 1, be.s = .1
                                   , al.aK0 = 1, be.aK0 = .5, al.n = 1
                                   , be.n = .1, al.H = 1, be.H = .1
                                   , mean.s = mean.s, mean.b=mean.b,mean.aK0 = 0
                                   , mean.H = mean.H, Harv.data = (Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, null = F, estaK0 = T,homo = T)
mean.invK0 = apply(SIMULATION_RESDD$invK0.mcmc,2,mean)
hist(SIMULATION_RESDD$invK0.mcmc)

