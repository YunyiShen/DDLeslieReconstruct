source('DDLislie_harvest.R')

Survival = matrix(c(0.5,0.4,0.3),nrow = 3,ncol = 1)
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
mean.K0 = as.matrix( K0)
mean.H = Harvpar
Harv.data = floor(data_inhomo) # this harv have density dependency
prop.var = data.frame() # not sure the col names
proj.periods = (ncol(Harv.data)-1)
