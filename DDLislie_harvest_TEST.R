source('DDLislie_harvest.R')

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.2,0.15,0.15))
non0ferc = c(2,3)
nage = 3
K0 = 800
period = 10
Survival_inhomo = matrix( rep((Survival),period),nrow = nage,ncol = period )
Ferc_inhomo = matrix( rep((Ferc),period),nrow = nage,ncol = period )

### Test of helper functions
Lislie = getLislie(Survival,Ferc,T)
H = GetHarvest(Harvpar,nage)
H_prime = GetHarvest(0.05,nage) # age imspecific harvest rate
D = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),K0,T)
I = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),K0,T)
ProjectHarvest_helper(c(1,1,1), Lislie, H, global=T, E0=NULL, K0=K0, null = F)
data_homo = ProjectHarvest_homo(Survival=Survival,Harvpar=Harvpar,Ferc=Ferc,E0=NULL,K0=K0,global = T,null=F,bl=c(10,10,10),period = period,nage)
data_inhomo = ProjectHarvest_inhomo(Survival=Survival_inhomo,Harvpar=Harvpar,Ferc=Ferc_inhomo,E0=NULL,K0=K0,global = T,null=F,bl=as.matrix(c(1,1,1)),period = period,nage)
living_homo = getLivingIdividuals(H,data_homo)
log.lhood(log(data_homo),log(data_inhomo+matrix(rnorm(24,0,0.05),3,11)),1)

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

### SIMULATION START HERE ###

##without assume DD

SIMULATION_RES = HDDLislie.sampler(n.iter = 5000, burn.in = 300, mean.f = mean.f
                                   ,al.f = 1, be.f = .5, al.s = 1, be.s = .5
                                   , al.K0 = 1, be.K0 = .01, al.n = 1
                                   , be.n = .01, al.H = 1, be.H = .01
                                  , mean.s = mean.s, mean.b=mean.b,mean.K0 = mean.K0
                                  , mean.H = mean.H, Harv.data = (Harv.data+1e-4 * (Harv.data==0))
                                  , prop.vars = prop.vars)
mean.surv = apply(SIMULATION_RES$surv.prop.mcmc,2,mean)
mean.surv.2 = mean.surv[(1:21)%%3==2]
plot(mean.surv.2)

mean.harv = apply(SIMULATION_RES$lx.mcmc,2,mean)
mean.harv.2 = mean.harv[(1:21)%%3==2]
plot(mean.harv.2,mean.surv.2)

data.age.2 = data.frame(mean.harv.2,mean.surv.2)
DDsurv = lm(mean.surv.2~.,data.age.2)

mean.f.est = apply(SIMULATION_RES$fert.rate.mcmc,2,mean)
mean.f.2 = mean.f[(1:21)%%3==2]
plot(mean.f.2)

require(ggplot2)

ggplot(data.age.2,aes(x=mean.harv.2,y=mean.surv.2) )+
  geom_point()+
  geom_line()+
  stat_smooth(method = "lm")

data.harvest = data.frame(t=1:11,harv = apply(Harv.data,2,sum))
ggplot(data.harvest,aes(x=t,y=harv))+
  geom_point()+
  geom_line()
  #stat_smooth(method = "lm")

## asume DD
SIMULATION_RESDD = HDDLislie.sampler(n.iter = 5000, burn.in = 300, mean.f = mean.f
                                   ,al.f = 1, be.f = .1, al.s = 1, be.s = .1
                                   , al.K0 = 1, be.K0 = .005, al.n = 1
                                   , be.n = .1, al.H = 1, be.H = .1
                                   , mean.s = mean.s, mean.b=mean.b,mean.K0 = mean.K0
                                   , mean.H = mean.H, Harv.data = (Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, null = F, estK0 = T)
mean.survDD = apply(SIMULATION_RESDD$surv.prop.mcmc,2,mean)
mean.survDD.2 = mean.survDD[(1:21)%%3==2]
plot(mean.surv.2)

mean.harvDD = apply(SIMULATION_RESDD$lx.mcmc,2,mean)
mean.harvDD.2 = mean.harvDD[(1:21)%%3==2]
plot(mean.harvDD.2,mean.survDD.2)

mean.K0.est = apply(SIMULATION_RESDD$K0.mcmc,2,mean)
mean.fDD = apply(SIMULATION_RESDD$fert.rate.mcmc,2,mean)
mean.fDD.2 = mean.fDD[(1:21)%%3==2]
mean.fDD.3 = mean.fDD[(1:21)%%3==0]

data.age.2 = data.frame(mean.harvDD.2,mean.survDD.2)

require(ggplot2)

ggplot(data.age.2,aes(x=mean.harvDD.2,y=mean.survDD.2) )+
  geom_point()+
  geom_line()+
  stat_smooth(method = "lm")
