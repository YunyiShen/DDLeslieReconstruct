source('./R/DDLeslie.R')
source('./R/misc.R')
nage = c(8,3) # nage is female first and then male, a vector with lenght usually 2
period = 14


mean.s = read.csv("./data/Survival_mean_Etter.csv",row.names = 1)
mean.f = read.csv("./data/Fecundity_mean.csv",row.names = 1)
mean.SRB = read.csv("./data/SRB_mean.csv",row.names = 1)
Harv.data = read.csv("./data/Culling.csv",row.names = 1)
Aeri.data = read.csv("./data/Aerial_count.csv",row.names = 1)
mean.b = Harv.data[,1]
mean.H = read.csv("./data/Harvest_rate_prior.csv",row.names = 1)
Harv_assump = read.csv("./data/Harv_assump.csv",header = F)
Harv_assump = as.matrix(Harv_assump) # this is the assumption matrix for specific harvest!
#Surv_assump = read.csv("./data/Surv_assump_age.csv",header = F)

Assumptions = list()

Assumptions$Fec = list(time = eyes(period),age = as.matrix(eyes(nage[1])))
Assumptions$Surv = list(time = eyes(period),age = as.matrix(eyes(sum(nage))))

# for direct inference:
#Assumptions$Fec = list(time = matrix(1,1,period),age = eyes((nage[1])))
#Assumptions$Surv = list(time = matrix(1,1,period),age = eyes(sum(nage)))
# end direct inference

Assumptions$SRB = list(time = eyes(period),age = eyes(1))
Assumptions$AerialDet  = list(time = eyes(period+1),age = eyes(1))
Assumptions$Harv = list(time = eyes(period+1),age = Harv_assump) # tons of assumptions on vital rates
Assumptions$aK0 = list(eyes(nage[1]),eyes(sum(nage)),matrix(1,1,1))
# full matrix for e.g. Harvest will be:
#  Assumptions$Harv$age %*% as.matrix(mean.H) %*% Assumptions$Harv$time
#  It is a good idea to try the command above to see how to use assumption matrices.

mean.A = matrix(0.7,1,period+1)
mean.aK0 = list(matrix(0,nage[1],1),matrix(0,sum(nage),1))
prop.vars = list(fert.rate = matrix(1,nrow = nage[1],ncol = period),
                 surv.prop = matrix(1,nrow = sum(nage), ncol = period),
                 SRB = matrix(.1,nage[1],period), # vital rates has period cols
                 A = matrix(1,1,period+1),
                 H = matrix(1,nrow = 4,ncol = period+1),
                 aK0=list(5e-8,5e-8,50),
                 baseline.pop.count = matrix(1,nrow = sum(nage),ncol = 1))

set.seed(42)

Chicago_RES = HDDLislie.sampler( n.iter = 50000, burn.in = 5000,thin.by = 25, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = 1e-2, al.s = 1, be.s = .05
                                   , al.SRB = 1, be.SRB = .05
                                   , al.aK0 = list(matrix(-.001,nage[1],1),matrix(-.001,sum(nage),1),0)
                                   , be.aK0 = list(matrix(.001,nage[1],1),matrix(.001,sum(nage),1),500)
                                   , al.H = 1, be.H = .05
                                   , al.A = 1, be.A = .05
                                   , mean.s = as.matrix(mean.s)
                                   , mean.b= as.matrix(mean.b)
                                   , mean.H = as.matrix(mean.H)
                                   , mean.SRB = as.matrix( mean.SRB)
                                   , mean.A = as.matrix( mean.A)
                                   , Assumptions = Assumptions
                                   , start.sigmasq.f = .05, start.sigmasq.s = .05, start.sigmasq.SRB = .05
                                   , start.sigmasq.H = .05
                                   , start.sigmasq.A = .05
                                   , Harv.data = as.matrix(Harv.data)
                                   , Aerial.data = as.matrix( Aeri.data)
                                   , prop.vars = prop.vars, estFer = T,nage = nage,estaK0 = F,null = T)


