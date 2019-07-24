source('DDLeslie.r')
nage = c(8,3) # nage is female first and then male, a vector with lenght usually 2
period = 14


mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,-(1)]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,-(1)]
mean.SRB = read.csv("./data/SRB_mean.csv")
mean.SRB = mean.SRB[,-1]
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,-(1)]
Aeri.data = read.csv("./data/Aerial_count.csv")
Aeri.data = Aeri.data[,-1]
mean.b = Harv.data[,1]
mean.H = read.csv("./data/Harvest_rate_prior.csv",row.names = 1)


mean.A = matrix(0.7,1,period+1)

prop.vars = list(fert.rate = matrix(1,nrow = nage[1],ncol = period),
                 surv.prop = matrix(1,nrow = sum(nage), ncol = period),
                 SRB = matrix(.1,nage[1],period), # vital rates has period cols
                 A = matrix(1,1,period+1),
                 H = matrix(1,nrow = 4,ncol = period+1),
                 aK0=1e-3,
                 baseline.pop.count = matrix(1,nrow = sum(nage),ncol = 1))

set.seed(42)

Chicago_RES = HDDLislie.sampler( n.iter = 50000, burn.in = 5000,thin.by = 25, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = 1e-2, al.s = 1, be.s = .05
                                   , al.SRB = 1, be.SRB = .05
                                   , al.aK0 = 1, be.aK0 = 1e-1
                                   , al.H = 1, be.H = .05
                                   , al.A = 1, be.A = .05
                                   , mean.s = as.matrix(mean.s)
                                   , mean.b= as.matrix(mean.b)
                                   , mean.aK0 = matrix(0,1,2)
                                   , mean.H = as.matrix(mean.H)
                                   , mean.SRB = as.matrix( mean.SRB)
                                   , mean.A = as.matrix( mean.A)
                                   , start.sigmasq.f = .05, start.sigmasq.s = .05, start.sigmasq.SRB = .05
                                   , start.sigmasq.aK0 = .05, start.sigmasq.H = .05
                                   , start.sigmasq.A = .05
                                   , Harv.data = as.matrix(Harv.data)
                                   , Aerial.data = as.matrix( Aeri.data)
                                   , prop.vars = prop.vars, estFer = T,nage = nage,homo = F,estaK0 = F)


