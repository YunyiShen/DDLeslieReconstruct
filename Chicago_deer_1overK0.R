source('DDLislie_harvest_1overK0.r')
nage = 8
period = 14


prop.vars = list(fert.rate = matrix(.01,nrow = nage,ncol = period),
                 surv.prop = matrix(.01,nrow = nage, ncol = period),
                 H = matrix(.01,nrow = nage,ncol = 1),
                 aK0=.0001,
                 baseline.pop.count = matrix(.01,nrow = nage,ncol = 1))

mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,2]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,2]
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,-1]
mean.b = (584 * Harv.data[,1]/sum(Harv.data[,1]))/0.7

Chicago_RES = HDDLislie.sampler( n.iter = 50000, burn.in = 300, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = .005, al.s = 1, be.s = .01
                                   , al.aK0 = 1, be.aK0 = .1, al.n = 1
                                   , be.n = .01, al.H = 1, be.H = .01
                                   , mean.s = as.matrix(mean.s), mean.b= as.matrix(mean.b),mean.aK0 = matrix(0,1,1)
                                   , mean.H = matrix(0.7,nage,1), Harv.data = as.matrix(Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, estFer = T,nage = nage,homo = T)


hist(Chicago_RES$invK0.mcmc)

mean.harv = apply(Chicago_RES$lx.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
mean.total.harv = apply(mean.harv.matrix,2,sum)
plot(mean.total.harv)

mean.ferc = apply(Chicago_RES$fert.rate.mcmc,2,mean)
