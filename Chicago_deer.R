source('DDLislie_harvest.R')
nage = 8
period = 14


prop.vars = list(fert.rate = matrix(.1,nrow = nage,ncol = period),
                 surv.prop = matrix(1,nrow = nage, ncol = period),
                 H = matrix(1,nrow = nage,ncol = period),
                 K0=.1,
                 baseline.pop.count = matrix(1,nrow = nage,ncol = 1))

mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,-1]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,-1]
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,-1]
mean.b = (584 * Harv.data[,1]/sum(Harv.data[,1]))/0.7

Chicago_RES = HDDLislie.sampler( n.iter = 15000, burn.in = 500, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = .01, al.s = 1, be.s = .01
                                   , al.K0 = 1, be.K0 = .01, al.n = 1
                                   , be.n = .0001, al.H = 1, be.H = .01
                                   , mean.s = as.matrix(mean.s), mean.b= as.matrix(mean.b),mean.K0 = matrix(700,1,1)
                                   , mean.H = matrix(0.7,nage,period), Harv.data = as.matrix(Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, estFer = T,nage = nage)

mean.surv = apply(Chicago_RES$surv.prop.mcmc,2,mean)
mean.ferc = apply(Chicago_RES$fert.rate.mcmc,2,mean)

mean.harv = apply(Chicago_RES$lx.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
mean.harv.matrix = cbind(mean.harv.matrix)
mean.total.harv = apply(mean.harv.matrix,2,sum)
mean.harv.por = apply(Chicago_RES$H.mcmc,2,mean)
mean.harv.por_matrix = matrix(mean.harv.por,nrow = nage,ncol = period)

mean.living = (mean.harv.matrix/mean.harv.por_matrix)*(1-mean.harv.por_matrix)
mean.living.total = apply(mean.living,2,sum)

mean.age.surv = list()
mean.age.ferc = list()
mean.age.harv = list()
DDsurv = list()
DDferc = list()
harv_surv = list()
harv_ferc = list()

for(i in 1:8){
  
  mean.age.surv[[i]] = mean.surv[(1:(nage*period))%%nage==(i%%nage)]
  mean.age.ferc[[i]] = mean.ferc[(1:(nage*period))%%nage==(i%%nage)]
  mean.age.harv[[i]] = mean.harv.por[(1:(nage*period))%%nage==(i%%nage)]
  data.temp = data.frame(mean.age.surv = mean.age.surv[[i]]
                         ,mean.age.ferc = mean.age.ferc[[i]]
                         ,mean.total.harv
                         ,mean.age.harv = mean.age.harv[[i]])
  
  DDsurv[[i]] = lm(mean.age.surv~mean.total.harv,data = data.temp)
  DDferc[[i]] = lm(mean.age.ferc~mean.total.harv,data = data.temp)
  harv_surv[[i]] = lm(mean.age.surv~mean.age.harv,data = data.temp)
  harv_ferc[[i]] = lm(mean.age.ferc~mean.age.harv,data = data.temp)
  
}


#mean.surv.1 = mean.surv[(1:(nage*period))%%nage==1]

lapply(DDsurv, summary)
lapply(DDferc, summary)
lapply(harv_surv, summary)
lapply(harv_ferc, summary)


#mean.harv.2 = mean.harv[(1:(nage*period))%%nage==2]
plot(mean.total.harv,mean.age.surv[[7]])




require(ggplot2)

ggplot(DDsurv[[8]]$model,aes(x=mean.total.harv,y=mean.age.surv) )+
  geom_point()+
  geom_line()+
  stat_smooth(method = "lm")


mean.ferc.2 = mean.ferc[(1:(nage*period))%%nage==2]
plot(mean.ferc.2)

plot(mean.total.harv,mean.ferc.2)
