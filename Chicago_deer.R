source('DDLislie_harvest.R')
nage = 8
period = 14


prop.vars = list(fert.rate = matrix(.1,nrow = nage,ncol = period),
                 surv.prop = matrix(1,nrow = nage, ncol = period),
                 H = matrix(1,nrow = nage,ncol = 1),
                 K0=.1,
                 baseline.pop.count = matrix(1,nrow = nage,ncol = 1))

mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,-1]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,-1]
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,c(-1,-2)]
mean.b = (584 * Harv.data[,1]/sum(Harv.data[,1]))/0.7

Chicago_RES = HDDLislie.sampler( n.iter = 5000, burn.in = 300, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = .001, al.s = 1, be.s = .5
                                   , al.K0 = 1, be.K0 = .001, al.n = 1
                                   , be.n = .001, al.H = 1, be.H = .01
                                   , mean.s = as.matrix(mean.s), mean.b= as.matrix(mean.b),mean.K0 = matrix(700,1,1)
                                   , mean.H = matrix(0.7,nage,1), Harv.data = as.matrix(Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, estFer = T,nage = nage)

mean.surv = apply(Chicago_RES$surv.prop.mcmc,2,mean)
mean.ferc = apply(Chicago_RES$fert.rate.mcmc,2,mean)

mean.harv = apply(Chicago_RES$lx.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
mean.total.harv = apply(mean.harv.matrix,2,sum)

mean.age.surv = list()
mean.age.ferc = list()
DDsurv = list()
DDferc = list()

for(i in 1:7){
  
  mean.age.surv[[i]] = mean.surv[(1:(nage*period))%%nage==i]
  mean.age.ferc[[i]] = mean.ferc[(1:(nage*period))%%nage==i]
  data.temp = data.frame(mean.age.surv = mean.age.surv[[i]]
                         ,mean.age.ferc = mean.age.ferc[[i]]
                         ,mean.total.harv)
  DDsurv[[i]] = lm(mean.age.surv~mean.total.harv,data = data.temp)
  DDferc[[i]] = lm(mean.age.ferc~mean.total.harv,data = data.temp)
  
}

mean.age.surv[[8]] = mean.surv[(1:(nage*period))%%nage==0]
mean.age.ferc[[8]] = mean.ferc[(1:(nage*period))%%nage==0]
data.temp = data.frame(mean.age.surv = mean.age.surv[[8]]
                       ,mean.age.ferc = mean.age.ferc[[8]]
                       ,mean.total.harv)
DDsurv[[8]] = lm(mean.age.surv~mean.total.harv,data = data.temp)
DDferc[[8]] = lm(mean.age.ferc~mean.total.harv,data = data.temp)
#mean.surv.1 = mean.surv[(1:(nage*period))%%nage==1]

lapply(DDsurv, summary)
lapply(DDferc, summary)


#mean.harv.2 = mean.harv[(1:(nage*period))%%nage==2]
plot(mean.total.harv,mean.age.surv[[8]])

data.age.2 = data.frame(mean.total.harv,mean.surv.1)

for(i in 1:8){
  formulea = as.formula( paste0("mean.surv.",i,"~mean.total.harv"))
  DDsurv[[i]] = lm(formulea,data=data.frame(mean.total.harv))
  
}
summaries = lapply(DDsurv,summary)


require(ggplot2)

ggplot(data.age.2,aes(x=mean.total.harv,y=mean.surv.2) )+
  geom_point()+
  geom_line()+
  stat_smooth(method = "lm")


mean.ferc.2 = mean.ferc[(1:(nage*period))%%nage==2]
plot(mean.ferc.2)

plot(mean.total.harv,mean.ferc.2)
