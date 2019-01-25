source('DDLislie_harvest_1overK0.r')
nage = 8
period = 14


prop.vars = list(fert.rate = matrix(.001,nrow = nage,ncol = period),
                 surv.prop = matrix(.001,nrow = nage, ncol = period),
                 H = matrix(.001,nrow = nage,ncol = 1),
                 aK0=.0001,
                 baseline.pop.count = matrix(.001,nrow = nage,ncol = 1))

mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,2]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,2]
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,-1]
mean.b = (584 * Harv.data[,1]/sum(Harv.data[,1]))/0.7

prop.vars = list(fert.rate = matrix(.001,nrow = nage,ncol = period),
                 surv.prop = matrix(.001,nrow = nage, ncol = period),
                 H = matrix(.001,nrow = nage,ncol = 1),
                 aK0=.0001,
                 baseline.pop.count = matrix(.001,nrow = nage,ncol = 1))

set.seed(42)
Chicago_RES = HDDLislie.sampler( n.iter = 50000, burn.in = 300, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = .01, al.s = 1, be.s = .01
                                   , al.aK0 = 1, be.aK0 = 4e-4, al.n = 1
                                   , be.n = .1, al.H = 1, be.H = .1
                                   , mean.s = as.matrix(mean.s), mean.b= as.matrix(mean.b),mean.aK0 = matrix(0,1,1)
                                   , mean.H = matrix(0.8,nage,1), Harv.data = as.matrix(Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, estFer = T,nage = nage,homo = T)


hist(Chicago_RES$invK0.mcmc)
invK0_post = data.frame(invK0 = Chicago_RES$invK0.mcmc)

ggplot(invK0_post,aes(x=invK0)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.05) +
  geom_vline(xintercept = 0)
ggsave("./figs/Harv_reconstruct/invK0post.jpg")

mean.harv = apply(Chicago_RES$lx.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
#harv_mean = data.frame(age = 1:8,mean.harv.matrix)
mean.total.harv = apply(mean.harv.matrix,2,sum)
plot(mean.total.harv)

BI.low.harv = apply(Chicago_RES$lx.mcmc,2,quantile,probs = .05)
BI.low.harv.matrix = matrix(BI.low.harv,nrow = nage,ncol = period)
BI_harv_low = data.frame(age = 1:8,BI.low.harv.matrix)


BI.high.harv = apply(Chicago_RES$lx.mcmc,2,quantile,probs = .95)
BI.high.harv.matrix = matrix(BI.high.harv,nrow = nage,ncol = period)
BI_harv_high = data.frame(age = 1:8,BI.high.harv.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","mean","low","high","time")
har_data = har_data[-1,]

for(i in 1:8){
  temp = data.frame(age = i,mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1:14)
  har_data = rbind(har_data,temp)
}


for(i in 1:8){
  temp1 = data.frame(point = "model predict",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1:14)
  temp2 = data.frame(point = "data",mean = t(Harv.data[i,2:15]),low =t( Harv.data[i,2:15]),high = t(Harv.data[i,2:15]),time = 1:14)
  colnames(temp2) = colnames(temp1)
  temp = rbind(temp1,temp2)
  write.csv(temp,paste0("./figs/Harv_reconstruct/age",i,".csv"))
  rm(temp1)
  rm(temp2)
  filename = paste0("./figs/Harv_reconstruct/age",i,".jpg")
  #jpeg(filename)
  ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
    geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
    #geom_line() +
    geom_point()
  #Sys.sleep(5)
  ggsave(filename, plot = last_plot())
  #dev.off()
}




mean.ferc = apply(Chicago_RES$fert.rate.mcmc,2,mean)
BI.low.ferc = apply(Chicago_RES$fert.rate.mcmc,2,quantile,probs = .05)
BI.high.ferc = apply(Chicago_RES$fert.rate.mcmc,2,quantile,probs = .95)


require(ggplot2)
ferc_post = data.frame(age = 1:8,mean_ferc = mean.ferc,BI_low = BI.low.ferc,BI_high = BI.high.ferc)
write.csv(ferc_post,paste0("./figs/Harv_reconstruct/ferc_post",".csv"))
ggplot(ferc_post, aes(x=age, y=mean_ferc)) + 
  geom_errorbar(aes(ymin=BI_low, ymax=BI_high), width=.1) +
  geom_line() +
  geom_point()
ggsave("./figs/Harv_reconstruct/ferc.jpg")

mean.surv = apply(Chicago_RES$surv.prop.mcmc,2,mean)
BI.low.surv = apply(Chicago_RES$surv.prop.mcmc,2,quantile,probs = .05)
BI.high.surv = apply(Chicago_RES$surv.prop.mcmc,2,quantile,probs = .95)


require(ggplot2)
surv_post = data.frame(age = 1:8,mean_surv = mean.surv,BI_low = BI.low.surv,BI_high = BI.high.surv)
write.csv(surv_post,"./figs/Harv_reconstruct/surv_post.csv")
ggplot(surv_post, aes(x=age, y=mean_surv)) + 
  geom_errorbar(aes(ymin=BI_low, ymax=BI_high), width=.1) +
  geom_line() +
  geom_point()
ggsave("./figs/Harv_reconstruct/surv_post.jpg")

mean.harv = apply(Chicago_RES$H.mcmc,2,mean)
BI.low.harv = apply(Chicago_RES$H.mcmc,2,quantile,probs = .05)
BI.high.harv = apply(Chicago_RES$H.mcmc,2,quantile,probs = .95)


require(ggplot2)
harv_post = data.frame(age = 1:8,mean_Harv = mean.harv,BI_low = BI.low.harv,BI_high = BI.high.harv)
write.csv(harv_post,"./figs/Harv_reconstruct/harv_post.csv")
ggplot(harv_post, aes(x=age, y=mean_Harv)) + 
  geom_errorbar(aes(ymin=BI_low, ymax=BI_high), width=.1) +
  geom_line() +
  geom_point()
ggsave("./figs/Harv_reconstruct/harv_post.jpg")
