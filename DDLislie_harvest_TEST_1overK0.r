source('DDLislie_harvest_1overK0.r')

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.2,0.15,0.15))
non0ferc = c(2,3)
nage = 3
aK0 = c(-1/1000,-1/1000)
period = 10
#Survival_inhomo = matrix( rep((Survival),period),nrow = nage,ncol = period )
#Ferc_inhomo = matrix( rep((Ferc),period),nrow = nage,ncol = period )

### Test of helper functions, DONE 0122
Lislie = getLislie(Survival,Ferc,F)
H = GetHarvest(Harvpar,nage)
H_prime = GetHarvest(0.05,nage) # age imspecific harvest rate
D = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),aK0,F)
I = DensityDependcy(T,c(1,1,1),E0=c(0.8,.1,0.1),aK0,T)
ProjectHarvest_helper(c(5,5,5), Lislie, Harvpar, global=T, E0=NULL, aK0=aK0, null = F)
data_homo = 
  ProjectHarvest_homo(Survival=Survival,Harvpar=Harvpar,Ferc=Ferc,E0=NULL,aK0=aK0,global = T,null=F,bl=c(10,10,10),period = period,nage)
data_homo = data_homo + matrix(rnorm(nrow(data_homo)*ncol(data_homo)),nrow = nrow(data_homo),ncol(data_homo))
data_homo = round(data_homo)
#data_inhomo = ProjectHarvest_inhomo(Survival=Survival_inhomo,Harvpar=Harvpar,Ferc=Ferc_inhomo,E0=NULL,K0=K0,global = T,null=F,bl=as.matrix(c(10,10,10)),period = period,nage)
#living_homo = getLivingIdividuals(Harvpar,data_homo)
#log.lhood(log(data_homo),log(data_homo+matrix(rnorm(24,0,0.05),3,11)),1)

# skip the simplest functions

### Test of sampler function

mean.f = Ferc
mean.s = Survival
mean.b = as.matrix( c(10,10,10)/Harvpar)
mean.aK0 = 0
mean.H = Harvpar
Harv.data = floor(data_homo) # this harv have density dependency
prop.vars = list(fert.rate = matrix(.01,nrow = length(non0ferc),ncol = period),
                 surv.prop = matrix(.01,nrow = nage, ncol = period),
                 H = matrix(.01,nrow = nage,ncol = 1),
                 aK0=.001,
                 baseline.pop.count = matrix(.01,nrow = nage,ncol = 1)) # should be a list, fert.rate for all age classes can give birth to new 
proj.periods = (ncol(Harv.data)-1)
ptm = proc.time()

### SIMULATION START HERE ###

 ## statistical inference of 1/K0
SIMULATION_RESDD = HDDLislie.sampler(n.iter = 50000, burn.in = 300, mean.f = mean.f
                                   ,al.f = 1, be.f = .1, al.s = 1, be.s = .1
                                   , al.aK0 = 1, be.aK0 = 4e-6, al.n = 1
                                   , be.n = .1, al.H = 1, be.H = .1
                                   , mean.s = mean.s, mean.b=mean.b,mean.aK0 = 0
                                   , mean.H = mean.H, Harv.data = (Harv.data+1e-4 * (Harv.data==0))
                                   , prop.vars = prop.vars, null = F, estaK0 = T,homo = T)
mean.invK0 = apply(SIMULATION_RESDD$invK0.mcmc,2,mean)
plot(SIMULATION_RESDD$invK0.mvmv[,1])
hist(SIMULATION_RESDD$invK0.mcmc)
sum(SIMULATION_RESDD$invK0.mcmc>0)/50000

invK0_post = data.frame(invK0 = SIMULATION_RESDD$invK0.mcmc)

ggplot(invK0_post,aes(x=invK0)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.9) +
  geom_vline(xintercept = 0)
ggsave("./figs/Simu/invK0post.jpg", plot = last_plot())
mean.harv = apply(SIMULATION_RESDD$lx.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
#harv_mean = data.frame(age = 1:8,mean.harv.matrix)
mean.total.harv = apply(mean.harv.matrix,2,sum)
plot(mean.total.harv)

BI.low.harv = apply(SIMULATION_RESDD$lx.mcmc,2,quantile,probs = .025)
BI.low.harv.matrix = matrix(BI.low.harv,nrow = nage,ncol = period)
BI_harv_low = data.frame(age = 1:3,BI.low.harv.matrix)


BI.high.harv = apply(SIMULATION_RESDD$lx.mcmc,2,quantile,probs = .975)
BI.high.harv.matrix = matrix(BI.high.harv,nrow = nage,ncol = period)
BI_harv_high = data.frame(age = 1:3,BI.high.harv.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","mean","low","high","time")
har_data = har_data[-1,]



for(i in 1:3){
  temp1 = data.frame(point = "model predict",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1:10)
  temp2 = data.frame(point = "data",mean = matrix(Harv.data[i,-1]),low =matrix( Harv.data[i,-1]),high = matrix(Harv.data[i,-1]),time = 1:10)
  colnames(temp2) = colnames(temp1)
  temp = rbind(temp1,temp2)
  write.csv(temp,paste0("./figs/Simu/age",i,".csv"))
  rm(temp1)
  rm(temp2)
  filename = paste0("./figs/Simu/age",i,".jpg")
  #jpeg(filename)
  ggplot(temp,aes(x=time, y=mean, colour = point)) + 
    geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
    #geom_line() +
    geom_point()
  #Sys.sleep(5)
  ggsave(filename, plot = last_plot())
  #dev.off()
}




mean.ferc = apply(SIMULATION_RESDD$fert.rate.mcmc,2,mean)
BI.low.ferc = apply(SIMULATION_RESDD$fert.rate.mcmc,2,quantile,probs = .025)
BI.high.ferc = apply(SIMULATION_RESDD$fert.rate.mcmc,2,quantile,probs = .975)


require(ggplot2)
ferc_post = data.frame(age = 2:3,mean_ferc = mean.ferc,BI_low = BI.low.ferc,BI_high = BI.high.ferc)
write.csv(ferc_post,paste0("./figs/Simu/ferc_post",".csv"))
ggplot(ferc_post, aes(x=age, y=mean_ferc)) + 
  geom_errorbar(aes(ymin=BI_low, ymax=BI_high), width=.1) +
  geom_line() +
  geom_point()
ggsave("./ferc.jpg")

mean.surv = apply(SIMULATION_RESDD$surv.prop.mcmc,2,mean)
BI.low.surv = apply(SIMULATION_RESDD$surv.prop.mcmc,2,quantile,probs = .025)
BI.high.surv = apply(SIMULATION_RESDD$surv.prop.mcmc,2,quantile,probs = .975)


require(ggplot2)
surv_post = data.frame(age = 1:3,mean_surv = mean.surv,BI_low = BI.low.surv,BI_high = BI.high.surv)
write.csv(surv_post,"./figs/Simu/surv_post.csv")
ggplot(surv_post, aes(x=age, y=mean_surv)) + 
  geom_errorbar(aes(ymin=BI_low, ymax=BI_high), width=.1) +
  geom_line() +
  geom_point()
ggsave("./figs/Simu/surv_post.jpg")

