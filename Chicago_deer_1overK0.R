source('DDLeslie_harvest_1overK0.r')
nage = c(8,3) # nage is female first and then male
period = 14


mean.s = read.csv("./data/Survival_mean_Etter.csv")
mean.s = mean.s[,-(1)]
mean.f = read.csv("./data/Fecundity_mean.csv")
mean.f = mean.f[,-(1)]
mean.SRB = read.csv("./data/SRB_mean.csv")
mean.SRB = mean.SRB[,-1]
#mean.f[1,] = mean.f[1,]+1e-6
Harv.data = read.csv("./data/Culling.csv")
Harv.data = Harv.data[,-(1)]
Aeri.data = read.csv("./data/Aerial_count.csv")
Aeri.data = Aeri.data[,-1]
mean.b = (977 * Harv.data[,1]/sum(Harv.data[,1]))



#mean.b[7]=10
mean.H = 0.7 * matrix(c(.3,.5,.4,.5,0.5,0.6,0.4,0.5,0.4,0.3,0.4,0.5,0.5,0.5,0.5))
mean.H = matrix(rep(mean.H,3),3,period+1,byrow = T)
mean.H[1,]=0.5*mean.H[1,]

mean.A = matrix(0.7,1,period+1)

prop.vars = list(fert.rate = matrix(.1,nrow = nage[1],ncol = period),
                 surv.prop = matrix(.1,nrow = sum(nage), ncol = period),
                 SRB = matrix(.1,nage[1],period), # vital rates has period cols
                 A = matrix(.1,1,period+1),
                 H = matrix(.1,nrow = 3,ncol = period+1),
                 aK0=1e-3,
                 baseline.pop.count = matrix(.1,nrow = sum(nage),ncol = 1))

set.seed(42)

Chicago_RES = HDDLislie.sampler( n.iter = 5000, burn.in = 1000, mean.f = as.matrix( mean.f)
                                   ,al.f = 1, be.f = .08, al.s = 1, be.s = .1
                                   , al.SRB = 1, be.SRB = .01
                                   , al.aK0 = 1, be.aK0 = 1e-2
                                   , al.n = 1, be.n = .01
                                   , al.ae = 1, be.ae = .01
                                   , al.H = 1, be.H = .1
                                   , al.A = 1, be.A = .1
                                   , mean.s = as.matrix(mean.s)
                                   , mean.b= as.matrix(mean.b)
                                   , mean.aK0 = matrix(0,1,2)
                                   , mean.H = (mean.H)
                                   , mean.SRB = as.matrix( mean.SRB)
                                   , mean.A = as.matrix( mean.A)
                                   #, mean.H = 0.6
                                   , Harv.data = as.matrix(Harv.data+1e-4 * (Harv.data==0))
                                   , Aerial.data = as.matrix( Aeri.data)
                                   , prop.vars = prop.vars, estFer = T,nage = nage,homo = F,estaK0 = F)


br_K0 = Chicago_RES$invK0.mcmc
plot(Chicago_RES$invK0.mcmc[,1])

hist(Chicago_RES$invK0.mcmc)
invK0_post = data.frame(invK0 = Chicago_RES$invK0.mcmc)

ggplot(invK0_post,aes(x=invK0)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.05) +
  geom_vline(xintercept = 0)
ggsave("./figs/temp/invK0post.jpg")

mean.harv = apply(Chicago_RES$harvest.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = sum(nage),ncol = period)
#harv_mean = data.frame(age = 1:8,mean.harv.matrix)
mean.total.harv = apply(mean.harv.matrix,2,sum)
plot(mean.total.harv)

BI.low.harv = apply(Chicago_RES$harvest.mcmc,2,quantile,probs = .025)
BI.low.harv.matrix = matrix(BI.low.harv,nrow = sum(nage),ncol = period)
BI_harv_low = data.frame(age = c(paste0("F",1:8),paste0("M",1:3)),BI.low.harv.matrix)


BI.high.harv = apply(Chicago_RES$harvest.mcmc,2,quantile,probs = .975)
BI.high.harv.matrix = matrix(BI.high.harv,nrow = sum(nage),ncol = period)
BI_harv_high = data.frame(age = c(paste0("F",1:8),paste0("M",1:3)),BI.high.harv.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","mean","low","high","time")
har_data = har_data[-1,]

for(i in 1:8){
  temp = data.frame(age = i,mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1993:2006)
  har_data = rbind(har_data,temp)
}
require(ggplot2)

for(i in 1:11){
  temp1 = data.frame(point = "model predict (95% CI)",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1993:2006)
  temp2 = data.frame(point = "data",mean = t(Harv.data[i,2:15]),low =t( Harv.data[i,2:15]),high = t(Harv.data[i,2:15]),time = 1993:2006)
  colnames(temp2) = colnames(temp1)
  temp = rbind(temp1,temp2)
  write.csv(temp,paste0("./figs/temp/",BI_harv_high$age[i],".csv"))
  rm(temp1)
  rm(temp2)
  filename = paste0("./figs/temp/",BI_harv_high$age[i],".jpg")
  #jpeg(filename)
  ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
    geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
    #geom_line() +
    geom_point()
  #Sys.sleep(5)
  ggsave(filename, plot = last_plot())
  #dev.off()
}
H_full = 0 * Chicago_RES$lx.mcmc
H_full = cbind(H_full[,1:nage],H_full)

for(i in 0:period+1){
  for(j in 1:nage){
    ind = (i-1)*nage + j
    plc = (i-1)*2 + 2 - (j==1)
    H_full[,ind]=Chicago_RES$H.mcmc[,plc]
    
  }
}



living_inid = (Chicago_RES$lx.mcmc/H_full[,-(1:nage)])*(1-H_full[,-(1:nage)])

plotthings(YD_obj=living_inid,pathsave="./figs/temp/living_af_culling_age",nage,period,1993:2006)

total_living = matrix(0,nrow(living_inid),ncol = period)

for(i in 1:period){
  total_living[,i] = rowSums(living_inid[,1:nage+(i-1)*nage])
  
}
bl = rowSums(Chicago_RES$baseline.count.mcmc*(1-H_full[,(1:nage)]))
bl_mean = colMeans(Chicago_RES$baseline.count.mcmc*(1-H_full[,(1:nage)]))
total_living_bl = cbind(bl,total_living)

plotthings(YD_obj=total_living_bl,pathsave="./figs/temp/living_af_culling_all",1,period+1,1992:2006)

plotthings(YD_obj=(Chicago_RES$surv.prop.mcmc),pathsave="./figs/temp/survival_age",nage=8,period,1993:2006)
plotthings(Chicago_RES$fert.rate.mcmc,pathsave="./figs/temp/fec_age",nage=7,period,1993:2006)
plotthings(YD_obj=Chicago_RES$H.mcmc,pathsave="./figs/temp/Harvpor",2,period+1,1992:2006)

