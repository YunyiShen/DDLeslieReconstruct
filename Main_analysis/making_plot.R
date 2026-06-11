require(ggplot2)
source('./R/misc.R')
# Just for plotting things


## Harvest prediction
mean.harv = apply(Chicago_RES$mcmc.objs$harvest.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,nrow = sum(nage),ncol = period)
mean.total.harv = apply(mean.harv.matrix,2,sum)
plot(mean.total.harv)

BI.low.harv = apply(Chicago_RES$mcmc.objs$harvest.mcmc,2,quantile,probs = .025)
BI.low.harv.matrix = matrix(BI.low.harv,nrow = sum(nage),ncol = period)
BI_harv_low = data.frame(age = c(paste0("F",1:8),paste0("M",1:3)),BI.low.harv.matrix)


BI.high.harv = apply(Chicago_RES$mcmc.objs$harvest.mcmc,2,quantile,probs = .975)
BI.high.harv.matrix = matrix(BI.high.harv,nrow = sum(nage),ncol = period)
BI_harv_high = data.frame(age = c(paste0("F",1:8),paste0("M",1:3)),BI.high.harv.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","mean","low","high","time")
har_data = har_data[-1,]

require(ggplot2)

for(i in 1:11){
  temp1 = data.frame(point = "model predict (95% CI)",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1993:2006)
  temp2 = data.frame(point = "data",mean = t(Harv.data[i,2:15]),low =t( Harv.data[i,2:15]),high = t(Harv.data[i,2:15]),time = 1993:2006)
  colnames(temp2) = colnames(temp1)
  temp = rbind(temp1,temp2)
  write.csv(temp,paste0("./Main_analysis/figs/temp/",BI_harv_high$age[i],".csv"))
  rm(temp1)
  rm(temp2)
  filename = paste0("./Main_analysis/figs/temp/",BI_harv_high$age[i],".jpg")
  ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
    geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
    geom_point()
  ggsave(filename, plot = last_plot())
}


## for ploting living inividuals
H_full = 0 * Chicago_RES$mcmc.objs$harvest.mcmc
H_full = cbind(H_full[,1:sum(nage)],H_full)

for(i in 0:period+1){
  for(j in 1:sum(nage)){
    ind = (i-1)*sum(nage) + j
    plc = (i-1)*4 + (j==1) +  3*(j==9) + 2 * (j>=2 & j<=8) + 4 * (j>=10)
    H_full[,ind]=Chicago_RES$mcmc.objs$H.mcmc[,plc]
    
  }
}



living_inid = (Chicago_RES$mcmc.objs$harvest.mcmc/H_full[,-(1:sum(nage))])*(1-H_full[,-(1:sum(nage))])
plotthings(YD_obj=living_inid,pathsave="./Main_analysis/figs/temp/living_af_culling_age",sum(nage),period,1993:2006)

total_living = matrix(0,nrow(living_inid),ncol = period)

for(i in 1:period){
  total_living[,i] = rowSums(living_inid[,1:sum(nage)+(i-1)*sum(nage)])
  
}
bl = rowSums(Chicago_RES$mcmc.objs$baseline.count.mcmc/(H_full[,(1:sum(nage))])*(1-H_full[,(1:sum(nage))]))
bl_mean = colMeans(Chicago_RES$mcmc.objs$baseline.count.mcmc/(H_full[,(1:sum(nage))])*(1-H_full[,(1:sum(nage))]))
total_living_bl = cbind(bl,total_living)

plotthings(YD_obj=total_living_bl,pathsave="./Main_analysis/figs/temp/living_af_culling_all",1,period+1,1992:2006)

mean.total.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,median)
plot(mean.total.aeri)

BI.low.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,quantile,probs = .025)


BI.high.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,quantile,probs = .975)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","median","low","high","time")
har_data = har_data[-1,]

temp1 = data.frame(point = "model predict (95% CI)",mean = mean.total.aeri,low = BI.low.aeri,high = BI.high.aeri,time = 1992:2006)
temp2 = data.frame(point = "data",mean =t( Aeri.data),low =t(Aeri.data),high = t(Aeri.data),time = 1992:2006)
colnames(temp2) = colnames(temp1)
temp = rbind(temp1,temp2)
write.csv(temp,paste0("./Main_analysis/figs/temp/","Aerial_count.csv"))
rm(temp1)
rm(temp2)
filename = paste0("./Main_analysis/figs/temp/","Aerial_count.jpg")
ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  geom_point()
ggsave(filename, plot = last_plot())

## vital rates
plotthings(YD_obj=(Chicago_RES$mcmc.objs$surv.prop.mcmc),pathsave="./Main_analysis/figs/temp/survival_age",nage=11,period,1993:2006)
plotthings(Chicago_RES$mcmc.objs$fert.rate.mcmc,pathsave="./Main_analysis/figs/temp/fec_age",nage=8,period,1993:2006)
plotthings(YD_obj=Chicago_RES$mcmc.objs$H.mcmc,pathsave="./Main_analysis/figs/temp/Harvpor",4,period+1,1992:2006)
plotthings(YD_obj=Chicago_RES$mcmc.objs$SRB.mcmc,pathsave="./Main_analysis/figs/temp/SRB",1,period,1993:2006)
plotthings(YD_obj=Chicago_RES$mcmc.objs$aerial.detection.mcmc,pathsave="./Main_analysis/figs/temp/Aerial_detectionrate",1,period+1,1992:2006)


## Aerial counts prediction
mean.total.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,median)
plot(mean.total.aeri)

BI.low.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,quantile,probs = .025)
BI.high.aeri = apply(Chicago_RES$mcmc.objs$aerial.count.mcmc,2,quantile,probs = .975)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","median","low","high","time")
har_data = har_data[-1,]

temp1 = data.frame(point = "model predict (95% CI)",mean = mean.total.aeri,low = BI.low.aeri,high = BI.high.aeri,time = 1992:2006)
temp2 = data.frame(point = "data",mean =t( Aeri.data),low =t(Aeri.data),high = t(Aeri.data),time = 1992:2006)
colnames(temp2) = colnames(temp1)
temp = rbind(temp1,temp2)
write.csv(temp,paste0("./Main_analysis/figs/temp/","Aerial_count.csv"))
rm(temp1)
rm(temp2)
filename = paste0("./Main_analysis/figs/temp/","Aerial_count.jpg")
ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  geom_point()
ggsave(filename, plot = last_plot())


sd_est = lapply(Chicago_RES$mcmc.objs,apply,2,sd)

sd_fec = matrix(sd_est$fert.rate.mcmc,8,14)
sd_surv = matrix(sd_est$surv.prop.mcmc,11,14)

row.names(sd_fec) = 0:7+0.5
colnames(sd_fec) = 1992+1:14

row.names(sd_surv) = row.names(mean.s)
colnames(sd_surv) = 1:14+1992

write.csv(sd_fec,"fecundity_sd.csv")
write.csv(sd_surv,"survival_sd.csv")


last_year_inid = living_inid[,1:11 + 13 * 11]

laster_year_female_YA = last_year_inid[,2:8]

laster_year_male_YA = last_year_inid[,c(10,11)]

laster_year_female_YA_sum = rowSums(laster_year_female_YA)
laster_year_male_YA_sum = rowSums(laster_year_male_YA)
