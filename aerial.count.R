mean.total.aeri = apply(Chicago_RES$aerial.count.mcmc,2,median)
#mean.aeri.matrix = matrix(mean.aeri,nrow = sum(nage),ncol = period)
#aeri_mean = data.frame(age = 1:8,mean.aeri.matrix)
#mean.total.aeri = apply(mean.aeri.matrix,2,sum)
plot(mean.total.aeri)

BI.low.aeri = apply(Chicago_RES$aerial.count.mcmc,2,quantile,probs = .025)


BI.high.aeri = apply(Chicago_RES$aerial.count.mcmc,2,quantile,probs = .975)
#BI.high.aeri.matrix = matrix(BI.high.aeri,nrow = sum(nage),ncol = period)
#BI_aeri_high = data.frame(age = c(paste0("F",1:8),paste0("M",1:3)),BI.high.aeri.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","median","low","high","time")
har_data = har_data[-1,]

temp1 = data.frame(point = "model predict (95% CI)",mean = mean.total.aeri,low = BI.low.aeri,high = BI.high.aeri,time = 1992:2006)
temp2 = data.frame(point = "data",mean =t( Aeri.data),low =t(Aeri.data),high = t(Aeri.data),time = 1992:2006)
colnames(temp2) = colnames(temp1)
temp = rbind(temp1,temp2)
write.csv(temp,paste0("./figs/temp/","Aerial_count.csv"))
rm(temp1)
rm(temp2)
filename = paste0("./figs/temp/","Aerial_count.jpg")
#jpeg(filename)
ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  #geom_line() +
  geom_point()
#Sys.sleep(5)
ggsave(filename, plot = last_plot())
