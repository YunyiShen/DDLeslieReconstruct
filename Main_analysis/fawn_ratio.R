##female fawn:Y+A
require(coda)
ratio_mcmc = mcmc(matrix(nrow = nrow(living_inid),ncol = (period)))
fawn_mcmc = ratio_mcmc
YpA_mcmc = ratio_mcmc

for (i in 1:(period)){
  temp1 = living_inid[,(i-1)*sum(nage)+1:nage[1]]
  temp = apply(temp1,1,function(w){w[1]/(sum(w)-w[1])})
  ratio_mcmc[,i] = temp
  temp = apply(temp1,1,function(w){w[1]})
  fawn_mcmc[,i] = temp
  temp = apply(temp1,1,function(w){sum(w)-w[1]})
  YpA_mcmc[,i] = temp
}

plotthings(YD_obj=ratio_mcmc,pathsave="./Main_analysis/figs/temp/female_fawn_ratio",1,period,1993:2006)
plotthings(YD_obj=fawn_mcmc,pathsave="./Main_analysis/figs/temp/female_fawn_count",1,period,1993:2006)
plotthings(YD_obj=YpA_mcmc,pathsave="./Main_analysis/figs/temp/female_YpA_count",1,period,1993:2006)

YpA_count = read.csv("./Main_analysis/figs/temp/female_YpA_count1.csv",row.names = 1,as.is=T)
fawn_count = read.csv("./Main_analysis/figs/temp/female_fawn_count1.csv",row.names = 1,as.is=T)

YpA_count$point=sub("model predict","Y+A",YpA_count$point)
fawn_count$point=sub("model predict","fawn",fawn_count$point)

plot_data = rbind(YpA_count,fawn_count)
ggplot(data = plot_data,aes(x=time,y=mean,color=point))+
  geom_line()+
  geom_point() + 
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) 

ggsave("./Main_analysis/figs/temp/fawn_YpA_count.jpg", plot = last_plot(),dpi = 500)
