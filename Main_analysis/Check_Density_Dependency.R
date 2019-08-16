## Density Dependency, run after making_plot
mean.surv = apply(Chicago_RES$mcmc.objs$surv.prop.mcmc,2,median)
mean.ferc = apply(Chicago_RES$mcmc.objs$fert.rate.mcmc,2,median)

mean.living.total =c(median(bl), apply( total_living,2,median))

mean.age.surv = list()
mean.age.ferc = list()
DDsurv = list()
DDfec = list()

for(i in 1:sum(nage)){
  
  mean.age.surv[[i]] = mean.surv[(1:((sum(nage))*period))%%sum(nage)==(i%%(sum(nage)))]
  data.temp = data.frame(mean.age.surv = mean.age.surv[[i]]
                         ,mean.living.total = mean.living.total[1:period])
  DDsurv[[i]] = lm(mean.age.surv~mean.living.total,data = data.temp)
}

for(i in 1:nage[1]){
  
  mean.age.ferc[[i]] = mean.ferc[(1:(nage[1]*period))%%nage[1]==(i%%nage[1])]
  data.temp = data.frame( mean.age.fec = mean.age.ferc[[i]]
                         ,mean.living.total = mean.living.total[1:period])
  DDfec[[i]] = lm(mean.age.fec~mean.living.total,data = data.temp)
}



## check results
summary_DDsurv = lapply(DDsurv, summary)
summary_DDfec = lapply(DDfec, summary)

pvalues_DDsurv = Reduce(rbind, lapply(summary_DDsurv,function(w){w$coefficients[2,]}))
adjRsqr_DDsurv = unlist( lapply(summary_DDsurv,function(w){w$adj.r.squared}))

res_DDsurv = data.frame( pvalues_DDsurv,adjRsqr = adjRsqr_DDsurv)
write.csv(res_DDsurv,"./Main_analysis/figs/temp/DDsurv.csv")

pvalues_DDfec = Reduce(rbind, lapply(summary_DDfec,function(w){w$coefficients[2,]}))
adjRsqr_DDfec = unlist( lapply(summary_DDfec,function(w){w$adj.r.squared}))

res_DDfec = data.frame(pvalues_DDfec,adjRsqr = adjRsqr_DDfec)
write.csv(res_DDfec,"./Main_analysis/figs/temp/DDfec.csv")

require(ggplot2)

for(i in 1:nage[1]){
  ggplot(DDfec[[i]]$model,aes(x=mean.living.total,y=mean.age.fec) )+
    geom_point()+
    stat_smooth(method = "lm")
  filename = paste0("./Main_analysis/figs/temp/DDfec_age",i,".jpg")
  ggsave(filename, plot = last_plot())
}

for(i in 1:sum(nage)){
  ggplot(DDsurv[[i]]$model,aes(x=mean.living.total,y=mean.age.surv) )+
    geom_point()+
    stat_smooth(method = "lm")
  filename = paste0("./Main_analysis/figs/temp/DDsurv_age",i,".jpg")
  ggsave(filename, plot = last_plot())
}
