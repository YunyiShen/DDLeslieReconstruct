##female fawn:Y+A
require(coda)
ratio_mcmc = mcmc(matrix(nrow = nrow(living_inid),ncol = (period)))
for (i in 1:(period)){
  temp = living_inid[,(i-1)*sum(nage)+1:nage[1]]
  temp = apply(temp,1,function(w){w[1]/(sum(w)-w[1])})
  ratio_mcmc[,i] = temp
}

plotthings(YD_obj=ratio_mcmc,pathsave="./figs/temp/female_fawn_ratio",1,period,1993:2006)
