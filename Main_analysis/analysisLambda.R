analyisoflambda =
  analysisLambda(Chicago_RES$mcmc.objs,Assumptions,nage,14)

mean_lambda = Reduce("+",analyisoflambda)/nrow(Chicago_RES$mcmc.objs$surv.prop.mcmc)

list_each_lambda = lapply(1:length(mean_lambda),function(i,analyisoflambda){
  temp = lapply(analyisoflambda,function(ana,i){ana[i]},i)
  Reduce(rbind,temp)
},analyisoflambda)

lower_025 = sapply(list_each_lambda,quantile,probs = .025)
lower_025 = matrix(lower_025,ncol=14)

higher_975 = sapply(list_each_lambda,quantile,probs = .975)
higher_975 = matrix(higher_975,ncol=14)

year = 1:14+1992

observed = data.frame(point = "observed (w/ culling)"
                      ,lambda=mean_lambda[6,]
                      ,low = lower_025[6,]
                      ,high = higher_975[6,]
                      ,year = year)
even = data.frame(point = "uniform age structure"
                      ,lambda=mean_lambda[2,]
                      ,low = lower_025[2,]
                      ,high = higher_975[2,]
                      ,year = year)
nocull = data.frame(point = "skip culling"
                      ,lambda=mean_lambda[4,]
                      ,low = lower_025[4,]
                      ,high = higher_975[4,]
                      ,year = year)
stable = data.frame(point = "stable age structure"
                    ,lambda=mean_lambda[3,]
                    ,low = lower_025[3,]
                    ,high = higher_975[3,]
                    ,year = year)

plot_data = rbind(observed,even,nocull,stable)

require(ggplot2)

ggplot(data = plot_data,aes(x=year,y=lambda,color=point))+
  geom_line()+
  geom_point() +
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  labs(y = "Lambda   X(t+1)/X(t)")

ggsave("./lambda.jpg", plot = last_plot(),dpi = 500)
