source("DDtest.R")

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.2,0.2,0.2))
nage = 3
aK0 = 1/800
period = 20
bl=c(5,5,5)
nonspcharv = T
nonspcDD=T
ZeroFerc = 1

data_homo_lambda = ProjectHarvest_homo(Survival=Survival
                                ,Harvpar=Harvpar
                                ,Ferc=Ferc
                                ,aK0=aK0
                                ,E0 = c(1,1,1)
                                ,global = T
                                ,null=F
                                ,bl=bl
                                ,period = period
                                ,nage=nage)



data_homo = apply(data_homo_lambda,c(1,2),rpois,n=1)
data_homo = data_homo[,-1]
plot(data_homo[1,])

p = 3*nage + 1 + (!nonspcharv)*(nage-1) + (!nonspcDD) * nage + 1 # length of pars
initial = rnorm(p,0,.01)
initial = c(rnorm(p-1,0,.01),1/400)
initial[1:nage] = runif(nage)
initial[1:nage + 2*nage] = log( data_homo[,1] )

real = c(log(Ferc),log(Survival/(1-Survival)),log(bl),log(.2/(1-.2)),1/800)
noisy = real + rnorm(p,0,0.1)
noisy[p] = abs(noisy[p])

logLiklihood_fun(par=real,Harvdata=data_homo,nage,period, nonspcharv, nonspcDD,ZeroFerc)

kk = DDLislie_fit(data_homo,nage,period, nonspcharv, nonspcDD,ZeroFerc,initial=real,method = "BFGS",control = list(maxit = 1000))
