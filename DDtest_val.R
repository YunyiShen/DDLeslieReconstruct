source("DDtest.R")

Survival = matrix(c(0.8,0.6,0.4),nrow = 3,ncol = 1)
Ferc = matrix(c(0,1.3,1.5),nrow = 3,ncol = 1)
Harvpar = as.matrix(c(0.2,0.2,0.2))
nage = 3
aK0 = 1/800
period = 20
bl=c(15,15,15)
nonspcharv = T
nonspcDD=T
ZeroFerc = 1

Lislie = getLislie(Survival,Ferc)

data_homo_lambda = ProjectHarvest_homo(Survival=Survival
                                ,Harvpar=Harvpar
                                ,Ferc=Ferc
                                ,aK0=aK0
                                ,period = period
                                ,nage=nage
                                ,bl = bl)



data_homo = apply(data_homo_lambda,c(1,2),rpois,n=1)
plot(data_homo[1,])

p = 2*nage + 1 + (!nonspcharv)*(nage-1) + 1 # length of pars
initial = rnorm(p,0,.01)
initial = c(rnorm(p-1,0,.01),1/400)
initial[1:nage] = runif(nage)
#initial[1:nage + 2*nage] = log( data_homo[,1] )

real = c(log(Ferc),log(Survival/(1-Survival)),log(.2/(1-.2)),1/800)
noisy = real + rnorm(p,0,0.1)
noisy[p] = abs(noisy[p])

Lislie_noisy = getLislie(invlogit(noisy[1:nage + nage]),exp(noisy[1:nage]))

logLiklihood_fun(par=noisy,Harvdata=data_homo,nage,period, nonspcharv,ZeroFerc)

kk = DDLislie_fit(data_homo,nage,period, nonspcharv,ZeroFerc,initial=real,method = "BFGS",control = list(maxit = 1000))
