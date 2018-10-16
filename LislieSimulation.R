source('DDLislie_harvest.R')

Survival = matrix(c(0.5,-.1,-.2),nrow = 3,ncol = 1)
Ferc = matrix(c(0,3,1.5),nrow = 3,ncol = 1)
nage = 3
Lislie = getLislie(Survival,F,Ferc,T)
Havpar = as.matrix(c(-1.50,-1.90,-1.50))
H = GetHarvest(Havpar,3)
K0 = 800

nperiod = 100

Population = matrix(0,ncol = nperiod + 1,nrow = 3)
Population[,1] = c(100,100,100)
data = Population
I = matrix(0,nage,nage)
diag(I) = 1

for(i in 1:nperiod + 1){
  Repoduce = (Lislie %*% DensityDependcy(global = T, Xn=Population[,i-1], E0=c(.5,.3,.2), K0=K0) %*% Population[,i-1] + Population[,i-1])
  Population[,i] = (I-H) %*% Repoduce
  data[,i-1] = ( H %*% Repoduce + rnorm(3,0,.1))
}

data = data[,-(nperiod + 1)]
pars = rbind(as.matrix(Survival),Havpar)

sqrResidue(pars,Fercest=F,Ferc=Ferc, estK0 = F, E0=c(.3,.5,.2), K0 = K0, global = T,data)

opt = optim(pars + rnorm(length(pars),0,0.01),
            sqrResidue,NULL, 
            Fercest=F,Ferc=Ferc, estK0 = F, E0=c(.3,.5,.2), K0 = K0, global = T,data=data)

pars_K0 = c(pars,K0)
opt_estK0 = optim(pars_K0 ,
                  sqrResidue,NULL, 
                  Fercest=F,Ferc=Ferc, estK0 = T, E0=c(.3,.5,.2), K0 = NULL, global = T,data=data)


