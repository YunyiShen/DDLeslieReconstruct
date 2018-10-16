# get Lislie Matrix using specific parameters, can also return L-I if needed
getLislie = function(Survival,Fercest=T,Ferc=NULL,minus1=F){
  nage = length(Survival)/(1 + Fercest)
  Tr = matrix(0,nage,nage)
  if(!Fercest){
    Tr[1,] = Ferc
  }
  else{
    Tr[1,] = Survival[(nage+1):(2*nage)]
  }
  Nextage = Survival[1:(nage-1)]
  Nextage = exp(Nextage)/(1+exp(Nextage))
  diag(Tr[2:(nage),1:(nage-1)]) = Nextage
  Tr[nage,nage] = exp( Survival[nage])/(1+exp(Survival[nage]))
  if(minus1){
    diag(Tr)=diag(Tr)-1
  }
  return(Tr)
}

# Linear density dependency matrix, can be global of age specific
DensityDependcy = function(global = F, Xn, E0, K0){
  nage = length(Xn)
  D = matrix(0,ncol = nage,nrow = nage)
  if(global){
    den = as.matrix( 1-sum(Xn)/K0)
    diag(D) = den
  }
  else{
    diag(D) = 1-Xn/(K0 * E0)
  }
  return(D)
}

# get age specific harvest proportion matrix
GetHarvest = function(Harvpar,nage){
  Harvest = matrix(0,nage,nage)
  diag(Harvest) =exp( Harvpar)/(1+exp(Harvpar))
  return(Harvest)
}

# Project the harvest-after-reproducing model from year i to i+1, given harvest # and return harvest # of year i+1 
ProjectHarvest_helper = function(data_n, Lislie, H, global, E0, K0){
	I = matrix(0,nrow(data),nrow(data))	
	diag(I) = 1
    X_n1 = (I-H) %*% solve(H,data_n)
    D = DensityDependcy(global = global, Xn=X_n1, E0=E0, K0=K0)
	data_n1 = H %*% (Lislie %*% D %*% X_n1 + X_n1)
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  }

# Project harvest model from a initial harvest
ProjectHarvest = function(pars,Fercest=T,Ferc=NULL, estK0 = T, E0, K0 = NULL, global = F,data){
  period = ncol(data)
  nage = length(E0)
  if(estK0){
    K0 = pars[length(pars)]
    pars = pars[-length(pars)]
  }
  Survival = pars[1:(nage * (1 + Fercest))]
  Harvpar = pars[-(1:(nage * (1 + Fercest)))]
  Lislie = getLislie(Survival,Fercest=Fercest,Ferc=Ferc,minus1=T)
  H = GetHarvest(Harvpar,nage)
  Harvest[,1] = data[,1] 
  for(i in 2 : period){
	Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Lislie = Lislie, E0=E0, K0=K0,H=H) # project harvest from very initial
  }              
  return(Harvest)
}

# Given harvest data, calculate the living individual in certain year
getLivingIdividuals = function(H,data){
  I = matrix(0,nrow(data),nrow(data))	
  diag(I) = 1
  IminusH = I-H
  LivingIdividuals = apply(data,2,function(ww,H,IminusH){IminusH%*%solve(H,ww)},H=H,IminusH=IminusH)
}

# Problem is how to deal with the noise in the time series, 
#  least square given here will cause error to accumulate
sqrResidue = function(pars,Fercest=T,Ferc=NULL, estK0 = T, E0, K0 = NULL, global = F,data){
  period = ncol(data)
  nage = length(E0)
  if(estK0){
    K0 = pars[length(pars)]
    pars = pars[-length(pars)]
  }
  Survival = pars[1:(nage * (1 + Fercest))]
  Harvpar = pars[-(1:(nage * (1 + Fercest)))]
  Lislie = getLislie(Survival,Fercest=Fercest,Ferc=Ferc)
  H = GetHarvest(Harvpar,nage)
  Harest = apply(data[,-period],2,
                 function(data_n,
                          Lislie,
                          H,
                          global, 
                          E0, 
                          K0){
	I = matrix(0,nrow(data),nrow(data))	
	diag(I) = 1
    X_n1 = (I-H) %*% solve(H,data_n)
    D = DensityDependcy(global = global, Xn=X_n1, E0=E0, K0=K0)
	data_n1 = H %*% (Lislie %*% D %*% X_n1 + X_n1)
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  },global = global, Lislie = Lislie, E0=E0, K0=K0,H=H)
  Residue = sum((data[,-1]-Harest)^2)
  return(Residue)
}

########
# Coming functions are
# modified from popReconstruct package
########
### ----------------------- HELPER FUNCTIONS ---------------------- ###
### --------------------------------------------------------------- ###


## ........... Misc Functions .......... ##
## ..................................... ##

logitf = function(p) log(p / (1 - p))
invlogit = function(x){
    if(any(is.infinite(exp(x)))) {
        y = x
        y[is.infinite(exp(x))] = 1
        y[!is.infinite(exp(x))] =
            invlogit(y[!is.infinite(exp(x))])
        return(y)
    }
    else return(exp(x) / (1 + exp(x)))
}

##--- Generates random draws from inverse gamma ---##
rinvGamma = function(n, shape, scale){
    return(1/rgamma(n, shape = shape, rate = scale))
}

##--- Returns value of inverse gamma pdf ---##
dinvGamma = function(x, shape, scale, log = FALSE){
    if(log) d <-
        shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
    else d <- scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
    return(d)
}

## ............. Likelihood ............ ##
## ..................................... ##

log.lhood =function(log.n.census, log.n.hat, ll.var){
    ##.. log.n.census and log.n.hat should already be logged

    ##.. log.n.hat should be log projected counts for census years
    ##   interpolated if necessary


    ##-- value of log likelihoods --##

    density = dnorm(log.n.census,
                     mean = log.n.hat,
                     sd = sqrt(ll.var),
                     log = TRUE
                     )

    ##-- joint log likelihood --##

    return(sum(density))
}

