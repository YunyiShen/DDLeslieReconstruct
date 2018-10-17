# get Lislie Matrix using specific parameters, can also return L-I if needed, first nage entries should be survival(inv logit transformed) and next fertility rate, 
getLislie = function(Survival,estFer=T,Ferc=NULL,minus1=F){
  nage = length(Survival)/(1 + estFer)
  Tr = matrix(0,nage,nage)
  if(!estFer){
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

# Linear density dependency matrix, can be global of age specific, or null = T for no density dependency
DensityDependcy = function(global = F, Xn, E0, K0, null = F){
  nage = length(Xn)
  D = matrix(0,ncol = nage,nrow = nage)
  if(global){
    den = as.matrix( 1-sum(Xn)/K0)
    diag(D) = (1-null)*den + null
  }
  else{
    diag(D) = 1-(1-null)*Xn/(K0 * E0)
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
ProjectHarvest_helper = function(data_n, Lislie, H, global, E0, K0, null = F){
	I = matrix(0,nrow(data),nrow(data))	
	diag(I) = 1
    X_n1 = (I-H) %*% solve(H,data_n)
    D = DensityDependcy(global = global, Xn=X_n1, E0=E0, K0=K0, null = null)
	data_n1 = H %*% (Lislie %*% D %*% X_n1 + X_n1)
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  }

# Project harvest model from a initial harvest
ProjectHarvest = function(pars,estFer=T,Ferc=NULL, estK0 = T, E0, K0 = NULL, global = F, null = F,data){
  period = ncol(data)
  nage = length(E0)
  if(estK0){
    K0 = pars[length(pars)]
    pars = pars[-length(pars)]
  }
  Survival = pars[1:(nage * (1 + estFer))]
  Harvpar = pars[-(1:(nage * (1 + estFer)))]
  Lislie = getLislie(Survival,estFer=estFer,Ferc=Ferc,minus1=T)
  H = GetHarvest(Harvpar,nage)
  Harvest[,1] = data[,1] 
  for(i in 2 : period){
	Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Lislie = Lislie, E0=E0, K0=K0,H=H,null = null) # project harvest from very initial
  }              
  return(Harvest)
}

# Given harvest data, calculate the living individual in certain year
getLivingIdividuals = function(H,data){
  I = matrix(0,nrow(data),nrow(data))	
  diag(I) = 1
  IminusH = I-H
  LivingIdividuals = apply(data,2,function(ww,H,IminusH){IminusH%*%solve(H,ww)},H=H,IminusH=IminusH)
  return(LivingIdividuals)
}

# Problem is how to deal with the noise in the time series, 
#  least square given here will cause error to accumulate
sqrResidue = function(pars,estFer=T,Ferc=NULL, estK0 = T, E0, K0 = NULL, global = F,data){
  period = ncol(data)
  nage = length(E0)
  if(estK0){
    K0 = pars[length(pars)]
    pars = pars[-length(pars)]
  }
  Survival = pars[1:(nage * (1 + estFer))]
  Harvpar = pars[-(1:(nage * (1 + estFer)))]
  Lislie = getLislie(Survival,estFer=estFer,Ferc=Ferc)
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

logitf = function(p){
  log(p / (1 - p))
 }

## logit function  
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

#log likelihood function of gaussian distributed
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

## .............. Posterior ............ ##
## ..................................... ##  keep it, change in sampler function, but no migration here, should add estK0

log.post = function(## estimated vitals
                           f, s, baseline.n, K0, H
						   , estFer, estK0
                           ## fixed prior means on vitals
                           ,prior.mean.f, prior.mean.s
                           ,prior.mean.b, prior.mean.K0
                           ## fixed prior parameters on variance distns
                           ,alpha.f, beta.f, alpha.s, beta.s
                           ,alpha.n, beta.n, alpha.K0, beta.K0
						   ,alpha.H, beta.H
                           ## updated variances on prior distns
                           ,sigmasq.f, sigmasq.s, sigmasq.n, sigmasq.K0
                           ## value of the log likelihood
                           ,log.like
                           ## non zero rows of fertility matrix
                           ,non.zero.fert
                           ){

    ##-- Values of prior densities for vitals --##

    ##.. Note that log densities are calculated for numerical stability.
    ##     f, baseline.n, prior.mean.f, prior.mean.b are logged coming
    ##     in, s, prior.mean.s is logit transformed coming in, g and
    ##     prior.mean.g are not transformed coming in.
	##-- prior for f and baseline K0 if needed to be estimatedd --##
	if(estFer){
	  log.f.prior = dnorm(as.vector(f[non.zero.fert,])
                         ,mean = as.vector(prior.mean.f[non.zero.fert,])
                         ,sd = sqrt(sigmasq.f)
                         ,log = TRUE)
	  log.sigmasq.f.prior =
        log(dinvGamma(sigmasq.f, alpha.f, beta.f))
	}
	else {
	  log.f.prior = 0
	  log.sigmasq.f.prior = 0
	}
    
	if(estK0){
	  log.K0.prior = dnorm(K0, mean = prior.mean.K0
						   , sd = sqrt(sigmasq.K0)
						   , log = T )
	  log.sigmasq.K0.prior =
        log(dinvGamma(sigmasq.K0, alpha.K0, beta.K0))
	}
	else {
	  log.K0.prior = 0
	  log.sigmasq.K0.prior = 0
	}

    ##-- prior for s and baseline n, and hunting rate H --##
	log.s.prior = dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                         ,log = TRUE)
    log.b.prior = dnorm(baseline.n, mean = prior.mean.b
                         ,sd = sqrt(sigmasq.n)
                         ,log = TRUE)
	log.H.prior = dnorm(H, mean = prior.mean.b
                         ,sd = sqrt(sigmasq.n)
                         ,log = TRUE)

    log.sigmasq.s.prior =
        log(dinvGamma(sigmasq.s, alpha.s, beta.s))
    log.sigmasq.n.prior =
        log(dinvGamma(sigmasq.n, alpha.n, beta.n))
	log.sigmasq.H.prior =
        log(dinvGamma(sigmasq.H, alpha.n, beta.n))
	
    ##-- The log posterior is the SUM of these with the log.like --##

    return(sum(log.f.prior, log.s.prior, log.b.prior, log.K0.prior, log.H.prior
               ,log.sigmasq.f.prior
               ,log.sigmasq.s.prior
               ,log.sigmasq.n.prior
			   ,log.sigmasq.K0.prior
			   ,log.sigmasq.H.prior
               ,log.like))

}

## ......... Acceptance Ratio .......... ##
## ..................................... ##
acc.ra = function(log.prop, log.current){
    min(1, exp(log.prop - log.current))
}

acc.ra.var = function(log.prop.post, log.curr.post, log.prop.var, log.curr.var){
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post))
}

