getLislie = function(Survival,Ferc,minus1=F){
  nage = length(Survival)
  Tr = matrix(0,nage,nage)
  Tr[1,] = Ferc
  diag(Tr[2:(nage),1:(nage-1)]) = Survival[-nage]
  Tr[nage,nage] = Survival[nage]
  if(minus1){
    diag(Tr)=diag(Tr)-1
  }
  return(Tr)
} #checked 10/24/2018

# Linear density dependency matrix, can be global of age specific, or null = T for no density dependency. in estimaton process, if E0 is not fixed, it should be estimated as kk = solve(H, data[,1]) E0 = kk/sum(kk), by assuming that first year is at equalibrium
# get age specific harvest proportion matrix
GetHarvest = function(Harvpar,nage){
  Harvest = matrix(0,nage,nage)
  diag(Harvest) = Harvpar
  return(Harvest)
} # checked 10/24/2018, if not age-specific, give a single Harvpar

# Project the (density dependent) harvest-after-reproducing model from year i to i+1, given harvest # and return harvest # of year i+1 
ProjectHarvest_helper = function(data_n, Lislie, Harvpar, aK0, null = F){
	# I = matrix(0,length(data_n),length(data_n))	
	# diag(I) = 1
    X_n1 = (1-Harvpar) * (data_n/Harvpar) # population remained to reproduce
	data_n1 = Harvpar *  ((1 - aK0 * sum(X_n1)) * Lislie %*% X_n1 + X_n1)
	data_n1= data_n1 * (data_n1>0) + (1e-15 * (data_n1<=0)) # avoid NAs 
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  } # checked 10/24/2018

# Project harvest model from a initial harvest, survival is col vector with all survival rate of all age class, this nrow(Survival)=nage
ProjectHarvest_homo = function(Survival, Harvpar,Ferc, aK0 , period, nage,bl){
  Lislie = getLislie(Survival,Ferc = Ferc,minus1=T)
  H = GetHarvest(Harvpar,nage)
  Harvest = matrix(0,nage,period)
  Harvest[,1] = bl 
  for(i in 2 : period ){
	Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1], Lislie = Lislie, aK0=aK0, Harvpar=Harvpar) # project harvest from very initial
  }              
  return(Harvest)
} #checked 10/24/2018

logLiklihood_fun = function(par,Harvdata,nage,period, nonspcharv = T, ZeroFerc){
	p = length(par)
	Ferc = exp( par[1:nage])
	Ferc[ZeroFerc]=0
	Surv = invlogit( par[1:nage + nage] )
	#bl = Harvdata[,1]
	bl = exp(par[1:nage + 2*nage])
	if(nonspcharv) Harvpar = invlogit( par[3*nage + 1])
	else Harvpar = invlogit( par[1:nage + 3*nage])
	aK0 = par[p]
	#E0 = E0/sum(E0)
	Proj = ProjectHarvest_homo(Surv, Harvpar,Ferc, aK0 = aK0 , period, nage,bl)
	# Proj = Proj[,-1]
	logll = sum(log(dpois(Harvdata,Proj)+1e-15))
	#if(min(Proj)<=0 | sum(is.na(Proj))>0 | logll==-Inf) return(1e50)
	#cat(par,'\n')
	return(-logll)
}

invlogit = function(x){
    if(any(is.infinite(exp(x)))) {
        y = x
        y[is.infinite(exp(x))] = 1
        y[!is.infinite(exp(x))] =
            invlogit(y[!is.infinite(exp(x))])
        return(y)
    }
    else return(exp(x) / (1 + exp(x)))
} # checked 10/24/2018 



DDLislie_fit = function(Harvest,nage,period, nonspcharv = T, nonspcDD = T,ZeroFerc,initial,...){ # period is years of data (w/o baseline year)
	p = 3*nage + 1 + (!nonspcharv)*(nage-1) + nonspcDD * nage + 1 # length of pars
	if(missing(initial)){
	  initial = c(rnorm(p-1,0,.1),1/400)
	  initial[1:nage] = runif(nage)
	  initial[1:nage + 2*nage] = log( Harvest[,1] )
	}
	opt = optim(initial,fn = logLiklihood_fun,Harvdata = Harvest,nage=nage,period=period, nonspcharv=nonspcharv, ZeroFerc = ZeroFerc,...)
  opt$par[ZeroFerc]=0
	return(opt)
}