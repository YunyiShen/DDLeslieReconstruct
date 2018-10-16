getLislie = function(Survival,Fercest=T,Ferc=NULL){
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
  return(Tr)
}

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

GetHarvest = function(Harvpar,nage){
  Harvest = matrix(0,nage,nage)
  diag(Harvest) =exp( Harvpar)/(1+exp(Harvpar))
  return(Harvest)
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

