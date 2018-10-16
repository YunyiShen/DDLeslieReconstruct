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
  Tr[nage,nage] = Survival[nage]
  return(Tr)
}

DensityDependcy = function(global = F, Xn, E0, K0){
  if(global){
    D = 1-sum(Xn)/K0
  }
  else{
    nage = length(Xn)
    D = matrix(0,ncol = nage,nrow = nage)
    diag(D) = 1-Xn/(K0 * E0)
  }
  return(D)
}

GetHarvest = function(Harvpar,nage){
  Harvest = matrix(0,nage,nage)
  diag(Harvest) =exp( Harvpar)/(1+exp(Harvpar))
  return(Harvest)
}

sqrResidue = function(pars,Fercest=T,Ferc=NULL, estK0 = T, E0, K0 = NULL, global = F,data){
  period = ncol(data)
  nage = length(E0)
  if(estK0){
    K0 = pars[length(pars)]
    pars = pars[-length(pars)]
  }
  Survival = pars[1:(nage * (1 + Fercest))]
  Harvpar = pars[-1:(nage * (1 + Fercest))]
  Lislie = getLislie(Survival,Fercest=T,Ferc=NULL)
  H = GetHarvest(Harvpar,nage)
  Harest = apply(data,2,
                 function(Xn,
                          Lislie,
                          H,
                          global, 
                          E0, 
                          K0){
    Popu_before_harvest = solve(H,Xn)
    eyes = diag(rep(1,length(Xn)))
    D = DensityDependcy(global = global, Xn=Xn, E0=E0, K0=K0)
    Popu_after = (eyes-H)Popu_before_harvest
    Xn1 = H %*% (D %*% Popu_after)
    return(Xn1)
  },global = global, Xn=Xn, E0=E0, K0=K0)
  Residue = sum((data[,-1]-Harest)^2)
  return(Residue)
}

