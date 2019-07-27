HarvestSen = function(Fec,Surv,SRB,Harvpar,nage,Harv_assump){
    source("DDLeslie.r")
    L = getLeslie(Surv,Fec,SRB) # intrinsic growth
    H = matrix(0,sum(nage),sum(nage))
    diag(H) = Harv_assump %*% Harvpar # propotional harvest
    
    EL = H %*% L # effactive growth
    
    left_eig = eigen(t(EL))$vectors
    right_eig = eigen(EL)$vectors
    E_Sen = (left_eig[,1]%*%t(right_eig[,1]))/
      as.numeric(t(left_eig[,1])%*%EL%*%right_eig[,1]) # sensitivity analysis of the effective growth
    # apply chain rule:
    H_Sen = rowSums( t(Harv_assump) %*% (E_Sen*L))
    return(H_Sen)
} 

