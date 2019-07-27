HarvestSen = function(Fec,Surv,SRB,Harvpar,nage,Harv_assump){
    source("DDLeslie.r")
    L = getLeslie(Surv,Fec,SRB) # intrinsic growth
    H = matrix(0,sum(nage),sum(nage))
    diag(H) = Harv_assump %*% Harvpar # propotional harvest
    
    EL = H %*% L # effactive growth
    
    ev = eigen(EL)
    lmax = which(Re(ev$values) == max(Re(ev$values)))
    lambda = Re(ev$values[lmax])
    W = ev$vectors
    w = abs(Re(W[, lmax]))
    V = Conj(solve(W))
    v = abs(Re(V[lmax, ]))
    E_Sen = v %o% w # sensitivity analysis of the effective growth
    
    # apply chain rule:
    H_Sen = rowSums( t(Harv_assump) %*% (E_Sen*L)) # sensitivity of harvest
    return(H_Sen)
} 

