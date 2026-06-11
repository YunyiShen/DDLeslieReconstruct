########
### ----------------------- HELPER FUNCTIONS ---------------------- ###
### --------------------------------------------------------------- ###

# get Leslie Matrix using specific parameters, can also return L-I if needed, first nage entries should be survival(no transiform, 0<S<1,transformation occor in sampler function) and next fertility rate, 


plotthings = function(YD_obj,pathsave="./figs/temp/age",nage,period,years,ppt=F){ # YD_obj should be a mcmc object, with vital rates in it
    if(ppt){require(export)}
    mean.harv = apply(YD_obj,2,mean)
    mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
    
    BI.low.harv = apply(YD_obj,2,quantile,probs = .025)
    BI.low.harv.matrix = matrix(BI.low.harv,nrow = nage,ncol = period)
    BI_harv_low = data.frame(age = 1:nage,BI.low.harv.matrix)
    
    
    BI.high.harv = apply(YD_obj,2,quantile,probs = .975)
    BI.high.harv.matrix = matrix(BI.high.harv,nrow = nage,ncol = period)
    BI_harv_high = data.frame(age = 1:nage,BI.high.harv.matrix)
    
    har_data = data.frame(matrix(nrow = 1,ncol = 5))
    colnames(har_data) = c("age","mean","low","high","time")
    har_data = har_data[-1,]
    
    for(i in 1:nage){
        temp = data.frame(age = i,mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = years)
        har_data = rbind(har_data,temp)
    }
    require(ggplot2)
    
    for(i in 1:nage){
        temp = data.frame(point = "model predict (95% CI)",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = years)
        write.csv(temp,paste0(pathsave,i,".csv"))
        filename = paste0(pathsave,i,".jpg")
        ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
            geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
            geom_point()
        if(ppt)    graph2ppt(file=sub("jpg","pptx",filename))
        else ggsave(filename, plot = last_plot())
    }
}


########
# Coming functions are
# modified from popReconstruct package

## ........... Misc Functions .......... ##
## ..................................... ##

logitf = function(p){
    log(p / (1 - p))
 } # checked 10/24/2018 

## eyes function in matlab
eyes = function(n){
        re = matrix(0,n,n)
        diag(re) = 1
        return(re)
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
} # checked 10/24/2018 

##--- Generates random draws from inverse gamma ---##
rinvGamma = function(n, shape, scale){
        return(1/(rgamma(n, shape = shape, rate = scale)))
} # checked 10/24/2018 # debug mode with print

##--- Returns value of inverse gamma pdf ---##
dinvGamma = function(x, shape, scale, log = FALSE){
        if(log) d =
                shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
        else d = scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
        return(d)
} # checked 10/24/2018

## ............. Likelihood ............ ##
## ..................................... ##

#log likelihood function of gaussian distributed
log.lhood =function(n.census, n.hat){
    ##-- value of log likelihoods --##
    density = dpois(x=n.census,lambda = n.hat,log = TRUE) # Use Poisson instead in Poisson_likelihood
    ##-- joint log likelihood --##

        return(sum(density))
} # checked 10/24/2018

## .............. Posterior ............ ##
## ..................................... ##    keep it, change in sampler function, but no migration here, should add estaK0

log.post = function(## estimated vitals
                    f, s, SRB,baseline.n, aK0, A, H # A for Aerial count detection probability
                    , estFer, estaK0
                    ## fixed prior means on vitals
                    , prior.mean.f, prior.mean.s, prior.mean.SRB
                    , prior.mean.b #, prior.mean.aK0
                    , prior.mean.A, prior.mean.H
                    ## fixed prior parameters on variance distns
                    , alpha.f, beta.f, alpha.s, beta.s, alpha.SRB, beta.SRB
                    , alpha.aK0, beta.aK0
                    , alpha.A, beta.A, alpha.H, beta.H
                    ## updated variances on prior distns
                    , sigmasq.f, sigmasq.s, sigmasq.SRB,sigmasq.n#, sigmasq.aK0
                    , sigmasq.A ,sigmasq.H
                    ## value of the log likelihood
                    , log.like
                    ## non zero rows of fertility matrix
                    , non.zero.fert){

        ##-- Values of prior densities for vitals --##

        ##.. Note that log densities are calculated for numerical stability.
        ##         f, baseline.n, prior.mean.f, prior.mean.b are logged coming
        ##         in, s, prior.mean.s is logit transformed coming in, g and
        ##         prior.mean.g are not transformed coming in.
        ##-- prior for f and baseline K0 if needed to be estimatedd --##
        
        if(estFer){
            log.f.prior = dnorm(as.vector(f[non.zero.fert,])
                                , mean = as.vector(prior.mean.f[non.zero.fert,])
                                , sd = sqrt(sigmasq.f)
                                , log = TRUE)
            log.sigmasq.f.prior =
                log(dinvGamma(sigmasq.f, alpha.f, beta.f))
        }
        else {
            log.f.prior = 0
            log.sigmasq.f.prior = 0
        }
        
        if(estaK0){
            log.aK0.prior =#sum( dunif(aK0[[1]],alpha.aK0[[1]],beta.aK0[[1]],T) , dunif(aK0[[2]],alpha.aK0[[2]],beta.aK0[[2]],T))
                Reduce(sum, 
                       lapply(1:length(aK0),
                              function(i,aK0,al.aK0,be.aK0){
                                  dunif(aK0[[i]],al.aK0[[i]],be.aK0[[i]],T)}
                              ,aK0,alpha.aK0,beta.aK0))
            
        }
        else {
            log.aK0.prior = 0
            log.sigmasq.aK0.prior = 0
        }

        ##-- prior for s and Sex Ratio at Birth (SRB), Aerial count detection rate, and hunting rate H --##
        log.s.prior = dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                            ,log = TRUE)
        log.SRB.prior = dnorm(SRB, mean = prior.mean.SRB, sd = sqrt(sigmasq.SRB)
                              ,log = TRUE)
        log.H.prior = dnorm(H, mean = prior.mean.H
                            ,sd = sqrt(sigmasq.H)
                            ,log = TRUE)
        log.A.prior = dnorm(A, mean = prior.mean.A
                            ,sd = sqrt(sigmasq.A)
                            ,log = TRUE)
        log.sigmasq.s.prior =
                log(dinvGamma(sigmasq.s, alpha.s, beta.s))
        log.sigmasq.SRB.prior =
                log(dinvGamma(sigmasq.SRB, alpha.SRB, beta.SRB))
        log.sigmasq.A.prior =
                log(dinvGamma(sigmasq.A, alpha.A, beta.A))
        log.sigmasq.H.prior =
                log(dinvGamma(sigmasq.H, alpha.H, beta.H))
        
        ##-- The log posterior is the SUM of these with the log.like --##

    return(sum(log.f.prior, log.s.prior, log.SRB.prior#, log.b.prior
    , log.aK0.prior, log.H.prior,log.A.prior,
                             log.sigmasq.f.prior
                             ,log.sigmasq.s.prior
                             ,log.sigmasq.SRB.prior                             
                             #,log.sigmasq.aK0.prior
                             ,log.sigmasq.H.prior
                             ,log.sigmasq.A.prior
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

## sensitivity analysis
HarvestSen = function(Fec,Surv,SRB,Harvpar,nage,Harv_assump){
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