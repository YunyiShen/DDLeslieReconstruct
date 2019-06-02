########
### ----------------------- HELPER FUNCTIONS ---------------------- ###
### --------------------------------------------------------------- ###

# get Leslie Matrix using specific parameters, can also return L-I if needed, first nage entries should be survival(no transiform, 0<S<1,transformation occor in sampler function) and next fertility rate, 

getLeslie = function(Survival,Fec,SRB,minus1=F){ # SRB, sax ratio at birth
  nage_female = length(Fec)
  nage_male = length(Survival)-length(Fec)
  Tr = matrix(0,nage_female+nage_male,nage_female+nage_male)
  Tr[1,1:nage_female] = Fec*SRB
  Tr[1+nage_female,1:nage_female] = Fec * (1-SRB)
  diag(Tr[2:(nage_female),1:(nage_female-1)]) = Survival[2:nage_female - 1]
  Tr[nage_female,nage_female] = Survival[nage_female]
  diag(Tr[2:(nage_male) + nage_female,1:(nage_male-1) + nage_female]) = Survival[2:nage_male - 1 + nage_female]
  Tr[nage_male+nage_female,nage_male+nage_female] = Survival[nage_male+nage_female]
  if(minus1){
    diag(Tr)=diag(Tr)-1
  }
  return(Tr)
} #checked 10/24/2018

# Linear density dependency matrix, can be global of age specific, or null = T for no density dependency. in estimaton process, if E0 is not fixed, it should be estimated as kk = solve(H, data[,1]) E0 = kk/sum(kk), by assuming that first year is at equalibrium
DensityDependcy = function(global = F, Xn, E0, aK0, null = F){
  E0 = E0/sum(E0)
  nage = length(Xn)
  #D = matrix(0,ncol = nage,nrow = nage)
  if(global){
    den = as.matrix( 1 + (aK0) * sum(Xn) )
    D = (1-null)*den + null
  }
  else{
    D = 1-(1-null)* aK0 * Xn
  }
  return(as.numeric(D))
} #checked 10/24/2018


# Project the (density dependent) harvest-after-reproducing model from year i to i+1, given harvest # and return harvest # of year i+1 
ProjectHarvest_helper = function(data_n, Leslie, H_n,H_np1, global, E0, aK0, null = F,nage_female){
    nage_tol = ncol(Leslie)
    nage_male = nage_tol-nage_female
  #I = matrix(0,length(data_n),length(data_n))    
  #H_n_temp = as.numeric(H_n)
  #H_np1_temp = as.numeric(H_np1)
  #H_n[2:nage_female]=H_n_temp[2] # for doe Y+A
  #H_n[c(1,nage_female + 1)] = H_n_temp[1] # for fawns
  #H_n[2:nage_male + nage_female] = H_n_temp[3] # for buck Y+A
  #H_np1[2:nage_female]=H_np1_temp[2]
  #H_np1[c(1,nage_female + 1)] = H_np1_temp[1]
  #H_np1[2:nage_male + nage_female] = H_np1_temp[3]
  #H_logit = logitf(H)
  
  
    #diag(I) = 1
    X_n1 = (1-H_n) * (data_n/H_n)
    #DDH = (1+aK0[3]*sum(X_n1))*H_logit
    #H = invlogit(DDH)
    D_bir = DensityDependcy(global = global, Xn=X_n1, E0=E0, aK0=aK0[1], null = null)
    D_dea = DensityDependcy(global = global, Xn=X_n1, E0=E0, aK0=aK0[2], null = null)
    Leslie[c(1,nage_female+1),] =  D_bir *(Leslie[c(1,nage_female+1),])
    Leslie[-c(1,nage_female+1),] = D_dea * Leslie[-c(1,nage_female+1),]
    #Leslie = as.numeric(1+t(aK0[1:nage,]) %*% ( X_n1)) * Leslie
  #ata_n1 = H * (Leslie %*% (D_bir * X_n1) - D_dea * X_n1 +X_n1)
    data_n1 = H_np1 * (Leslie %*% ( X_n1))
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  } # checked 10/24/2018

# Project harvest model from a initial harvest, survival is col vector with all survival rate of all age class, this nrow(Survival)=nage
ProjectHarvest_homo = function(Survival, Harvpar,Fec, SRB,E0=NULL, aK0 = NULL, global = F, null = F,bl , period, nage){
  
  
  Leslie = getLeslie(Survival,Fec=Fec,SRB = SRB,minus1=F)
  # H = GetHarvest(Harvpar,nage)
  Harvest = matrix(0,nage,period + 1)
  Harvest[,1] = bl 
  if(length(E0)==0){
    E0 = (bl/as.numeric( Harvpar))
    E0 = E0/(sum(E0))
  }
  else E0 = E0/(sum(E0))
  for(i in 1 : period + 1){
    Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Leslie = Leslie, E0=E0, aK0=aK0, H_n=Harvpar,H_np1 = Harvpar,null = null,nage_female = nage[1]) # project harvest from very initial
  }              
  return(Harvest)
} #checked 10/24/2018

# project when time inhomo survival and ferc, survival should be a matrix with nage x period entries, so do ferc
ProjectHarvest_inhomo = function(Survival, Harvpar,Fec, SRB,E0=NULL, aK0 = NULL, global = F, null = F,bl , period, nage){
  Harvest = matrix(0,sum(nage),period+1) # nage is a vector with female first
  Harvest[,1] = bl
  # H = GetHarvest(Harvpar,nage)
  Harvpar = apply(Harvpar,2,getfullHarvpar,nage = nage)
  if(length(E0)==0){
    E0 = (bl/( Harvpar[,1]))
    E0 = E0/(sum(E0))
  }
  else E0 = E0/(sum(E0))
  for(i in 1 : period + 1){
    Leslie = getLeslie(Survival[,i-1],Fec[,i-1],SRB[i-1],minus1=F) # inhomo, Survival rows are age structure, cols are time, SRB has only one row or col, for 
    Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Leslie = Leslie, E0=E0, aK0=aK0, H_n=Harvpar[,i-1],H_np1 = Harvpar[,i],null = null,nage_female = nage[1])
  }
  return(Harvest)
} # checked 10/24/2018

getfullHarvpar = function(H_n,nage){
  H_n_temp = as.numeric(H_n)
  H_n[2:nage[1]]=H_n_temp[2] # for doe Y+A
  H_n[1] = H_n_temp[1] # for female fawns
  H_n[nage[1]+1] = H_n_temp[3] # male fawns
  H_n[2:nage[2] + nage[1]] = H_n_temp[4] # for buck Y+A
  return(H_n)


}


# Given harvest data, calculate the living individual in certain year
getLivingIdividuals = function(H,data){

  IminusH = 1-H
  LivingIdividuals = apply(data,2,function(ww,H,IminusH){IminusH*(ww/H)},H=H,IminusH=IminusH)
  return(LivingIdividuals)
} # checked 10/24/2018

getLivingIdividuals_preharv = function(H,data){

  #IminusH = 1-H
  LivingIdividuals = apply(data,2,function(ww,H){(ww/H)},H=H)
  return(LivingIdividuals)
} # checked 10/24/2018

getLivingIdividuals_inhomo = function(Harv,H,nage){
    #nage = .5 * nrow(Harv)
    H = apply(H,2,getfullHarvpar,nage)
    Harv*(1/H-1)
}

getAerialCount = function(Harv,H,A,nage){
    Xt = getLivingIdividuals_inhomo(Harv,H,nage)
    Ae = colSums(Xt) * A # aerial counts
    return(Ae)
}

plotthings = function(YD_obj,pathsave="./figs/temp/age",nage,period,years){
  mean.harv = apply(YD_obj,2,mean)
  mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)
  #harv_mean = data.frame(age = 1:8,mean.harv.matrix)
  #mean.total.harv = apply(mean.harv.matrix,2,sum)
  #plot(mean.total.harv)
  
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
    #temp2 = data.frame(point = "data",mean = t(Harv.data[i,2:12]),low =t( Harv.data[i,2:12]),high = t(Harv.data[i,2:12]),time = 1996:2006)
    #colnames(temp2) = colnames(temp1)
    #temp = rbind(temp1)
    write.csv(temp,paste0(pathsave,i,".csv"))
    #rm(temp1)
    #rm(temp2)
    filename = paste0(pathsave,i,".jpg")
    #jpeg(filename)
    ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) + 
      geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
      #geom_line() +
      geom_point()
    #Sys.sleep(5)
    ggsave(filename, plot = last_plot())
    #dev.off()
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
  #print(shape) # debug mode
  #print(scale)
    return(1/(rgamma(n, shape = shape, rate = scale)))
} # checked 10/24/2018 # debug mode with print

##--- Returns value of inverse gamma pdf ---##
dinvGamma = function(x, shape, scale, log = FALSE){
    if(log) d <-
        shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
    else d <- scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
    return(d)
} # checked 10/24/2018

## ............. Likelihood ............ ##
## ..................................... ##

#log likelihood function of gaussian distributed
log.lhood =function(log.n.census, log.n.hat, ll.var){
    ##.. log.n.census and log.n.hat should already be logged

    ##.. log.n.hat should be log projected counts for census years
    ##   interpolated if necessary


    ##-- value of log likelihoods --##

    density = dnorm(log.n.hat,
                     mean = log.n.census,
                     sd = sqrt(ll.var),
                     log = TRUE
                     )
    ##-- joint log likelihood --##

    return(sum(density))
} # checked 10/24/2018

## .............. Posterior ............ ##
## ..................................... ##  keep it, change in sampler function, but no migration here, should add estaK0

log.post = function(## estimated vitals
                           f, s, SRB,baseline.n, aK0, A, H # A for Aerial count detection probability
                                       ,estFer, estaK0
                           ## fixed prior means on vitals
                           ,prior.mean.f, prior.mean.s, prior.mean.SRB
                           ,prior.mean.b, prior.mean.aK0
                           ,prior.mean.A, prior.mean.H
                           ## fixed prior parameters on variance distns
                           ,alpha.f, beta.f, alpha.s, beta.s, alpha.SRB, beta.SRB
                           ,alpha.n, beta.n, alpha.aK0, beta.aK0
                            ,alpha.A, beta.A, alpha.H, beta.H, alpha.ae, beta.ae
                           ## updated variances on prior distns
                           ,sigmasq.f, sigmasq.s, sigmasq.SRB,sigmasq.n, sigmasq.aK0
                            ,sigmasq.A ,sigmasq.H,sigmasq.ae
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
    
    if(estaK0){
      log.aK0.prior = dnorm(aK0, mean = prior.mean.aK0
                           , sd = sqrt(sigmasq.aK0)
                           , log = T )
      log.sigmasq.aK0.prior =
        log(dinvGamma(sigmasq.aK0, alpha.aK0, beta.aK0))
    }
    else {
      log.aK0.prior = 0
      log.sigmasq.aK0.prior = 0
    }

    ##-- prior for s and baseline n, and hunting rate H --##
    log.s.prior = dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                         ,log = TRUE)
    log.SRB.prior = dnorm(SRB, mean = prior.mean.SRB, sd = sqrt(sigmasq.SRB)
                         ,log = TRUE)
    log.b.prior = dnorm(baseline.n, mean = prior.mean.b
                         ,sd = sqrt(sigmasq.n)
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
    log.sigmasq.n.prior =
        log(dinvGamma(sigmasq.n, alpha.n, beta.n))
    log.sigmasq.A.prior =
        log(dinvGamma(sigmasq.A, alpha.A, beta.A))
    log.sigmasq.H.prior =
        log(dinvGamma(sigmasq.H, alpha.H, beta.H))
    log.sigmasq.ae.prior =
        log(dinvGamma(sigmasq.ae, alpha.ae, beta.ae))
    
    ##-- The log posterior is the SUM of these with the log.like --##

  return(sum(log.f.prior, log.s.prior, log.SRB.prior, log.b.prior, log.aK0.prior, log.H.prior,log.A.prior,
               log.sigmasq.f.prior
               ,log.sigmasq.s.prior
               ,log.sigmasq.SRB.prior
               ,log.sigmasq.n.prior
               ,log.sigmasq.aK0.prior
               ,log.sigmasq.H.prior
               ,log.sigmasq.A.prior
               ,log.sigmasq.ae.prior
               ,log.like))

} # not checked yet, assume to be right

## ......... Acceptance Ratio .......... ##
## ..................................... ##
acc.ra = function(log.prop, log.current){
    min(1, exp(log.prop - log.current))
} # assume to be right

acc.ra.var = function(log.prop.post, log.curr.post, log.prop.var, log.curr.var){
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post))
} # assume to be right


### --------------------------- SAMPLER --------------------------- ###
### --------------------------------------------------------------- ### hardest part is coming

HDDLislie.sampler <-
    function(#.. number of iterations and burn-in (not saved)
             n.iter, burn.in = 0, thin.by = 1

             #.. fixed variance hyper-parameters
             ,al.f = 1, be.f = 0.0109, al.s = 1, be.s = 0.0109,al.SRB = 1,be.SRB = 0.0109
             , al.aK0 = 1, be.aK0 = 0.0436
             , al.n = 1, be.n = 0.0109, al.ae = 1, be.ae = 0.0109
             , al.H = 1, be.H = 0.0436, al.A = 1, be.A = 0.0436

             #.. fixed prior means
             ,mean.f, mean.s, mean.SRB, mean.b, mean.aK0, mean.H, mean.A

             #.. inital values for vitals and variances
             #   *vitals not transformed coming in* all not transfer, will 
             ,start.f = mean.f, start.s = mean.s, start.SRB = mean.SRB
             ,start.b = mean.b, start.aK0 = mean.aK0, start.H = mean.H
             ,start.A = mean.A
             ,start.sigmasq.f = 5, start.sigmasq.s = 5, start.sigmasq.SRB = 5
             ,start.sigmasq.n = 5,start.sigmasq.ae = 5, start.sigmasq.aK0 = 5, start.sigmasq.H = 5
             ,start.sigmasq.A = 5

             #.. census data
             #   *not transformed coming in*
             ,Harv.data
             ,Aerial.data
             #.. **variances** for proposal distributions used in M-H
             #   steps which update vital rates.
             ,prop.vars # col names should be "fert.rate", 
                   # "fert.rate", "surv.prop", "SRB","H", "A","aK0"
                   # ,"baseline.count"

             #.. number of periods to project forward over (e.g.,
             #     number of five-year steps)
             ,proj.periods = (ncol(Harv.data)-1)

             #.. age group width
             ,nage = 8 
             
             ,estFer=T, Fec=rep(1,nage), estaK0 = T
             ,E0=NULL , aK0 = 0, global = T, null = T # control parameters for the model, global is whether density dependency is global rather than age specific, null is whether exist density dependency.
             ,homo = F # whether assume time homogeneous 

             #.. print algorithm progress
             ,verb = FALSE

             #.. tolerance defining allowable survival probabilities
             ,s.tol = 10^(-10)
             )
{
    require(coda)
    ## .............. Sampler .............. ##
    ## ..................................... ##

    mean.aK0 = as.matrix(mean.aK0)
    ## -------- Begin timing ------- ##

    ptm <- proc.time()


    ## ------ Match functions ------ ##


    ## ------- Check input dimensions ------- ##
    ## from popReconstruct, not needed yet.
    # input.dims <- sapply(list(start.s = start.s[1:(nrow(start.s)-1),]
                              # , start.g = start.g
             # ,mean.f = mean.f, mean.s = mean.s[1:(nrow(mean.s)-1),]
             # ,mean.g = mean.g
             # )
           # ,"dim")
    # mismatch.dims <- apply(input.dims, 2, "identical", dim(start.f))
    # if(!all(mismatch.dims))
        # stop("Dims of these inputs do not match 'dim(start.f)'", "\n"
             # ,paste(names(mismatch.dims)[!mismatch.dims], collapse = "  "))


    # ## ------- Check Years --------- ##

    # all.vr.years <-
        # list(start.s = colnames(start.s), start.g = colnames(start.g)
             # ,mean.f = colnames(mean.f), mean.s = colnames(mean.s)
             # ,mean.g = colnames(mean.g)
             # )
    # mismatch.yrs <- sapply(all.vr.years, "identical", colnames(start.f))
    # if(!all(mismatch.dims))
        # stop("Years of these inputs do not match years of 'start.f'", "\n"
             # ,paste(names(mismatch.yrs)[!mismatch.yrs], collapse = "  "))

    # all.vr.years.eq <-
        # sapply(all.vr.years, FUN = function(z) all.equal(colnames(start.f), z))
    # if(!all(all.vr.years.eq)) {
        # stop(paste("colnames(", names(all.vr.years)[min(which(!all.vr.years.eq))], ")"
                   # ," != colnames(start.f). There may be more...", sep = "")
             # )
    # }

    # proj.years <-
        # seq(from = as.numeric(colnames(start.f)[1]), by = age.size
            # ,length = proj.periods + 1)

    # if(!all.equal(proj.years[1:ncol(start.f)], as.numeric(colnames(start.f))))
        # stop("colnames(start.f) !=  seq(from = as.numeric(colnames(start.f)[1]), by = age.size, length = ncol(start.f))")

    # vr.years <- as.numeric(colnames(start.f))
    # baseline.year <- as.numeric(colnames(start.b))
    # census.years <- as.numeric(colnames(pop.data))


    # ## ------- Determine fert.rows --------- ##

    zero.elements <- mean.f == 0
    fert.rows <- as.logical(apply(zero.elements, 1, function(z) !all(z)))

    # ## ------- Type of migration data ------- ##

    # ## No reason for this anymore; remove later
    # mig.string <- "prop"


    ## ---------- Storage ---------- ##

    #.. MCMC objects for posterior samples
    # Samples are stored as 2D arrays for compatibility with coda's
    # mcmc format with iterations as rows, year*age.group as columns.
    # Age.group cycles fastest across columns, e.g.,
    # _____________________________________________________
    #   1960  | 1960  | 1960  | ... | 1965  | 1965  | ...
    #   15.19 | 20.24 | 25.29 | ... | 15.19 | 20.24 | ...
    # 1  --   |  --   |  --   | ... |  --   |  --   | ...
    # 2  --   |  --   |  --   | ... |  --   |  --   | ...
    #   etc.
    # _____________________________________________________

    ## How many (samples) stored?
    cat("Allocating for RAMs...\n")
    n.stored = ceiling(n.iter / thin.by)
    ntimes = (!homo) * ncol(start.s) + (homo) # whether assume time homogeneous of survival etc.
      # Fertility
        
      if(estFer){
      fert.rate.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = sum(fert.rows) * ntimes)
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(fert.rate.mcmc) = NULL
          
    } 
    else{fert.rate.mcmc = NULL}
      # Survival proportions
      surv.prop.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.s) * ntimes)
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(surv.prop.mcmc) = NULL
       # Sex Ratio at Birth
      SRB.mcmc <-
          mcmc(matrix(nrow = n.stored
                      ,ncol = ntimes)
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(surv.prop.mcmc) = NULL

      # lx
      # this is for current population, for us, it is current culling data
      lx.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.b) * (proj.periods))
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(lx.mcmc) = NULL
      
      ae.mcmc = 
          mcmc(matrix(nrow = n.stored
                      ,ncol = proj.periods+1)
                      ,start = burn.in + 1
                      ,thin = thin.by)      # this is for aerial count
      
      # carrying capacity assumed to be time homogeneous
      if(estaK0){
        aK0.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = length(mean.aK0))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(aK0.mcmc) = NULL
      } 
      else{
        aK0.mcmc = NULL
      }

      # Harvest proportion, can be either time homo or not
      H.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.H) * (ntimes+!homo))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(H.mcmc) = NULL
      # Aerial counts
      A.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.A) * (ntimes+!homo))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(H.mcmc) = NULL
      
      # baseline counts
      baseline.count.mcmc =
          mcmc(matrix(nrow = n.stored, ncol = nrow(start.b))
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(baseline.count.mcmc) = NULL

      # variances
      variances.mcmc =
          mcmc(matrix(nrow = n.stored, ncol = 8)
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(variances.mcmc) =
          c("fert.rate.var", "surv.prop.var", "SRB.var","H.var", "A.var","aK0.var"
            ,"population.count.var","aerial.count.var") # may need check here since K0 and fert can be known (though I prefer estimating it)

    #.. Record acceptance rate

    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(as.matrix( mean.f[fert.rows,]))
             ,ncol = ncol(as.matrix( mean.f[fert.rows,]))
             ,dimnames = dimnames(( mean.f[fert.rows,]))
             )
             ,surv.prop = matrix(0, nrow = nrow(as.matrix( mean.s))
              ,ncol = ncol(as.matrix( mean.s))
              ,dimnames = dimnames(mean.s)
              )
             ,SRB = matrix(0, nrow = nrow(as.matrix( mean.SRB))
              ,ncol = ncol(as.matrix( mean.SRB))
              ,dimnames = dimnames(mean.SRB))
             ,A = matrix(0, nrow = nrow(as.matrix(mean.A)), ncol = ncol(as.matrix( mean.A))
              ,dimnames = dimnames(mean.A)
              )
             ,H = matrix(0, nrow = nrow(as.matrix(mean.H)), ncol = ncol(as.matrix( mean.H))
              ,dimnames = dimnames(mean.H)
              )
             ,aK0 = matrix(0, nrow = nrow(as.matrix( mean.aK0)), ncol = ncol( as.matrix(mean.aK0))
              ,dimnames = dimnames(mean.aK0)
              )
             ,baseline.count = matrix(0, nrow = nrow(as.matrix( mean.b))
              ,dimnames = dimnames( mean.b)
              )
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.SRB = 0
             ,sigmasq.A = 0
             ,sigmasq.H = 0
             ,sigmasq.aK0 = 0
             ,sigmasq.n = 0
             ,sigmasq.ae = 0 # aerial counts
             )


    #.. Count how often acceptance ratio missing or na

    ar.na <- acc.count


    #.. Count how often projection gives negative population

    pop.negative <-
        list(fert.rate = matrix(0, nrow = nrow(as.matrix( mean.f[fert.rows,]))
             ,ncol = ncol(as.matrix( mean.f[fert.rows,]))
             ,dimnames = dimnames(( mean.f[fert.rows,]))
             )
             ,surv.prop = matrix(0, nrow = nrow(as.matrix( mean.s))
              ,ncol = ncol(as.matrix( mean.s))
              ,dimnames = dimnames(mean.s)
              )
             ,SRB = matrix(0, nrow = nrow(as.matrix( mean.SRB))
              ,ncol = ncol(as.matrix( mean.SRB)))
              ,dimnames = dimnames(mean.SRB)
             ,A = matrix(0, nrow = nrow(as.matrix(mean.A)), ncol = ncol(as.matrix( mean.A))
              ,dimnames = dimnames(mean.A)
              )
             ,H = matrix(0, nrow = nrow(as.matrix(mean.H)), ncol = ncol(as.matrix( mean.H))
              ,dimnames = dimnames(mean.H)
              )
             ,aK0 = matrix(0, nrow = nrow(as.matrix( mean.aK0)), ncol = ncol( as.matrix(mean.aK0))
              ,dimnames = dimnames(mean.aK0)
              )
             ,baseline.count = matrix(0, nrow = nrow(as.matrix( mean.b))
              ,dimnames = dimnames( mean.b)
              )
             )


    #.. Count how often surv probs are outside tolerance

    s.out.tol <- matrix(0, nrow = nrow(mean.s), ncol = ncol(mean.s)
                        ,dimnames = dimnames(mean.s))


    ## -------- Initialize -------- ## Restart here in 10/19/2018
    cat("Initializing...")
    #.. Set current vitals and variances to inital values
    #   Take logs/logits here where required
    if(estFer){
      log.curr.f = log(start.f)  #<-- log(0) stored as "-Inf". Gets
      log.prop.f = log(start.f)  #    converted to 0 under exponentiation
      
    }
    else{
      log.curr.f =  (!estFer)*log(start.f) #<-- log(0) stored as "-Inf". Gets
      log.prop.f =  (!estFer)*log(start.f) #    converted to 0 under exponentiation
    }
    logit.curr.s = logitf(start.s)
    logit.curr.SRB = logitf(start.SRB)
    logit.curr.H = logitf(start.H)
    logit.curr.A = logitf(start.A)
    curr.aK0=(start.aK0)
    #curr.aK0=(aK0)}
    #curr.aK0 = estaK0 * log(start.aK0) + (!estaK0) * log(K0)
    log.curr.b = log(start.b)

    curr.sigmasq.f = start.sigmasq.f
    curr.sigmasq.s = start.sigmasq.s
    curr.sigmasq.SRB = start.sigmasq.SRB
    curr.sigmasq.A = start.sigmasq.A
    curr.sigmasq.H = start.sigmasq.H
    curr.sigmasq.aK0 = start.sigmasq.aK0
    curr.sigmasq.n = start.sigmasq.n
    curr.sigmasq.ae = start.sigmasq.ae


    #.. Fixed means for vitals and baseline
    #   Set these to inputs, take logs where required.

    log.mean.f = log(mean.f)
    logit.mean.s = logitf(mean.s)
    logit.mean.SRB = logitf(mean.SRB)
    logit.mean.A = logitf(mean.A)
    logit.mean.H = logitf(mean.H)
    log.mean.b = log(mean.b)


    #.. Fixed Harvest data
    #   Take logs here

    log.Harv.mat = log(Harv.data)
    log.Aeri.mat = log(Aerial.data)


    #.. Set current projection: base on initial values # homo or not is important, determin it use homo = T
    if(homo){
      log.curr.proj =
        log(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB),E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
    }
    else{
      log.curr.proj =
        log(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB),E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
    }
    log.curr.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))
## stop here 10/19/2018   
## restart here 10/22/2018

    #.. Current log posterior
    ## shit so many parameters to pass here... really hard to read...
    
# checked here,  solved, -Inf return 10/25/2018 14:31, Done, because of test data has zero likelihood (changed parameter)
    log.curr.posterior =
        log.post(f = log.curr.f
                       ,s = logit.curr.s
                       ,SRB = logit.curr.SRB
                       ,A = logit.curr.A
                       ,H = logit.curr.H
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat
                                ,log.n.hat = log.curr.proj,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri
                                ,ll.var = curr.sigmasq.ae) #　both harvest and aerial count, potentially reconstruction?
                       ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                       )

# stop checking here 10/25/2018
    ## -------- Begin loop ------- ##
    #...............................#

    if(verb) {
        cat("\n\ntotal iterations = ", n.iter+burn.in
            ,"\nburn in = ", burn.in
            ,"\nthin = ", thin.by
            ,",\nnumber stored = ", n.stored, sep = "")
        cat("\n\nfert.rows = ", which(fert.rows)
            )
        cat("\n\n"
            ,"iter ", " quantity\n", "---- ", " --------"
            ,sep = "")
    }
    cat("done\n")
    cat("Start sampling...\n")
    for(i in 1:(n.iter + burn.in)) {
       #if(i %% 100==0) cat("Iteration #",i,"out of",n.iter + burn.in,"total\n")
      svMisc::progress(((i-1)/(n.iter + burn.in))*100,progress.bar = T)
      # k is the index into the storage objects
      k <- (i - burn.in - 1) / thin.by + 1


      ## -------- Vital Rate M-H Steps ------- ##
    if(estFer){
      ##...... Fertility .....##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Fertility")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(log.curr.f[fert.rows,])) {

        #.. make a matrix conformable w fertility rate matrix
        log.prop.f.mat <-
            matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
        log.prop.f.mat[fert.rows,][j] <-
            rnorm(1, 0, sqrt(prop.vars$fert.rate[j])) #pop vars col names 

        #.. make proposal
        log.prop.f <- log.curr.f + log.prop.f.mat
        # - Run CCMP (project on the original scale)
        #   ** Don't allow negative population
        if(homo){
            full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.prop.f) #<-- use proposal
                , SRB = invlogit(logit.curr.SRB)
                , E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * as.numeric( invlogit(logit.curr.H)) , period = proj.periods, nage = nage))
        }
        else{
            full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.prop.f)#<-- use proposal
                , SRB = invlogit(logit.curr.SRB)
                , E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
        }


        if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
           || is.nan(sum(full.proj))) {
            if(i > burn.in) {
                pop.negative$fert.rate[j] =
                    pop.negative$fert.rate[j] + 1/n.iter
            }
        } else {
            log.prop.proj = log(full.proj)
            log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))
          # - Calculate log posterior of proposed vital under projection          
          ## shit again so many parameters to pass here... really hard to read...

          log.prop.posterior =
              log.post(f = log.prop.f #<-- use proposal
                       ,s = logit.curr.s
                       ,SRB = logit.curr.SRB
                       ,A = logit.curr.A
                       ,H = logit.curr.H
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )
          #- Acceptance ratio
          ar <- acc.ra(log.prop = log.prop.posterior,
                             log.current = log.curr.posterior)

          # - Move or stay
          #.. stay if acceptance ratio 0, missing, infinity, etc.
          if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$fert.rate[j] <-
                ar.na$fert.rate[j] + 1/n.iter
          } else {
            #.. if accept, update current fert rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
              if(i > burn.in) acc.count$fert.rate[j] <-
                  acc.count$fert.rate[j] + 1/n.iter
              log.curr.f <- log.prop.f
              log.curr.proj <- log.prop.proj
              log.curr.aeri <- log.prop.aeri
              log.curr.posterior <- log.prop.posterior # change log curr
            }
            #.. if reject, leave current fert rates and projections
            #   alone

          } # close else after checking for ar=0, missing, inf

        } # close else after checking neg or zero population

      } # close loop over all age-spec fertility rates

      #.. Store proposed fertility rate matrix
      if(k %% 1 == 0 && k > 0) fert.rate.mcmc[k,] <-
          as.vector(exp(log.curr.f[fert.rows,]))
    }
    # pause 0519
      ##...... Survival ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Survival")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.s)) {

        #.. make a matrix conformable w rate matrix
        ## TEST 10/26/2018, propble occur, logit.prop.s.mat is almost all 0
        ##   this is a strange structure inherite from popReconstruct 
        ##   the logit.prop.s.mat renew every loop, do not know why they do this.
        logit.prop.s.mat <-
            matrix(0, nrow = nrow(logit.curr.s)
                   ,ncol = ncol(logit.curr.s)) # this result depends on whether time-homo assumed.
        logit.prop.s.mat[j] <- rnorm(1, 0, sqrt(prop.vars$surv.prop[j]))

        #.. make proposal
        logit.prop.s <- logit.curr.s + logit.prop.s.mat

        #.. If proposal resulted in back-transformed s = 0 or 1, do
        #   nothing
        if(invlogit(logit.prop.s[j]) > 1 - s.tol ||
           invlogit(logit.prop.s[j]) < s.tol) {
          #.. leave current surv rates and projections
          #   alone (simply do not propose
          #   extreme survival probabilities)
          s.out.tol[j] <- s.out.tol[j] + 1/n.iter
        } else {

          # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
            if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.prop.s) #<-- use proposal
                , Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.prop.s)#<-- use proposal
                , Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
                log.prop.proj = log(full.proj)
                log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

            # - Calculate log posterior of proposed vital under projection
            log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.prop.s #<-- use proposal
                       ,SRB = logit.curr.SRB
                       ,A = logit.curr.A
                       ,H = logit.curr.H
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )

            #- Acceptance ratio
            ar <- acc.ra(log.prop = log.prop.posterior,
                              log.current = log.curr.posterior)

            # - Move or stay
            #.. stay if acceptance ratio 0, missing, infinity, etc.
            if(is.na(ar) || is.nan(ar) || ar < 0) {
              if(i > burn.in) ar.na$surv.prop[j] <-
                  ar.na$surv.prop[j] + 1/n.iter
            } else {
              #.. if accept, update current surv rates,
              #   update current projection and count acceptance
              if(runif(1) <= ar) {
                if(i > burn.in) acc.count$surv.prop[j] <-
                    acc.count$surv.prop[j] + 1/n.iter
                logit.curr.s <- logit.prop.s
                log.curr.proj <- log.prop.proj
                log.curr.aeri <- log.prop.aeri
                log.curr.posterior <- log.prop.posterior
              } 

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(k %% 1 == 0 && k > 0) surv.prop.mcmc[k,] <-
        as.vector(invlogit(logit.curr.s))
        

        
        
        
          ##...... SRB ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " SRB")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.SRB)) {

        #.. make a matrix conformable w rate matrix
        ## TEST 10/26/2018, propble occur, logit.prop.s.mat is almost all 0
        ##   this is a strange structure inherite from popReconstruct 
        ##   the logit.prop.s.mat renew every loop, do not know why they do this.
        logit.prop.SRB.mat <-
            matrix(0, nrow = nrow(logit.curr.SRB)
                   ,ncol = ncol(logit.curr.SRB)) # this result depends on whether time-homo assumed.
        logit.prop.SRB.mat[j] <- rnorm(1, 0, sqrt(prop.vars$SRB[j]))

        #.. make proposal
        logit.prop.SRB <- logit.curr.SRB + logit.prop.SRB.mat

        #.. If proposal resulted in back-transformed s = 0 or 1, do
        #   nothing
        if(invlogit(logit.prop.SRB[j]) > 1 - s.tol ||
           invlogit(logit.prop.SRB[j]) < s.tol) {
          #.. leave current surv rates and projections
          #   alone (simply do not propose
          #   extreme survival probabilities)
          s.out.tol[j] <- s.out.tol[j] + 1/n.iter
        } else {

          # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
            if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s) , Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.prop.SRB)#<-- use proposal
                , E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s)
                , Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.prop.SRB)#<-- use proposal
                , E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
                log.prop.proj = log(full.proj)
                log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

            # - Calculate log posterior of proposed vital under projection
            log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.prop.SRB #<-- use proposal
                       ,A = logit.curr.A
                       ,H = logit.curr.H
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )

            #- Acceptance ratio
            ar <- acc.ra(log.prop = log.prop.posterior,
                              log.current = log.curr.posterior)

            # - Move or stay
            #.. stay if acceptance ratio 0, missing, infinity, etc.
            if(is.na(ar) || is.nan(ar) || ar < 0) {
              if(i > burn.in) ar.na$SRB[j] <-
                  ar.na$SRB[j] + 1/n.iter
            } else {
              #.. if accept, update current surv rates,
              #   update current projection and count acceptance
              if(runif(1) <= ar) {
                if(i > burn.in) acc.count$SRB[j] <-
                    acc.count$SRB[j] + 1/n.iter
                    logit.curr.SRB <- logit.prop.SRB
                    log.curr.proj <- log.prop.proj
                    log.curr.aeri <- log.prop.aeri
                    log.curr.posterior <- log.prop.posterior
              } 

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(k %% 1 == 0 && k > 0) SRB.mcmc[k,] <-
        as.vector(invlogit(logit.curr.SRB))    
        

        
      ##...... Harvesting ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Harvesting")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.H)) {

        #.. make a matrix conformable w rate matrix
        prop.H.mat <-
            0*logit.curr.H
        prop.H.mat[j] <- rnorm(1, 0, sqrt(prop.vars$H[j])) # if need age-imspecific harvest, simply give a 1 by 1 start.H

        #.. make proposal
        logit.prop.H <- logit.curr.H + prop.H.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
        # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
      if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = ( invlogit(logit.prop.H))#<-- use proposal
                ,Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * as.numeric( invlogit(logit.prop.H)) , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar =( invlogit(logit.prop.H))#<-- use proposal
                ,Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar( invlogit(logit.prop.H[,1]),nage) , period = proj.periods, nage = nage))
            }

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
                log.prop.proj = log(full.proj)
                log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.prop.H #<-- use proposal
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )

        #- Acceptance ratio
        ar <- acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$H[j] <-
                ar.na$H[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$H[j] <-
                    acc.count$H[j] + 1/n.iter
                logit.curr.H = logit.prop.H
                log.curr.proj <- log.prop.proj
                log.curr.aeri <- log.prop.aeri
                log.curr.posterior <- log.prop.posterior
            } 

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

    } # close loop over all age-specific Harvest proportions

      #.. Store proposed Harvest proportion matrix
      if(k %% 1 == 0 && k > 0) H.mcmc[k,] <- as.vector(invlogit(logit.curr.H))
        
        
        # add Aerial count here 
    if(verb && identical(i%%1000, 0)) cat("\n", i, " Aerial Counts Detection")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.A)) {

        #.. make a matrix conformable w rate matrix
        prop.A.mat <-
            0*logit.curr.A
        prop.A.mat[j] <- rnorm(1, 0, sqrt(prop.vars$A[j])) # if need age-imspecific harvest, simply give a 1 by 1 start.H

        #.. make proposal
        logit.prop.A <- logit.curr.A + prop.A.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
        # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
      if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H)
                ,Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * as.numeric( invlogit(logit.prop.H)) , period = proj.periods, nage = nage))# no influence on projection (i.e. it is not a vital rate)
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H)
                ,Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar( invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$A[j] =
                        pop.negative$A[j] + 1/n.iter
                }
            } else {
                log.prop.proj = log(full.proj)
                log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.prop.A))) #<-- use proposal

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.prop.A #<-- use proposal
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )

        #- Acceptance ratio
        ar <- acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$A[j] <-
                ar.na$A[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$A[j] <-
                    acc.count$A[j] + 1/n.iter
                logit.curr.A = logit.prop.A
                log.curr.proj <- log.prop.proj
                log.curr.aeri <- log.prop.aeri
                log.curr.posterior <- log.prop.posterior
            } 

        } # close else after checking for ar=na, nan, zero

      } # close else after checking for negative population

    } # close loop over all age-specific Harvest proportions

      #.. Store proposed Harvest proportion matrix
      if(k %% 1 == 0 && k > 0) A.mcmc[k,] <- as.vector(invlogit(logit.curr.A))
          

      ##...... Carrying Capacity ......##
    if(estaK0){
      prop.aK0 = curr.aK0 + rnorm(length(curr.aK0), 0, sqrt(prop.vars$aK0))
      
            if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = prop.aK0 #<-- use proposal
                , global = global, null = null, bl = exp(log.curr.b) * as.numeric( invlogit(logit.curr.H)) , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = prop.aK0#<-- use proposal
                , global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
                log.prop.proj = log(full.proj)
                log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = prop.aK0 #<-- use proposal
                       ,baseline.n = log.curr.b
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )

        #- Acceptance ratio
        ar <- acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$aK0[j] <-
                ar.na$aK0[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$K0[j] <-
                    acc.count$aK0[j] + 1/n.iter
                curr.aK0 = prop.aK0
                log.curr.proj = log.prop.proj
                log.curr.posterior = log.prop.posterior
            } 

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population
    }
      #.. Store proposed K0 matrix
      if(k %% 1 == 0 && k > 0 && estaK0)  aK0.mcmc[k,] <- as.vector((curr.aK0))    
      

      ##...... Baseline population ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

      # - Proposal

      #.. cycle through components (never update last
      #   value as this affects years beyond the estimation period)
      for(j in 1:length(log.curr.b)) {

      #.. make a matrix conformable w rate matrix
      #log.prop.b.mat <- matrix(0, nrow = nrow(log.curr.b), ncol = 1)
        log.prop.b.mat = 0 * log.curr.b
      log.prop.b.mat[j] <- rnorm(1, 0, sqrt(prop.vars$baseline.pop.count[j]))

      #.. make proposal
      log.prop.b <- log.curr.b + log.prop.b.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
            if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.prop.b) * as.numeric( invlogit(logit.curr.H)) #<-- use proposal
                , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H), SRB = invlogit(logit.curr.SRB),Fec=exp(log.curr.f), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.prop.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) #<-- use proposal
                , period = proj.periods, nage = nage))
            }


      if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
         || is.nan(sum(full.proj))) {
        if(i > burn.in) {
          pop.negative$baseline.count[j] <-
              pop.negative$baseline.count[j] + 1/n.iter
        }
      } else {
        log.prop.proj <- log(full.proj)
        log.prop.aeri = log( getAerialCount(nage = nage, Harv = exp( log.curr.proj),H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

        # - Calculate log posterior of proposed vital under projection
         log.prop.posterior =
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.prop.b #<-- use proposal
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.prop.proj#<-- use proposal
                                ,ll.var = curr.sigmasq.n) +
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.prop.aeri#<-- use proposal
                                ,ll.var = curr.sigmasq.ae)
                       ,non.zero.fert = fert.rows 
                       )


        #- Acceptance ratio
        ar = acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$baseline.count[j] <-
                ar.na$baseline.count[j] + 1/n.iter
        } else {
            #.. if accept, update current mig rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$baseline.count[j] <-
                    acc.count$baseline.count[j] + 1/n.iter
                    log.curr.b = log.prop.b
                    log.curr.proj = log.prop.proj
                    log.curr.posterior = log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

  } # close loop over all age-specific baseline counts

      #.. Store proposed baseline count matrix
      if(k %% 1 == 0 && k > 0) baseline.count.mcmc[k,] <-
          as.vector(exp(log.curr.b))

# TEST stop here 10/26/2018
      ## ------- Variance Updates ------- ##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Variances")

      ##...... Fertility rate ......##
    if(estFer){ # if not est Fer, this is not needed
      prop.sigmasq.f <-
        rinvGamma(1, al.f +
                         length(mean.f[fert.rows,])/2,
                  be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                  log.mean.f[fert.rows,])^2)
                  )

        # - Calculate log posterior of proposed vital under projection
        
        log.prop.posterior =
                      log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB
                       ,A = logit.curr.A 
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(    
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.f <-
                ar.na$sigmasq.f + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.f <-
                    acc.count$sigmasq.f + 1/n.iter
                curr.sigmasq.f <- prop.sigmasq.f
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"fert.rate.var"] <- curr.sigmasq.f
    }
      ##...... Survival Proportion ......##

      prop.sigmasq.s <-
        rinvGamma(1, al.s + length(mean.s)/2,
                  be.s +
                    0.5*sum((logit.curr.s - logit.mean.s)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = prop.sigmasq.s #<-- use proposal
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )
        #pause here 20190520 during adding aerial count

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.s <-
                ar.na$sigmasq.s + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.s <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.s <- prop.sigmasq.s
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"surv.prop.var"] <- curr.sigmasq.s

      
      ##...... Sex Ratio at Birth ......## 
      prop.sigmasq.SRB <-
        rinvGamma(1, al.SRB + length(mean.SRB)/2,
                  be.SRB +
                    0.5*sum((logit.curr.SRB - logit.mean.SRB)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.SRB = prop.sigmasq.SRB #<-- use proposal
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.SRB
                              ,al.SRB + length(mean.SRB)/2
                              ,be.SRB + 0.5*sum((logit.curr.SRB -
                                               logit.mean.SRB)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.SRB
                              ,al.SRB + length(mean.SRB)/2
                              ,be.SRB + 0.5*sum((logit.curr.SRB -
                                               logit.mean.SRB)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.SRB <-
                ar.na$sigmasq.SRB + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.SRB <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.SRB <- prop.sigmasq.SRB
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"SRB.var"] <- curr.sigmasq.SRB  
      
      
      ##...... Aerial Count Detection ......## not changed yet
      prop.sigmasq.A <-
        rinvGamma(1, al.A + length(mean.A)/2,
                  be.A +
                    0.5*sum((logit.curr.A - logit.mean.A)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.SRB = curr.sigmasq.SRB
                       ,sigmasq.A = prop.sigmasq.A #<-- use proposal
                       ,sigmasq.H = curr.sigmasq.H
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )
        #pause here 20190520 during adding aerial count

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.A
                              ,al.A + length(mean.A)/2
                              ,be.A + 0.5*sum((logit.curr.A -
                                               logit.mean.A)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.A
                              ,al.A + length(mean.A)/2
                              ,be.A + 0.5*sum((logit.curr.A -
                                               logit.mean.A)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.A <-
                ar.na$sigmasq.A + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.A <-
                    acc.count$sigmasq.A + 1/n.iter
                curr.sigmasq.A <- prop.sigmasq.A
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"A.var"] <- curr.sigmasq.A      
      
      ##...... Harvest Proportion ......##

      prop.sigmasq.H <-
        rinvGamma(1, al.H + length(logit.mean.H)/2,
                  be.H +
                    0.5*sum((logit.curr.H - logit.mean.H)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.SRB = curr.sigmasq.SRB 
                       ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = prop.sigmasq.H #<-- use proposal
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.H
                              ,al.H + length(mean.H)/2
                              ,be.H + 0.5*sum((logit.curr.H -
                                               logit.mean.H)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.H
                              ,al.H + length(mean.H)/2
                              ,be.H + 0.5*sum((logit.curr.H -
                                               logit.mean.H)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.H <-
                ar.na$sigmasq.H + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.H <-
                    acc.count$sigmasq.H + 1/n.iter
                curr.sigmasq.H <- prop.sigmasq.H
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"H.var"] <- curr.sigmasq.H
    
      ##...... Carrying Capacity ......##    # Feel it will not be used  
      if(estaK0){
        prop.sigmasq.aK0 = rinvGamma(1,al.aK0+0.5,be.aK0 + 0.5*sum((curr.aK0-mean.aK0)^2)) #bug here give NA
            log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
					   ,sigmasq.SRB = curr.sigmasq.SRB 
					   ,sigmasq.A = curr.sigmasq.A
                       ,sigmasq.H = curr.sigmasq.H 
                       ,sigmasq.aK0 = prop.sigmasq.aK0 #<-- use proposal
                       ,sigmasq.n = curr.sigmasq.n
                       ,sigmasq.ae = curr.sigmasq.ae
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )

        #- Acceptance ratio
      ar = acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.aK0
                              ,al.aK0 + length(mean.aK0)/2
                              ,be.aK0 + 0.5*sum((curr.aK0 -
                                               mean.aK0)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.aK0
                              ,al.aK0 + length(mean.aK0)/2
                              ,be.aK0 + 0.5*sum((curr.aK0 -
                                               mean.aK0)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.aK0 <-
                ar.na$sigmasq.aK0 + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.aK0 <-
                    acc.count$sigmasq.aK0 + 1/n.iter
                curr.sigmasq.aK0 <- prop.sigmasq.aK0
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"aK0.var"] <- curr.sigmasq.aK0

      
      
      
      }
        
        
      ##...... Harvest Count ......##

      prop.sigmasq.n <-
        rinvGamma(1, al.n + (length(mean.b) +
                                    length(log.Harv.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2) +
                  sum((log.Harv.mat - log.curr.proj)^2) # be careful if we add reconstruction, things may go into here
                  )
                              )

        # - Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A 
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.SRB = curr.sigmasq.SRB                     
                       ,sigmasq.A = curr.sigmasq.A 
                       ,sigmasq.H = curr.sigmasq.H 
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                       ,sigmasq.ae = curr.sigmasq.ae 
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = prop.sigmasq.n) +#<-- use proposal
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = curr.sigmasq.ae)#<-- use current
                       ,non.zero.fert = fert.rows 
                       )

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.Harv.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.Harv.mat - log.curr.proj)^2))
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.Harv.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.Harv.mat - log.curr.proj)^2))
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.n <-
                ar.na$sigmasq.n + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.n <-
                    acc.count$sigmasq.n + 1/n.iter
                curr.sigmasq.n <- prop.sigmasq.n
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) {
        variances.mcmc[k,"population.count.var"] <- curr.sigmasq.n
      }





      ##...... Aerial Count ......##

      prop.sigmasq.ae <-
        rinvGamma(1, al.ae + (length(log.Aeri.mat))/2,
                be.ae + 0.5 * (sum((log.Aeri.mat - log.curr.aeri)^2)))

        # - Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,SRB = logit.curr.SRB 
                       ,A = logit.curr.A 
                       ,H = logit.curr.H 
                       ,aK0 = curr.aK0 
                       ,baseline.n = log.curr.b 
                       ,estFer=estFer, estaK0=estaK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.SRB = logit.mean.SRB
                       ,prior.mean.A = logit.mean.A
                       ,prior.mean.H = logit.mean.H
                       ,prior.mean.aK0 = mean.aK0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                       ,alpha.A = al.A, beta.A = be.A
                       ,alpha.H = al.H, beta.H = be.H
                       ,alpha.aK0 = al.aK0, beta.aK0 = be.aK0
                       ,alpha.n = al.n, beta.n = be.n
                       ,alpha.ae = al.ae, beta.ae = be.ae
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.SRB = curr.sigmasq.SRB                     
                       ,sigmasq.A = curr.sigmasq.A 
                       ,sigmasq.H = curr.sigmasq.H 
                       ,sigmasq.aK0 = curr.sigmasq.aK0
                       ,sigmasq.n = curr.sigmasq.n 
                       ,sigmasq.ae = prop.sigmasq.ae #<-- use proposal
                       ,log.like = 
                            log.lhood(
                                log.n.census = log.Harv.mat 
                                ,log.n.hat = log.curr.proj#<-- use current
                                ,ll.var = curr.sigmasq.n) +#<-- use current
                            log.lhood(
                                log.n.census = log.Aeri.mat
                                ,log.n.hat = log.curr.aeri#<-- use current
                                ,ll.var = prop.sigmasq.ae)#<-- use proposal
                       ,non.zero.fert = fert.rows 
                       )

      #- Acceptance ratio
      ar <- acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.n
                              ,al.ae + ( length(log.Aeri.mat))/2
                              ,be.ae + 0.5 * (sum((log.Aeri.mat - log.curr.aeri)^2))
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.ae
                              ,al.ae + (
                                       length(log.Aeri.mat))/2
                              ,be.ae + 0.5 * ( sum((log.Aeri.mat - log.curr.aeri)^2))
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.ae <-
                ar.na$sigmasq.ae + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.ae <-
                    acc.count$sigmasq.ae + 1/n.iter
                curr.sigmasq.ae <- prop.sigmasq.ae
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) {
        variances.mcmc[k,"aerial.count.var"] <- curr.sigmasq.ae
      }









      ## ------- Store current population ------- ##
            if(homo){
                full.proj =
                (ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * as.numeric( invlogit(logit.curr.H)) , period = proj.periods, nage = nage))
            }
            else{
                full.proj =
                (ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Fec=exp(log.curr.f), SRB = invlogit(logit.curr.SRB), E0=E0, aK0 = (curr.aK0), global = global, null = null, bl = exp(log.curr.b) * getfullHarvpar(invlogit(logit.curr.H[,1]),nage) , period = proj.periods, nage = nage))
            }

      lx.mcmc[k,] =
          as.vector(full.proj)[-(1:ncol(baseline.count.mcmc))] # to delete all age class' baseline count, because of as.vector,thus need to do like this
      ae.mcmc[k,] = 
          as.vector(getAerialCount(nage = nage, Harv = full.proj,H = invlogit(logit.curr.H),A = invlogit(logit.curr.A)))

      if(verb && identical(i%%1000, 0)) cat("\n\n")

  
      #print(i)
      #warnings() # debug mode 
      #Sys.sleep(5)
      } # Ends outer-most loop

    ## ......... End Loop ........ ##
    #...............................#


    ## ---------- Output --------- ##

    #cat("inital values", "\n\n")
    #.. initial values
    start.vals <- list(fert.rate = start.f
                      ,surv.prop = start.s
                      ,SRB = start.SRB
                      ,H = start.H
                      ,Aerial.detection = start.A
                      ,K0 = start.aK0
                      ,baseline.count = start.b
                      ,start.sigmasq.f = start.sigmasq.f
                      ,start.sigmasq.s = start.sigmasq.s
                      ,start.sigmasq.SRB = start.sigmasq.SRB
                      ,start.sigmasq.A = start.sigmasq.A
                      ,start.sigmasq.H = start.sigmasq.H
                      ,start.sigmasq.aK0 = start.sigmasq.aK0
                      ,start.sigmasq.n = start.sigmasq.n
                      ,start.sigmasq.ae = start.sigmasq.ae
                      ,Harv.data = Harv.data
                      ,Aerial.data = Aerial.data
                      )

    #.. fixed parameters
    fixed.params <- list(alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.SRB = al.SRB
                         ,beta.SRB = be.SRB
                         ,alpha.aerial.det = al.A
                         ,beta.aerial.det = be.A
                         ,alpha.Harvest = al.H
                         ,beta.Hervest = be.H
                         ,alpha.1overK = al.aK0
                         ,beta.1overK = be.aK0
                         ,alpha.harvest.count = al.n
                         ,beta.harvest.count = be.n
                         ,alpha.aerial.count = al.ae
                         ,beta.aerial.count = be.ae
                         ,mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.SRB = mean.SRB
                         ,mean.aerial.detection = mean.A
                         ,mean.Harvest.proportion = mean.H
                         ,mean.1overK = mean.aK0
                         ,mean.baseline.count = mean.b
                         ,mean.Harv.data = Harv.data
                         ,Aerial.data = Aerial.data
                         )



    #cat("algorithm statistics", "\n\n")
    #.. algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             #,acc.prop.adj4neg = mapply(FUN = function(a, b, n) {
            #     (a * n) / (n - b)
             #},
              #acc.count[1:5], pop.negative, MoreArgs = list(n = n.iter)
              #)
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    #cat("algorithm parameters", "\n\n")
    #.. algorithm parameters
    alg.params <- list(prop.vars = prop.vars
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit",SRB = "logit",aerial.detection="logit", H = "logit", invK0="identical"
                        ,baseline.count = "log"
                        ,harvest.count = "log"
                        ,aerial.count = "log")
                       ,projection.periods = proj.periods
                       ,age.gp.size = nage
                       ,non.zero.fert.rows = fert.rows
                       ,surv.tolerance = s.tol
                       ,burn.in = burn.in
                       ,iters = n.iter
                        )
                       

    #.. results
    ret.list <- list(fert.rate.mcmc = fert.rate.mcmc
                  ,surv.prop.mcmc = surv.prop.mcmc
                  ,SRB.mcmc = SRB.mcmc
                  ,aerial.detection.mcmc = A.mcmc
                  ,H.mcmc = H.mcmc
                  ,invK0.mcmc = aK0.mcmc
                  ,baseline.count.mcmc = baseline.count.mcmc
                  ,harvest.mcmc = lx.mcmc
                  ,aerial.count.mcmc = ae.mcmc
                  ,variances.mcmc = variances.mcmc
                  ,alg.stats = alg.stats
                  ,fixed.params = fixed.params
                  ,start.vals = start.vals
                  ,alg.params = alg.params
                  )
    cat("\n done \n")
    return(ret.list)
# stop here 10/23/2018 15:51 start to test tomorrow.
# TEST of code done here, start SIMULATION
}
