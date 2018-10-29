########
### ----------------------- HELPER FUNCTIONS ---------------------- ###
### --------------------------------------------------------------- ###

# get Lislie Matrix using specific parameters, can also return L-I if needed, first nage entries should be survival(no transiform, 0<S<1,transformation occor in sampler function) and next fertility rate, 
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
DensityDependcy = function(global = F, Xn, E0, K0, null = F){
  E0 = E0/sum(E0)
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
} #checked 10/24/2018

# get age specific harvest proportion matrix
GetHarvest = function(Harvpar,nage){
  Harvest = matrix(0,nage,nage)
  diag(Harvest) = Harvpar
  return(Harvest)
} # checked 10/24/2018, if not age-specific, give a single Harvpar

# Project the (density dependent) harvest-after-reproducing model from year i to i+1, given harvest # and return harvest # of year i+1 
ProjectHarvest_helper = function(data_n, Lislie, H, global, E0, K0, null = F){
	I = matrix(0,length(data_n),length(data_n))	
	diag(I) = 1
    X_n1 = (I-H) %*% solve(H,data_n)
    D = DensityDependcy(global = global, Xn=X_n1, E0=E0, K0=K0, null = null)
	data_n1 = H %*% (Lislie %*% D %*% X_n1 + X_n1)
    #Popu_after = (eyes-H)%*%Popu_before_harvest 
    return(data_n1)
  } # checked 10/24/2018

# Project harvest model from a initial harvest, survival is col vector with all survival rate of all age class, this nrow(Survival)=nage
ProjectHarvest_homo = function(Survival, Harvpar,Ferc, E0=NULL, K0 = NULL, global = F, null = F,bl , period, nage){
  Lislie = getLislie(Survival,Ferc=Ferc,minus1=T)
  H = GetHarvest(Harvpar,nage)
  Harvest = matrix(0,nage,period + 1)
  Harvest[,1] = bl 
  if(length(E0)==0){
	E0 = solve(H,bl)
	E0 = E0/(sum(E0))
  }
  else E0 = E0/(sum(E0))
  for(i in 1 : period + 1){
	Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Lislie = Lislie, E0=E0, K0=K0, H=H,null = null) # project harvest from very initial
  }              
  return(Harvest)
} #checked 10/24/2018

# project when time inhomo survival and ferc, survival should be a matrix with nage x period entries, so do ferc
ProjectHarvest_inhomo = function(Survival, Harvpar,Ferc, E0=NULL, K0 = NULL, global = F, null = F,bl , period, nage){
  Harvest = matrix(0,nage,period+1)
  Harvest[,1] = bl
  H = GetHarvest(Harvpar,nage)
  if(length(E0)==0){
	E0 = solve(H,bl)
	E0 = E0/(sum(E0))
  }
  else E0 = E0/(sum(E0))
  for(i in 1 : period + 1){
    Lislie = getLislie(Survival[,i-1],Ferc[,i-1],minus1=T) # inhomo, Survival rows are age structure, cols are time
    Harvest[,i] = ProjectHarvest_helper(Harvest[,i-1],global = global, Lislie = Lislie, E0=E0, K0=K0, H=H,null = null)
  }
  return(Harvest)
} # checked 10/24/2018

# Given harvest data, calculate the living individual in certain year
getLivingIdividuals = function(H,data){
  I = matrix(0,nrow(data),nrow(data))	
  diag(I) = 1
  IminusH = I-H
  LivingIdividuals = apply(data,2,function(ww,H,IminusH){IminusH%*%solve(H,ww)},H=H,IminusH=IminusH)
  return(LivingIdividuals)
} # checked 10/24/2018


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
    return(1/rgamma(n, shape = shape, rate = scale))
} # checked 10/24/2018 

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

    density = dnorm(log.n.census,
                     mean = log.n.hat,
                     sd = sqrt(ll.var),
                     log = TRUE
                     )

    ##-- joint log likelihood --##

    return(sum(density))
} # checked 10/24/2018

## .............. Posterior ............ ##
## ..................................... ##  keep it, change in sampler function, but no migration here, should add estK0

log.post = function(## estimated vitals
                           f, s, baseline.n, K0, H
						               ,estFer, estK0
                           ## fixed prior means on vitals
                           ,prior.mean.f, prior.mean.s
                           ,prior.mean.b, prior.mean.K0
						               ,prior.mean.H
                           ## fixed prior parameters on variance distns
                           ,alpha.f, beta.f, alpha.s, beta.s
                           ,alpha.n, beta.n, alpha.K0, beta.K0
						               ,alpha.H, beta.H
                           ## updated variances on prior distns
                           ,sigmasq.f, sigmasq.s, sigmasq.n, sigmasq.K0
						               ,sigmasq.H
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
	log.H.prior = dnorm(H, mean = prior.mean.H
                         ,sd = sqrt(sigmasq.H)
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
             ,al.f = 1, be.f = 0.0109, al.s = 1, be.s = 0.0109
             , al.K0 = 1, be.K0 = 0.0436, al.n = 1
			       , be.n = 0.0109, al.H = 1, be.H = 0.0436

             #.. fixed prior means
             ,mean.f, mean.s, mean.b, mean.K0, mean.H

             #.. inital values for vitals and variances
             #   *vitals not transformed coming in* all not transfer, will 
             ,start.f = mean.f, start.s = mean.s
			       ,start.b = mean.b, start.K0 = mean.K0, start.H = mean.H
             ,start.sigmasq.f = 5, start.sigmasq.s = 5
             ,start.sigmasq.n = 5, start.sigmasq.K0 = 5, start.sigmasq.H = 5

             #.. census data
             #   *not transformed coming in*
             ,Harv.data

             #.. **variances** for proposal distributions used in M-H
             #   steps which update vital rates.
             ,prop.vars # col names should be "fert.rate", 
			       # "fert.rate.var", "surv.prop.var", "H.var", "K0.var"
			       # ,"population.count.var"
             #.. ccmp function

             #.. number of periods to project forward over (e.g.,
             #     number of five-year steps)
             ,proj.periods = (ncol(Harv.data)-1)

             #.. age group width
             ,nage = 3 
			 
			       ,estFer=T, Ferc=rep(1,nage), estK0 = T
			       ,E0=NULL , K0 = 0, global = T, null = T # control parameters for the model, global is whether density dependency is global rather than age specific, null is whether exist density dependency.
			       ,timehomo = F # whether assume time homogeneous 

             #.. print algorithm progress
             ,verb = FALSE

             #.. tolerance defining allowable survival probabilities
             ,s.tol = 10^(-10)
             )
{
	require(coda)
    ## .............. Sampler .............. ##
    ## ..................................... ##


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
	
    n.stored = ceiling(n.iter / thin.by)
	ntimes = (!timehomo) * ncol(start.s) + (timehomo) # whether assume time homogeneous of survival etc.
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

      # lx
	  # this is for current population, for us, it is current culling data
      lx.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.b) * (proj.periods))
               ,start = burn.in + 1
               ,thin = thin.by
               )
      colnames(lx.mcmc) = NULL
	  
	  # carrying capacity assumed to be time homogeneous
	  if(estK0){
		K0.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = 1)
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(K0.mcmc) = NULL
	  } 
	  else{
	    K0.mcmc = NULL
	  }

      # Harvest proportion, can be either time homo or not
      H.mcmc =
          mcmc(matrix(nrow = n.stored
                      ,ncol = nrow(start.H) * ntimes)
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
          mcmc(matrix(nrow = n.stored, ncol = 5)
               ,start = burn.in + 1
               ,thin = thin.by)
      colnames(variances.mcmc) =
          c("fert.rate.var", "surv.prop.var", "H.var", "K0.var"
            ,"population.count.var") # may need check here since K0 and fert can be known (though I prefer estimating it)

    #.. Record acceptance rate

    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,H = matrix(0, nrow = nrow(mean.H), ncol = ncol(mean.H)
              ,dimnames = dimnames(mean.H)
              )
			 ,K0 = matrix(0, nrow = nrow(mean.K0), ncol = ncol(mean.K0)
              ,dimnames = dimnames(mean.K0)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.H = 0
			 ,sigmasq.K0 = 0
             ,sigmasq.n = 0
             )


    #.. Count how often acceptance ratio missing or na

    ar.na <- acc.count


    #.. Count how often projection gives negative population

    pop.negative <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,H = matrix(0, nrow = nrow(mean.H), ncol = ncol(mean.H)
              ,dimnames = dimnames(mean.H)
              )
			 ,K0 = matrix(0, nrow = nrow(mean.K0), ncol = ncol(mean.K0)
              ,dimnames = dimnames(mean.K0))
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             )


    #.. Count how often surv probs are outside tolerance

    s.out.tol <- matrix(0, nrow = nrow(mean.s), ncol = ncol(mean.s)
                        ,dimnames = dimnames(mean.s))


    ## -------- Initialize -------- ## Restart here in 10/19/2018

    #.. Set current vitals and variances to inital values
    #   Take logs/logits here where required
    if(estFer){
      log.curr.f = log(start.f)  #<-- log(0) stored as "-Inf". Gets
      log.prop.f = log(start.f)  #    converted to 0 under exponentiation
      
    }
    else{
      log.curr.f =  (!estFer)*log(Ferc) #<-- log(0) stored as "-Inf". Gets
      log.prop.f =  (!estFer)*log(Ferc) #    converted to 0 under exponentiation
    }
    logit.curr.s = logitf(start.s)
	  logit.curr.H = logitf(start.H)
    log.curr.K0 = estK0 * log(start.K0) + (!estK0) * log(K0)
    log.curr.b = log(start.b)

    curr.sigmasq.f = start.sigmasq.f
    curr.sigmasq.s = start.sigmasq.s
    curr.sigmasq.H = start.sigmasq.H
	curr.sigmasq.K0 = start.sigmasq.K0
    curr.sigmasq.n = start.sigmasq.n


    #.. Fixed means for vitals and baseline
    #   Set these to inputs, take logs where required.

    log.mean.f = log(mean.f)
    logit.mean.s = logitf(mean.s)
    logit.mean.H = logitf(mean.H)
	log.mean.K0 = log(mean.K0)
    log.mean.b = log(mean.b)


    #.. Fixed Harvest data
    #   Take logs here

    log.Harv.mat = log(Harv.data)


    #.. Set current projection: base on initial values # timehomo or not is important, determin it use homo = T
	if(homo){
	  log.curr.proj =
        log(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
	}
	else{
	  log.curr.proj =
        log(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
	}
## stop here 10/19/2018   
## restart here 10/22/2018

    #.. Current log posterior
	## shit so many parameters to pass here... really hard to read...
    
# checked here,  solved, -Inf return 10/25/2018 14:31, Done, because of test data has zero likelihood (changed parameter)
    log.curr.posterior =
        log.post(f = log.curr.f
                       ,s = logit.curr.s
                       ,H = logit.curr.H
					   ,K0 = log.curr.K0
                       ,baseline.n = log.curr.b
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj
                        ,ll.var = curr.sigmasq.n)
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

    for(i in 1:(n.iter + burn.in)) {

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
				(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.prop.f) #<-- use proposal
				, E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
		}
		else{
			full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.prop.f)#<-- use proposal
				, E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
		}


        if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
           || is.nan(sum(full.proj))) {
            if(i > burn.in) {
                pop.negative$fert.rate[j] =
                    pop.negative$fert.rate[j] + 1/n.iter
            }
        } else {
            log.prop.proj = log(full.proj)

          # - Calculate log posterior of proposed vital under projection		  
		  ## shit again so many parameters to pass here... really hard to read...

		  log.prop.posterior =
			  log.post(f = log.prop.f #<-- use proposal
                       ,s = logit.curr.s
                       ,H = logit.curr.H
					   ,K0 = log.curr.K0
                       ,baseline.n = log.curr.b
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.prop.proj #<-- use proposal
                        ,ll.var = curr.sigmasq.n)
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
				, Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}
			else{
				full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.prop.s)#<-- use proposal
				, Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
				log.prop.proj = log(full.proj)

            # - Calculate log posterior of proposed vital under projection
			log.prop.posterior =
			  log.post(f = log.curr.f 
                       ,s = logit.prop.s #<-- use proposal
                       ,H = logit.curr.H
					   ,K0 = log.curr.K0
                       ,baseline.n = log.curr.b
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.prop.proj #<-- use proposal
                        ,ll.var = curr.sigmasq.n)
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
                log.curr.posterior <- log.prop.posterior
              } 

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(k %% 1 == 0 && k > 0) surv.prop.mcmc[k,] <-
        as.vector(invlogit(logit.curr.s))


      ##...... Harvesting ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Harvesting")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.H)) {

        #.. make a matrix conformable w rate matrix
        prop.H.mat <-
            matrix(0, nrow = nrow(logit.curr.H), ncol = ncol(logit.curr.H))
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
				(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.prop.H)#<-- use proposal
				,Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}
			else{
				full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.prop.H)#<-- use proposal
				,Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
				log.prop.proj = log(full.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
			  log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.prop.H #<-- use proposal
					   ,K0 = log.curr.K0
                       ,baseline.n = log.curr.b
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.prop.proj #<-- use proposal
                        ,ll.var = curr.sigmasq.n)
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
                log.curr.posterior <- log.prop.posterior
            } 

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

    } # close loop over all age-specific Harvest proportions

      #.. Store proposed Harvest proportion matrix
      if(k %% 1 == 0 && k > 0) H.mcmc[k,] <- as.vector(invlogit(logit.curr.H))


      ##...... Carrying Capacity ......##
	if(estK0){
	  log.prop.K0 = log.curr.K0 + rnorm(1, 0, sqrt(prop.vars$K0))
	  
	        if(homo){
				full.proj =
				(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.prop.K0) #<-- use proposal
				, global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}
			else{
				full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.prop.K0)#<-- use proposal
				, global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}

            if(sum(full.proj < 0) > 0 || is.na(sum(full.proj))
               || is.nan(sum(full.proj))) {
                if(i > burn.in) {
                    pop.negative$surv.prop[j] =
                        pop.negative$surv.prop[j] + 1/n.iter
                }
            } else {
				log.prop.proj = log(full.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
			  log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.curr.H 
					   ,K0 = log.prop.K0 #<-- use proposal
                       ,baseline.n = log.curr.b
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.prop.proj #<-- use proposal
                        ,ll.var = curr.sigmasq.n)
                       ,non.zero.fert = fert.rows 
                       )

        #- Acceptance ratio
        ar <- acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$K0[j] <-
                ar.na$K0[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$K0[j] <-
                    acc.count$K0[j] + 1/n.iter
                curr.K0 = prop.K0
                log.curr.proj = log.prop.proj
                log.curr.posterior = log.prop.posterior
            } 

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population
	}
      #.. Store proposed K0 matrix
      if(k %% 1 == 0 && k > 0) K0.mcmc[k,] <- as.vector(exp(log.curr.K0))    
      
# stop here for flu shot 10/22/2018	  17:07
# restart here 10/22/2018 20:45. Left arm is So0o0o0o0o painful after the flu shot!!!!! 
	  ##...... Baseline population ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

      # - Proposal

      #.. cycle through components (never update last
      #   value as this affects years beyond the estimation period)
      for(j in 1:length(log.curr.b)) {

      #.. make a matrix conformable w rate matrix
      log.prop.b.mat <- matrix(0, nrow = nrow(log.curr.b), ncol = 1)
      log.prop.b.mat[j] <- rnorm(1, 0, sqrt(prop.vars$baseline.pop.count[j]))

      #.. make proposal
      log.prop.b <- log.curr.b + log.prop.b.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
	        if(homo){
				full.proj =
				(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.prop.b) #<-- use proposal
				, period = proj.periods, nage = nage))
			}
			else{
				full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.prop.b) #<-- use proposal
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

        # - Calculate log posterior of proposed vital under projection
         log.prop.posterior =
			  log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.prop.b #<-- use proposal
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.prop.proj #<-- use proposal
                        ,ll.var = curr.sigmasq.n)
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
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.curr.b 
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj 
                        ,ll.var = curr.sigmasq.n) #<-- use current
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
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.curr.b 
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = prop.sigmasq.s #<-- use proposal
                       ,sigmasq.H = curr.sigmasq.H
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj 
                        ,ll.var = curr.sigmasq.n) #<-- use current
                       ,non.zero.fert = fert.rows 
                       )

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


      ##...... Harvest Proportion ......##

      prop.sigmasq.H <-
        rinvGamma(1, al.H + length(logit.mean.H)/2,
                  be.H +
                    0.5*sum((logit.curr.H - logit.mean.H)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.curr.b 
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.H = prop.sigmasq.H #<-- use proposal
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj 
                        ,ll.var = curr.sigmasq.n) #<-- use current
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

# stop here 10/22/2018 21:36 Left arm is too pain to continue.
# restart here 10/23/2018 14:34	
      ##...... Carrying Capacity ......##	
	  if(estK0){
		prop.sigmasq.K0 = rinvGamma(1,al.K0+0.5,be.K0 + 0.5*sum((log.curr.K0-log.mean.K0)^2))
		    log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.curr.b 
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.H = curr.sigmasq.H 
					   ,sigmasq.K0 = prop.sigmasq.K0 #<-- use proposal
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj 
                        ,ll.var = curr.sigmasq.n) #<-- use current
                       ,non.zero.fert = fert.rows 
                       )

	    #- Acceptance ratio
      ar = acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = dinvGamma(prop.sigmasq.K0
                              ,al.K0 + length(mean.H)/2
                              ,be.K0 + 0.5*sum((logit.curr.H -
                                               logit.mean.H)^2)
                              ,log = TRUE)
                             ,log.curr.var = dinvGamma(curr.sigmasq.K0
                              ,al.K0 + length(mean.K0)/2
                              ,be.K0 + 0.5*sum((log.curr.K0 -
                                               log.mean.K0)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.K0 <-
                ar.na$sigmasq.K0 + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.K0 <-
                    acc.count$sigmasq.K0 + 1/n.iter
                curr.sigmasq.K0 <- prop.sigmasq.K0
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(k %% 1 == 0 && k > 0) variances.mcmc[k,"K0.var"] <- curr.sigmasq.K0

	  
	  
	  
	  }
		





      ##...... Population Count ......##

      prop.sigmasq.n <-
        rinvGamma(1, al.n + (length(mean.b) +
                                    length(log.Harv.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2) +
                  sum((log.Harv.mat - log.curr.proj)^2)
                  )
                              )

        # - Calculate log posterior of proposed vital under projection
		        log.prop.posterior <-
              log.post(f = log.curr.f 
                       ,s = logit.curr.s 
                       ,H = logit.curr.H 
					   ,K0 = log.curr.K0 
                       ,baseline.n = log.curr.b 
					   ,estFer=estFer, estK0=estK0
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.H = logit.mean.H
					   ,prior.mean.K0 = log.mean.K0
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.H = al.H, beta.H = be.H
					   ,alpha.K0 = al.K0, beta.K0 = be.K0
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s 
                       ,sigmasq.H = curr.sigmasq.H 
					   ,sigmasq.K0 = curr.sigmasq.K0
                       ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                       ,log.like = log.lhood(
                        log.n.census = log.Harv.mat
                        ,log.n.hat = log.curr.proj 
                        ,ll.var = prop.sigmasq.n) #<-- use proposal
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


      ## ------- Store current population ------- ##
			if(homo){
				full.proj =
				(ProjectHarvest_homo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}
			else{
				full.proj =
				(ProjectHarvest_inhomo(Survival = invlogit(logit.curr.s), Harvpar = invlogit(logit.curr.H),Ferc=exp(log.curr.f), E0=E0, K0 = exp(log.curr.K0), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))
			}

      lx.mcmc[k,] =
          as.vector(full.proj)[-(1:ncol(baseline.count.mcmc))] # to delete all age class' baseline count, because of as.vector,thus need to do like this


      if(verb && identical(i%%1000, 0)) cat("\n\n")

  } # Ends outer-most loop

    ## ......... End Loop ........ ##
    #...............................#


    ## ---------- Output --------- ##

    #cat("inital values", "\n\n")
    #.. initial values
    start.vals <- list(fert.rate = start.f
                      ,surv.prop = start.s
                      ,H = start.H
					  ,K0 = start.K0
                      ,baseline.count = start.b
                      ,start.sigmasq.f = start.sigmasq.f
                      ,start.sigmasq.s = start.sigmasq.s
                      ,start.sigmasq.H = start.sigmasq.H
					  ,start.sigmasq.K0 = start.sigmasq.K0
                      ,start.sigmasq.n = start.sigmasq.n
                      ,Harv.data = Harv.data
                      )

    #.. fixed parameters
    fixed.params <- list(alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.Harvest = al.H
                         ,beta.Hervest = be.H
						 ,alpha.K = al.K0
						 ,beta.K = be.K0
                         ,alpha.population.count = al.n
                         ,beta.population.count = be.n
                         ,mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.Harvest.proportion = mean.H
						 ,mean.K = mean.K0
                         ,mean.baseline.count = mean.b
                         ,mean.Harv.data = Harv.data
                         )



    #cat("algorithm statistics", "\n\n")
    #.. algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             ,acc.prop.adj4neg = mapply(FUN = function(a, b, n) {
                 (a * n) / (n - b)
             },
              acc.count[1:5], pop.negative, MoreArgs = list(n = n.iter)
              )
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    #cat("algorithm parameters", "\n\n")
    #.. algorithm parameters
    alg.params <- list(prop.vars = prop.vars
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit", H = "logit", K0="log"
                        ,baseline.count = "log"
                        ,population.count = "log")
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
                  ,H.mcmc = H.mcmc
				  ,K0.mcmc = K0.mcmc
                  ,baseline.count.mcmc = baseline.count.mcmc
                  ,lx.mcmc = lx.mcmc
                  ,variances.mcmc = variances.mcmc
                  ,alg.stats = alg.stats
                  ,fixed.params = fixed.params
                  ,start.vals = start.vals
                  ,alg.params = alg.params
                  )

    return(ret.list)
# stop here 10/23/2018 15:51 start to test tomorrow.
# TEST of code done here, start SIMULATION
}
