sim_jm.PRO <- function(n, t.max, minobs, betas, phi, ntrial,
                       alpha, nu, gamma, D, type = NULL){
  
  if(is.null(type) == TRUE){
    type <- "p"
  }
  
  if(!type %in% c("p", "mp")){
    stop("relationship between survival and longitudinal data should be p or mp")
  }
  
  #simulate random effects
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  
  # survival time
  if(type == "mp"){
    invS <- function(t, u, i) {
      h <- function(s) {
        eta.s <- betas[1] + betas[2]*s + b[i,1] + b[i,2]*s
        exp(log(nu) + (nu - 1) * log(s) + gamma[1] + (ntrial*exp(eta.s)/(1+exp(eta.s))) * alpha)
      }
      integrate(h, lower = 0, upper = t)$value + log(u)
    }
  } else {
    invS <- function(t, u, i) {
      h <- function(s) {
        eta.s <- betas[1] + betas[2]*s + b[i,1] + b[i,2]*s
        exp(log(nu) + (nu - 1) * log(s) + gamma[1] + (exp(eta.s)/(1+exp(eta.s))) * alpha)
      }
      integrate(h, lower = 0, upper = t)$value + log(u)
    }
  }
  
  u <- runif(n)
  
  trueTimes <- numeric(n)
  for (i in 1:n) {
    Up <- 50
    tries <- 10
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 250
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  
  trueTimes <- ifelse(is.na(trueTimes),t.max,trueTimes)
  
  # simulate censoring times from a uniform distribution
  
  Ctimes <- numeric(n)
  Ctimes <- runif(n, 0,t.max)
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  
  
  # at which time points longitudinal measurements are supposed to be taken
  
  # Longitudinal time points
  nmaxobs <- ceiling(t.max) + minobs - 1
  nobs <- ceiling(Time) + minobs
  times <- NULL
  for(i in 1:n){
    times <- c(times,c(seq(0,Time[i],len=nobs[i])[-nobs[i]],rep(NA,nmaxobs-nobs[i]+1)))
  }
  
  # Simulating longitudinal data
  
  repb0 <- rep(b[,1],times=rep(nmaxobs,n))
  repb1 <- rep(b[,2],times=rep(nmaxobs,n))
  
  eta.y <- betas[1] + repb0 +  (betas[2]+repb1)*times 
  mu.y <- as.vector(exp(eta.y)/(1+exp(eta.y)))
  
  Y<- c()
  for (i in 1:length(mu.y)) {
    if(is.na(mu.y[i])){
      Y[i]<-NA
    } else{
      Y[i]<- PROreg::rBB(1,ntrial ,mu.y[i],phi)
    }
    
  }
  
  
  id <- rep(1:n,times=rep(nmaxobs,n))
  dta.long <- data.frame(id = id, y = Y, time = times)
  
  # Longitudinal and survival data
  dat.long <- dta.long[which(times <= rep(Time,rep(nmaxobs,n))),]
  dat.surv <- data.frame(id = unique(id), Time = Time, status = event)
  
  
  out <- list(longitudinal=dat.long, survival=dat.surv)
  out
}
