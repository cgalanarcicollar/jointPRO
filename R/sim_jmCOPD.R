sim_jm.COPD <- function(n, t.max, betas, phi, ntrial, alpha, 
                        nu, gamma, D, cens, type = NULL){
  
  if( is.null(type) == TRUE){
    type <- "p"
  }

  if(!type %in% c("p","mp")){
    stop("relationship between survival and longitudinal data should be p or mp")
  }
  
  #simulate random effects
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  
  # time measurements 
  dat=data.frame(id = as.factor(1:n),t0 = runif(n,0,1.5))
  dat$t1<-dat$t0+1
  dat$t2<-dat$t1+1
  dat$t3<-dat$t2+3
  
  
  # generate longitudinal data set in long format
  dat.long<- tidyr::pivot_longer(dat, cols=2:5, names_to = "meas", values_to = "time" )
  
  # model matrices
  X <- model.matrix(~time, data=dat.long) 
  Z.i <- model.matrix(~id-1, data=dat.long)
  Z.t <- model.matrix(~id:time-1, data=dat.long)
  Z <- cbind(Z.i, Z.t)
  
  #generate longitudinal PRO
  eta.y <- as.vector(X%*%betas + Z%*%c(b[,1 ],b[,2]))
  mu.y <- exp(eta.y)/(1+exp(eta.y))
  y<-PROreg::rBB(n*4, ntrial, mu.y, phi)
  dat.long$y <- y 
  
  #generate survival data
  if(type == "mp"){
  Surv<-function(t,j,u){
    h <-function (s) { 
      eta.s <-betas[1]+betas[2]*s+b[i,1]+b[i,2]*s
      nu*s^(nu-1)*exp(gamma[1]+alpha*(ntrial*exp(eta.s)/(1+exp(eta.s))))
    }
    integrate(h, lower = dat$t0[j], upper = t)$value + log(u)
  }
  }else{
    Surv<-function(t,j,u){
      h <-function (s) { 
        eta.s <-betas[1]+betas[2]*s+b[i,1]+b[i,2]*s
        nu*s^(nu-1)*exp(gamma[1]+alpha*(exp(eta.s)/(1+exp(eta.s))))
      }
      integrate(h, lower = dat$t0[j], upper = t)$value + log(u)
    }
  }
  
  u <- runif(n)
  e.time<-c()
  
  for (i in 1:n) {
    Up <- 50
    tries <- 10
    Root <- try(uniroot(Surv, interval = c(1e-05, Up), u = u[i], j = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 500
      Root <- try(uniroot(Surv, interval = c(1e-05, Up), u = u[i], j = i)$root, TRUE)
    }
    e.time[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  
  
  e.time <- ifelse(is.na(e.time), t.max, e.time) #Observed event time
  
  # Random censoring times 
  c.time <- e.time
  censor <- sample(1:n, floor(cens*n)) # random censored id 
  c.time[censor] <- runif(floor(cens*n), dat[censor,]$t0, pmin(e.time,dat[,ncol(dat)-1])[censor])
  
  #True observed event time, including administrative censoring
  Time <- pmin(e.time,c.time,dat$t3) 
  status <-as.numeric(Time==e.time) # event indicator 
  
  
  # final generated data set 
  
  data_long <- dat.long[which(dat.long$time <= rep(Time,rep(4,n))),]  #deleted not observed data
  
  
  # Longitudinal and survival data
  dat.jm <- data.frame(id = unique(dat$id), Time = Time, status = status)
  
  
  out <- list(longitudinal=data_long[,-2],survival=dat.jm)
  out
}

