MAP_fc <- function(x){
  
  lim.inf <- min(x)-1
  lim.sup <- max(x)+1
  s <- density(x,from=lim.inf,to=lim.sup,bw=0.2)
  n <- length(s$y)
  v1 <- s$y[1:(n-2)]
  v2 <- s$y[2:(n-1)]
  v3 <- s$y[3:n]
  ix <- 1+which((v1<v2)&(v2>v3))
  out <- s$x[which(s$y==max(s$y))]
  
  return( out )
  
}

update_bi <- function(data,m, fit, k, iter=1000){
  
  y <- data$longitudinal$y #PRO

  
  times <- data$longitudinal$time #measurement times
  Time <- data$survival$Time #survival time
  status <- data$survival$status #survival status
  
  # Parameters
  # LONGITUDINAL
  betas <- rstan::extract(fit, "betas1")$betas #longitudinal fix eff
  Sigma <- rstan::extract(fit, "Sigma")$Sigma #Variance matrix 
  phi <- rstan::extract(fit, "phi1")$phi #dispersion parameter 
  
  # SURVIVAL
  Alpha <- rstan::extract(fit, "Alpha")$Alpha # asociation parameter
  nu <- rstan::extract(fit, "nu")$nu     # weibul parameter 1
  gamma <- rstan::extract(fit, "gamma")$gamma    # weibul parameter 2
  
  # STEP 1: SIMULATE THE FIXED PARAMETERS ACCORDING TO MODEL ESTIMATIONS 
  prop_betas <- matrix(NA, nrow = iter, ncol = 2)
  prop_Sigma <- array(NA, dim = c(iter, 2, 2))
  prop_phi <- c()
  prop_Alpha <- c()
  prop_nu <- c()
  prop_gamma <- c()
  
  MAP_Sigma <- diag(apply(Sigma, c(2,3), MAP_fc))
  bi <- matrix(0, nrow=iter, ncol=2)
  bi_prev <- rep(0, 2)
  
  #para cada uno de los theta simulados, simulamos el bi correspondiente
  for(j in 1:iter){
    pos <- sample(1: 2000, size=1)
    
    prop_betas[j,] <- betas[pos,] 
    

    prop_Sigma[j,,] <- Sigma[pos,,]
    prop_phi[j] <- phi[pos]
    prop_Alpha[j] <- Alpha[pos]
    prop_nu[j] <- nu[pos]
    prop_gamma[j] <- gamma[pos]
    # full conditional posterior function
    p.log <- function(b){
      
      # Longitudinal
      lp_y <- (prop_betas[j,1] + b[1]) + (prop_betas[j,2] + b[2])*times 
      mu_y <- exp(lp_y)/(1+exp(lp_y))
      log_like_yj <- c()
      for (i in 1:length(mu_y)) {
        # log_like_yj[i] <- extraDistr::dbbinom(y[i], m, mu_y[i]*prop_phi[j], prop_phi[j]*(1-mu_y[i]), log = TRUE)
        log_like_yj[i] <- log(choose( m, y[i])) + lgamma(1/prop_phi[j]) - lgamma(1/prop_phi[j]+m)+ lgamma(mu_y[i]/prop_phi[j]+y[i]) - lgamma(mu_y[i]/prop_phi[j]) + lgamma((1-mu_y[i])/prop_phi[j]+m-y[i]) - lgamma((1-mu_y[i])/prop_phi[j]) ;
      }
      log_like_y <- sum(log_like_yj)
      
      log_re_y <- emdbook::dmvnorm(b[1:2], mu=c(0,0), Sigma=prop_Sigma[j,,], log=TRUE)
      
      # Survival
      glq <- statmod::gauss.quad(k, kind = "legendre")
      xk <- glq$nodes   # nodes
      wk <- glq$weights # weights
      K <- length(xk)   # K-points
      
      loghaz =  prop_nu[j] + (prop_nu[j]-1) * log(Time) + prop_gamma[j] + 
        prop_Alpha[j] *   exp((prop_betas[j,1] + b[1]) + (prop_betas[j,2] + b[2])*Time)/(1+exp((prop_betas[j,1] + b[1]) + (prop_betas[j,2] + b[2])*Time));
      
      cumHazK<-c()
      for(K in 1:k){
        cumHazK[K] =  prop_nu[j] * (Time/2*(xk[K]+1))^(prop_nu[j]-1) *
          exp( prop_gamma[j] + prop_Alpha[j] *   exp((prop_betas[j,1] + b[1]) + (prop_betas[j,2] + b[2])*Time/2*(xk[K]+1))/(1+exp((prop_betas[j,1] + b[1]) + (prop_betas[j,2] + b[2])*Time/2*(xk[K]+1))));
      }
      
      cumHaz = Time/ 2 * pracma::dot(wk, cumHazK);
      log_like_surv <- status*loghaz - cumHaz
      
      
      # Joint likelihood
      return( log_like_y + 
                log_re_y +
                log_like_surv)
    }
    
    # Adaptive MCMC 
    invisible(capture.output( bb <- adaptMCMC::MCMC(p.log, n=1000, init=bi_prev, scale=MAP_Sigma,
                                         adapt=TRUE, acc.rate=0.234, showProgressBar=T) )) # try with optim 
    
    
    bi[j,] <- apply(bb$samples[-c(1:500),], 2, median) # comentar si optim
    
    bi_prev <- as.vector(bi[j,] )
  }
  
  
  out <- data.frame(betas=prop_betas, Sigma=prop_Sigma, bi=bi[,1:2], phi=prop_phi,
                    Alpha=prop_Alpha, gamma=prop_gamma, nu = prop_nu)
  
  return( out )
  
}


