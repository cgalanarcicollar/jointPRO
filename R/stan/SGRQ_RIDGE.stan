functions{
  // ------------------------------------------------------
    //     LINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
  // ------------------------------------------------------ 
    vector linear_predictor( vector times, int[] ID, vector beta, matrix bi){
      int N = num_elements(times);
      vector[N] out;
      
      out = beta[1] + bi[ID,1] + beta[2]*times + rows_dot_product(bi[ID,2],times);
      
      return out;
    } 
  // ------------------------------------------------------ 
    
    
    // ------------------------------------------------------
    //    BETA-BINOMIAL LOG-LIKELIHOOD               
  // ------------------------------------------------------ 
    real logBB(int y, real mu, real phi, int m){
      real lpdf;
      lpdf = log(choose( m, y)) + lgamma(1/phi) - lgamma(1/phi+m)+ lgamma(mu/phi+y) - lgamma(mu/phi) + lgamma((1-mu)/phi+m-y) - lgamma((1-mu)/phi) ; 
      
      return lpdf;
    } 
  // ------------------------------------------------------ 
    
}


data{
  int N; // number of observations
  int n; // number of ids
  int<lower=0> m1; // the maximum score 
  int<lower=0> m2; // the maximum score 
  int<lower=0> m3; // the maximum score 
  int<lower=0,upper=m1> y1[N]; // longitudinal response 1
  int<lower=0,upper=m2> y2[N]; // longitudinal response 2
  int<lower=0,upper=m3> y3[N]; // longitudinal response 2
  vector[N] times; // measurement times 
  
  int<lower=1,upper=n> ID[N]; // id indicator 
  int<lower=1,upper=N> start[n]; // row for first medition of id=n
  int<lower=1,upper=N> stop[n]; // row for last medition
  
  vector[n] Time; // survival time 
  vector[n] status; // status time
  int K; // GLQ number of points
  vector[K] xk;
  vector[K] wk;
}


parameters{
  vector[2] betas1;
  vector[2] betas2;
  vector[2] betas3;
  vector[3] Alpha;
  real<lower=0> phi1;
  real<lower=0> phi2;
  real<lower=0> phi3;
  cov_matrix[6] Sigma;
  matrix[n,6] bi;
  
  // Weibull baseline hazard parameters  
  real<lower=0> nu;
  real gamma;
  
  // RIDGE REGULARIZATION
  vector<lower=0>[3] tau2;
  real v;
  real s;
  
}


model{
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  1              
  // ------------------------------------------------------
    {
      vector[N] invlogitmu1; 
      vector[N] mu1; 
      
      vector[num_elements(y1)] lBB1;
      
      // Linear predictor
      invlogitmu1 = linear_predictor( times, ID, betas1, bi[,1:2]);
      mu1 = inv_logit(invlogitmu1);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y1) ){
        lBB1[i] = logBB( y1[i], mu1[i], phi1, m1);
      }
      target +=  sum(lBB1);
    }
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  2              
  // ------------------------------------------------------
    {
      vector[N] invlogitmu2; 
      vector[N] mu2; 
      
      vector[num_elements(y2)] lBB2;
      
      // Linear predictor
      invlogitmu2 = linear_predictor( times, ID, betas2, bi[,3:4]);
      mu2 = inv_logit(invlogitmu2);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y2) ){
        lBB2[i] = logBB( y2[i], mu2[i], phi2, m2);
      }
      target +=  sum(lBB2);
    }
  
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  3              
  // ------------------------------------------------------
    {
      vector[N] invlogitmu3; 
      vector[N] mu3; 
      
      vector[num_elements(y3)] lBB3;
      
      // Linear predictor
      invlogitmu3 = linear_predictor( times, ID, betas3, bi[,5:6]);
      mu3 = inv_logit(invlogitmu3);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y3) ){
        lBB3[i] = logBB( y3[i], mu3[i], phi3, m3);
      }
      target +=  sum(lBB3);
    }
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
  // ------------------------------------------------------
    {
      vector[n] haz;
      matrix[n,K] cumHazK;
      vector[n] cumHaz;
      
      for(i in 1:n){
        // Hazard function
        haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha[1] *  inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i])  + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]) );
        
        // Hazard function evaluated at Gauss-Legendre quadrature integration points
        for(j in 1:K){
          cumHazK[i,j] =  nu * (Time[i]/2*(xk[j]+1))^(nu-1) * exp( gamma + Alpha[1] * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[j]+1) ) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i]/2*(xk[j]+1) ) + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]/2*(xk[j]+1) ) );
        }
        
        // Cumulative hazard function with Gauss-Legendre quadrature
        cumHaz[i] = Time[i] / 2 * dot_product(wk, cumHazK[i,]);
        
        target += status[i]*log(haz[i]) - cumHaz[i];
      }
    }   
  // ------------------------------------------------------
    //                       LOG-PRIORS                       
  // ------------------------------------------------------
    // Longitudinal fixed effects
  betas1 ~ normal(0,10);
  betas2 ~ normal(0,10);
  betas3 ~ normal(0,10);
  
  
  // Dispersion beta-binomial parameter
  phi1 ~ cauchy(0, 2);
  phi2 ~ cauchy(0, 2);
  phi3 ~ cauchy(0, 2);
  
  // Baseline hazard Weibull parameters 
  target += inv_gamma_lpdf(nu | 0.1, 0.1);
  target += normal_lpdf(gamma | 0, 10);
  
  // Random-effects variance-covariance matrix
  Sigma ~ inv_wishart(6, diag_matrix(rep_vector(1.0,6)));
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:6] | rep_vector(0.0,6), Sigma); }
  
  // RIDGE REGULARIZATION
  // Association parameter
  Alpha ~ normal(0,tau2);
  tau2 ~ inv_gamma(v/2,(v/2)*s^2);
  1/v ~ uniform(0,1);
  s ~ uniform(0,100);
}

generated quantities{
  vector[n] log_lik;
  vector[N] longit;
  vector[N] linpred1 = linear_predictor( times, ID, betas1, bi[,1:2]);
  vector[N] linpred2 = linear_predictor( times, ID, betas2, bi[,3:4]);
  vector[N] linpred3 = linear_predictor( times, ID, betas3, bi[,5:6]);
  vector[n] haz;
  matrix[n, K] cumHazK;
  vector[n] cumHaz;
  
  for(j in 1:N){
    longit[j] = logBB( y1[j], inv_logit(linpred1[j]), phi1, m1)+ 
                logBB( y2[j], inv_logit(linpred2[j]), phi2, m2)+
                logBB( y3[j], inv_logit(linpred3[j]), phi3, m3);
  }
  
  // SURVIVAL PART
  for(i in 1:n){
    haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha[1] *  inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i])  + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]) );
    
    // Hazard function at integration points
    for(k in 1:K){
      cumHazK[i, k] = nu * pow(Time[i] / 2 * (xk[k]+1), nu-1) * exp( gamma + 
                                                                       Alpha[1] * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[k]+1) ) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i]/2*(xk[k]+1) ) + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]/2*(xk[k]+1) ) );
    }
    
    cumHaz[i] = Time[i] / 2 * dot_product(wk, cumHazK[i,]);
    
    log_lik[i] = sum(longit[start[i]:stop[i]]) + status[i]*log(haz[i]) - cumHaz[i];
  }
}