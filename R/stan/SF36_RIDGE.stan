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
  int<lower=0> m4; // the maximum score 
  int<lower=0> m5; // the maximum score 
  int<lower=0> m6; // the maximum score 
  int<lower=0> m7; // the maximum score 
  int<lower=0> m8; // the maximum score 
  
  
  int<lower=0,upper=m1> y1[N]; // longitudinal response 1
  int<lower=0,upper=m2> y2[N]; // longitudinal response 2
  int<lower=0,upper=m3> y3[N]; // longitudinal response 3
  int<lower=0,upper=m4> y4[N]; // longitudinal response 4
  int<lower=0,upper=m5> y5[N]; // longitudinal response 5
  int<lower=0,upper=m6> y6[N]; // longitudinal response 6
  int<lower=0,upper=m7> y7[N]; // longitudinal response 7
  int<lower=0,upper=m8> y8[N]; // longitudinal response 8
  
  int<lower=1,upper=N> start[n]; // row for first medition of id=n
  int<lower=1,upper=N> stop[n]; // row for last medition  
  vector[N] times; // measurement times 
  int<lower=1,upper=n> ID[N]; // id indicator 
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
  vector[2] betas4;
  vector[2] betas5;
  vector[2] betas6;
  vector[2] betas7;
  vector[2] betas8;
  
  vector[8] Alpha;
  
  real<lower=0> phi1;
  real<lower=0> phi2;
  real<lower=0> phi3;
  real<lower=0> phi4;
  real<lower=0> phi5;
  real<lower=0> phi6;
  real<lower=0> phi7;
  real<lower=0> phi8;
  
  cov_matrix[16] Sigma;
  matrix[n,16] bi;
  
  // Weibull baseline hazard parameters  
  real<lower=0> nu;
  real gamma;
  
  // RIDGE REGULARIZATION
  vector<lower=0>[8] tau2;
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
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  4              
  // ------------------------------------------------------
    {
      vector[N] invlogitmu4; 
      vector[N] mu4; 
      
      vector[num_elements(y4)] lBB4;
      
      // Linear predictor
      invlogitmu4 = linear_predictor( times, ID, betas4, bi[,7:8]);
      mu4 = inv_logit(invlogitmu4);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y4) ){
        lBB4[i] = logBB( y4[i], mu4[i], phi4, m4);
      }
      target +=  sum(lBB4);
    }
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  5             
  // ------------------------------------------------------
    {
      vector[N] invlogitmu5; 
      vector[N] mu5; 
      
      vector[num_elements(y5)] lBB5;
      
      // Linear predictor
      invlogitmu5 = linear_predictor( times, ID, betas5, bi[,9:10]);
      mu5 = inv_logit(invlogitmu5);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y5) ){
        lBB5[i] = logBB( y5[i], mu5[i], phi5, m5);
      }
      target +=  sum(lBB5);
    }
  
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  6             
  // ------------------------------------------------------
    {
      vector[N] invlogitmu6; 
      vector[N] mu6; 
      
      vector[num_elements(y6)] lBB6;
      
      // Linear predictor
      invlogitmu6 = linear_predictor( times, ID, betas6, bi[,11:12]);
      mu6 = inv_logit(invlogitmu6);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y6) ){
        lBB6[i] = logBB( y6[i], mu6[i], phi6, m6);
      }
      target +=  sum(lBB6);
    }
  
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  7            
  // ------------------------------------------------------
    {
      vector[N] invlogitmu7; 
      vector[N] mu7; 
      
      vector[num_elements(y7)] lBB7;
      
      // Linear predictor
      invlogitmu7 = linear_predictor( times, ID, betas7, bi[,13:14]);
      mu7 = inv_logit(invlogitmu7);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y7) ){
        lBB7[i] = logBB( y7[i], mu7[i], phi7, m7);
      }
      target +=  sum(lBB7);
    }
  
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL  8             
  // ------------------------------------------------------
    {
      vector[N] invlogitmu8; 
      vector[N] mu8; 
      
      vector[num_elements(y8)] lBB8;
      
      // Linear predictor
      invlogitmu8 = linear_predictor( times, ID, betas8, bi[,15:16]);
      mu8 = inv_logit(invlogitmu8);
      
      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:num_elements(y8) ){
        lBB8[i] = logBB( y8[i], mu8[i], phi8, m8);
      }
      target +=  sum(lBB8);
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
        haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha[1] *  inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i])  + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]) + Alpha[4] * inv_logit(betas4[1] + bi[i,7] + (betas4[2] + bi[i,8])*Time[i]) + Alpha[5] * inv_logit(betas5[1] + bi[i,9] + (betas5[2] + bi[i,10])*Time[i])+ Alpha[6] * inv_logit(betas6[1] + bi[i,11] + (betas6[2] + bi[i,12])*Time[i])+ Alpha[7] * inv_logit(betas7[1] + bi[i,13] + (betas7[2] + bi[i,14])*Time[i])+ Alpha[8] * inv_logit(betas8[1] + bi[i,15] + (betas8[2] + bi[i,16])*Time[i]));
        
        // Hazard function evaluated at Gauss-Legendre quadrature integration points
        for(j in 1:K){
          cumHazK[i,j] =  nu * (Time[i]/2*(xk[j]+1))^(nu-1) * exp( gamma + Alpha[1] * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[j]+1) ) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i]/2*(xk[j]+1) ) + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]/2*(xk[j]+1) ) + Alpha[4] * inv_logit(betas4[1] + bi[i,7] + (betas4[2] + bi[i,8])*Time[i]/2*(xk[j]+1) ) + Alpha[5] * inv_logit(betas5[1] + bi[i,9] + (betas5[2] + bi[i,10])*Time[i]/2*(xk[j]+1) )+ Alpha[6] * inv_logit(betas6[1] + bi[i,11] + (betas6[2] + bi[i,12])*Time[i]/2*(xk[j]+1) ) + Alpha[7] * inv_logit(betas7[1] + bi[i,13] + (betas7[2] + bi[i,14])*Time[i]/2*(xk[j]+1) )+ Alpha[8] * inv_logit(betas8[1] + bi[i,15] + (betas8[2] + bi[i,16])*Time[i]/2*(xk[j]+1) ));
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
  betas4 ~ normal(0,10);
  betas5 ~ normal(0,10);
  betas6 ~ normal(0,10);
  betas7 ~ normal(0,10);
  betas8 ~ normal(0,10);
  
  
  // Dispersion beta-binomial parameter
  phi1 ~ cauchy(0, 2);
  phi2 ~ cauchy(0, 2);
  phi3 ~ cauchy(0, 2);
  phi4 ~ cauchy(0, 2);
  phi5 ~ cauchy(0, 2);
  phi6 ~ cauchy(0, 2);
  phi7 ~ cauchy(0, 2);
  phi8 ~ cauchy(0, 2);
  
  // Baseline hazard Weibull parameters 
  target += inv_gamma_lpdf(nu | 0.1, 0.1);
  target += normal_lpdf(gamma | 0, 10);
  
  // Random-effects variance-covariance matrix
  Sigma ~ inv_wishart(16, diag_matrix(rep_vector(1.0,16)));
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:16] | rep_vector(0.0,16), Sigma); }
  
  // RIDGE REGULARIZATION
  // Association parameter
  Alpha ~ normal(0,tau2);
  tau2 ~ inv_gamma(v/2,(v/2)*s^2);
  1/v ~ uniform(0,1);
  s ~ uniform(0,100);
}

generated quantities{
  vector[N] linpred1 = linear_predictor( times, ID, betas1, bi[,1:2]);
  vector[N] linpred2 = linear_predictor( times, ID, betas2, bi[,3:4]);
  vector[N] linpred3 = linear_predictor( times, ID, betas3, bi[,5:6]);   
  vector[N] linpred4 = linear_predictor( times, ID, betas4, bi[,7:8]);
  vector[N] linpred5 = linear_predictor( times, ID, betas5, bi[,9:10]);
  vector[N] linpred6 = linear_predictor( times, ID, betas6, bi[,11:12]);
  vector[N] linpred7 = linear_predictor( times, ID, betas7, bi[,13:14]);
  vector[N] linpred8 = linear_predictor( times, ID, betas8, bi[,15:16]);
  
  vector[n] haz;
  matrix[n,K] cumHazK;
  vector[n] cumHaz;
  
  vector[n] log_lik;
  vector[N] longit;

  
  
  
  for(j in 1:N){
    longit[j] = logBB( y1[j], inv_logit(linpred1[j]), phi1, m1)+ 
                logBB( y2[j], inv_logit(linpred2[j]), phi2, m2)+
                logBB( y3[j], inv_logit(linpred3[j]), phi3, m3)+
                logBB( y4[j], inv_logit(linpred4[j]), phi4, m4)+ 
                logBB( y5[j], inv_logit(linpred5[j]), phi5, m5)+
                logBB( y6[j], inv_logit(linpred6[j]), phi6, m6)+
                logBB( y7[j], inv_logit(linpred7[j]), phi7, m7)+
                logBB( y8[j], inv_logit(linpred8[j]), phi8, m8);
  }

  
  
  for(i in 1:n){
    // Hazard function
    haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha[1] *  inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i])  + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]) + Alpha[4] * inv_logit(betas4[1] + bi[i,7] + (betas4[2] + bi[i,8])*Time[i]) + Alpha[5] * inv_logit(betas5[1] + bi[i,9] + (betas5[2] + bi[i,10])*Time[i])+ Alpha[6] * inv_logit(betas6[1] + bi[i,11] + (betas6[2] + bi[i,12])*Time[i])+ Alpha[7] * inv_logit(betas7[1] + bi[i,13] + (betas7[2] + bi[i,14])*Time[i])+ Alpha[8] * inv_logit(betas8[1] + bi[i,15] + (betas8[2] + bi[i,16])*Time[i]));
    
    // Hazard function evaluated at Gauss-Legendre quadrature integration points
    for(j in 1:K){
      cumHazK[i,j] =  nu * (Time[i]/2*(xk[j]+1))^(nu-1) * exp( gamma + Alpha[1] * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[j]+1) ) + Alpha[2] * inv_logit(betas2[1] + bi[i,3] + (betas2[2] + bi[i,4])*Time[i]/2*(xk[j]+1) ) + Alpha[3] * inv_logit(betas3[1] + bi[i,5] + (betas3[2] + bi[i,6])*Time[i]/2*(xk[j]+1) ) + Alpha[4] * inv_logit(betas4[1] + bi[i,7] + (betas4[2] + bi[i,8])*Time[i]/2*(xk[j]+1) ) + Alpha[5] * inv_logit(betas5[1] + bi[i,9] + (betas5[2] + bi[i,10])*Time[i]/2*(xk[j]+1) )+ Alpha[6] * inv_logit(betas6[1] + bi[i,11] + (betas6[2] + bi[i,12])*Time[i]/2*(xk[j]+1) ) + Alpha[7] * inv_logit(betas7[1] + bi[i,13] + (betas7[2] + bi[i,14])*Time[i]/2*(xk[j]+1) )+ Alpha[8] * inv_logit(betas8[1] + bi[i,15] + (betas8[2] + bi[i,16])*Time[i]/2*(xk[j]+1) ));
    }
    
    
    // Cumulative hazard function with Gauss-Legendre quadrature
    cumHaz[i] = Time[i] / 2 * dot_product(wk, cumHazK[i,]);
    
    // Add the survival part for each id to their long recorded previously
    
    log_lik[i] = sum(longit[start[i]:stop[i]]) + status[i]*log(haz[i]) - cumHaz[i];
    
    
  }
  
}