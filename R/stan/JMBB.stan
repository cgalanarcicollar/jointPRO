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
  int<lower=0,upper=m1> y1[N]; // longitudinal response
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
  real Alpha;
  real<lower=0> phi1;
  cov_matrix[2] Sigma;
  matrix[n,2] bi;
  
  // Weibull baseline hazard parameters  
  real<lower=0> nu;
  real gamma;
}


model{
  // ------------------------------------------------------
    //        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
  // ------------------------------------------------------
    {
      vector[N] invlogitmu; 
      vector[N] mu; 
      
      vector[N] lBB;
      
      // Linear predictor
      invlogitmu = linear_predictor( times, ID, betas1, bi);
      mu = inv_logit(invlogitmu);
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:N){
        lBB[i] = logBB( y1[i], mu[i], phi1, m1);
      }
      target +=  sum(lBB);
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
        haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha * m1 * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) );
        
        // Hazard function evaluated at Gauss-Legendre quadrature integration points
        for(j in 1:K){
          cumHazK[i,j] =  nu * (Time[i]/2*(xk[j]+1))^(nu-1) * exp( gamma + Alpha * m1 * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[j]+1) ) );
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
  
  
  // Association parameter
  target += normal_lpdf(Alpha | 0, 10);
  
  // Dispersion beta-binomial parameter
  phi1 ~ cauchy(0, 2);
  
  // Baseline hazard Weibull parameters 
  target += inv_gamma_lpdf(nu | 0.1, 0.1);
  target += normal_lpdf(gamma | 0, 10);
  
  // Random-effects variance-covariance matrix
  Sigma ~ inv_wishart(2, diag_matrix(rep_vector(1.0,2)));
  
  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:2] | rep_vector(0.0,2), Sigma); }
  
}

generated quantities{
  vector[N] linpred1 = linear_predictor( times, ID, betas1, bi[,1:2]);
  vector[N]  mu = inv_logit(linpred1);
      
  matrix[n,K] cumHazK;
  vector[n] cumHaz;
   vector[n] haz;
  vector[n] log_lik;
  vector[N] longit;


   // Linear predictor

      
      
      // Longitudinal Beta-Binomial  log-likelihood
      for(i in 1:N ){
        longit[i] = logBB( y1[i], mu[i], phi1, m1);
      }
      
      
     for(i in 1:n){
        // Hazard function
        haz[i] =  nu * Time[i]^( nu-1) * exp(gamma + Alpha * m1 * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]) );
        
        // Hazard function evaluated at Gauss-Legendre quadrature integration points
        for(j in 1:K){
          cumHazK[i,j] =  nu * (Time[i]/2*(xk[j]+1))^(nu-1) * exp( gamma + Alpha * m1 * inv_logit(betas1[1] + bi[i,1] + (betas1[2] + bi[i,2])*Time[i]/2*(xk[j]+1) ) );
        }
        
        // Cumulative hazard function with Gauss-Legendre quadrature
        cumHaz[i] = Time[i] / 2 * dot_product(wk, cumHazK[i,]);
     }
      for(j in 1:n){ log_lik[j] = sum(longit[start[j]:stop[j]]) + status[j]*log(haz[j]) - cumHaz[j]; }  


}