functions {
  matrix factor_loadings(vector lambda_lower, vector lambda_diag) {
    int K = num_elements(lambda_diag);
    int M = (num_elements(lambda_lower) + choose(K,2))/K + 1;
    matrix[M,K] Lambda;
    int i=1;
    for (m in 1:M) {
      for (n in 1:K) {
        if (m == n) {
          Lambda[m, n] = lambda_diag[m];
        } else if (m > n) {
          Lambda[m, n] = lambda_lower[i];
          i += 1;
        } else if (m < n) {
          Lambda[m, n] = 0;
        }
      }
    }
    
    return(Lambda);
  }
  
}

data {
  int<lower=0> N;
  int<lower=2> M;
  int<lower=1> K;
  int<lower=1> P;
  int<lower=1> R;
  matrix[M,N] Y;
  matrix[P*R,N] X;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix[P,R] beta_raw[M];
  matrix[P,R] mu_beta;
  matrix<lower=0>[P,R] sigma_beta;
  cholesky_factor_corr[P] L_p;
  cholesky_factor_corr[R] L_r;
  
  real<lower=0> sigma_lambda;
  vector[(M-1)*K - choose(K,2)] lambda_lower;
  vector<lower=0>[K] lambda_diag;
  matrix[K,N] u;
  
  real mu_sigma;
  real sigma_sigma;
  vector[M] sigma_raw;
}

transformed parameters {
  matrix[P*R,M] beta;
  vector[M] sigma = exp(mu_sigma + sigma_sigma * sigma_raw);
  matrix[M,K] Lambda = factor_loadings(sigma_lambda*lambda_lower, sigma_lambda*lambda_diag);
  
  for (i in 1:M) { 
    beta[,i] = to_vector(mu_beta + L_p * beta_raw[i] * L_r');
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[M,N] eta = X' * beta + Lambda * u;
  matrix[M,N] Sigma = rep_matrix(sigma,N);
  
  to_vector(Y) ~ normal(to_vector(eta) , to_vector(Sigma));
  to_vector(u) ~ normal(0,1);
  for (i in 1:M) to_vector(beta_raw[i]) ~ normal(0,1);
  lambda_lower ~ normal(0,1);
  lambda_diag ~ normal(0,1);
  
  to_vector(mu_beta) ~ normal(0,1);
  to_vector(sigma_beta) ~ normal(0,1);
  L_p ~ lkj_corr_cholesky(2);
  L_r ~ lkj_corr_cholesky(2);
  
  sigma_lambda ~ normal(0,1);
  
  sigma_raw ~ normal(0,1);
  mu_sigma ~ normal(0,1);
  sigma_sigma ~ normal(0,1);
}

generated quantities {
  real log_lik[M*N];
  vector[M*N] y = to_vector(Y);
  vector[M*N] Sigma = to_vector(rep_matrix(sigma,N));
  vector[M*N] eta = to_vector(X' * beta + Lambda * u);
  
  for (i in 1:(M*N)) {
    log_lik[i] = normal_lpdf(y[i] | eta[i],Sigma[i]);
  }
}
