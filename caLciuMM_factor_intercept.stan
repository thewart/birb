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
  // int<lower=1> P;
  matrix[M,N] Y;
  // matrix[P,N] X;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[M] alpha_raw;
  real mu_alpha;
  real<lower=0> sigma_alpha;
  
  real<lower=0> sigma_lambda;
  vector[(M-1)*K - choose(K,2)] lambda_lower;
  vector<lower=0>[K] lambda_diag;
  matrix[K,N] u;
  
  real mu_sigma;
  real sigma_sigma;
  vector[M] sigma_raw;
}

transformed parameters {
  vector[M] alpha = mu_alpha + sigma_alpha * alpha_raw;
  vector[M] sigma = exp(mu_sigma + sigma_sigma * sigma_raw);
  matrix[M,K] Lambda = factor_loadings(sigma_lambda*lambda_lower, sigma_lambda*lambda_diag);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[M,N] eta = rep_matrix(alpha,N) + Lambda * u;
  matrix[M,N] Sigma = rep_matrix(sigma,N);
  
  to_vector(Y) ~ normal(to_vector(eta) , to_vector(Sigma));
  to_vector(u) ~ normal(0,1);
  alpha_raw ~ normal(0,1);
  lambda_lower ~ normal(0,1);
  lambda_diag ~ normal(0,1);
  
  mu_alpha ~ normal(0,1);
  sigma_alpha ~ normal(0,1);
  sigma_lambda ~ normal(0,1);
  
  sigma_raw ~ normal(0,1);
  mu_sigma ~ normal(0,1);
  sigma_sigma ~ normal(0,1);
}
