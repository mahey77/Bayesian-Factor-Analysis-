functions {
  real log_density_Lambda_diag(vector x, real C_0) {
    int k = num_elements(x);
    real ldens = 0.0;
    for(i in 1:k) {
      ldens += (k - i) * log(x[i]) - square(x[i]) / ((2 * C_0));
    }
    return ldens;
  }
}
data {
  int<lower=1> n; // number of observations
  int<lower=1> p; // number of variables
  int<lower=1, upper=p> k; // number of factors
  matrix[n, p] y; // data
  real<lower=0> a_idio; // shape parameter 1 in prior for idiosyncratic variances
  real<lower=0> b_idio; // shape parameter 2 in prior for idiosyncratic variances
  real<lower=0> C_Lambda; // scale parameter in prior for factor loadings
  int<lower=1> n_star;
  
}
parameters {
  vector<lower=0>[p] diag_Sigma; // idiosyncratic variances
  vector<lower=0>[k] PLT_Lambda_diag; // diagonal elements of Lambda
  vector[p*k-k*(k+1)/2] PLT_Lambda_lowertri;   // below-diagonal elements of Lambda
  matrix[n, k] f; // factors
  matrix[k, k] Gamma; // 
}
transformed parameters {
  matrix[p, k] PLT_Lambda = rep_matrix(0.0, p, k); // PLT Lambda
  {
    int l = 1;
    for(i in 1:p) {
      for(j in 1:min(i, k)) {
        if(i==j) PLT_Lambda[i, j] = PLT_Lambda_diag[i];
        else {
          PLT_Lambda[i, j] = PLT_Lambda_lowertri[l];
          l += 1;
        }
      }
    }
  }
}
model {
  // Likelihood:
  f[1,] ~ std_normal();
  y[1,] ~ normal(PLT_Lambda * f[1,]', sqrt(diag_Sigma));
  for(i in 2:n) {
    f[i,] ~ normal(Gamma*f[i-1,]', rep_vector(1.0, k));
    y[i,] ~ normal(PLT_Lambda * f[i,]', sqrt(diag_Sigma));
  }
  // Prior:
  diag_Sigma ~ inv_gamma((a_idio*0.5), (a_idio*square(b_idio)*0.5));
  PLT_Lambda_lowertri ~ normal(0.0, C_Lambda);
  for(i in 1:k) Gamma[,i] ~ normal(0, 0.5);
  target += log_density_Lambda_diag(PLT_Lambda_diag, C_Lambda);
}

generated quantities {
  matrix[(n_star+1), k] f_star; 
  matrix[n_star, p] y_star;
  f_star[1,] = f[n,];
  for(i in 1:n_star) {
    f_star[i+1,] = multi_normal_rng(Gamma * f_star[i,]', 
    diag_matrix(rep_vector(1.0, k)))';
    y_star[i,] = multi_normal_rng(PLT_Lambda * f_star[i+1,]', 
    diag_matrix(diag_Sigma))';
  }
}



