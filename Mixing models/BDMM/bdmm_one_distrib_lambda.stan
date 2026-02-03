functions {
  // OPTIMIZATION 1: Efficient Polynomial Basis
  // Avoids repeated expensive pow() calls by using recursive multiplication
  vector P_poly(real time_raw, int K, int D, matrix beta, real t_mean, real t_sd) {
    real z = (time_raw - t_mean) / t_sd;
    vector[D + 1] basis;
    basis[1] = 1; 
    for (d in 2:(D + 1)) {
      basis[d] = basis[d-1] * z;
    }

    vector[K] f = beta * basis; 

    // REMOVED: f -= mean(f); <--- Not needed with reference category
    return softmax(f);
  }
  
  // Interpolation (Keep as is, but ensure t is sorted in R)
  real interpolate(vector t, vector y, real t_interp) {
    int N = num_elements(t);
    int i = 1;
    
    // Optimization: If t_interp is exactly the last element, return it to avoid out-of-bounds
    if(t_interp >= t[N]) return y[N];
    if(t_interp <= t[1]) return y[1];

    for (n in 1:(N - 1)) {
      if (t[n] <= t_interp && t_interp <= t[n + 1]) {
        i = n;
        break;
      }
    }
    real weight = (t_interp - t[i]) / (t[i + 1] - t[i]);
    return y[i] + weight * (y[i + 1] - y[i]);
  }

  // Faster Mix calculation
  real Mix(int K, real time, matrix S, vector t_S, vector q, int D, matrix beta, real t_mean, real t_sd) {
    vector[K] p = P_poly(time, K, D, beta, t_mean, t_sd);
    vector[K] numerator;
    // We can extract the specific column we need for interpolation here to avoid matrix slicing overhead
    for (k in 1:K) {
       numerator[k] = p[k] * q[k] * interpolate(t_S, col(S, k), time);
    }
    // Dot product is faster than sum(elementwise_product)
    real denom = dot_product(p, q);
    return sum(numerator) / denom;
  }

  // ODE System
  vector bdmm(real time, vector y, 
              // Parameters packed for ODE solver
              real lambda, int K, matrix s1, vector t_S1, matrix s2, vector t_S2, 
              vector q1, vector q2, int D, matrix beta, real t_mean, real t_sd) {
    
    vector[2] derivative;
    // Calculate mixes once
    real mix1 = Mix(K, time, s1, t_S1, q1, D, beta, t_mean, t_sd);
    real mix2 = Mix(K, time, s2, t_S2, q2, D, beta, t_mean, t_sd);
    
    // dy/dt = lambda * (Equilibrium - Current)
    derivative[1] = lambda * (mix1 - y[1]);
    derivative[2] = lambda * (mix2 - y[2]);
    
    return derivative;
  }
}

data {
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=0> max_num_steps;
  int<lower=1> n_obs;        
  int<lower=1> n_fit;        
  int<lower=1> K;            
  matrix[n_obs, 2] y_obs;    
  vector[n_obs] t_obs;       
  array[n_fit] real t_fit;   
  real t0;                   
  real<lower=0> lambda;      

  // Priors for random effects:
  matrix[n_obs, K] mu_s1; 
  matrix[n_obs, K] mu_s2; 
  matrix<lower=0>[n_obs, K] sigma_s1; 
  matrix<lower=0>[n_obs, K] sigma_s2; 
  
  int<lower=0> D;            
  real t_mean;               
  real t_sd;                 

  vector<lower=0, upper=1>[K] q1; 
  vector<lower=0, upper=1>[K] q2; 
}

parameters {
  vector[2] y0;                   
  vector<lower=0.01>[2] sigma;        
  matrix[K - 1, D + 1] beta_raw; // Estimate K-1 sources
  
  // STRATEGY 2: Non-Centered Parameterization (NCP)
  // We sample from a standard normal (z-scores) and transform later
  matrix[n_obs, K] s1_raw;        
  matrix[n_obs, K] s2_raw;        
}

transformed parameters {
  // Reconstruct the actual s1 and s2 values
  // s = mu + raw * sigma
  matrix[n_obs, K] s1 = mu_s1 + s1_raw .* sigma_s1;
  matrix[n_obs, K] s2 = mu_s2 + s2_raw .* sigma_s2;
  matrix[K, D + 1] beta;
  for (d in 1:(D + 1)) {
    beta[1, d] = 0; // Pin the first source to zero
    for (k in 2:K) {
      beta[k, d] = beta_raw[k - 1, d];
    }
  }
}

model {
  // Priors
  sigma ~ normal(0, 1); // Truncation handled by lower bound in parameters block
  y0 ~ normal(0, 10);
  
  // Regularization for Polynomial coefficients
  for (d in 1:(D + 1)) {
    real scale = 2.0 / pow(2.0, d - 1);
    col(beta_raw, d) ~ normal(0, scale);
  }
  
  // NCP Priors: These must be Standard Normal
  to_vector(s1_raw) ~ std_normal();
  to_vector(s2_raw) ~ std_normal();

  // STRATEGY 1: Vectorized ODE Solver
  array[n_obs] vector[2] mu_arr;
  
  // FIX: wrap t_obs in to_array_1d() for the 4th argument
  mu_arr = ode_bdf_tol(bdmm, y0, t0, to_array_1d(t_obs), 
                       rel_tol, abs_tol, max_num_steps, 
                       lambda, K, s1, t_obs, s2, t_obs, q1, q2, D, beta, t_mean, t_sd);
                    
                    
  // Convert array of vectors to matrix for vectorised likelihood
  matrix[n_obs, 2] mu_vec;
  for(i in 1:n_obs) {
    mu_vec[i] = mu_arr[i]';
  }

  // Vectorized Likelihood
  y_obs[, 1] ~ normal(mu_vec[, 1], sigma[1]);
  y_obs[, 2] ~ normal(mu_vec[, 2], sigma[2]);
}

generated quantities {
  // Predictions work exactly the same way, solving the ODE once
  array[n_fit] vector[2] y_fit;
  matrix[n_fit, K] f_fit;
  matrix[n_fit, K] p_fit;

  // Use the transformed parameters s1 and s2 (or means if you prefer mean predictions)
  // To match your original code which used means:
  matrix[n_obs,K] s1_pred = mu_s1;
  matrix[n_obs,K] s2_pred = mu_s2;

  y_fit = ode_bdf_tol(bdmm, y0, t0, to_array_1d(t_fit), 
                       rel_tol, abs_tol, max_num_steps, 
                       lambda, K, s1, t_obs, s2, t_obs, q1, q2, D, beta, t_mean, t_sd);


  // Basis calculation optimization in GQ
  matrix[n_fit, D + 1] B;
  for (i in 1:n_fit) {
    real z = (t_fit[i] - t_mean) / t_sd;
    B[i, 1] = 1;
    for (d in 2:(D + 1)) B[i, d] = B[i, d-1] * z;
  }
  
  f_fit = B * beta'; 

  // Efficient softmax calculation
  for (i in 1:n_fit) {
    vector[K] r = f_fit[i]';
    f_fit[i] = (r - mean(r))'; // Center
    p_fit[i, ] = to_row_vector(softmax(f_fit[i]'));
  }
}
