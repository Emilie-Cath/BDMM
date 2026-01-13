functions {
  // Softmax of polynomial logits in standardised time
  vector P_poly(real time_raw,
                int K,
                int D,                   // polynomial degree
                matrix beta,             // K x (D+1), includes intercept
                real t_mean,
                real t_sd) {
    real z = (time_raw - t_mean) / t_sd;  // standardised time
    vector[D + 1] basis;
    basis[1] = 1;                          // intercept
    for (d in 2:(D + 1)) basis[d] = pow(z, d - 1);

    vector[K] f;
    for (k in 1:K) {
      f[k] = dot_product(beta[k], basis);
    }

    // centre for softmax identifiability
    f -= mean(f);
    return softmax(f);
  }
  
  
  // interpolation function
  real interpolate(vector t, vector y, real t_interp) {
    int N = num_elements(t);
    int i = 1;

    // Find interval [t[i], t[i+1]] such that t[i] <= t_interp <= t[i+1]
    for (n in 1:(N - 1)) {
      if (t[n] <= t_interp && t_interp <= t[n + 1]) {
        i = n;
        break;
      }
    }
    // Compute linear interpolation
    real weight = (t_interp - t[i]) / (t[i + 1] - t[i]);
    return y[i] + weight * (y[i + 1] - y[i]);
  }

  real Mix(int K, real time, matrix S,vector t_S, vector q, int D, matrix beta, real t_mean, real t_sd) {
    vector[K] p = P_poly(time, K, D, beta, t_mean, t_sd);
    vector[K] numerator;
    for (i in 1:K) {numerator[i]= p[i] * q[i] * interpolate(t_S,S[,i],time);}
    vector[K] denominator = p .* q;
    return sum(numerator) / sum(denominator);
  }

  vector bdmm(real time,
              vector y,
              real lambda,
              int K,
              matrix s1,
              vector t_S1,
              matrix s2,
              vector t_S2,
              vector q1,
              vector q2,
              int D,
              matrix beta, 
              real t_mean, 
              real t_sd
              ) {
    vector[2] derivative;
    derivative[1] = lambda * (-y[1] + Mix(K, time, s1,t_S1, q1, D, beta, t_mean, t_sd));
    derivative[2] = lambda * (-y[2] + Mix(K, time, s2,t_S2, q2, D, beta, t_mean, t_sd));
    return derivative;
  }
}

data {
  int<lower=1> n_obs;        // Number of observations
  int<lower=1> n_fit;        // Number of fitted values
  int<lower=1> K;            // Number of sources
  matrix[n_obs, 2] y_obs;    // Observed isotopes
  vector[n_obs] t_obs;   // Observed times
  array[n_fit] real t_fit;   // Times for predictions
  real t0;                   // Initial time

  // Priors for random effects:
  matrix[n_obs, K] mu_s1; // Observation level source mean values of isotope 1
  matrix[n_obs, K] mu_s2; // Observation level source mean values of isotope 2
  matrix<lower=0>[n_obs, K] sigma_s1; // Obs level SD of sources for isotope 1
  matrix<lower=0>[n_obs, K] sigma_s2; // Obs level SD of sources for isotope 2
  
  // Hyperparameters for lambda prior (Gamma shape–rate)
  real<lower=0> alpha_lambda;
  real<lower=0> beta_lambda;
  
  // Polynomial regression details
  int<lower=0> D;            // Degree e.g. 2 or 3
  real t_mean;               // from R - for standardisation
  real t_sd;                 // from R

  // Concentration dependence
  vector<lower=0, upper=1>[K] q1; // CD for isotope 1
  vector<lower=0, upper=1>[K] q2; // CD for isotope 2
}

parameters {
  vector[2] y0;                      // Initial isotope measurement
  vector<lower=0.01>[2] sigma;       // Measurement error
  matrix[K, D + 1] beta;             // per-class coefficients (incl. intercept)
  real<lower=0> lambda;              // Turnover rate
  matrix[n_obs, K] s1;               // Random effects isotope 1 (obs-level)
  matrix[n_obs, K] s2;               // Random effects isotope 2 (obs-level)
}

model {
  array[n_obs] vector[2] mu_arr;
  matrix[n_obs, 2] mu_vec;

  // Priors:
  sigma ~ normal(0, 1) T[0.01, ];
  y0 ~ normal(0, 10);
  for (d in 1:(D + 1)) {
    real scale = 2.0 / pow(2.0, d - 1);
    col(beta, d) ~ normal(0, scale);
  }

  // Prior on turnover rate
  lambda ~ gamma(alpha_lambda, beta_lambda);

  // Random effects priors (hierarchical):
  for (k in 1:K){
    s1[, k] ~ normal(mu_s1[,k], sigma_s1[,k]);
    s2[, k] ~ normal(mu_s2[,k], sigma_s2[,k]);
  }

  // Likelihood (solve ODE separately for each observation interval):
  for (i in 1:n_obs) {
    
    mu_arr[i] = ode_rk45(bdmm, (i == 1 ? y0 : mu_arr[i - 1]), 
                         (i == 1 ? t0 : t_obs[i - 1]), 
                         {t_obs[i]}, lambda, K, 
                         s1,t_obs, s2,t_obs, q1, q2, D, beta, t_mean, t_sd)[1];

    mu_vec[i, 1] = mu_arr[i, 1];
    mu_vec[i, 2] = mu_arr[i, 2];
  }

  // Observations:
  y_obs[, 1] ~ normal(mu_vec[, 1], sigma[1]);
  y_obs[, 2] ~ normal(mu_vec[, 2], sigma[2]);
}

generated quantities {
  array[n_fit] vector[2] y_fit;
  matrix[n_fit, K] f_fit;
  matrix[n_fit, K] p_fit;

  // For predictions use FIRST source values (mu_s1 and mu_s2):
  // NOTE THIS WILL NEED TO BE CHANGED FOR FINAL MODEL RUN
  matrix[n_obs,K] s1_mean = mu_s1;
  matrix[n_obs,K] s2_mean = mu_s2;

  y_fit = ode_rk45(bdmm, y0, t0, t_fit,
                   lambda, K, s1_mean, t_obs, s2_mean,t_obs, q1, q2,
                   D, beta, t_mean, t_sd);

  // ---- logits under polynomial model ----
  {
    matrix[n_fit, D + 1] B;  // basis: [1, z, z^2, ..., z^D]
    for (i in 1:n_fit) {
      real z = (t_fit[i] - t_mean) / t_sd;
      B[i, 1] = 1;
      for (d in 2:(D + 1)) B[i, d] = pow(z, d - 1);
    }
    f_fit = B * beta';  // n_fit x K

    // If P_poly centres logits before softmax, mirror that here:
    for (i in 1:n_fit) {
      row_vector[K] r = f_fit[i];
      real m = mean(to_vector(r));
      f_fit[i] = r - rep_row_vector(m, K);
    }
  }

  // probabilities
  for (i in 1:n_fit) {
    p_fit[i, ] = to_row_vector(softmax(f_fit[i]'));
  }
}
