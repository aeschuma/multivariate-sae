data {
    int<lower=1> R; // number of observations
    real<lower = 0> Sigma[R]; // fixed block variances for each observation 
    vector[R] y; // outcome
    int<lower=0> N_edges;
    int<lower=1, upper=R> node1[N_edges];
    int<lower=1, upper=R> node2[N_edges];
    real<lower=0> scaling_factor; // scales the variance of the spatial effects
    real<lower=0> rho_beta_1; // beta dist hyperpar 1 for % spatial paramter
    real<lower=0> rho_beta_2; // beta dist hyperpar 2 for % spatial paramter
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    real beta; // FE
    vector[R] v; // RE on region, centered on beta1
    vector[R] u;
    real<lower=0> sigma;        // overall standard deviation
    real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
}
transformed parameters {
  vector[R] convolved_re;
  
  // variance of each component should be approximately equal to 1
  convolved_re =  (sqrt(1 - rho) * v) + (sqrt(rho / scaling_factor) * u);
}
model {

    y ~ normal(beta + convolved_re * sigma, Sigma); // bivariate normal observations
    
    // the following computes the prior on phi on the unit scale with sd = 1
    target += -0.5 * dot_self(u[node1] - u[node2]);
    
    // soft sum-to-zero constraint on phi)
    sum(u) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    
    rho ~ beta(rho_beta_1, rho_beta_2);
    sigma ~ normal(0, sigma_normal_sd);
    
    v ~ normal(0, 1); // IID normal REs on region
    
    // beta ~ normal(0,5);
}
generated quantities {
    real log_sigma;
    vector[R] preds;
    real log_lik[R]; // log likelihood for use with WAIC or PSIS-LOO
   
    log_sigma = log(sigma);
    preds = beta + convolved_re * sigma;
    
    for (r in 1:R) {
        log_lik[r] = normal_lpdf(y[r] | beta + convolved_re[r] * sigma, Sigma[r]);
    }
}
