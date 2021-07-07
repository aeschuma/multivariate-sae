data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> C; // number of causes
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[C, C] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[C] y[N]; // outcome
}
parameters {
    real<lower=0> sigma_gamma; // standard deviation of the IID REs on region
    vector[C] gamma[R]; // REs on region X cause
    vector[C] beta; // FEs on cause
}
transformed parameters {
    vector[C] gamma_obs[N]; // values of REs for each bivariate normal observation
    matrix[R, C] gamma_tx; // matrix form of gamma for sum to zero constraint specification
    
    for (i in 1:N) {
        gamma_obs[i] = gamma[regions[i]]; // repeat the correct RE for both causes
    }
    for (r in 1:R) {
        gamma_tx[r,] = to_row_vector(gamma[r]);
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(beta + gamma_obs[i], Sigma[i]); // bivariate normal observations
    }
    for (r in 1:R) {
        gamma[r] ~ normal(0, sigma_gamma); // IID normal REs on region
    }
    sum(to_vector(gamma_tx)) ~ normal(0, 0.001 * R * C);  // equivalent to mean(gamma_tx) ~ normal(0,0.001)
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
