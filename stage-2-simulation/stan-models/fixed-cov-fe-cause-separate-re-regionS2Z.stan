data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> C; // number of causes
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[C, C] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[C] y[N]; // outcome
}
parameters {
    real<lower=0> sigma_gamma[C]; // standard deviation of the IID REs on region
    matrix[R, C] gamma; // REs on region X cause
    vector[C] beta; // FEs on cause
}
transformed parameters {
    matrix[N, C] gamma_obs; // values of REs for each bivariate normal observation

    for (i in 1:N) {
        for (c in 1:C) {
             gamma_obs[i, c] = gamma[regions[i], c]; // repeat the correct RE for both causes
        }
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(beta + to_vector(gamma_obs[i, ]), Sigma[i]); // bivariate normal observations
    }
    for (c in 1:C) {
        gamma[, c] ~ normal(0, sigma_gamma[c]); // IID normal REs on region
        sum(gamma[, c]) ~ normal(0, 0.001 * R);  // equivalent to mean(gamma) ~ normal(0,0.001)
    }
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
