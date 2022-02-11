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
    matrix[R, C] alpha; // REs on region X cause
    vector[C] beta; // FEs on cause
}
transformed parameters {
    matrix[N, C] alpha_obs; // values of REs for each bivariate normal observation

    for (i in 1:N) {
        for (c in 1:C) {
             alpha_obs[i, c] = alpha[regions[i], c]; // repeat the correct RE for both causes
        }
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(alpha_obs[i, ]), Sigma[i]); // bivariate normal observations
    }
    for (c in 1:C) {
        alpha[, c] ~ normal(beta[c], sigma_gamma[c]); // IID normal REs on region
    }
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
