data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> C; // number of causes
    int regions[N]; // index for which region is for which bivariate observation
    matrix[C, C] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[C] y[N]; // outcome
}
parameters {
    real<lower=0> sigma_gamma; // standard deviation of the IID REs on region
    real gamma[R]; // REs on region
    vector[C] beta; // FEs on cause
}
transformed parameters {
    vector[N] gamma_obs; // values of REs for each bivariate normal observation

    for (i in 1:N) {
        gamma_obs[i] = gamma[regions[i]]; // repeat the correct RE for both causes
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(beta + gamma_obs[i], Sigma[i]); // bivariate normal observations
    }
    gamma ~ normal(0, sigma_gamma); // IID normal REs on region
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior on the standard deviation (because sigma is lower bounded by 0)
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas
}
