data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int regions[N]; // index for which region is for which bivariate observation
    vector[N] sigma; // fixed block diagonal matrix of covariances for each observed pair
    vector[N] y; // outcome
}
parameters {
    real<lower=0> sigma_gamma; // standard deviation of the IID REs on region
    real gamma[R]; // REs on region
    real beta; // FEs on cause
}
transformed parameters {
    vector[N] gamma_obs; // values of REs for each bivariate normal observation

    for (i in 1:N) {
        gamma_obs[i] = gamma[regions[i]]; // repeat the correct RE for both causes
    }
}
model {
    for (i in 1:N) {
        y[i] ~ normal(beta + gamma_obs[i], sigma); // bivariate normal observations
    }
    gamma ~ normal(0, sigma_gamma); // IID normal REs on region
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior on the standard deviation (because sigma is lower bounded by 0)
    sum(gamma) ~ normal(0, 0.001 * R);  // equivalent to mean(gamma) ~ normal(0,0.001)
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas
}
