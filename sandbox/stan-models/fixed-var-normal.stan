data {
    int<lower=1> N; // number of observations
    int<lower=1> nreg; // number of regions
    int<lower=1> ncause; // number of causes
    vector[N] Sigma; // fixed block diagonal matrix of covariances for each observed pair
    vector[N] y; // outcome
    matrix[N, ncause] X; // model matrix for cause FEs
    matrix[N, nreg] Z; // model matrix for region REs
}
parameters {
    real<lower=0> sigma; // standard deviation of the IID REs on region
    vector[nreg] gamma; // REs on region
    vector[ncause] beta; // FEs on cause
}
transformed parameters {
    vector[N] beta_obs; // values of REs for each bivariate normal observation
    vector[N] gamma_obs; // values of REs for each bivariate normal observation

    beta_obs = X * beta;
    gamma_obs = Z * gamma;
}
model {
    y ~ normal(beta_obs + gamma_obs, Sigma); // normal observations
    gamma ~ normal(0, sigma); // IID normal REs on region
    sum(gamma) ~ normal(0, 0.00001 * nreg);  // equivalent to mean(gamma_mat[l,]) ~ normal(0,0.00001)
    sigma ~ student_t(3,0,1); // leads to a half t prior on the standard deviation (because sigma is lower bounded by 0)
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas
}
