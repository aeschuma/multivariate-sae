data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
}
parameters {
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    vector[R] alpha1; // REs on region for cause 1 (shared)
    vector[R] alpha2; // RE on region for cause 2
    vector[R] delta; // shared RE on region
    real<lower=0> sigma_delta; // standard deviation of the shared IID REs on region
    vector[2] beta; // FEs on cause
    real<lower=0> lambda; // scaling coefficient on shared component
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs

    for (i in 1:N) {
         mu[i, 1] = alpha1[regions[i]] + (lambda * delta[regions[i]]); //
         mu[i, 2] = alpha2[regions[i]] + (lambda * delta[regions[i]]); //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i]), Sigma[i]); // bivariate normal observations
    }
    alpha1 ~ normal(beta[1], sigma_gamma[1]); // IID normal REs on region
    sum(alpha1) ~ normal(R * beta[1], 0.001 * R);  // equivalent to mean(alpha1) ~ normal(0,0.001)
    alpha2 ~ normal(beta[2], sigma_gamma[2]); // IID normal REs on region
    sum(alpha2) ~ normal(R * beta[2], 0.001 * R);  // equivalent to mean(alpha2) ~ normal(0,0.001)
    delta ~ normal(0, sigma_gamma[2]); // shared IID normal REs on region
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
