data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
}
parameters {
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    vector[R] alpha1star; // RE on region for cause 1, centered on beta1
    vector[R] alpha2star; // RE on region for cause 2, centered on beta2
    vector[R] delta; // shared RE on region
    real<lower=0> sigma_delta; // standard deviation of the shared IID REs on region
    // vector[2] beta; // FEs on cause
    real beta; // overall intercept
    real<lower=0> lambda; // scaling coefficient on shared component
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs
    // vector[R] gamma1; // RE on region for cause 1 (mean 0)
    // vector[R] gamma2; // RE on region for cause 2 (mean 0)
    vector[R] alpha1; // RE on region for cause 1, centered on beta1
    vector[R] alpha2; // RE on region for cause 2, centered on beta2

    for (rr in 1:R) {
        alpha1[rr]  = alpha1star[rr]  * sigma_gamma[1];
        alpha2[rr]  = alpha2star[rr]  * sigma_gamma[2];
    }
    
    for (i in 1:N) {
         mu[i, 1] = alpha1[regions[i]] + (delta[regions[i]] * lambda); //
         mu[i, 2] = alpha2[regions[i]] + (delta[regions[i]] / lambda); //
    }

    // gamma1 = alpha1 - (beta[1]);
    // gamma2 = alpha2 - (beta[2]);
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i]), Sigma[i]); // bivariate normal observations
    }
    alpha1star ~ normal(0, 1); // IID normal REs on region
    // sum(alpha1) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha1) ~ normal(0,0.001); jacobian of transformation is 1
    alpha2star ~ normal(0, 1); // IID normal REs on region
    // sum(alpha2) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha2) ~ normal(0,0.001); jacobian of transformation is 1
    delta ~ normal(beta, sigma_delta); // shared IID normal REs on region
    // sum(delta) ~ normal(0, 0.001 * R); // equivalent to mean(delta) ~ normal(0,0.001); jacobian of transformation is 1
    sigma_gamma ~ student_t(3,0,5); // leads to a half t prior
    sigma_delta ~ student_t(3,0,5); // leads to a half t prior
    lambda ~ lognormal(0,5); // leads to a lognormal prior
    beta ~ normal(0,5); 
}
