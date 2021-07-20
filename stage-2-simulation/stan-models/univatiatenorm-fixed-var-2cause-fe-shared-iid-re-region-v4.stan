data {
    int<lower=1> N; // number of normal observations (=I*R*2)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which observation
    int<lower=1> causes[N]; // index for which cause is for which observation
    real<lower=0> Sigma[N]; // fixed vector of variances for each observed pair
    vector[N] y; // outcome
}
parameters {
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    vector[R] alpha1; // RE on region for cause 1, centered on beta1
    vector[R] alpha2; // RE on region for cause 2, centered on beta2
    vector[R] delta; // shared RE on region
    real<lower=0> sigma_delta; // standard deviation of the shared IID REs on region
    vector[2] beta; // FEs on cause
    real<lower=0> lambda; // scaling coefficient on shared component
}
transformed parameters {
    vector[N] mu; // means of normal obs

    for (i in 1:N) {
         mu[i] = alpha1[regions[i]] + (delta[regions[i]] * lambda); //
         mu[i] = alpha2[regions[i]] + (delta[regions[i]] / lambda); //
    }

    // gamma1 = alpha1 - (beta[1]);
    // gamma2 = alpha2 - (beta[2]);
}
model {
    for (i in 1:N) {
        y[i] ~ normal(mu[i], Sigma[i]); // normal observations
    }
    alpha1 ~ normal(beta[1], sigma_gamma[1]); // IID normal REs on region
    // sum(alpha1) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha1) ~ normal(0,0.001); jacobian of transformation is 1
    alpha2 ~ normal(beta[2], sigma_gamma[2]); // IID normal REs on region
    // sum(alpha2) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha2) ~ normal(0,0.001); jacobian of transformation is 1
    delta ~ normal(0, sigma_delta); // shared IID normal REs on region
    // sum(delta) ~ normal(0, 0.001 * R); // equivalent to mean(delta) ~ normal(0,0.001); jacobian of transformation is 1
    // for (rr in 1:R) {
    //     (alpha1[rr] - beta[1] + (delta[rr] * lambda)) ~ normal(0, 0.001 * 2);
    //     (alpha2[rr] - beta[2] + (delta[rr] / lambda)) ~ normal(0, 0.001 * 2);
    // }
    sigma_gamma ~ student_t(3,0,5); // leads to a half t prior
    sigma_delta ~ student_t(3,0,5); // leads to a half t prior
    lambda ~ lognormal(log(2),0.001); // leads to a lognormal prior
    beta ~ normal(0,5); 
}
