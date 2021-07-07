data {
    int<lower=1> N; // number of observations (=I*R)
    int<lower=1> R; // number of regions
    int regions[N]; // index for which region is for which observation
    vector[N] sigma; // fixed variances for each observation
    vector[N] y; // outcome
}
parameters {
    real<lower=0> sigma_gamma; // standard deviation of the IID REs on region
    real alpha[R]; // RE on region
    real beta; // fixed intercept
}
transformed parameters {
    vector[N] mu; // values of REs for each bivariate normal observation

    for (i in 1:N) {
        mu[i] = alpha[regions[i]]; // repeat the correct RE for both causes
    }
}
model {
    for (i in 1:N) {
        y[i] ~ normal(mu[i], sigma); // bivariate normal observations
    }
    alpha ~ normal(beta, sigma_gamma); // IID normal REs on region
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5);
}
