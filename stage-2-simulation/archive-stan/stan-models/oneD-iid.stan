data {
    int<lower=1> R; // number of observations
    real<lower = 0> Sigma[R]; // fixed block variances for each observation 
    vector[R] y; // outcome
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    real beta; // FE
    vector[R] gamma; // RE on region, centered on beta1
    real<lower=0> sigma;        // overall standard deviation
}
model {
    for (i in 1:R) {
        y[i] ~ normal(gamma[i], Sigma[i]); // bivariate normal observations
    }
    
    sigma ~ normal(0, sigma_normal_sd);
    gamma ~ normal(beta, sigma); // IID normal REs on region
    // beta ~ normal(0,5);
}
