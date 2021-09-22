data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    vector[2] beta; // FEs on cause
    vector[R] gamma_1; // IID RE on region for cause 1, centered on beta1
    vector[R] gamma_2; // IID RE on region for cause 2, centered on beta2
    vector[R] gamma_3; // shared IID RE on region
    real<lower=0> sigma[3]; // overall standard deviations
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs
    
    for (i in 1:N) {
         mu[i, 1] = beta[1] + gamma_1[regions[i]] * sigma[1] + (gamma_3[regions[i]] * sigma[3]); //
         mu[i, 2] = beta[2] + gamma_2[regions[i]] * sigma[2] + (gamma_3[regions[i]] * sigma[3]); //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i,]), Sigma[i]); // bivariate normal observations
    }
    
    sigma ~ normal(0, sigma_normal_sd);
    
    gamma_1 ~ normal(0, 1); // IID normal REs on region
    gamma_2 ~ normal(0, 1); // IID normal REs on region
    gamma_3 ~ normal(0, 1); // IID normal REs on region
    
    beta ~ normal(0,5);
}
