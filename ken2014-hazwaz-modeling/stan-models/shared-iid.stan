data {
    int<lower=1> R; // number of regions
    int<lower=1> regions[R]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[R]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[R]; // outcome
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    vector[2] beta; // FEs on cause
    vector[R] v_1; // IID RE on region for cause 1, centered on beta1
    vector[R] v_2; // IID RE on region for cause 2, centered on beta2
    real<lower=0> sigma[2]; // IID standard deviations
    real<lower=0> lambda; // proportion scale for shared RE
}
transformed parameters {
    matrix[R, 2] mu; // means of bivariate normal obs
    
    for (i in 1:R) {
        mu[i, 1] = v_1[i] + lambda * v_2[i]; //
        mu[i, 2] = v_2[i]; //
    }
}
model {
    for (i in 1:R) {
        y[i] ~ multi_normal(to_vector(mu[i,]), Sigma[i]); // bivariate normal observations
    }
    
    sigma ~ normal(0, sigma_normal_sd);
    
    v_1 ~ normal(0, 1); // IID normal REs on region
    v_2 ~ normal(0, 1); // IID normal REs on region

    lambda ~ normal(0, 1);
    
    beta ~ normal(0,5);
}
generated quantities {
    real log_sigma[2];
    matrix[R, 2] preds;
    real log_lik[R]; // log likelihood for use with WAIC or PSIS-LOO
    
    log_sigma = log(sigma);
    preds = mu;
    
    for (r in 1:R) {
        log_lik[r] = multi_normal_lpdf(y[r] | to_vector(mu[r,]), Sigma[r]);
    }
}
