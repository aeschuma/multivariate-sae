data {
    int<lower=1> R; // number of regions
    int<lower=1> regions[R]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[R]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[R]; // outcome
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    vector[2] beta; // FEs on cause
    vector[R] gamma_1; // RE on region for cause 1, centered on beta1
    vector[R] gamma_2; // RE on region for cause 2, centered on beta2
    real<lower=0> sigma[2]; // overall standard deviations
}
transformed parameters {
    matrix[R, 2] mu; // means of bivariate normal obs

    for (i in 1:R) {
         mu[i, 1] = beta[1] + gamma_1[i]; //
         mu[i, 2] = beta[2] + gamma_2[i]; //
    }
}
model {
    for (i in 1:R) {
        y[i] ~ multi_normal(to_vector(mu[i,]), Sigma[i]); // bivariate normal observations
    }

    sigma ~ gamma(5, 0.00005);
    gamma_1 ~ normal(0, sigma[1]); // IID normal REs on region
    gamma_2 ~ normal(0, sigma[2]); // IID normal REs on region
    // beta ~ normal(0, 5);
}
generated quantities {
    real log_sigma[2];
    vector[R] v_1;
    vector[R] v_2;
    matrix[R, 2] preds;
    real log_lik[R]; // log likelihood for use with WAIC or PSIS-LOO
    
    log_sigma = log(sigma);
    preds = mu;
    
    for (r in 1:R) {
        log_lik[r] = multi_normal_lpdf(y[r] | to_vector(mu[r,]), Sigma[r]);
    }
}
