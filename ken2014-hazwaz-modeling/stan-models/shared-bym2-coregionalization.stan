data {
    int<lower=1> R; // number of regions
    int<lower=1> regions[R]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[R]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[R]; // outcome
    int<lower=0> N_edges;
    int<lower=1, upper=R> node1[N_edges];
    int<lower=1, upper=R> node2[N_edges];
    real<lower=0> scaling_factor; // scales the variances of the spatial effects
    real<lower=0> rho_beta_1; // beta dist hyperpar 1 for % spatial paramter
    real<lower=0> rho_beta_2; // beta dist hyperpar 2 for % spatial paramter
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    vector[2] beta; // FEs on cause
    vector[R] v_1; // IID RE on region for cause 1, centered on beta1
    vector[R] v_2; // IID RE on region for cause 2, centered on beta2
    vector[R] u_1; // ICAR RE on region for cause 1
    vector[R] u_2; // ICAR RE on region for cause 2, shared with cause 1
    real<lower=0> sigma[2]; // overall standard deviations
    real<lower=0, upper=1> rho[2]; // proportion unstructured vs. spatially structured variances
    real<lower=0> lambda; // proportion scale for shared RE
}
transformed parameters {
    matrix[R, 2] mu; // means of bivariate normal obs
    vector[R] convolved_re_1;
    vector[R] convolved_re_2;
    
    // variance of each component should be approximately equal to 1
    convolved_re_1 =  (sqrt(1 - rho[1]) * v_1) + (sqrt(rho[1] / scaling_factor) * u_1);
    convolved_re_2 =  (sqrt(1 - rho[2]) * v_2) + (sqrt(rho[2] / scaling_factor) * u_2);
    
    for (i in 1:R) {
         mu[i, 1] = beta[1] + convolved_re_1[regions[i]] * sigma[1] + (lambda * u_2[regions[i]]); //
         mu[i, 2] = beta[2] + convolved_re_2[regions[i]] * sigma[2]; // recall that the convolved RE includes z = u_2 / delta
    }
}
model {
    for (i in 1:R) {
        y[i] ~ multi_normal(to_vector(mu[i,]), Sigma[i]); // bivariate normal observations
    }
    
    // the following computes the prior on phi on the unit scale with sd = 1
    target += -0.5 * dot_self(u_1[node1] - u_1[node2]);
    target += -0.5 * dot_self(u_2[node1] - u_2[node2]);

    // soft sum-to-zero constraint on phi)
    sum(u_1) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    sum(u_2) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)

    rho ~ beta(rho_beta_1, rho_beta_2);
    sigma ~ normal(0, sigma_normal_sd);
    
    v_1 ~ normal(0, 1); // IID normal REs on region
    v_2 ~ normal(0, 1); // IID normal REs on region

    lambda ~ normal(0, 100);
    
    // beta ~ normal(0,5);
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
