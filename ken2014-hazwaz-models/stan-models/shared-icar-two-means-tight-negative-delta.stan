data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];
    int<lower=1, upper=N> node2[N_edges];
    real<lower=0> scaling_factor; // scales the variances of the spatial effects
    real<lower=0> rho_beta_1; // beta dist hyperpar 1 for % spatial paramter
    real<lower=0> rho_beta_2; // beta dist hyperpar 2 for % spatial paramter
    real<lower=0> sigma_normal_sd; // half normal SD for total var parameter
}
parameters {
    vector[2] beta; // FEs on cause
    vector[R] gamma_1; // IID RE on region for cause 1, centered on beta1
    vector[R] gamma_2; // IID RE on region for cause 2, centered on beta2
    vector[R] gamma_3; // shared IID RE on region
    vector[R] phi_1; // ICAR RE on region for cause 1
    vector[R] phi_2; // ICAR RE on region for cause 2
    vector[R] phi_3; // shared ICAR RE
    real<lower=0> sigma[3]; // overall standard deviations
    real<lower=0, upper=1> rho[3]; // proportion unstructured vs. spatially structured variances
    real<lower=0> delta; // proportion scale for shared RE
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs
    vector[R] convolved_re_1;
    vector[R] convolved_re_2;
    vector[R] convolved_re_3;
    
    // variance of each component should be approximately equal to 1
    convolved_re_1 =  sqrt(1 - rho[1]) * gamma_1 + sqrt(rho[1] / scaling_factor) * phi_1;
    convolved_re_2 =  sqrt(1 - rho[2]) * gamma_2 + sqrt(rho[2] / scaling_factor) * phi_2;
    convolved_re_3 =  sqrt(1 - rho[3]) * gamma_3 + sqrt(rho[3] / scaling_factor) * phi_3;
    
    for (i in 1:N) {
         mu[i, 1] = beta[1] + convolved_re_1[regions[i]] * sigma[1] + (convolved_re_3[regions[i]] * sigma[3]) * delta; //
         mu[i, 2] = beta[2] + convolved_re_2[regions[i]] * sigma[2] + (convolved_re_3[regions[i]] * sigma[3]) / delta; //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i,]), Sigma[i]); // bivariate normal observations
    }
    
    // the following computes the prior on phi on the unit scale with sd = 1
    target += -0.5 * dot_self(phi_1[node1] - phi_1[node2]);
    target += -0.5 * dot_self(phi_2[node1] - phi_2[node2]);
    target += -0.5 * dot_self(phi_3[node1] - phi_3[node2]);
    
    // soft sum-to-zero constraint on phi)
    sum(phi_1) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    sum(phi_2) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    sum(phi_3) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)

    rho ~ beta(rho_beta_1, rho_beta_2);
    sigma ~ normal(0, sigma_normal_sd);
    
    gamma_1 ~ normal(0, 1); // IID normal REs on region
    gamma_2 ~ normal(0, 1); // IID normal REs on region
    gamma_3 ~ normal(0, 1); // IID normal REs on region
    
    delta ~ lognormal(log(0.85), 0.001);
    
    beta ~ normal(0,5);
}
