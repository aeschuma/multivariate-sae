data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];
    int<lower=1, upper=N> node2[N_edges];
    real<lower = 0> tau_gamma_hyper_alpha;
    real<lower = 0> tau_gamma_hyper_beta;
    real<lower = 0> tau_phi_hyper_alpha;
    real<lower = 0> tau_phi_hyper_beta;
}
parameters {
    real<lower=0> tau_gamma[2]; // IID precision
    real<lower=0> tau_theta;    
    vector[R] gamma1; // RE on region for cause 1, centered on beta1
    vector[R] gamma2; // RE on region for cause 2, centered on beta2
    vector[R] theta; // shared RE on region
    vector[2] beta; // FEs on cause
    real<lower=0> delta; // scaling coefficient on shared component
    vector[R] phi_1;
    vector[R] phi_2;
    vector[R] phi_3;
    real<lower=0> tau_phi[3]; // ICAR precision
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    real<lower=0> sigma_theta;
    real<lower=0> sigma_phi[3]; // standard deviation of the spatial REs on region
    
    for (j in 1:2) {
        sigma_gamma[j] = inv(sqrt(tau_gamma[j]));
        sigma_phi[j] = inv(sqrt(tau_phi[j]));
    }
    sigma_theta = inv(sqrt(tau_theta));
    sigma_phi[3] = inv(sqrt(tau_phi[3]));
    
    for (i in 1:N) {
         mu[i, 1] = beta[1] + (gamma1[regions[i]] * sigma_gamma[1]) + (phi_1[regions[i]] * sigma_phi[1]) + (delta * ((theta[regions[i]] * sigma_theta) + (phi_3[regions[i]] * sigma_phi[3]))); //
         mu[i, 2] = beta[2] + (gamma2[regions[i]] * sigma_gamma[2]) + (phi_2[regions[i]] * sigma_phi[2]) + (((theta[regions[i]] * sigma_theta) + (phi_3[regions[i]] * sigma_phi[3])) / delta); //
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
    
    gamma1 ~ normal(0, 1); // IID normal REs on region
    gamma2 ~ normal(0, 1); // IID normal REs on region
    theta ~ normal(0, 1); // shared IID normal REs on region
    tau_gamma ~ gamma(tau_gamma_hyper_alpha, tau_gamma_hyper_beta);
    tau_theta ~ gamma(tau_gamma_hyper_alpha, tau_gamma_hyper_beta);
    tau_phi ~ gamma(tau_phi_hyper_alpha, tau_phi_hyper_beta);
    delta ~ lognormal(log(1.5), 0.001); // leads to a lognormal prior
    beta ~ normal(0, 5);
}
