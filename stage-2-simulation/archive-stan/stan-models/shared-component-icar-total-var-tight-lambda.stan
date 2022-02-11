data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];
    int<lower=1, upper=N> node2[N_edges];
}
parameters {
    real<lower=0> sigma_total; // total variation of REs
    simplex[3] kappa; //  percent variability for each RE
    vector[R] alpha1; // RE on region for cause 1, centered on beta1
    vector[R] alpha2; // RE on region for cause 2, centered on beta2
    vector[R] theta; // shared RE on region
    vector[2] beta; // FEs on cause
    real<lower=0> delta; // scaling coefficient on shared component
    vector[R] phi_1;
    vector[R] phi_2;
    vector[R] phi_3;
    real<lower=0> sigma_phi[3]; // total variation of REs
}
transformed parameters {
    matrix[N, 2] mu; // means of bivariate normal obs
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    real<lower=0> sigma_theta;
    vector[R] gamma1; // RE on region for cause 1, centered on beta1
    vector[R] gamma2; // RE on region for cause 2, centered on beta2
    
    sigma_gamma[1] = sqrt(kappa[1] * (sigma_total^2)); 
    sigma_gamma[2] = sqrt(kappa[2] * (sigma_total^2)); 
    sigma_theta = sqrt(kappa[3] * (sigma_total^2)); 
    
    for (i in 1:N) {
         mu[i, 1] = alpha1[regions[i]] + (phi_1[regions[i]] * sigma_phi[1]) + (delta * (theta[regions[i]] + (phi_3[regions[i]] * sigma_phi[3]))); //
         mu[i, 2] = alpha2[regions[i]] + (phi_2[regions[i]] * sigma_phi[2]) + ((theta[regions[i]] + (phi_3[regions[i]] * sigma_phi[3])) / delta); //
    }
    
    gamma1 = alpha1 - beta[1];
    gamma2 = alpha2 - beta[2];
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i]), Sigma[i]); // bivariate normal observations
    }
    
    // the following computes the prior on phi on the unit scale with sd = 1
    target += -0.5 * dot_self(phi_1[node1] - phi_1[node2]);
    target += -0.5 * dot_self(phi_2[node1] - phi_2[node2]);
    target += -0.5 * dot_self(phi_3[node1] - phi_3[node2]);
    // soft sum-to-zero constraint on phi)
    sum(phi_1) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    sum(phi_2) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
    sum(phi_3) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)

    sigma_phi ~ student_t(3,0,1);
    alpha1 ~ normal(beta[1], sigma_gamma[1]); // IID normal REs on region
    alpha2 ~ normal(beta[2], sigma_gamma[2]); // IID normal REs on region
    theta ~ normal(0, sigma_theta); // shared IID normal REs on region
    sigma_total ~ student_t(3,0,1); // leads to a half t prior
    kappa ~ dirichlet(rep_vector(1, 3));
    delta ~ lognormal(log(1.5),0.001); // leads to a lognormal prior
    // beta ~ normal(0,5); 
}
