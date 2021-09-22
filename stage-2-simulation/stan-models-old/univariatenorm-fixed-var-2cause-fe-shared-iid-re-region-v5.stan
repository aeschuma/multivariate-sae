data {
    int<lower=1> N; // number of normal observations (=I*R*C)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which observation
    int<lower=0, upper=1> cause1ind[N]; // index for which cause is for which observation
    int<lower=0, upper=1> cause2ind[N]; // index for which cause is for which observation
    real<lower=0> Sigma[N]; // fixed vector of variances for each observed pair
    vector[N] y; // outcome
}
parameters {
    real<lower=0> sigma_total; // total variation of REs
    simplex[3] theta; //  percent variability for each RE
    vector[R] alpha1; // RE on region for cause 1, centered on beta1
    vector[R] alpha2; // RE on region for cause 2, centered on beta2
    vector[R] delta; // shared RE on region
    vector[2] beta; // FEs on cause
    real<lower=0> lambda; // scaling coefficient on shared component
}
transformed parameters {
    vector[N] mu; // means of normal obs
    vector[R] gamma1; // RE on region for cause 1, centered on beta1
    vector[R] gamma2; // RE on region for cause 2, centered on beta2
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    real<lower=0> sigma_delta;

    sigma_gamma[1] = sqrt(theta[1] * (sigma_total^2)); 
    sigma_gamma[2] = sqrt(theta[2] * (sigma_total^2)); 
    sigma_delta = sqrt(theta[3] * (sigma_total^2)); 
    
    for (i in 1:N) {
         mu[i] = cause1ind[i] * (beta[1] + alpha1[regions[i]] + (delta[regions[i]] * lambda)) + cause2ind[i] * (beta[2] + alpha2[regions[i]] + (delta[regions[i]] / lambda)); //
    }

    gamma1 = alpha1 - (beta[1]);
    gamma2 = alpha2 - (beta[2]);
}
model {
    for (i in 1:N) {
        y[i] ~ normal(mu[i], Sigma[i]); // normal observations
    }
    alpha1 ~ normal(0, sigma_gamma[1]); // IID normal REs on region
    // sum(alpha1) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha1) ~ normal(0,0.001); jacobian of transformation is 1
    alpha2 ~ normal(0, sigma_gamma[2]); // IID normal REs on region
    // sum(alpha2) ~ normal(0, 0.001 * R);  // equivalent to mean(alpha2) ~ normal(0,0.001); jacobian of transformation is 1
    delta ~ normal(0, sigma_delta); // shared IID normal REs on region
    // sum(delta) ~ normal(0, 0.001 * R); // equivalent to mean(delta) ~ normal(0,0.001); jacobian of transformation is 1

    sigma_total ~ student_t(7,0,2); // leads to a half t prior
    theta ~ dirichlet(rep_vector(1, 3));
    
    lambda ~ lognormal(log(1.5),0.01); // leads to a lognormal prior
    
    beta ~ normal(0,5); 
}
