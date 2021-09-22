data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
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
    matrix[N, 2] mu; // means of bivariate normal obs
    real<lower=0> sigma_gamma[2]; // standard deviation of the IID REs on region
    real<lower=0> sigma_delta;
    
    sigma_gamma[1] = sqrt(theta[1] * (sigma_total^2)); 
    sigma_gamma[2] = sqrt(theta[2] * (sigma_total^2)); 
    sigma_delta = sqrt(theta[3] * (sigma_total^2)); 

    for (i in 1:N) {
         mu[i, 1] = alpha1[regions[i]] + (delta[regions[i]] * lambda); //
         mu[i, 2] = alpha2[regions[i]] + (delta[regions[i]] / lambda); //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i]), Sigma[i]); // bivariate normal observations
    }
    alpha1 ~ normal(beta[1], sigma_gamma[1]); // IID normal REs on region
    alpha2 ~ normal(beta[2], sigma_gamma[2]); // IID normal REs on region
    delta ~ normal(0, sigma_delta); // shared IID normal REs on region
    sigma_total ~ student_t(7,0,1); // leads to a half t prior
    theta ~ dirichlet(rep_vector(1, 3));
    lambda ~ lognormal(log(1.5), 0.001); // leads to a lognormal prior
    beta ~ normal(0,5); 
}
