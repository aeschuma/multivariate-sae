data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
}
parameters {
    vector<lower=0>[2] sigma_gamma; // standard deviations of the IID REs on region
    vector[2] alpha[R]; // REs on region for cause 1 
    vector[2] beta; // FEs on cause
    cholesky_factor_corr[2] Lcorr;// cholesky factor (L_u matrix for Omega)
}
transformed parameters {
    vector[2] mu[N]; // means of bivariate normal obs
    corr_matrix[2] Omega; // correlation matrix for REs
    cov_matrix[2] V; // covaraince matrix for REs
    
    Omega = multiply_lower_tri_self_transpose(Lcorr);
    V = quad_form_diag(Omega, sigma_gamma); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
    for (i in 1:N) {
         mu[i] = alpha[regions[i]]; //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(to_vector(mu[i]), Sigma[i]); // bivariate normal observations
    }
    alpha ~ multi_normal(beta, V); // IID normal REs on region
    Lcorr ~ lkj_corr_cholesky(1); 
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
