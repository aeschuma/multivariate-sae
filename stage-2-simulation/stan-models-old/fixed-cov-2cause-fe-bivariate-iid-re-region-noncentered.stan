data {
    int<lower=1> N; // number of bivariate normal observations (=I*R)
    int<lower=1> R; // number of regions
    int<lower=1> regions[N]; // index for which region is for which bivariate observation
    matrix[2, 2] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[2] y[N]; // outcome
}
parameters {
    vector<lower=0>[2] sigma_gamma; // standard deviations of the IID REs on region
    matrix[R, 2] gamma_star; // REs on region (bivariate), on normal(0,1) scale to be transformed
    vector[2] beta; // FEs on cause
    cholesky_factor_corr[2] Lcorr;// cholesky factor (L_u matrix for Omega)
}
transformed parameters {
    vector[2] mu[N]; // means of bivariate normal obs
    corr_matrix[2] Omega; // correlation matrix for REs
    cov_matrix[2] V; // covaraince matrix for REs
    matrix[R, 2] gamma_mat; // alphas in matrix form
    matrix[2, 2] L; // cholesky decomposition of V
    
    Omega = multiply_lower_tri_self_transpose(Lcorr);
    V = quad_form_diag(Omega, sigma_gamma); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
    L = cholesky_decompose(V);
    gamma_mat = gamma_star * L'; // b in matrix form
    for (i in 1:N) {
         mu[i] = beta + to_vector(gamma_mat[regions[i], ]); //
    }
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(mu[i], Sigma[i]); // bivariate normal observations
    }
    to_vector(gamma_star) ~ normal(0, 1); // IID normal REs on region
    Lcorr ~ lkj_corr_cholesky(1); 
    sigma_gamma ~ student_t(3,0,1); // leads to a half t prior
    beta ~ normal(0,5); 
}
generated quantities {
    real rho_gamma;
    
    rho_gamma = Omega[1, 2];
}
