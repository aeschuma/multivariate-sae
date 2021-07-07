data {
    int<lower=1> N; // number of bivariate normal observations
    int<lower=1> C; // number of causes
    matrix[C, C] Sigma[N]; // fixed block diagonal matrix of covariances for each observed pair
    vector[C] y[N]; // outcome
}
parameters {
    vector[C] beta; // FEs on cause
}
model {
    for (i in 1:N) {
        y[i] ~ multi_normal(beta, Sigma[i]); // bivariate normal observations
    }
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas
}
