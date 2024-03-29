data {
int<lower=0> C; // number of cells
int<lower=0> G; // number of genes
int<lower=0> K; // number of celltypes
row_vector<lower=0>[G] meanExprR[K]; // predictor matrix
int<lower=0> countsR[C,G];
real<lower=0> prior[K];
vector<lower=0>[2] eta_prior_params;
}

parameters {
real<lower = 0> eta;
//real<lower = 0> alpha;
}

transformed parameters {
matrix[K,C] lp;
row_vector<lower=0>[G] meanExprT[K];
for (k in 1:K){
{
meanExprT[k] = meanExprR[k] * eta;
for (c in 1:C){
lp[k,c] = log(prior[k]) + neg_binomial_2_lpmf(countsR[c] | meanExprT[k],50);
}
}
}
}

model {
// sensitivity prior
//alpha ~ gamma(10,10);
eta ~ gamma(eta_prior_params[1], eta_prior_params[2]);
// add to the log density
target += log_sum_exp(to_vector(lp));
}


