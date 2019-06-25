data {
int<lower=0> C; // number of cells
int<lower=0> G; // number of genes
int<lower=0> K; // number of celltypes
row_vector<lower=0>[G] meanExprR[K]; // predictor matrix
int<lower=0> countsR[C,G];
real<lower=0> prior[K];
}

// parameters {
// vector<lower = 0>[G] eta;
// //real<lower = 0> alpha;
// }

transformed parameters {
matrix[K,C] lp;
for (k in 1:K){
{
for (c in 1:C){
lp[k,c] = log(prior[k]) + neg_binomial_2_lpmf(countsR[c]| meanExprR[k],2);
}
}
}
}

model {
// sensitivity prior
//alpha ~ gamma(10,2);
//eta ~ gamma(alpha, 10);
//eta ~ gamma(2, 1);
// add to the log density
target += log_sum_exp(to_vector(lp));
}


