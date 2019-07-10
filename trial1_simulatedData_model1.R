### Simulate expected data to check if model works as expected
require(myUtils)
require(rstan)
require(ComplexHeatmap)
source('/home/jovyan/RNAscopeAnalysis/RNAscopeAnalysis.R')

### Simulate data:

RNAscopeGenes1 = c("Slc1a3", "Gad1", "Plp1", "Pdgfra")

# Mean expression profiles in scRNAseq data:
options(stringsAsFactors = FALSE)
counts = readRDS('/home/jovyan/data/Allen/mouse_VISp_2018-06-14_exon+intron_counts-matrix.rds')
rowdata = read.delim('/home/jovyan/data/Allen/mouse_ALM_2018-06-14_genes-rows.csv', sep = ',')
rownames(counts) = rowdata[,1]
counts = counts[,2:dim(counts)[2]]
coldata = read.delim('/home/jovyan/data/Allen/mouse_VISp_2018-06-14_samples-columns.csv', sep = ',')
specific_type = as.character(coldata[,'subclass'])
names(specific_type) = colnames(counts)
classes = coldata[,'class']
names(classes) = colnames(counts)
keep = (classes %in% c('GABAergic', 'Endothelial', 'Glutamatergic', 'Non-Neuronal'))
classes = classes[keep]
specific_type = specific_type[keep]
counts = counts[RNAscopeGenes1,keep]
classes[specific_type == 'Astro'] = 'Astro'
classes[specific_type == 'Endo'] = 'Endo'
classes[specific_type == 'Macrophage'] = 'Macrophage'
classes[specific_type == 'Oligo'] = 'Oligo'
keep = (classes %in% c('Astro', 'Endothelial', 'GABAergic', 'Glutamatergic', 'Macrophage', 'Oligo'))
classes = classes[keep]
specific_type = specific_type[keep]
counts = counts[RNAscopeGenes1,keep]
meanExpr = do.call("cbind", tapply(names(classes), classes, function(x) rowMeans(counts[RNAscopeGenes1,x])))

set.seed(999) # Random number seed
n1 = 10000 # Number of cells to simulate
eta = c(0.5,1,1.5,1)
size = c(2,2,2,2)
names(eta) = names(size) = RNAscopeGenes1

# RNAscope data:
meanExpr_scope = apply(meanExpr,2, function(x) x*eta+0.001)
celltypes_prob = table(classes)/sum(table(classes))
celltypes = sample(x = names(celltypes_prob), size = n1, prob = celltypes_prob, replace = TRUE)
counts_RNAscope = sapply(celltypes, function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = size[x], mu = meanExpr_scope[x,y])))
# effect = rnorm(length(celltypes_prob), mean = 1, sd = 0.05)
# celltypes_prob1 = effect*table(specific_type)/sum(table(specific_type))
# celltypes_prob1 = celltypes_prob1/sum(celltypes_prob1)
# celltypes1 = sample(x = names(celltypes_prob1), size = n2, prob = celltypes_prob1, replace = TRUE)
# counts_RNAscope1 = sapply(celltypes1, function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = size[x], mu = meanExpr_scope[x,y])))

### Perform inference on data:
classes = colnames(counts_RNAscope) #c(colnames(counts_RNAscope), colnames(counts_RNAscope1))
countsR = t(counts_RNAscope) #t(cbind(counts_RNAscope, counts_RNAscope1))
meanExprR = t(meanExpr) + 0.001
prior = celltypes_prob
prior = prior[rownames(meanExprR)]
C = dim(countsR)[1]
K = dim(meanExprR)[1]
G = dim(countsR)[2]
stan_rdump(c('prior', 'countsR','meanExprR', 'C', 'K', 'G'), file = "trial1_model2_data.txt", append = FALSE)
trial1_model2_data <- read_rdump("trial1_model2_data.txt") 

fit1 <- stan(
  file = "/home/jovyan/RNAscopeAnalysis/smFISH_model2.stan",       # Stan program
  data = trial1_model2_data,         # named list of data
  chains = 4,            # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 10000,             # total number of iterations per chain
  cores = 4,             # number of cores (using 2 just for the vignette)
  refresh = 50,
  algorithm = "Fixed_param"
)

### Analyse the results:
fit_summary <- summary(fit1)
metrics = accuracyMetrics(fit1, classes, meanExprR, prior, subset = 1:n1)

# Total cell type numbers:
samples = 1000
k_matrix = t(metrics[[1]])
kTotal_matrix = matrix(0,dim(k_matrix)[1], samples)
rownames(kTotal_matrix) = rownames(k_matrix)
for (i in 1:samples){
  print(i)
  tab = table(unlist(lapply(1:n1, function(x) sample(rownames(k_matrix), size = 1, prob = k_matrix[,x]))))
  kTotal_matrix[names(tab),i] = tab
}

## Run vb algorithm:

model = stan_model(file = '/home/jovyan/RNAscopeAnalysis/smFISH_model2.stan')
fitVB = vb(model, data = trial1_model2_data)

fit_summaryVB <- summary(fitVB)
