### Simulate expected data to check if model works as expected
require(myUtils)
require(rstan)
require(ComplexHeatmap)

### Simulate data:

RNAscopeGenes1 = c("Slc1a3", "Gad1", "Plp1", "Pdgfra")
cellGroupsUp = c('Neuron.Gad1Gad2.Drd1-Nefm', 'Neuron.Gad1Gad2.Adora2a-Nefm', 'Neuron.Gad1Gad2.Drd1-Fos')
cellGroupsDown = c('Neuron.Gad1Gad2.Pvalb-Rgs12','Neuron.Gad1Gad2.Pnoc')

# Mean expression profiles in scRNAseq data:
counts = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_counts.rds')
coldata = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_coldata.rds')
specific_type = as.character(coldata[,2])
names(specific_type) = colnames(counts)
meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(counts[RNAscopeGenes1,x])))

set.seed(999) # Random number seed
n1 = 1000 # Number of cells to simulate
n2 = 1000
eta = c(10,10,10,10)
size = c(2,2,2,2)
names(eta) = names(size) = RNAscopeGenes1

# RNAscope data:
meanExpr_scope = meanExpr*eta+0.001
celltypes_prob = table(specific_type)/sum(table(specific_type))
celltypes = sample(x = names(celltypes_prob), size = n1, prob = celltypes_prob, replace = TRUE)
counts_RNAscope = sapply(celltypes, function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = size[x], mu = meanExpr_scope[x,y])))
effect = rnorm(length(celltypes_prob), mean = 1, sd = 0.05)
celltypes_prob1 = effect*table(specific_type)/sum(table(specific_type))
celltypes_prob1 = celltypes_prob1/sum(celltypes_prob1)
celltypes1 = sample(x = names(celltypes_prob1), size = n2, prob = celltypes_prob1, replace = TRUE)
counts_RNAscope1 = sapply(celltypes1, function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = size[x], mu = meanExpr_scope[x,y])))

### Perform inference on data:
classes = c(colnames(counts_RNAscope), colnames(counts_RNAscope1))
countsR = t(cbind(counts_RNAscope, counts_RNAscope1))
meanExprR = t(meanExpr_scope) + 0.001
prior = celltypes_prob
prior = prior[rownames(meanExprR)]
C = dim(countsR)[1]
K = dim(meanExprR)[1]
G = dim(countsR)[2]
stan_rdump(c('prior', 'countsR','meanExprR', 'C', 'K', 'G'), file = "trial1_model1_data.txt", append = FALSE)
trial1_model1_data <- read_rdump("trial1_model1_data.txt") 

fit1 <- stan(
  file = "/home/jovyan/RNAscopeAnalysis/smFISH_model1.stan",       # Stan program
  data = trial1_model1_data,         # named list of data
  chains = 4,            # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
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
