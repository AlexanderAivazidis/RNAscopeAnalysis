### Simulate expected data to check if model works as expected

require(rstan)

setwd('/home/jovyan/RNAscopeAnalysis/')

### Simulate data:

RNAscopeGenes1 = c("Hdac6", "Myh8", "Plcxd3", "Klf5")
RNAscopeGenes2 = c("Slc1a3", "Gad1", "Plp1", "Pdgfra")

# Mean expression profiles in scRNAseq data:
counts = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_counts.rds')
coldata = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_coldata.rds')
keep = which(rownames(counts) %in% RNAscopeGenes1)
counts = counts[keep,]
specific_type = as.character(coldata[,2])
names(specific_type) = colnames(counts)
meanExpr = do.call("cbind", tapply(names(specific_type), specific_type, function(x) rowMeans(counts[,x])))

set.seed(999) # Random number seed
n = 1000 # Number of cells to simulate
eta = c(100,100,100,100)
size = c(2,2,2,2)
names(eta) = names(size) = RNAscopeGenes1

# RNAscope data:
meanExpr_scope = meanExpr*eta
celltypes_prob = table(specific_type)/sum(table(specific_type))
celltypes = sample(x = names(celltypes_prob), size = 1000, prob = celltypes_prob, replace = TRUE)
counts_RNAscope = sapply(celltypes, function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = size[x], mu = meanExpr_scope[x,y])))

### Perform inference on data:

classes = colnames(counts_RNAscope)
countsR = t(counts_RNAscope)
meanExprR = t(meanExpr_scope) + 0.001
prior = rep(1,length(celltypes_prob))
names(prior) = rownames(meanExprR)
prior = prior/sum(prior)
#prior = prior[rownames(meanExprR)]
C = dim(countsR)[1]
K = dim(meanExprR)[1]
G = dim(countsR)[2]
stan_rdump(c('prior', 'countsR','meanExprR', 'C', 'K', 'G'), file = "trial1_model1_data.txt", append = FALSE)
trial1_model1_data <- read_rdump("trial1_model1_data.txt") 

fit <- stan(
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

fit_summary <- summary(fit)



