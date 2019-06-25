### Test RNAscope analysis algorithms 

require(AllenData)

# # Load Allen Data data
# readRDS(reads, file = '/home/jovyan/data/Allen/mouse_VISp_2018-06-14_exon+intron_counts-matrix.rds')
# reads = data[[1]]
# coldata1 = read.delim('/home/jovyan/data/Allen/mouse_VISp_2018-06-14_samples-columns.csv', sep = ',')
# rowdata1 = read.delim('/home/jovyan/data/Allen/mouse_VISp_2018-06-14_genes-rows.csv', sep = ',')
# 
# # Remove cells that occure less than 5 times:
# reads = reads[,table(coldata1[,'cluster'])[coldata1[,'cluster']] > 5]
# coldata1 = coldata1[table(coldata1[,'cluster'])[coldata1[,'cluster']] > 5,]
# rownames(reads) = rowdata1[,1]
# reads = reads[,2:dim(reads)[2]]
# 
# # Calculate mean expression:
# classes = as.character(coldata1[,'class'])
# keep = (classes %in% c('GABAergic', 'Endothelial', 'Glutamatergic', 'Non-Neuronal'))
# reads = reads[,keep]
# coldata1 = coldata1[keep,]
# classes = paste(as.character(coldata1[,'class']), as.character(coldata1[,'subclass']), sep = '_')
# # classes[coldata1[,'class'] == 'Endothelial'] = 'Endothelial'
# # classes[coldata1[,'class'] == 'GABAergic'] = 'GABAergic'
# # classes[coldata1[,'class'] == 'Glutamatergic'] = 'Glutamatergic'

# Load DropSeq data:

counts = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_counts.rds')
coldata = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_coldata.rds')
rowdata = readRDS('/home/jovyan/data/Saunders/Saunders_Mouse_Striatum_rowdata.rds')

subclasses = as.character(coldata$subclass)
classes = unlist(lapply(1:length(subclasses), function(x) strsplit(subclasses[x], split = '\\.')[[1]][1]))
names(classes) = colnames(counts)
meanExpr = do.call("cbind", tapply(names(classes), classes, function(x) rowMeans(counts[,x])))
runGenes = c("Hdac6", "Myh8", "Plcxd3", "Klf5")
runGenes = c("Slc1a3", "Gad1", "Plp1", "Pdgfra")

## Calculate probabilities for each cell class:

setwd('/home/jovyan/')

require(rstan)

keep = sample(1:dim(counts)[2], 10000)
countsR = t(counts[runGenes,keep])
meanExprR = t(meanExpr[runGenes,]) + 0.001
#meanExprR = t(apply(meanExprR, 1, function(x) x*c(0.5,0.5,0.5,0.5)))
prior = table(as.character(classes))/length(classes)
prior = prior[rownames(meanExprR)]
C = dim(countsR)[1]
K = dim(meanExprR)[1]
G = length(runGenes)

stan_rdump(c('prior', 'countsR','meanExprR', 'C', 'K', 'G'), file = "QC_data.txt", append = FALSE)
QC_data <- read_rdump("QC_data.txt") 

## working seed: 1002649637

fit <- stan(
  file = "QC1.stan",       # Stan program
  data = QC_data,         # named list of data
  chains = 3,            # number of Markov chains
  warmup = 1000,           # number of warmup iterations per chain
  iter = 2000,             # total number of iterations per chain
  cores = 3,             # number of cores (using 2 just for the vignette)
  refresh = 50,
  algorithm = "Fixed_param"
)

fit2 = vb(stan_model("QC1.stan"), data = QC_data, iter = 10^4)

saveRDS(fit, file = 'fit1.rds')

## Analyse and plot fit:

# Probability for each cell class:

list_of_draws <- extract(fit)
k_prob = lapply(1:dim(list_of_draws[['lp']])[3], function(x) colMeans(exp(list_of_draws[['lp']][,,x])))
k_matrix = matrix(0,length(k_prob), length(k_prob[[1]]))
colnames(k_matrix) = rownames(meanExprR)
predictions = rep('',length(k_prob))
predictions_Confidence = rep('',length(k_prob))
for (i in 1:length(k_prob))
{
  k_prob[[i]] = k_prob[[i]]/sum(k_prob[[i]])
  k_matrix[i,] = k_prob[[i]]
  names(k_prob[[i]]) = colnames(meanExpr)
  k_prob[[i]] = sort(k_prob[[i]])
  predictions[i] = names(k_prob[[i]])[length(k_prob[[i]])]
  predictions_Confidence[i] = k_prob[[i]][length(k_prob[[i]])]
}

true_classes = matrix(0,C,K)
colnames(true_classes) = rownames(meanExprR)
for (i in 1:C){
  true_classes[i,classes[keep][i]] = 1
}

N = 10
accuracy_binned = rep(0,N)
accuracy_sampleSize = rep(0,N)
accuracy_sd = rep(0,N)
for (i in 1:N){
  relevant = (k_matrix > (i-1)*(1/N) & k_matrix < i*(1/N))
  accuracy_binned[i] = sum(true_classes[relevant])/sum(relevant)  
  accuracy_sampleSize[i] = sum(relevant)
  accuracy_sd[i] = sd(true_classes[relevant])
}

plot(seq(0.05,0.95,0.1), accuracy_binned, pch = 16, cex = 1, ylab = 'prediction_accuracy', xlab = 'prediction_confidence', ylim = c(0,1))
lines(seq(0,1,0.1), seq(0,1,0.1), lty = 'dashed')

d = data.frame(x = seq(0.05,0.95,0.1), y = accuracy_binned, se = accuracy_sd/sqrt(accuracy_sampleSize))

require(Hmisc)
plot(d$x, d$y,pch = 16, cex = 1, ylab = 'prediction_accuracy', xlab = 'prediction_confidence', ylim = c(0,1))
with (
  data = d
  , expr = errbar(x, y, y+se, y-se, add=T, pch=1, cap=.1)
)
lines(seq(0,1,0.1), seq(0,1,0.1), lty = 'dashed')

fit_summary <- summary(fit)
head(fit_summary$summary,5)

# Plot predicted cell type proportions:

# Plot accurcacy per cell type:

# Plot accuracy as function of confidence:





