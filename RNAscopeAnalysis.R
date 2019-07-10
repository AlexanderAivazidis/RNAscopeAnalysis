accuracyMetrics = function(fit, classes, meanExprR, prior, subset = NA){
  if (is.na(subset)){
    subset = 1:length(classes)
  }
  
  classes = classes[subset]
  C = length(classes)
  K = dim(meanExprR)[1]
  
  # Extract probability for each cell type and maxProp predictions
  list_of_draws <- extract(fit)
  k_prob = lapply(subset, function(x) colMeans(exp(list_of_draws[['lp']][,,x])))
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
  
  # Compare predictions to random choice:
  predicted_correct = sum(predictions == classes)/length(subset)
  random_correct = mean(unlist(lapply(1:100000, function(x) sum(sample(rownames(meanExprR),length(subset), replace = TRUE, prob = prior) == classes)/length(subset))))
  
  true_classes = matrix(0,C,K)
  colnames(true_classes) = rownames(meanExprR)
  for (i in 1:C){
    true_classes[i,classes[i]] = 1
  }
  
  # Get accuracy as function of confidence:
  N = 5
  accuracy_binned = rep(0,N)
  accuracy_sampleSize = rep(0,N)
  accuracy_mean = rep(0,N)
  for (i in 1:N){
    relevant = (k_matrix > (i-1)*(1/N) & k_matrix < i*(1/N))
    accuracy_mean[i] = mean(k_matrix[relevant])
    accuracy_binned[i] = sum(true_classes[relevant])/sum(relevant)  
    accuracy_sampleSize[i] = sum(relevant)
  }
  
  # Get accuracy as function of cell class:
  tab1 = table(classes[predictions == classes])
  tab2 = table(classes)
  tab1new = tab2
  tab1new[names(tab1)] = tab1
  tab1new[!names(tab1new) %in% names(tab1)] = 0
  tab1 = tab1new[names(tab2)]
  accuracy_cellClass = unlist(lapply(1:length(tab1), function(x) tab1[x]/tab2[x]))

  return(list(k_matrix, predicted_correct, random_correct, accuracy_cellClass))
}

plotAccuracyVsConfidence = function(){
  plot(accuracy_mean, accuracy_binned, pch = 16, cex = 1, ylab = 'prediction_accuracy', xlab = 'prediction_confidence', ylim = c(0,1))
  lines(seq(0,1,0.1), seq(0,1,0.1), lty = 'dashed')
}

# 
# # Efficiency of scRNAseq and RNAscope:
# eta1 = rgamma(4, shape = 10, rate = 1000)
# names(eta1) = rownames(counts)
# eta2 = rgamma(4, shape = 40, rate = 100)
# names(eta2) = RNAscopeGenes1
# 
# # Simulate true mean expression:
# meanExpr_true = apply(meanExpr,2, function(x) x/eta1)
# 
# # Calculate RNAscope mean expression:
# meanExpr_scope = apply(meanExpr_true[RNAscopeGenes1,], 2, function(x) x*eta2)
# 
# # Cell Area (normalized):
# # area = lapply(1:m, function(x) rtruncnorm(n, a=0, b=Inf, mean=1, sd=0.25))
# 
# # Simulate genotype effect size:
# celltype_prob0 = sapply(unique(specific_type), function(x) rexp(1,1))
# celltype_prob0 = celltype_prob0/sum(celltype_prob0)
# celltype_prob1 = sapply(unique(specific_type), function(x) rexp(1,1))
# celltype_prob1 = celltype_prob1/sum(celltype_prob1)
# 
# # Simulate RNAscope data (ignoring area for now):
# celltypes_prob = table(specific_type)/sum(table(specific_type))
# celltypes = lapply(1:6, function(x) sample(x = names(celltypes_prob), size = 3000, prob = (celltype_prob1+genotype[x])+(celltype_prob0 *(1-genotype[x])), replace = TRUE))
# counts_RNAscope = lapply(1:m, function(z) sapply(celltypes[[z]], function(y) sapply(RNAscopeGenes1, function(x) rnbinom(1, size = 2, mu = meanExpr_scope[x,y]))))
