# For methods

library(devtools)
install_github('k-perrakis/regjmix')
library(regjmix)
library(pROC)
library(MoEClust)
library(cluster)
library(SwarmSVM)
library(flexmix)

# For figures

library(ggplot2)
library(ggthemes)
library(paletteer)
library(vioplot)

############################### DATA ####################################

genes = read.table('genes.txt', header = TRUE)
head(genes)
dim(genes)

# 1st column has the true cancer type labels:
# 1 = BRCA, 2 = KIRC, 3 = LUAD, 4 = THCA
z = genes[,1]
table(z)
# 2nd column has the response (gene NAPSA) for Application 1 (NAPSA ~ X).
y = genes[,2]
# The rest of the columns are use as features X for Application 1 (NAPSA ~ X).
X = genes[,-c(1,2)]

######################## NAPSA as response ###############################

methods = c('k-means','k-medoids','fuzzy c-means','hclust','clustSVM',
            'mclust','MoEclust','RJM-NJ','RJM-Flasso','RJM-Rlasso')
K = 4
rand = rep(NA, length(methods))
names(rand) = methods
randkmeans = c()
randsvm = c()

### k-means ###

for (j in 1:100) {
  km = kmeans(cbind(y,X),centers=K)
  randkmeans[j] = adjustedRandIndex(z,km$cluster)
}
rand[1] = mean(randkmeans)

### k-medoids ###

kmed = pam(cbind(y,X),k = K)
rand[2] = adjustedRandIndex(z,kmed$clustering)

### fuzzy c means ###

fuzz = fanny(cbind(y,X), k = 4, memb.exp = 1.1)
rand[3]=adjustedRandIndex(z,fuzz$clustering)

### hclust ###

h.clust = hclust(dist(cbind(y,X), method = 'euclidean'),  method = 'ward.D')
z.est = cutree(h.clust, k=K)
rand[4]=adjustedRandIndex(z,z.est)


### cluster-SVM ###

for (j in 1:100) {
  svm = clusterSVM(X,y,centers = K)
  randsvm[j] = adjustedRandIndex(z,svm$label)
}
rand[5] = mean(randsvm)

### mclust ###

W = Mclust(cbind(y,X), G = K, verbose = TRUE)
W = Mclust(
  cbind(y,X),
  G = K,
  modelNames = W$modelName,
  verbose = TRUE
)
z.est = W$classification
rand[6] = adjustedRandIndex(z,z.est)

### moeclust ###

X1 = X
colnames(X1) = paste("V", 1:99, sep="")
mod1 = MoE_stepwise(y, X1, verbose=TRUE, initialG=K, stepG=FALSE, algo = 'CEM')
moe = MoE_compare(mod1,pick=1)$optimal
rand[7] = adjustedRandIndex(z,moe$classification)

## RJM methods ###

rjms=list()
length(rjms)=3
n.cluster = detectCores()
cl = n.cluster
cl = makeCluster(cl)
registerDoParallel(cl)
rjms[[1]] = rjm(y,X,K=4,method='nj',parallel=TRUE)
rjms[[2]] = rjm(y,X,K=4,method='lasso',lasso.pen='fixed',parallel=TRUE)
rjms[[3]] = rjm(y,X,K=4,method='lasso',lasso.pen='random',parallel=TRUE)
stopCluster(cl)

rand[8] = adjustedRandIndex(z,rjms[[1]]$z)
rand[9] = adjustedRandIndex(z,rjms[[2]]$z)
rand[10] = adjustedRandIndex(z,rjms[[3]]$z)

round(rand,2)

betasNJ = sort(as.vector(abs(coef(rjms[[1]])[-1,])),decreasing = TRUE)
betasFL = sort(as.vector(abs(coef(rjms[[1]])[-1,])),decreasing = TRUE)
betasRL = sort(as.vector(abs(coef(rjms[[1]])[-1,])),decreasing = TRUE)


################# All genes as responses - clustering ######################

# This will take a long time to run (about a day)

XX = cbind(y,X)
colnames(XX)[1] = 'NAPSA|9476'
XX1 = XX
p = ncol(XX1)
colnames(XX1) = paste("V", 1:p, sep="")
RAND = matrix(NA, p, 10)
colnames(RAND) = methods
randkmeans = c()
randsvm = c()
rjms=list()
length(rjms) = 3

n.cluster = detectCores()
cl = n.cluster
cl = makeCluster(cl)
registerDoParallel(cl)

for (j in 1:p) {

  ### k-means ###

  for (i in 1:100) {
    km=kmeans(XX[,-j],centers=K)
    randkmeans[i] = adjustedRandIndex(z,km$cluster)
  }
  RAND[j,1]=mean(randkmeans)

  ### k-medoids ###

  kmed = pam(XX[,-j],k = K)
  RAND[j,2] = adjustedRandIndex(z,kmed$clustering)

  ### fuzzy c means ###

  fuzz = fanny(XX[,-j], k = 4, memb.exp = 1.1)
  RAND[j,3] = adjustedRandIndex(z,fuzz$clustering)

  ### hclust ###

  h.clust = hclust(dist(XX[,-j], method = 'euclidean'),  method = 'ward.D')
  z.est = cutree(h.clust, k=K)
  RAND[j,4] = adjustedRandIndex(z,z.est)

  ### cluster-SVM ##

  for (i in 1:100) {
    svm = clusterSVM(XX[,-j],XX[,j],centers = K)
    randsvm[i] = adjustedRandIndex(z,svm$label)
  }
  RAND[j,5] = mean(randsvm)

  ### mclust ###

  W = Mclust(XX[,-j], G = K, verbose = FALSE)
  W = Mclust(
    XX[,-j],
    G = K,
    modelNames = W$modelName,
    verbose = TRUE
  )
  z.est = W$classification
  RAND[j,6] = adjustedRandIndex(z,z.est)

  ### moeclust ###

  mod1 = MoE_stepwise(XX1[,j], XX1[,-j], verbose=TRUE, initialG=K, stepG=FALSE, algo = 'CEM')
  moe  = MoE_compare(mod1,pick=1)$optimal
  RAND[j,7] = adjustedRandIndex(z,moe$classification)

  ## RJM methods ###

  rjms[[1]] = rjm(XX[,j],XX[,-j],K=4,method='nj',parallel=TRUE)
  rjms[[2]] = rjm(XX[,j],XX[,-j],K=4,method='lasso',lasso.pen='fixed',parallel=TRUE)
  rjms[[3]] = rjm(XX[,j],XX[,-j],K=4,method='lasso',lasso.pen='random',parallel=TRUE)

  if(is.null(rjms[[1]])==FALSE){
    RAND[j,8] = adjustedRandIndex(z,rjms[[1]]$z)
  }

  if(is.null(rjms[[2]])==FALSE){
    RAND[j,9] = adjustedRandIndex(z,rjms[[2]]$z)
  }

  if(is.null(rjms[[3]])==FALSE){
    RAND[j,10] = adjustedRandIndex(z,rjms[[3]]$z)
  }

  print(j)
}
stopCluster(cl)

RAND

################# All genes as responses - selection ######################

# This will take a long time to run (about a day)

model = matrix(NA,p,5)
colnames(model) = c('mclust','flexmix','NJbic','FLbic','RLbic')
n.cluster = detectCores()
cl = n.cluster
cl = makeCluster(cl)
registerDoParallel(cl)

for (j in 1:p) {

  ### mclust ###

  W = Mclust(XX[,-j], G = 2:4, verbose = FALSE)
  model[j,1] = W$G

  ### flexmix ###

  XX1 = XX
  XX1 = as.data.frame(XX1)
  colnames(XX1)[j] = "y"
  colnames(XX1)[-j] = paste("V", 1:(p-1), sep="")
  fmla1 = as.formula(paste("y ~ 1+", paste(colnames(XX1)[-j], collapse= "+")))
  moeflexmix = stepFlexmix(formula = fmla1,k = 2:4,
              model = FLXMRglmnet(family = 'gaussian', adaptive=FALSE),concomitant =  FLXPmultinom(formula = ~1), data = XX1, nrep=1)
  model[j,2] = length(getModel(moeflexmix, "BIC")@size)

  # ## RJM methods ###

  modelNJ2 = rjm(XX[,j],XX[,-j],K=2,method='nj',parallel=TRUE)
  modelNJ3 = rjm(XX[,j],XX[,-j],K=3,method='nj',parallel=TRUE)
  modelNJ4 = rjm(XX[,j],XX[,-j],K=4,method='nj',parallel=TRUE)
  modelFL2 = rjm(XX[,j],XX[,-j],K=2,method='lasso',lasso.pen='fixed',parallel=TRUE)
  modelFL3 = rjm(XX[,j],XX[,-j],K=3,method='lasso',lasso.pen='fixed',parallel=TRUE)
  modelFL4 = rjm(XX[,j],XX[,-j],K=4,method='lasso',lasso.pen='fixed',parallel=TRUE)
  modelRL2 = rjm(XX[,j],XX[,-j],K=2,method='lasso',lasso.pen='random',parallel=TRUE)
  modelRL3 = rjm(XX[,j],XX[,-j],K=3,method='lasso',lasso.pen='random',parallel=TRUE)
  modelRL4 = rjm(XX[,j],XX[,-j],K=4,method='lasso',lasso.pen='random',parallel=TRUE)

  model[j,5] = which.max(c(modelNJ2$BIC,modelNJ3$BIC,modelNJ4$BIC))+1
  model[j,6] = which.max(c(modelFL2$BIC,modelFL3$BIC,modelFL4$BIC))+1
  model[j,7] = which.max(c(modelRL2$BIC,modelRL3$BIC,modelRL4$BIC))+1

}
stopCluster(cl)
apply(model, 2, table)

############################### PLOTS ###################################


pdf('cancer_types.pdf',width = 9, height=4)
layout(matrix(c(1,2,3,1,4,5),2,3,byrow = TRUE))
plot(z,y, ylab='Response gene (NAPSA 9476)',xlab='',xaxt='n', bty='l', cex.lab = 1.4,
     col=c(rep(2,n.samples[1]),rep(7,n.samples[2]),rep(3,n.samples[3]),
     rep(4,n.samples[4])),pch=16)
legend('bottomright',legend = c('BRCA','KIRC','LUAD','THCA'),lty = 1,
       lwd=2,col = c(2,7,3,4),bty = 'n', cex = 1.2)
mtext("(a)",1,cex=0.8,font=2,outer=FALSE,line=0.5)
image(cov(X[z==1,]), main='BRCA', xaxt='n',yaxt='n')
mtext("(b)",1,cex=0.8,font=2,outer=FALSE,line=0.5)
image(cov(X[z==2,]), main='KIRC', xaxt='n',yaxt='n')
mtext("(c)",1,cex=0.8,font=2,outer=FALSE,line=0.5)
image(cov(X[z==3,]), main='LUAD', xaxt='n',yaxt='n')
mtext("(d)",1,cex=0.8,font=2,outer=FALSE,line=0.5)
image(cov(X[z==4,]), main='THCA', xaxt='n',yaxt='n')
mtext("(e)",1,cex=0.8,font=2,outer=FALSE,line=0.5)
dev.off()

MAX = max(abs(betasNJ),abs(betasFL),abs(betasRL))
pdf('coef.pdf',width = 8, height=4)
par(mfrow=c(1,4))
plot(1:(K*(p-1)),betasNJ, col = c(rep('green3', sum(betasNJ!=0)),rep('tomato1', sum(betasNJ==0))),
     pch = 19, ylim = c(0,MAX),bty='l',ylab='Absolute regression coefficients',xlab='',cex.lab=1.4,main='RJM-NJ',cex.main=1.4)
plot(1:(K*(p-1)),betasFL, col = c(rep('green3', sum(betasFL!=0)),rep('tomato1', sum(betasFL==0))),
     pch = 19, ylim = c(0,MAX),bty='l',ylab='Absolute regression coefficients',xlab='',cex.lab=1.4,main='RJM-FLasso',cex.main=1.4)
plot(1:(K*(p-1)),betasRL, col = c(rep('green3', sum(betasRL!=0)),rep('tomato1', sum(betasRL==0))),
     pch = 19, ylim = c(0,MAX),bty='l',ylab='Absolute regression coefficients',xlab='',cex.lab=1.4,main='RJM-RLasso',cex.main=1.4)
dev.off()

pdf('clustering_all.pdf',width = 16, height = 6)
par(bty='l',cex.lab=1.8,cex.axis=1.2)
vioplot(RAND[, c(1:2,4,3,7,5:6,8:10)],ylab='',col=paletteer_c("ggthemes::Red-Green-Gold Diverging", 10))
title(ylab="Adjusted Rand Index", line=2.8)
dev.off()
