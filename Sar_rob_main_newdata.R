# Name: Sar_rob_main_newdata.R
# Auth: u.niazi@imperial.ac.uk
# Date: 13 April 2015
# Desc: Variable selection for class separation using new dataset


source('Sar_rob_header.R')

## data loading
csFile = file.choose()
dfData = read.csv(csFile, header = T, row.names=1)
head(dfData)
# create factors
fGroups = gsub('(\\w+)\\.\\d+', '\\1', colnames(dfData))
fGroups = as.factor(fGroups)
table(fGroups)
fGroups.2 = as.character(fGroups)
fGroups.2[fGroups.2 != 'SA'] = 'Others'
fGroups.2 = factor(fGroups.2, levels=c('Others', 'SA'))

# remove NA data for the time being
dfData = na.omit(dfData)
# remove any infinite values
f = is.finite(rowSums(dfData))
dfData = dfData[f,]
# clean rownames for long gene names
rn = rownames(dfData)
rn = gsub('^(\\w+);+.+', replacement = '\\1', x = rn)
rownames(dfData) = rn
oDat.new = CData(dfData, fGroups, fGroups.2)

### load the old data set
csFile = file.choose()
dfData = read.csv(csFile, header = T, row.names=1)
head(dfData)
# create factors
fGroups = gsub('(\\w+)\\.\\d+', '\\1', colnames(dfData))
fGroups = as.factor(fGroups)
table(fGroups)
fGroups.2 = as.character(fGroups)
fGroups.2[fGroups.2 != 'SA'] = 'Others'
fGroups.2 = factor(fGroups.2, levels=c('Others', 'SA'))

# remove NA data for the time being
dfData = na.omit(dfData)
# remove any infinite values
f = is.finite(rowSums(dfData))
dfData = dfData[f,]
# clean rownames for long gene names
rn = rownames(dfData)
rn = gsub('^(\\w+);+.+', replacement = '\\1', x = rn)
rownames(dfData) = rn
oDat.old = CData(dfData, fGroups, fGroups.2)

## set variables choosing them from the class CData object
dfData.bk = oDat.new@df
fGroups = oDat.new@f1
fGroups.2 = oDat.new@f2
## OR
dfData.bk = oDat.old@df
fGroups = oDat.old@f1
fGroups.2 = oDat.old@f2
##
dfData = dfData.bk
# classify data test without scaling
mDat = t(dfData)
mDist = dist(mDat)
hc = hclust(mDist)
plot(hc)

# pr comp quality check
pr.out = prcomp(mDat, scale=T)
# check eigen values
plot(pr.out$sdev^2)
# plot the components
p.old = par(mfrow=c(2,2))
col.p = rainbow(length(unique(fGroups)))
col = col.p[as.numeric(fGroups)]
# plot the pca components
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
plot.new()
legend('center', legend = unique(fGroups), fill=col.p[as.numeric(unique(fGroups))])
par(p.old)

# does scaling the data make a difference
mDat = dfData
dim(mDat)
i = apply(mDat, 2, sd)
head(i)
mDat = sweep(mDat, 2, i, FUN = '/')
apply(mDat, 2, sd)
mDat = t(mDat)
# pr comp quality check
pr.out = prcomp(mDat, scale=T)
# check eigen values
plot(pr.out$sdev^2)
# plot the components
par(mfrow=c(2,2))
col.p = rainbow(length(unique(fGroups)))
col = col.p[as.numeric(fGroups)]
# plot the pca components
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
plot.new()
legend('center', legend = unique(fGroups), fill=col.p[as.numeric(unique(fGroups))])
par(p.old)
#dfData.bk = dfData
### it appears scaling does reduce the separation between the classes

### try variable selection without scaling
# dfData.bk = oDat.new@df
# fGroups = oDat.new@f1
# fGroups.2 = oDat.new@f2
set.seed(1)
dfData = t(dfData.bk)
dfData = data.frame(dfData, fGroups)
dim(dfData)
rf.fit = randomForest(fGroups ~., data=dfData, importance=T, ntree = 10000, mtry = 5)
varImpPlot(rf.fit)
dfImportance = as.data.frame(importance(rf.fit))
dfImportance = dfImportance[order(dfImportance$MeanDecreaseAccuracy, decreasing = T),]
hist(dfImportance$MeanDecreaseAccuracy)
dfImportance.SA = dfImportance[order(dfImportance$SA, decreasing = T),]

# select the top few genes
### TAG 1
# choose the top proteins for SA
hist(dfImportance.SA$SA)
f = which(dfImportance.SA$SA >= 4)
length(f)
cvTopGenes = rownames(dfImportance.SA)[f]

# subset the data based on these selected genes from old dataset
# dfData.bk = oDat.old@df
# fGroups = oDat.old@f1
# fGroups.2 = oDat.old@f2
dfData = dfData.bk[rownames(dfData.bk) %in% cvTopGenes,]
dfData = data.frame(t(dfData))
#cvTopGenes = colnames(dfData)
dfData$fGroups = fGroups.2

### Further variable classification check
### using CV and ROC
dfData.full = dfData
iBoot = 100

## as the 2 class proportions are not equal, fit random forest multiple times on random samples
## containing equal proportions of both classes and check variable importance measures
# fit random forests multiple times
# store results 
lVarImp = vector('list', iBoot)
for (i in 1:iBoot) {
  # get indices of the particular factors in data table
  ind.o = which(dfData.full$fGroups == 'Others')
  ind.p = which(dfData.full$fGroups == 'SA')
  # take sample of equal in size from group Others and SA 
  ind.o.s = sample(ind.o, size = length(ind.p), replace = F)
  # sample of SA groups, i.e. take everything as it is smaller group
  ind.p.s = sample(ind.p, size=length(ind.p), replace=F)
  # join the sample indices together
  ind = sample(c(ind.o.s, ind.p.s), replace=F)
  # take sample from the full dataset
  dfData = dfData.full[ind,]
  # fit model
  fit.rf = randomForest(fGroups ~., data=dfData, importance = TRUE, ntree = 5000)
  # get variables importance
  df = importance(fit.rf)
  df = df[order(df[,'MeanDecreaseAccuracy'], decreasing = T),]
  # put in list
  lVarImp[[i]] = df
} # for

## put data for each boot of each variable together in a dataframe
df = NULL
for (i in 1:iBoot) df = rbind(df, lVarImp[[i]])
# convert rownames i.e. gene names to factors
f = as.factor(rownames(df))
# calculate mean and sd for each gene
tapply(df[,'MeanDecreaseAccuracy'], f, mean)
tapply(df[,'MeanDecreaseAccuracy'], f, sd)
df = as.data.frame(df)
df$Symbol = rownames(df)
# boxplots 
par(mar=c(6,4,3,2)+0.1)
boxplot(df$MeanDecreaseAccuracy ~ df$Symbol, las=2)
rm(dfData)
par(p.old)

## stop and decide
## if we like the variable count here, then its fine, or else go back to TAG 1 and select more variables
# remove the variable of choice after testing
x = tapply(df[,'MeanDecreaseAccuracy'], f, mean)
sort(x)
summary(x)
i = which(x < 4)
n = names(x)[i]
i = which(cvTopGenes %in% n)
cvTopGenes = cvTopGenes[-i]

## look at the correlation of the genes
dfData = dfData.bk[rownames(dfData.bk) %in% cvTopGenes,]
mCor = cor(t(dfData))
i = findCorrelation(mCor, cutoff = 0.7)
n = colnames(mCor)[i]
# remove these correlated features 
i = which(cvTopGenes %in% n)
cvTopGenes = cvTopGenes[-i]


### cross validation with ROC
#### CV with ROC
dfData = dfData.bk[rownames(dfData.bk) %in% cvTopGenes,]
dfData = data.frame(t(dfData))
dfData$fGroups = fGroups.2
### Further variable classification check
### using CV and ROC
dfData.full = dfData
set.seed(1)
lPred = vector(mode = 'list', length = 20)
lLab = vector(mode = 'list', length=20)
iCv.error = NULL
for (oo in 1:20){
  t.lPred = NULL
  t.lLab = NULL
  # select a subset of equal numbers for others and SA
  ind.o = which(dfData.full$fGroups == 'Others')
  ind.p = which(dfData.full$fGroups == 'SA')
  ind.o.s = sample(ind.o, size = length(ind.p), replace = F)
  ind.p.s = sample(ind.p, size=length(ind.p), replace=F)
  ind = sample(c(ind.o.s, ind.p.s), replace=F)
  dfData = dfData.full[ind,]
  for (o in 1:1){
    # perform 10 fold cross validation
    k = 10
    folds = sample(1:k, nrow(dfData), replace = T, prob = rep(1/k, times=k))
    # choose the fold to fit and test the model
    for (i in 1:k){
      # check if selected fold leads to 0 for a class
      if ((length(unique(dfData$fGroups[folds != i])) < 2) || (length(unique(dfData$fGroups[folds == i])) < 2)) next
      # check if fold too small to fit model
      if (nrow(dfData[folds != i,]) < 3) next
      # fit model on data not in fold
      fit = lda(fGroups ~ ., data=dfData[folds != i,])
      # predict on data in fold
      pred = predict(fit, newdata = dfData[folds == i,])$posterior[,'SA']
      name = paste('pred',oo, o, i,sep='.' )
      t.lPred[[name]] = pred
      name = paste('label',oo,o, i,sep='.' )
      t.lLab[[name]] = dfData$fGroups[folds == i] == 'SA'
      pred = predict(fit, newdata = dfData[folds == i,])$class
      iCv.error = append(iCv.error, mean(pred != dfData$fGroups[folds == i]))
    }
  }
  t.lPred = unlist(t.lPred)
  t.lLab = unlist(t.lLab)
  lPred[[oo]] = t.lPred
  lLab[[oo]] = t.lLab
}

pred = prediction(lPred, lLab)
perf = performance(pred, 'tpr', 'fpr')
auc = performance(pred, 'auc')

plot(perf, main=paste('ROC Prediction of for', 'SA'),
     spread.estimate='stddev', avg='vertical', spread.scale=2)
auc = paste('auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
cv = paste('CV Error=', signif(mean(iCv.error), 3))
legend('bottomright', legend = c(auc, cv))
abline(0, 1, lty=2)
