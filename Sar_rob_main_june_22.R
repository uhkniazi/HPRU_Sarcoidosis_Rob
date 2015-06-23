# Name: Sar_rob_main_june_22.R
# Auth: u.niazi@imperial.ac.uk
# Date: 22 June 2015
# Desc: Variable selection for class separation using new dataset

library(MASS)
library(randomForest)
library(ROCR)
library(leaps)

csFile = file.choose()
dfData = read.csv(csFile, header = T)
head(dfData)

dfData = na.omit(dfData)
# create factors 
fGroups = as.character(colnames(dfData))
fGroups = gsub('(\\w+)\\.\\d+', replacement = '\\1', fGroups)
fGroups = fGroups[-length(fGroups)]
cGroups = fGroups
fGroups = factor(fGroups)

# get the last column which is protein names
cvProteins = as.character(dfData$X)
cvProteins = gsub('^(\\w+);+.+', replacement = '\\1', x = cvProteins)
cvProteins[length(cvProteins)] = 'Q9BYE2'

# create data frame to hold data
temp = dfData[,1:(ncol(dfData)-1)]
temp = t(temp)
colnames(temp) = cvProteins
dfData = data.frame(temp)

## import libraries to perform variable selection
source('../CCrossValidation/CCrossValidation.R')

# split the data into training and test sets
test = sample(1:nrow(dfData), size = 10, replace = F)
fGroups = rep(NA, length(cGroups))
i = grep('SA', cGroups)
fGroups[i] = 'SA'
fGroups[-i] = 'OD'
fGroups = factor(fGroups)

oVar.rf = CVariableSelection.RandomForest(data = dfData[-test,], groups = fGroups[-test], big.warn = F)
p.old = par()
plot.var.selection(oVar.rf)

dfVar = CVariableSelection.RandomForest.getVariables(ob = oVar.rf)

# select top 30 variables
cvProteins.top = rownames(dfVar)[1:30]

# perform subset selection
oVar.sub = CVariableSelection.ReduceModel(dfData[-test,cvProteins.top], groups = fGroups[-test], boot.num = 1000)
plot.var.selection(oVar.sub)

cvProteins.top = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 3)

oCV = CCrossValidation.LDA(dfData[test,cvProteins.top], dfData[-test,cvProteins.top], fGroups[test], fGroups[-test], 'SA')
plot.cv.performance(oCV)


# plot multivariable selection
par(mfrow=c(2,2))

for (i in 2:9){
  cvProteins.top = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = i)
  print(paste('Variable size', i))
  print(cvProteins.top)
  oCV = CCrossValidation.LDA(dfData[test,cvProteins.top], dfData[-test,cvProteins.top], fGroups[test], fGroups[-test], 'SA')
  plot.cv.performance(oCV)
}

par(p.old)

write.csv(dfVar, file='Results/june_22_variable_importance_2.csv')


