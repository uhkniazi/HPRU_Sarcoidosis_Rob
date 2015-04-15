# Name: Sar_rob_header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 13 April 2015
# Desc: header file with variables, libraries and functions


# libraries used
library(MASS)
library(randomForest)
library(ROCR)
library(leaps)
library(caret)


## class to hold the data
setClass('CData', slots=list(df='data.frame', f1='factor', f2='factor'))

CData = function(dfDat, fG1, fG2){
  new('CData', df=dfDat, f1=fG1, f2=fG2)
}