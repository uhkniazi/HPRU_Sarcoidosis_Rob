# Sarcoidosis_Rob
Biomarkers from proteomics data of Rob

# Sar_rob_header.R
header file with global variables, libraries and settings. the following libraries are loaded
library(MASS) - for using LDA classifier
library(randomForest) - Random Forest variable selection
library(ROCR) - plotting ROC curves and performance 
library(leaps) - for subset selection of variables
library(caret) - for findCorrelation function to remove colinear variables

# Sar_rob_main_newdata.R
this script starts with loading and formatting new data (training) and old data (test). the NA and INF values are removed for
this analysis, however in a future analysis we can think of replacing them with mean values for those genes. the protein ids
in the original result data included multiple names e.g. xx;xx;xx so we take the first name before the semicolon.
we create 3 objects of class CData which hold the data.frames for new, old, and full datasets along with 2 classification factors: TB, Healthy and Sarcoidosis; and Others (including TB and healthy) and Sarcoidosis. 
We perform the following analysis steps
1- TAG_2 is the place holder where we can choose to load a dataset from the object to use
2- some PCA and HCA plots to do some quality checks on data 
3- Random forest are used to select those proteins that have the most importance in class prediction,
this is done by looking at the distribution 
4- A resampling simulation is done using random forest to calculate the distribution of the error rates of the chosen proteins
5- Those proteins with a mean importance over zero are chosen
6- Colinear variables are removed
7- using training data, variable selection is performed using exhaustive search. The best number and combination of variables
is chosen using test error rate using test data
8- nested 10 fold cross validation is performed on (ALL data) to assess the performance of the LDA classifier on these variables.
9- LDA model is again fit on the training data and prediction done on the test data.
Both results from points 8 and 9 are plotted using ROC curves.
10- Boxplots of these proteins are produced.

# Sar_rob_main_newdata_2.R
similar to previous script but removing the healthy 
