
# v. 0.1
# 27/03/2018
#

# needed libraries
library('MASS')

#' Evaluate the applicability domain of a QSAR model for a new compound with the calculation of the leverage (hi). 
#' New compounds that are above the hi threshold are considered outside the AD. 
#' @reference A. Tropsha, P. Gramatica, V.K. Gombar, 
#' QSAR Comb. Sci. 22 (2003) 69-77 & A. Golbraikh, A. Tropsha, J. Mol. Graph. Model. 20 (2002) 269-276
#'
#' @param indicesTrain  matrix with descriptors of the Training Set
#' @param indicesTest Matrix with descriptors of the Test Set
#'
#' @return a vector woth the same number of elements of the test set 
#' @author Cosimo Toma cosimotoma88@gmail.com
#' assigning the label 'reliable' if the compound is under the hi threshold
#' and 'unreliable' if the compound is upper the hi treshold
#'
leverage.calculation = function ( indicesTrain, indicesTest){
require(MASS)
A<-as.matrix(indicesTrain)
B<-as.matrix(indicesTest)
temp = t(A)%*% A
temp = ginv(temp)
temp = B%*% temp
h1 = temp %*% t(B)
h = vector()
for (i in 1:ncol(h1)) {
  h[i] = h1[i,i]
}
n = nrow(A)
k = ncol(A)
limit = (3.0 * (k / n))
domain = vector()
for (i in 1:length(h)) {
  domain[i] = h[i]
}
Domain<- ifelse(domain> limit,"unreliable","reliable")
return(Domain)
}
