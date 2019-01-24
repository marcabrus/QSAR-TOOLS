# v. 0.1
# 10/01/2018
#
#' Calculates MCC also for Multiclass models. The formula is the same reported in the following article:
#' Jurman, Giuseppe, Samantha Riccadonna, and Cesare Furlanello. 
#' "A comparison of MCC and CEN error measures in multi-class prediction." PloS one 7.8 (2012): e41882.
#' 
#' It takes two factors, one for the observed values one for the predicted
#'
#' @param obs factor with observed values
#' @param pred factor with predicted values
#'
#' @return MCC
#'

multi_MCC <- function(obs,pred) {
  cont_matr<-as.data.frame.matrix(table(pred, obs))
  colnames(cont_matr)<-paste(colnames(cont_matr), "Ref", sep = "_")
  rownames(cont_matr)<-paste(rownames(cont_matr), "Pred", sep = "_")
  covxy = 0
  for (k in 1:length(cont_matr)[1]) {
    for (l in 1:length(cont_matr)[1])
      for (m in 1:length(cont_matr)[1])
        covxy <-
          ((cont_matr[k, k] * cont_matr[m, l]) - (cont_matr[l, k] * cont_matr[k, m])) + covxy
  }
  
  covxx = 0
  for (k in 1:length(cont_matr)[1]) {
    covxxa = 0
    for (l in 1:length(cont_matr)[1]) {
      covxxa <- covxxa + cont_matr[l, k]
    }
    covxxb = 0
    for (f in 1:length(cont_matr)[1]) {
      if (f == k) {
        next
      }
      
      for (g in 1:length(cont_matr)[1]) {
        covxxb <- covxxb + cont_matr[g, f]
      }
      
    }
    
    covxx <- covxx + (covxxa * covxxb)
  }
  
  covyy = 0
  for (k in 1:length(cont_matr)[1]) {
    covyya = 0
    for (l in 1:length(cont_matr)[1]) {
      covyya <- covyya + cont_matr[k, l]
    }
    covyyb = 0
    for (f in 1:length(cont_matr)[1]) {
      if (f == k) {
        next
      }
      
      for (g in 1:length(cont_matr)[1]) {
        covyyb <- covyyb + cont_matr[f, g]
      }
      
    }
    
    covyy <- covyy + (covyya * covyyb)
  }

  MCC <- covxy / (sqrt(covxx) * sqrt(covyy))
#print MCC  
  message(" ")
  cat(paste(" MCC : ",format(MCC,digits=2)))
return(MCC)

  }
aa<-multi_MCC(obs,pred)
