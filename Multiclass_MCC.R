# v. 0.1
# 10/01/2018
#
#' @reference Gorodkin, Jan. "Comparing two K-category assignments by a K-category correlation coefficient."
#' Computational biology and chemistry 28.5-6 (2004): 367-374.
#'
#' It takes two factors, one for the observed values one for the predicted
#'
#' @param obs factor with observed values
#' @param pred factor with predicted values
#' 
#' @author Cosimo Toma cosimotoma88@gmail.com
#'
#' @return MCC
#'

multi_MCC <- function(obs,pred) {
  # create the contingency table with the r base function table
  cont_matr<-as.data.frame.matrix(table(pred, obs))
  # not necessary, change names to contingency table
  colnames(cont_matr)<-paste(colnames(cont_matr), "Ref", sep = "_")
  rownames(cont_matr)<-paste(rownames(cont_matr), "Pred", sep = "_")

  #calculate covariance (obs,pred)
  covxy = 0
  for (k in 1:length(cont_matr)[1]) {
    for (l in 1:length(cont_matr)[1])
      for (m in 1:length(cont_matr)[1])
        covxy <-  ((cont_matr[k, k] * cont_matr[m, l]) - (cont_matr[l, k] * cont_matr[k, m])) + covxy
  }
  #calculate covariance (obs,pred)
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

  #calculate covariance (pred, pred)
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
