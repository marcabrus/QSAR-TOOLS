#' Several Regression metrics for QSAR model Validation
#'
#' @name Regression Metrics
#' @aliases rsquared
#'
#' @description this function This node computes certain statistics between
#' a numeric column's values (obs) and predicted (pred) values.
#'
#' @param obs a numeric vector of data
#' @param pred a numeric vector of data
#' @param SetColumn A character vector containing Set labels
#' @param SetLabel String identifying Training set in SetColumn
#'
#' @details
#'
#' @return
#' The function returns a named vector with the following parameters:
#' RMSE -Train R2 - Train MAE - Train CCC - Train k - Train ki - Train RMSE r2 MAE Q2-F1 Q2-F2 Q2-F3 CCC rsquared0 r2m r2mavg r2mDelta k ki F
#' If Setcolumn is not provided only these parameters are given:
#' RMSE r2(Pearson) r2 MAE Q2-F2 r20 r2m r2mavg r2mDelta k ki F CCC
#'
#' @references Gramatica, Paola, and Alessandro Sangion. "A historical excursus on the statistical validation parameters for QSAR models:
#' a clarification concerning metrics and terminology." Journal of chemical information and modeling 56.6 (2016): 1127-1131.
#'
#' @author Cosimo Toma cosimotoma88@gmail.com
#'
#'
#' @examples
#' a<-sample(seq(1,0.001,by = -0.001), 100)
#' b<-sample(seq(1,0.001,by = -0.001), 100)
#' Regression.metrics(a,b)
#'
#'
#' @keywords utilities
#'
#' @export Regression.metrics
#'
Regression.metrics <-
  function (obs,
            pred,
            SetColumn = NULL,
            SetLabel = "Train")  {
    isNA <- is.na(pred)
    pred <- pred[!isNA]
    obs <- obs[!isNA]
    if (is.null(SetColumn)) {
      cat("Set Column not provided, perrformances evaluated on entire set\n")
      if (length(obs) + length(pred) == 0) {
        out <- rep(NA, 3)
      }
      else {
        yPredYObs = sum((pred - mean(pred)) * (obs - mean(obs)))
        yPredYBar2 = sum((pred - mean(pred)) ^ 2)
        yObsYBar2 = sum((obs - mean(obs)) ^ 2)
        r2 <- (yPredYObs / sqrt(yPredYBar2 * yObsYBar2)) ^ 2
        rss <- sum((pred - obs) ^ 2)
        sst <- sum((obs - mean(obs)) ^ 2)
        sres <- rss / (length(pred) - 2)
        #   F<-sst/sres
        F <- (r2 * (length(obs) - 1 - 1)) / ((1 - r2) * 1)
        k <- sum(pred * obs) / sum(pred * pred)
        rss0 <- sum((obs - k * pred) ^ 2)
        sst0 <- sum((obs - mean(obs)) ^ 2)
        #compute with interchanged axes
        ki <- sum(pred * obs) / sum(obs * obs)
        rssi <- sum((pred - ki * obs) ^ 2)
        ssti <- sum((pred - mean(pred)) ^ 2)
        
        if (sst == 0 & rss == 0) {
          Q2F2 <- 1
          rsquared0 <- 1
          
        } else{
          Q2F2 <- 1 - (rss / sst)
          rsquared0 <- 1 - (rss0 / sst0)
          r2m <- r2 * (1 - sqrt(abs(r2 - rsquared0)))
          rsquaredi <- 1 - (rssi / ssti)
          r2mi <- r2 * (1 - sqrt(abs(r2 - rsquaredi)))
          r2mavg <- (r2m + r2mi) / 2
          r2mdelta <- abs(r2m - r2mi)
        }
        mse <- mean((pred - obs) ^ 2)
        mae <- mean(abs(pred - obs))
        CCC <- (2 * sum((obs - mean(obs)) * (pred - mean(pred)))) /
          (sum((obs - mean(obs)) ^ 2) + sum((pred - mean(pred)) ^ 2) +
             (length(obs) * (mean(obs) - mean(pred)) ^ 2))
        out <-
          c(
            sqrt(mse),
            r2,
            mae,
            Q2F2,
            rsquared0,
            r2m,
            r2mavg,
            r2mdelta,
            k,
            ki,
            F,
            CCC
          )
        names(out) <-
          c(
            "RMSE",
            "r2",
            "MAE",
            "Q2-F2",
            "rsquared0",
            "r2m",
            "r2mavg",
            "r2mDelta",
            "k",
            "ki",
            "F",
            "CCC"
          )
      }
    } else{
      cat(sprintf("Set Column provided, performances evaluated on entire set using %s as label\n",SetLabel))
      SetColumn <- SetColumn[!isNA]
      obsTrain = obs[which(SetColumn == SetLabel)]
      predTrain = pred[which(SetColumn == SetLabel)]
      obs = obs[-which(SetColumn == SetLabel)]
      pred = pred[-which(SetColumn == SetLabel)]
      if (length(obs) + length(pred) == 0| length(obsTrain) + length(predTrain) == 0) {
        out <- rep(NA, 3)
      }
      else {
      
      yPredYObs = sum((pred - mean(pred)) * (obs - mean(obs)))
      yPredYBar2 = sum((pred - mean(pred)) ^ 2)
      yObsYBar2 = sum((obs - mean(obs)) ^ 2)
      r2 <- (yPredYObs / sqrt(yPredYBar2 * yObsYBar2)) ^ 2
      rssTrain = sum((obsTrain - predTrain) ^ 2)
      sstTrain =  sum((obsTrain - mean(obsTrain)) ^ 2)
      r2Train <- 1 - rssTrain / sstTrain
      rss <- sum((pred - obs) ^ 2)
      sst <- sum((obs - mean(obs)) ^ 2)
      sres <- rss / (length(pred) - 2)
      F <- sst / sres
      #   F<-sst/sres
      F <- (r2 * (length(obs) - 1 - 1)) / ((1 - r2) * 1)
      k <- sum(pred * obs) / sum(pred * pred)
      kTrain <- sum(predTrain * obsTrain) / sum(predTrain * predTrain)
      rss0 <- sum((obs - k * pred) ^ 2)
      sst0 <- sum((obs - mean(obs)) ^ 2)
      #compute with interchanged axes
      ki <- sum(pred * obs) / sum(obs * obs)
      kiTrain <- sum(predTrain * obsTrain) / sum(obsTrain * obsTrain)
      rssi <- sum((pred - ki * obs) ^ 2)
      ssti <- sum((pred - mean(pred)) ^ 2)
      sst_ext_train <- sum((obs - mean(obsTrain)) ^ 2)
      sst_train <- sum((obsTrain - mean(obsTrain)) ^ 2)
      if (sst == 0 & rss == 0) {
        Q2F2 <- 1
        rsquared0 <- 1
      } else{
        Q2F1 <- 1 - (rss / sst_ext_train)
        Q2F2 <- 1 - (rss / sst)
        Q2F3 <- 1 - ((rss / length(obs)) / ((sst_train / length(obsTrain))))
        rsquared0 <- 1 - (rss0 / sst0)
        r2m <- r2 * (1 - sqrt(abs(r2 - rsquared0)))
        rsquaredi <- 1 - (rssi / ssti)
        r2mi <- r2 * (1 - sqrt(abs(r2 - rsquaredi)))
        r2mavg <- (r2m + r2mi) / 2
        r2mdelta <- abs(r2m - r2mi)
      }
      mse <- mean((pred - obs) ^ 2)
      mseTrain <- mean((predTrain - obsTrain) ^ 2)
      mae <- mean(abs(pred - obs))
      maeTrain <- mean(abs(predTrain - obsTrain))
      CCC <- (2 * sum((obs - mean(obs)) * (pred - mean(pred)))) /
        (sum((obs - mean(obs)) ^ 2) + sum((pred - mean(pred)) ^ 2) +
           (length(obs) * (mean(obs) - mean(pred)) ^ 2))
      CCCTrain <-
        (2 * sum((obsTrain - mean(obsTrain)) * (predTrain - mean(predTrain)))) /
        (sum((obsTrain - mean(obsTrain)) ^ 2) + sum((predTrain - mean(predTrain)) ^
                                                      2) +
           (length(obsTrain) * (mean(obsTrain) - mean(predTrain)) ^ 2))
      out <-
        c(
          sqrt(mseTrain),
          r2Train,
          maeTrain ,
          CCCTrain ,
          kTrain,
          kiTrain,
          sqrt(mse),
          r2,
          mae,
          Q2F1,
          Q2F2,
          Q2F3,
          CCC,
          rsquared0,
          r2m,
          r2mavg,
          r2mdelta,
          k,
          ki,
          F
        )
      
      names(out) <-
        c(
          "RMSE -Train",
          "R2 - Train",
          "MAE - Train",
          "CCC - Train",
          "k - Train",
          "ki - Train",
          "RMSE",
          "r2",
          "MAE",
          "Q2-F1",
          "Q2-F2",
          "Q2-F3",
          "CCC",
          "rsquared0",
          "r2m",
          "r2mavg",
          "r2mDelta",
          "k",
          "ki",
          "F"
        )
      }
    }
    if (any(is.nan(out)))
      out[is.nan(out)] <- NA
    return(out)
  }
