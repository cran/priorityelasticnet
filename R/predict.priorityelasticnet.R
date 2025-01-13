utils::globalVariables("missing_index_overview")

#' Predictions from priorityelasticnet
#'
#' Makes predictions for a \code{priorityelasticnet} object. It can be chosen between linear predictors or fitted values.
#'
#' \code{handle.missingtestdata} specifies how to deal with missing data.
#' The default \code{none} cannot handle missing data, \code{omit.prediction} does not make a prediction for observations with missing values and return \code{NA}. \code{set.zero} ignores
#' the missing data for the calculation of the prediction (the missing value is set to zero).
#' \code{impute.block} uses an imputation model to impute the offset of a missing block. This only works if the priorityelasticnet object was fitted with \code{handle.missingdata = "impute.offset"}.
#' If \code{impute.offset.cases = "complete.cases"} was used, then every observation can have only one missing block. For observations with more than one missing block, \code{NA} is returned.
#' If \code{impute.offset.cases = "available.cases"} was used, the missingness pattern in the test data has to be the same as in the train data. For observations with an unknown missingness pattern, \code{NA} is returned.
#'
#' @param object An object of class \code{priorityelasticnet}.
#' @param newdata (nnew \code{x} p) matrix or data frame with new values.
#' @param type Specifies the type of predictions. \code{link} gives the linear predictors for all types of response and \code{response} gives the fitted values.
#' @param handle.missingtestdata Specifies how to deal with missing data in the test data; possibilities are \code{none}, \code{omit.prediction}, \code{set.zero} and \code{impute.block}
#' @param include.allintercepts should the intercepts from all blocks included in the prediction? If \code{FALSE}, only the intercept from the first block is included (default in the past).
#' @param use.blocks determines which blocks are used for the prediction, the default is all. Otherwise one can specify the number of blocks which are used in a vector
#' @param alpha Elastic net mixing parameter used in the model fitting.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Predictions that depend on \code{type}.
#'
#' @importFrom checkmate assert_logical assert_numeric assert check_matrix check_data_frame
#' @importFrom stats predict
#' @import glmnet
#' @examples
#' \donttest{
#' pl_bin <- priorityelasticnet(X = matrix(rnorm(50*190),50,190), Y = rbinom(50,1,0.5),
#'                        family = "binomial", type.measure = "auc",
#'                        blocks = list(block1=1:13,block2=14:80, block3=81:190),
#'                        block1.penalization = TRUE, lambda.type = "lambda.min",
#'                        standardize = FALSE, nfolds = 3, alpha = 1)
#'
#' newdata_bin <- matrix(rnorm(10*190),10,190)
#'
#' predict(object = pl_bin, newdata = newdata_bin, type = "response", alpha = 1)
#' }
#' @export
predict.priorityelasticnet <- function(object,
                                       newdata = NULL,
                                       type = c("link", "response"),
                                       handle.missingtestdata = c("none",
                                                                  "omit.prediction",
                                                                  "set.zero",
                                                                  "impute.block"),
                                       include.allintercepts = FALSE,
                                       use.blocks = "all",
                                       alpha = 1,
                                       ...) {
  
  # Input checks and setup
  type <- match.arg(type)
  handle.missingtestdata <- match.arg(handle.missingtestdata)
  assert_logical(include.allintercepts)
  if (use.blocks[1] != "all" || length(use.blocks) != 1) {
    assert_numeric(use.blocks, min.len = 1, max.len = length(object$blocks),
                   unique = TRUE, lower = 1, upper = length(object$blocks))
  }
  
  if (is.null(newdata)) {
    if (!inherits(object$X, "matrix") || all(is.na(object$X))) {
      stop("No data provided by either the priorityelasticnet object or newdata.")
    } else {
      newdata <- as.matrix(object$X)
    }
  } else {
    assert(check_matrix(newdata), check_data_frame(newdata))
    newdata <- data.matrix(newdata)
  }
  
  if (ncol(newdata) != object$dim.x[2]) {
    stop("The newdata does not have the same number of covariates as the traindata in the priorityelasticnet object.")
  }
  
  if (handle.missingtestdata == "none" && sum(is.na(newdata)) > 0) {
    stop("X contains missing data. Please use another value than 'none' for handle.missingtestdata.")
  }
  
  # Check for single missing values within blocks
  lapply(object$blocks, function(block) {
    lapply(seq_len(nrow(newdata)), function(i) {
      if (sum(is.na(newdata[i, block])) != 0 &&
          sum(is.na(newdata[i, block])) != length(block)) {
        stop(paste0("Observation ", i, " contains a single missing value. This is not supported."))
      }
    })
  })
  
  # Check for observations with all missing values
  lapply(seq_len(nrow(newdata)), function(i) {
    if (sum(is.na(newdata[i, ])) == ncol(newdata)) {
      stop(paste0("Observation ", i, " only consists of missing values."))
    }
  })
  
  if (handle.missingtestdata != "none" &&
      sum(complete.cases(newdata)) == nrow(newdata)) {
    stop("The data consists only of complete cases. Please use handle.missingtestdata = 'none'.")
  }
  
  # Prepare coefficients and intercepts for different families
  if (object$family == "multinomial") {
    intercepts <- lapply(seq_along(object$glmnet.fit), function(i) {
      if (!is.null(object$glmnet.fit[[i]]$a0)) {
        object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]], , drop = FALSE]
      } else {
        matrix(0, nrow = 1, ncol = length(levels(object$Y)))
      }
    })
  } else if (object$family == "cox") {
    intercepts <- rep(0, times = length(object$blocks))
  } else {
    if (!is.null(object$block1unpen)) {
      intercepts <- sapply(1:length(object$glmnet.fit), function(i) {
        object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      })
    } else {
      intercepts <- sapply(1:length(object$glmnet.fit), function(i) {
        object$glmnet.fit[[i]]$a0[object$lambda.ind[[i]]]
      })
    }
  }
  
  # Initialize intercept model matrix
  intercept_model_matrix <- matrix(0, nrow = nrow(newdata),
                                   ncol = length(object$blocks))
  
  if (include.allintercepts) {
    if (use.blocks[1] == "all") {
      block_index <- 1:length(object$blocks)
    } else {
      block_index <- use.blocks
    }
    intercept_model_matrix[, block_index] <- 1
  } else {
    if (use.blocks[1] == "all" || 1 %in% use.blocks) {
      intercept_model_matrix[, 1] <- 1
    }
  }
  
  if (handle.missingtestdata == "set.zero") {
    for (i in seq_len(nrow(intercept_model_matrix))) {
      for (j in seq_len(ncol(intercept_model_matrix))) {
        if (missing_index_overview[i, j]) {
          intercept_model_matrix[i, j] <- 0
        }
      }
    }
  }
  
  # Determine the coefficients
  if (object$family == "multinomial") {
    pred <- matrix(0, nrow(newdata), length(levels(object$Y)))
    for (i in 1:length(object$glmnet.fit)) {
      if (use.blocks[1] == "all" || i %in% use.blocks) {
        block_pred <- newdata[, object$blocks[[i]], drop = FALSE] %*% object$glmnet.fit[[i]]$beta[, object$lambda.ind[[i]]]
        pred <- pred + block_pred + intercepts[[i]]
      }
    }
  } else if (object$family == "cox") {
    coeff <- object$coefficients
    pred <- newdata %*% coeff
  } else {
    if (is.null(object$block1unpen)) {
      coeff <- object$coefficients
    } else {
      coeff <- object$coefficients[-1]
    }
    pred <- newdata %*% coeff + intercept_model_matrix %*% intercepts
  }
  
  if (type == "response") {
    if (object$family == "binomial") {
      pred <- exp(pred) / (1 + exp(pred)) # fitted probabilities
    } else if (object$family == "multinomial") {
      pred <- exp(pred)
      pred <- pred / rowSums(pred)  # Normalize to obtain probabilities
    } else if (object$family == "cox") {
      pred <- exp(pred) # fitted relative risk (risk score (exp(lp)))
    }
  }
  
  if (object$family == "gaussian" && !is.null(object$y.scale.param)) {
    pred <- pred * object$y.scale.param$sd + object$y.scale.param$mean
  }
  
  return(pred)
}

