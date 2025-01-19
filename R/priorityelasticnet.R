#' Priority Elastic Net for High-Dimensional Data
#'
#' This function performs penalized regression analysis using the elastic net method, tailored for high-dimensional data with a known group structure. It also includes an optional feature to launch a Shiny application for model evaluation with weighted threshold optimization.
#'
#' @param X A numeric matrix of predictors.
#' @param Y A response vector. For family = "multinomial", Y should be a factor with more than two levels.
#' @param weights Optional observation weights. Default is NULL.
#' @param family A character string specifying the model type. Options are "gaussian", "binomial", "cox", and "multinomial". Default is "gaussian".
#' @param alpha The elastic net mixing parameter, with \eqn{0 \le \alpha \le 1}. The penalty is defined as \eqn{(1-\alpha)/2||\beta||_2^2 + \alpha||\beta||_1}. Default is 1.
#' @param type.measure Loss function for cross-validation. Options are "mse", "deviance", "class", "auc". Default depends on the family.
#' @param blocks A list where each element is a vector of indices indicating the predictors in that block.
#' @param max.coef A numeric vector specifying the maximum number of non-zero coefficients allowed in each block. Default is NULL, meaning no limit.
#' @param block1.penalization Logical. If FALSE, the first block will not be penalized. Default is TRUE.
#' @param lambda.type Type of lambda to select. Options are "lambda.min" or "lambda.1se". Default is "lambda.min".
#' @param standardize Logical flag for variable standardization, prior to fitting the model. Default is TRUE.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param foldid Optional vector of values between 1 and \code{nfolds} identifying what fold each observation is in. Default is NULL.
#' @param cvoffset Logical. If TRUE, a cross-validated offset is used. Default is FALSE.
#' @param cvoffsetnfolds Number of folds for cross-validation of the offset. Default is 10.
#' @param mcontrol Control parameters for handling missing data. Default is \code{missing.control()}.
#' @param scale.y Logical. If TRUE, the response variable Y is scaled. Default is FALSE.
#' @param return.x Logical. If TRUE, the function returns the input matrix X. Default is TRUE.
#' @param adaptive Logical. If \code{TRUE}, the adaptive elastic net is used, where penalties are adjusted based on the importance of the coefficients from an initial model fit. Default is \code{FALSE}.
#' @param initial_global_weight Logical. If TRUE (the default), global initial weights will be calculated based on all predictors. If FALSE, initial weights will be calculated separately for each block.
#' @param verbose Logical. If TRUE prints detailed logs of the process. Default is FALSE.
#' @param ... Additional arguments to be passed to \code{cv.glmnet}.
#' 
#' @return A list with the following components:
#' \item{lambda.ind}{Indices of the selected lambda values.}
#' \item{lambda.type}{Type of lambda used.}
#' \item{lambda.min}{Selected lambda values.}
#' \item{min.cvm}{Cross-validated mean squared error for each block.}
#' \item{nzero}{Number of non-zero coefficients for each block.}
#' \item{glmnet.fit}{Fitted \code{glmnet} objects for each block.}
#' \item{name}{Name of the model.}
#' \item{block1unpen}{Fitted model for the unpenalized first block, if applicable.}
#' \item{coefficients}{Coefficients of the fitted models.}
#' \item{call}{The function call.}
#' \item{X}{The input matrix X, if \code{return.x} is TRUE.}
#' \item{missing.data}{Logical vector indicating missing data.}
#' \item{imputation.models}{Imputation models used, if applicable.}
#' \item{blocks.used.for.imputation}{Blocks used for imputation, if applicable.}
#' \item{missingness.pattern}{Pattern of missing data, if applicable.}
#' \item{y.scale.param}{Parameters for scaling Y, if applicable.}
#' \item{blocks}{The input blocks.}
#' \item{mcontrol}{Control parameters for handling missing data.}
#' \item{family}{The model family.}
#' \item{dim.x}{Dimensions of the input matrix X.}
#'
#' @note Ensure that \code{glmnet} version >= 2.0.13 is installed. The function does not support single missing values within a block.
#'
#' @import glmnet
#' @import survival
#' @import glmSparseNet
#' @importFrom checkmate assert_logical
#' @importFrom stats complete.cases coef model.matrix
#' @importFrom utils packageVersion
#' @importFrom shiny shinyApp
#' @importFrom stats glm predict sd
#' @export
#'
#' @examples
#' \donttest{
#'
#'   # Simulation of multinomial data:
#'   set.seed(123)
#'   n <- 100
#'   p <- 50
#'   k <- 3
#'   x <- matrix(rnorm(n * p), n, p)
#'   y <- sample(1:k, n, replace = TRUE)
#'   y <- factor(y)
#'   blocks <- list(bp1 = 1:10, bp2 = 11:30, bp3 = 31:50)
#'   
#'   # Run priorityelasticnet:
#'   fit <- priorityelasticnet(x, y, family = "multinomial", alpha = 0.5, 
#'                      type.measure = "class", blocks = blocks,
#'                      block1.penalization = TRUE, lambda.type = "lambda.min", 
#'                      standardize = TRUE, nfolds = 5, 
#'                      adaptive = FALSE)
#'                      
#'    fit$coefficients
#' }


priorityelasticnet <- function(X,
                               Y,
                               weights = NULL,
                               family = c("gaussian", "binomial", "cox", "multinomial"),
                               alpha = 0.5,
                               type.measure,
                               blocks,
                               max.coef = NULL,
                               block1.penalization = TRUE,
                               lambda.type = "lambda.min",
                               standardize = TRUE,
                               nfolds = 10,
                               foldid = NULL,
                               cvoffset = FALSE,
                               cvoffsetnfolds = 10,
                               mcontrol = missing.control(),
                               scale.y = FALSE,
                               return.x = TRUE,
                               adaptive = FALSE,
                               initial_global_weight = TRUE,
                               verbose = FALSE, 
                               ...) {
  
  
  if(initial_global_weight){
    initial_weight_scope = "global"
  }else{
    initial_weight_scope = "block-wise"
  }
  
  if (packageVersion("glmnet") < "2.0.13") {
    stop("glmnet >= 2.0.13 needed for this function.", call. = FALSE)
  }
  
  if (verbose) {
    message("Starting priorityelasticnet with ", length(blocks), " blocks.")
  }
  
  if (is.null(max.coef)) {
    max.coef <- rep(+Inf, length(blocks))
  } else {
    if (min(max.coef) < +Inf && lambda.type == "lambda.1se") {
      warning("lambda.1se can only be chosen without restrictions of max.coef and is set to lambda.min.")
      lambda.type <- "lambda.min"
    }
    if (!setequal(length(blocks), length(max.coef))) {
      stop("The length of the entries of argument max.coef must equal the number of blocks.")
    }
  }
  
  if (sum(lapply(blocks, length) <= 1) != 0) {
    stop("A block has to contain at least two predictors.")
  }
  
  if (anyDuplicated(as.numeric(unlist(blocks))) != 0 || !setequal(as.numeric(unlist(blocks)), 1:ncol(X))) {
    stop("Each predictor should be included in exactly one block.")
  }
  
  if (!is.element(lambda.type, c("lambda.min", "lambda.1se"))) {
    stop("lambda.type must be either lambda.min or lambda.1se.")
  }
  
  if (verbose) {
    message("Checking family type and setting default type.measure if necessary...")
  }
  
  # Handle family-specific default settings for type.measure
  if (family == "gaussian") {
    if (type.measure != "mse")
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  if (family == "cox") {
    if (type.measure != "deviance")
      warning("type.measure is set to partial likelihood.")
    type.measure <- "deviance"
  }
  if (family == "multinomial") {
    if (type.measure != "class")
      warning("type.measure is set to class.")
    type.measure <- "class"
  }
  
  # Adaptive Elastic Net: Initial model to get adaptive weights
  if (adaptive) {
    if (verbose) {
      message("Calculating adaptive weights based on an initial model using ", initial_weight_scope, " approach...")
    }
    
    if (initial_weight_scope == "global") {
      # Global initial weights calculation
      initial_model <- glmnet(X, Y, family = family, alpha = alpha, standardize = standardize, ...)
      lambda_min_value <- min(initial_model$lambda)
      if(family == "cox"){
        initial_coeff <- abs(coef(initial_model, s = lambda_min_value))[,1]
       
      }
        
      else if(family == "multinomial"){
          
        initial_coeff <- abs(coef(initial_model, s = lambda_min_value)[[1]][,1])[-1]
          
      }else{
        
        initial_coeff <- abs(coef(initial_model, s = lambda_min_value))[-1] 
        
      }
      
      # Exclude intercept if present
      initial_coeff[initial_coeff == 0] <- 1e-6  # Avoid division by zero for zero coefficients
      adaptive_weights <- 1 / initial_coeff
      
    } else if (initial_weight_scope == "block-wise") {
      # Block-wise initial weights calculation
      adaptive_weights <- rep(NA, ncol(X))  # Preallocate
      for (i in seq_along(blocks)) {
        block_indices <- blocks[[i]]
        initial_model_block <- glmnet(X[, block_indices, drop = FALSE], Y, family = family, alpha = alpha, standardize = standardize, ...)
        lambda_min_value <- min(initial_model_block$lambda)
        if(family == "cox"){
          initial_coeff <- abs(coef(initial_model_block, s = lambda_min_value))
        }
        
        else if(family == "multinomial"){
          
          initial_coeff <- abs(coef(initial_model_block, s = lambda_min_value)[[1]][,1])[-1]
          
        }else{
          initial_coeff <- abs(coef(initial_model_block, s = lambda_min_value))[-1]
        }
        initial_coeff[initial_coeff == 0] <- 1e-6  # Avoid division by zero
        adaptive_weights[block_indices] <- 1 / initial_coeff
      }
    }
    
    if (verbose) {
      message("Adaptive weights calculated.")
    }
  } else {
    adaptive_weights <- rep(1, ncol(X))  # No adaptation, so weights are uniform
  }
  
  
  if (type.measure == "auc") {
    if (cvoffset) {
      if (nrow(X) * ((cvoffsetnfolds - 1) / cvoffsetnfolds) - nrow(X) * ((cvoffsetnfolds - 1) / cvoffsetnfolds) * (nfolds - 1) / nfolds < 10) {
        stop(paste("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; Use nfolds < ", floor(nrow(X) * ((cvoffsetnfolds - 1) / cvoffsetnfolds) / 10) + 1, ".", sep = ""))
      }
    } else {
      if (nrow(X) - nrow(X) * (nfolds - 1) / nfolds < 10)
        stop(paste("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; Use nfolds < ", floor(nrow(X) / 10) + 1, ".", sep = ""))
    }
  }
  
  if (is.null(weights)) {
    weights <- rep(1, nrow(X))
  } else {
    if (length(weights) != nrow(X))
      stop(paste("number of elements in weights (", length(weights),
                 ") not equal to the number of rows of X (", nrow(X),
                 ")", sep = ""))
  }
  
  if (!is.null(foldid)) {
    if (length(foldid) != nrow(X))
      stop(paste("number of elements in foldid (", length(foldid),
                 ") not equal to the number of rows of X (", nrow(X),
                 ")", sep = ""))
    else {
      if (nfolds != max(foldid)) {
        warning(paste("nfolds is set to", max(foldid)))
        nfolds <- max(foldid)
      }
    }
  } else {
    foldid <- sample(rep(seq(nfolds), length = nrow(X)))
  }
  
  if (verbose) {
    message("Handling missing data based on the provided mcontrol parameters...")
  }
  
  if (mcontrol$handle.missingdata == "none" && sum(is.na(X)) > 0) {
    stop("X contains missing data. Please use another value than 'none' for handle.missingdata.")
  }
  if (mcontrol$handle.missingdata != "none" && cvoffset) {
    stop("At the moment, a crossvalidated offset is only supported for complete data sets.")
  }
  
  if (mcontrol$handle.missingdata == "ignore" || mcontrol$handle.missingdata == "impute.offset") {
    foldid <- NULL
    warning(paste0("For handle.missingdata = ", mcontrol$handle.missingdata, ", the foldids of the observations are chosen individually for every block and not set globally. foldid is set to NULL"))
  }
  if (mcontrol$handle.missingdata == "impute.offset" && mcontrol$impute.offset.cases == "complete.cases") {
    perc_complete_cases <- sum(complete.cases(X)) / nrow(X)
    
    if (sum(complete.cases(X)) == 0) {
      stop("The dataset contains no complete cases (over all blocks). Imputation of the offsets not possible.")
    }
    if (perc_complete_cases < mcontrol$perc.comp.cases.warning) {
      warning(paste0("The fraction of complete cases only is ", round(perc_complete_cases, digits = 2)))
    }
    
    missing_index_overview <- matrix(FALSE, nrow = nrow(X), ncol = length(blocks))
    for (i in seq_along(blocks)) {
      missing_index_overview[, i] <- !complete.cases(X[, blocks[[i]]])
    }
    for (i in seq_len(nrow(missing_index_overview))) {
      if (sum(missing_index_overview[i, ]) > 1) {
        stop("For impute.offset.cases = 'complete.cases', every observation must only contain one missing block.")
      }
    }
  }
  
  if (mcontrol$handle.missingdata == "ignore" && mcontrol$offset.firstblock == "intercept" && family == "cox") {
    stop("offset.firstblock = 'intercept' can't be used with family = 'cox' as cox models don't fit an intercept")
  }
  
  if (scale.y && family != "gaussian") {
    stop("scale.y = TRUE can only be used with family = 'gaussian'")
  }
  
  if (scale.y) {
    y.scale.param <- list(mean = mean(Y), sd = sd(Y))
    Y <- scale(Y)
  } else {
    y.scale.param <- NULL
  }
  
  lapply(blocks, function(block) {
    lapply(seq_len(nrow(X)), function(i) {
      if (sum(is.na(X[i, block])) != 0 && sum(is.na(X[i, block])) != length(block)) {
        stop(paste0("Observation ", i, " contains a single missing value. This is not supported."))
      }
    })
  })
  
  observation_index <- lapply(blocks, function(block) {
    result <- which(complete.cases(X[, block]))
    if (length(result) == 0) {
      NULL
    } else {
      result
    }
  })
  missing_index <- lapply(blocks, function(block) {
    result <- which(!complete.cases(X[, block]))
    if (length(result) == 0) {
      NULL
    } else {
      result
    }
  })
  
  lambda.min <- list()
  lambda.ind <- list()
  min.cvm <- list()
  nzero <- list()
  glmnet.fit <- list()
  coeff <- list()
  lassoerg <- list()
  liste <- list(NULL)
  imputation_models <- list()
  blocks_used_for_imputation <- list()
  missingness_pattern <- list()
  missing.data <- list()
  start_block <- 1
  
  if (!block1.penalization) {
    if (length(blocks[[1]]) >= nrow(X)) {
      stop("An unpenalized block 1 is only possible if the number of predictors in this block is smaller than the number of observations.")
    }
    current_observations <- observation_index[[1]]
    current_missings <- missing_index[[1]]
    missing.data[[1]] <- seq(nrow(X)) %in% current_missings
    
    if (family != "cox") {
      block1erg <- glm(Y[current_observations] ~ X[current_observations, blocks[[1]]],
                       family = family,
                       weights = weights[current_observations])
      predict_type <- "link"
    } else {
      block1erg <- coxph(Y[current_observations, ] ~ X[current_observations, blocks[[1]]],
                         weights = weights[current_observations],
                         model = TRUE)
      predict_type <- "lp"
    }
    names(block1erg$coefficients) <- substr(names(block1erg$coefficients), start = 37, nchar(names(block1erg$coefficients)))
    
    if (cvoffset) {
      datablock1 <- data.frame(X[, blocks[[1]], drop = FALSE])
      datablock1$Y <- Y
      
      cvdiv <- makeCVdivision(n = nrow(X), K = cvoffsetnfolds, nrep = 1)[[1]]
      pred <- matrix(nrow = nrow(X), ncol = ifelse(family == "multinomial", length(levels(Y)), 1))
      for (count in seq(along = cvdiv)) {
        if (family != "cox") {
          block1ergtemp <- glm(Y ~ ., data = datablock1[cvdiv[[count]] == 1, ],
                               weights = weights[cvdiv[[count]] == 1],
                               family = family)
        } else {
          block1ergtemp <- coxph(Y ~ ., data = datablock1[cvdiv[[count]] == 1, ],
                                 weights = weights[cvdiv[[count]] == 1])
        }
        
        names(block1ergtemp$coefficients) <- substr(names(block1ergtemp$coefficients), start = 17, nchar(names(block1ergtemp$coefficients)))
        pred[cvdiv[[count]] == 0, ] <- as.matrix(predict(block1ergtemp, newdata = datablock1[cvdiv[[count]] == 0, ]), type = predict_type)
      }
    } else {
      pred <- as.matrix(predict(block1erg, type = predict_type))
    }
    
    start_block <- 2
    result_offsets <- calculate_offsets(current_missings = current_missings,
                                        current_observations = current_observations,
                                        mcontrol = mcontrol,
                                        current_block = 1,
                                        pred = pred,
                                        liste = liste,
                                        X = X,
                                        blocks = blocks,
                                        current_intercept = coef(block1erg)[1])
    liste[[2]] <- result_offsets[["new_offsets"]]
    imputation_models[[1]] <- result_offsets[["imputation_model"]]
    blocks_used_for_imputation[[1]] <- result_offsets[["blocks_used_for_imputation"]]
    missingness_pattern[[1]] <- result_offsets[["missingness_pattern"]]
    lassoerg <- list(block1erg)
    coeff[[1]] <- block1erg$coefficients
    
    if (verbose) {
      message("Finished processing unpenalized first block.")
    }
  } else {
    block1erg <- NULL
  }
  
  for (i in start_block:length(blocks)) {
    actual_block <- blocks[[i]]
    current_observations <- observation_index[[i]]
    current_missings <- missing_index[[i]]
    missing.data[[i]] <- seq(nrow(X)) %in% current_missings
    
    if (!is.null(liste[[i]]) && family == "multinomial") {
      offset_matrix <- matrix(liste[[i]][current_observations], nrow = length(current_observations), ncol = length(levels(Y)), byrow = TRUE)
    } else {
      offset_matrix <- liste[[i]][current_observations]
    }
    
    if (verbose) {
      message("Fitting model for block ", i, "...")
    }
    
    # Apply adaptive weights
    penalty.factor <- adaptive_weights[actual_block]
    
    lassoerg[[i]] <- cv.glmnet(X[current_observations, actual_block],
                               Y[current_observations],
                               weights[current_observations],
                               offset = offset_matrix,
                               family = family,
                               alpha = alpha,
                               type.measure = type.measure,
                               nfolds = nfolds,
                               foldid = foldid[current_observations],
                               standardize = standardize,
                               penalty.factor = penalty.factor,
                               ...)
    
    if (lambda.type == "lambda.1se") {
      lambda_to_use <- "lambda.1se"
      lambda.ind[i] <- which(lassoerg[[i]]$lambda == lassoerg[[i]][lambda.type])
      lambda.min[i] <- lassoerg[[i]][lambda.type]
    } else {
      which_lambda <- which(as.numeric(lassoerg[[i]]$nzero) <= max.coef[i])
      
      if (type.measure != "auc") {
        lcvmi <- lassoerg[[i]]$cvm
      } else {
        lcvmi <- -lassoerg[[i]]$cvm
      }
      
      lambda_to_use <- lassoerg[[i]]$lambda[which_lambda[which.min(lcvmi[which_lambda])[1]]]
      lambda.min[i] <- lambda_to_use
      lambda.ind[i] <- which(lassoerg[[i]]$lambda == lambda.min[i])
    }
    
    if (cvoffset) {
      cvdiv <- makeCVdivision(n = nrow(X), K = cvoffsetnfolds, nrep = 1)[[1]]
      pred <- matrix(nrow = nrow(X), ncol = ifelse(family == "multinomial", length(levels(Y)), 1))
      for (count in seq(along = cvdiv)) {
        if (!is.null(liste[[i]]) && family == "multinomial") {
          offset_matrix <- matrix(liste[[i]][cvdiv[[count]] == 1], nrow = sum(cvdiv[[count]] == 1), ncol = length(levels(Y)), byrow = TRUE)
        } else {
          offset_matrix <- liste[[i]][cvdiv[[count]] == 1]
        }
        
        lassoergtemp <- cv.glmnet(X[cvdiv[[count]] == 1, actual_block, drop = FALSE],
                                  Y[cvdiv[[count]] == 1],
                                  weights[cvdiv[[count]] == 1],
                                  offset = offset_matrix,
                                  family = family,
                                  alpha = alpha,
                                  type.measure = type.measure,
                                  nfolds = nfolds,
                                  foldid = foldid[cvdiv[[count]] == 1],
                                  standardize = standardize,
                                  penalty.factor = penalty.factor,
                                  ...)
        
        if (lambda.type == "lambda.1se") {
          lambda_to_use <- "lambda.1se"
        } else {
          which_lambdatemp <- which(as.numeric(lassoergtemp$nzero) <= max.coef[i])
          
          if (type.measure != "auc") {
            lcvmitemp <- lassoergtemp$cvm
          } else {
            lcvmitemp <- -lassoergtemp$cvm
          }
          
          lambda_to_use <- lassoergtemp$lambda[which_lambdatemp[which.min(lcvmitemp[which_lambdatemp])[1]]]
        }
        
        pred[cvdiv[[count]] == 0, ] <- predict(lassoergtemp,
                                               newx = X[cvdiv[[count]] == 0, actual_block],
                                               newoffset = offset_matrix,
                                               s = lambda_to_use,
                                               type = "link")
      }
    } else {
      pred <- predict(lassoerg[[i]],
                      newx = X[current_observations, actual_block],
                      newoffset = offset_matrix,
                      s = lambda_to_use,
                      type = "link")
    }
    
    result_offsets <- calculate_offsets(current_missings = current_missings,
                                        current_observations = current_observations,
                                        mcontrol = mcontrol,
                                        current_block = i,
                                        pred = pred,
                                        liste = liste,
                                        X = X,
                                        blocks = blocks,
                                        current_intercept = lassoerg[[i]]$glmnet.fit$a0[lambda.ind[[i]]])
    liste[[i + 1]] <- result_offsets[["new_offsets"]]
    imputation_models[[i]] <- result_offsets[["imputation_model"]]
    blocks_used_for_imputation[[i]] <- result_offsets[["blocks_used_for_imputation"]]
    missingness_pattern[[i]] <- result_offsets[["missingness_pattern"]]
    
    min.cvm[i] <- lassoerg[[i]]$cvm[lambda.ind[[i]]]
    nzero[i] <- lassoerg[[i]]$nzero[lambda.ind[[i]]]
    glmnet.fit[[i]] <- lassoerg[[i]]$glmnet.fit
    coeff[[i]] <- if (family == "multinomial") {
      do.call(cbind, lapply(1:length(glmnet.fit[[i]]$beta), function(j) glmnet.fit[[i]]$beta[[j]][, lambda.ind[[i]]]))
    } else {
      glmnet.fit[[i]]$beta[, lambda.ind[[i]]]
    }
    
    if (verbose) {
      message("Finished processing block ", i)
    }
  }
  
  name <- lassoerg[[i]]$name
  
  if (mcontrol$handle.missingdata != "impute.offset") {
    imputation_models <- NULL
  }
  
  if (return.x) {
    x_return_value <- X
  } else {
    x_return_value <- NA
  }
  
  finallist <- list(lambda.ind = lambda.ind,
                    lambda.type = lambda.type,
                    lambda.min = lambda.min,
                    min.cvm = min.cvm,
                    nzero = nzero,
                    glmnet.fit = glmnet.fit,
                    name = name,
                    block1unpen = block1erg,
                    coefficients = if (family == "multinomial") coeff else unlist(coeff),
                    call = match.call(),
                    X = x_return_value,
                    missing.data = missing.data,
                    imputation.models = imputation_models,
                    blocks.used.for.imputation = blocks_used_for_imputation,
                    missingness.pattern = missingness_pattern,
                    y.scale.param = y.scale.param,
                    blocks = blocks,
                    mcontrol = mcontrol,
                    family = family,
                    dim.x = dim(X),
                    pred = pred,
                    actuals = Y,
                    adaptive = adaptive,
                    adaptive_weights = if (adaptive) adaptive_weights else NULL,
                    initial_coeff = if (adaptive) initial_coeff else NULL,
                    initial_weight_scope = initial_weight_scope)
  
  class(finallist) <- c("priorityelasticnet", class(finallist))
  
  if (verbose) {
    message("priorityelasticnet completed successfully.")
  }
  
  return(finallist)
}

