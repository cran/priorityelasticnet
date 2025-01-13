#' priorityelasticnet with several block specifications
#'
#' Runs priorityelasticnet for a list of block specifications and gives the best results
#' in terms of cv error.
#'
#' @inheritParams priorityelasticnet
#' @param blocks.list list of the format \code{list(list(bp1=...,bp2=...,), list(bp1=,...,bp2=...,), ...)}. For the specification of the entries, see \code{\link[priorityelasticnet]{priorityelasticnet}}.
#' @param max.coef.list list of \code{max.coef} vectors. The first entries are omitted if \code{block1.penalization = FALSE}. Default is \code{NULL}.
#' @param alpha Elastic net mixing parameter. The elastic net penalty is defined as 
#' \deqn{(1 - \alpha)/2||\beta||_2^2 + \alpha||\beta||_1}{(1 - \alpha)/2||\beta||_2^2 + \alpha||\beta||_1}
#' Defaults to 1 (lasso penalty).
#' @param ... other arguments that can be passed to the function \code{priorityelasticnet}.
#'
#' @return object of class \code{cvm_priorityelasticnet} with the following elements. If these elements are lists, they contain the results for each penalized block of the best result.
#' \describe{
#' \item{\code{lambda.ind}}{list with indices of lambda for \code{lambda.type}.}
#' \item{\code{lambda.type}}{type of lambda which is used for the predictions.}
#' \item{\code{lambda.min}}{list with values of lambda for \code{lambda.type}.}
#' \item{\code{min.cvm}}{list with the mean cross-validated errors for \code{lambda.type}.}
#' \item{\code{nzero}}{list with numbers of non-zero coefficients for \code{lambda.type}.}
#' \item{\code{glmnet.fit}}{list of fitted \code{glmnet} objects.}
#' \item{\code{name}}{a text string indicating type of measure.}
#' \item{\code{block1unpen}}{if \code{block1.penalization = FALSE}, the results of either the fitted \code{glm} or \code{coxph} object.}
#' \item{\code{best.blocks}}{character vector with the indices of the best block specification.}
#' \item{\code{best.blocks.indices}}{list with the indices of the best block specification ordered by best to worst.}
#' \item{\code{best.max.coef}}{vector with the number of maximal coefficients corresponding to \code{best.blocks}.}
#' \item{\code{best.model}}{complete \code{priorityelasticnet} model of the best solution.}
#' \item{\code{coefficients}}{coefficients according to the results obtained with \code{best.blocks}.}
#' \item{\code{call}}{the function call.}
#' }
#'
#' @note The function description and the first example are based on the R package \code{ipflasso}.
#' @export
#'


cvm_priorityelasticnet <- function(X,
                                   Y,
                                   weights,
                                   family,
                                   type.measure,
                                   blocks.list,
                                   max.coef.list = NULL,
                                   block1.penalization = TRUE,
                                   lambda.type = "lambda.min",
                                   standardize = TRUE,
                                   nfolds = 10,
                                   foldid,
                                   cvoffset = FALSE,
                                   cvoffsetnfolds = 10,
                                   alpha = 1,
                                   ...){
  
  if(!is.null(max.coef.list)){
    if(length(blocks.list) != length(max.coef.list)){stop("blocks.list and max.coef.list must have the same length.")}
    if(any(unlist(lapply(blocks.list, length)) != unlist(lapply(max.coef.list, length)))){stop("blocks.list and the entries of max.coef.list must have the same length.")}
  }
  
  nlist <- length(blocks.list)
  all_res <- list()
  cvmin = +Inf
  
  for(j in 1:nlist){
    all_res[[j]] <- priorityelasticnet(X = X,
                                       Y = Y,
                                       weights = weights,
                                       family = family,
                                       type.measure = type.measure,
                                       blocks = blocks.list[[j]],
                                       max.coef = max.coef.list[[j]],
                                       block1.penalization = block1.penalization,
                                       lambda.type = lambda.type,
                                       standardize = standardize,
                                       nfolds = nfolds,
                                       foldid = foldid,
                                       cvoffset = cvoffset,
                                       cvoffsetnfolds = cvoffsetnfolds,
                                       alpha = alpha,
                                       ...)
    
    if (type.measure != "auc"){
      mincvmj <-  all_res[[j]]$min.cvm[[length(blocks.list[[j]])]]
    }
    if(type.measure == "auc"){
      mincvmj <-  -all_res[[j]]$min.cvm[[length(blocks.list[[j]])]]
    }
    
    if(mincvmj < cvmin){
      ind.best.pr <- j
      cvmin <- mincvmj
    }
  }
  
  best.var <- blocks.list[[ind.best.pr]]
  bfv <- lapply(best.var, '[[', 1)
  blv <- lapply(best.var, length)
  best.blocks <- vector()
  best.blocks.indices <- list()
  for(k in 1:length(bfv)){
    best.blocks <- c(best.blocks, paste("bp",k," = ", bfv[[k]],":",blv[[k]] + bfv[[k]] - 1, sep=""))
    best.blocks.indices[[k]] <- seq(from = bfv[[k]], to = blv[[k]] + bfv[[k]] - 1, by = 1)
  }
  
  finallist <- list(lambda.ind = all_res[[ind.best.pr]]$lambda.ind,
                    lambda.type = lambda.type,
                    lambda.min = all_res[[ind.best.pr]]$lambda.min,
                    min.cvm = all_res[[ind.best.pr]]$min.cvm,
                    nzero = all_res[[ind.best.pr]]$nzero,
                    glmnet.fit = all_res[[ind.best.pr]]$glmnet.fit,
                    name = all_res[[ind.best.pr]]$name,
                    block1unpen = all_res[[ind.best.pr]]$block1unpen,
                    best.blocks = best.blocks,
                    best.blocks.indices = best.blocks.indices,
                    best.max.coef = max.coef.list[[ind.best.pr]],
                    best.model = all_res[[ind.best.pr]],
                    coefficients = all_res[[ind.best.pr]]$coefficients,
                    call = match.call())
  
  class(finallist) <- c("cvm_priorityelasticnet", class(finallist))
  
  return(finallist)
}
