## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  class.output="scroll-100",
  cache.path = "cached/"
)
# Check if 'priorityelasticnet' is available
if (!requireNamespace("priorityelasticnet", quietly = TRUE)) {
  message("The 'priorityelasticnet' package is not installed. Please install it to fully reproduce this vignette.")
} else {
  library(priorityelasticnet)
}

## ----eval=FALSE---------------------------------------------------------------
# install.packages("priorityelasticnet")

## -----------------------------------------------------------------------------
# Simulate some data
set.seed(123)
n <- 100  # Number of observations
p <- 50   # Number of predictors

## -----------------------------------------------------------------------------
# Create a matrix of predictors
X <- matrix(rnorm(n * p), n, p)

## -----------------------------------------------------------------------------
# Generate a response vector based on a linear combination of some predictors
beta <- rnorm(10)  # Coefficients for the first 10 predictors
Y <- X[, 1:10] %*% beta + rnorm(n)  # Linear model with added noise

## -----------------------------------------------------------------------------
# Define predictor blocks
blocks <- list(
  block1 = 1:10,    # First block includes the first 10 predictors
  block2 = 11:30,   # Second block includes the next 20 predictors
  block3 = 31:50    # Third block includes the last 20 predictors
)

## -----------------------------------------------------------------------------
# Fit a priorityelasticnet model
fit <- priorityelasticnet(
  X = X, 
  Y = Y, 
  family = "gaussian", 
  blocks = blocks, 
  type.measure = "mse",
  alpha = 0.5
)

## -----------------------------------------------------------------------------
fit$lambda.ind

## -----------------------------------------------------------------------------
fit$lambda.type

## -----------------------------------------------------------------------------
fit$lambda.min

## -----------------------------------------------------------------------------
fit$min.cvm

## -----------------------------------------------------------------------------
fit$nzero

## -----------------------------------------------------------------------------
fit$glmnet.fit

## -----------------------------------------------------------------------------
fit$coefficients

## -----------------------------------------------------------------------------
head(cbind.data.frame(pred = fit$pred[,1], observed = fit$actuals))

## -----------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Number of observations and predictors
n <- 50  # Number of observations
p <- 300  # Number of predictors

# Number of non-zero coefficients
nzc <- trunc(p / 10)

# Simulate predictor matrix
x <- matrix(rnorm(n * p), n, p)

# Simulate regression coefficients for non-zero predictors
beta <- rnorm(nzc)

# Calculate linear predictor
fx <- x[, seq(nzc)] %*% beta / 3

# Calculate hazard function
hx <- exp(fx)

# Simulate survival times using exponential distribution
ty <- rexp(n, hx)

# Generate censoring indicator (30% censoring probability)
tcens <- rbinom(n = n, prob = .3, size = 1)

# Load survival library and create survival object
library(survival)
y <- Surv(ty, 1 - tcens)

## -----------------------------------------------------------------------------
blocks <- list(
  bp1 = 1:20,    # First block with predictors 1 to 20
  bp2 = 21:200,  # Second block with predictors 21 to 200
  bp3 = 201:300  # Third block with predictors 201 to 300
)

## -----------------------------------------------------------------------------
# Fit Cox model using priorityelasticnet
fit_cox <- priorityelasticnet(
  x, 
  y, 
  family = "cox", 
  alpha = 0.5, 
  type.measure = "deviance", 
  blocks = blocks,
  block1.penalization = TRUE,
  lambda.type = "lambda.min",
  standardize = TRUE,
  nfolds = 10,
  cvoffset = TRUE
  
)

## -----------------------------------------------------------------------------
fit_cox$min.cvm

## -----------------------------------------------------------------------------
fit_cox$coefficients

## -----------------------------------------------------------------------------
fit_cox$lambda.min

## ----eval = requireNamespace("glmSparseNet", quietly = TRUE)------------------

library(glmSparseNet)

# Extract coefficients from the fitted Cox model
chosen.btas <- fit_cox$coefficients
y <- data.frame(
  time = ty,          # Survival times
  status = 1 - tcens  # Event indicator
)

# Group patients and plot Kaplan-Meier survival curves
separate2GroupsCox(
  chosen.btas = chosen.btas,  # Coefficients from the model
  xdata = x,                  # Predictor matrix (xdata)
  ydata = y,                  # Survival data (ydata as Surv object)
  probs = c(0.4, 0.6),        # Median split (adjust if necessary)
  no.plot = FALSE,            # Plot the Kaplan-Meier curve
  plot.title = "Survival Curves",  # Plot title
  xlim = NULL,                # Automatic x-axis limits
  ylim = NULL,                # Automatic y-axis limits
  expand.yzero = FALSE,       # Don't force y-axis to start at zero
  legend.outside = FALSE      # Keep legend inside the plot
)


## -----------------------------------------------------------------------------

# Check if 'priorityelasticnet' is available
if (!requireNamespace("priorityelasticnet", quietly = TRUE)) {
  message("The 'priorityelasticnet' package is not installed. Please install it to fully reproduce this vignette.")
} else {
  library(priorityelasticnet)
  # Load the dataset only if the package is available
  data("Pen_Data", package = "priorityelasticnet")
  
}



## -----------------------------------------------------------------------------
dim(Pen_Data)



## -----------------------------------------------------------------------------
blocks <- list(
  block1 = 1:5,     # Block 1: First 5 predictors
  block2 = 6:179,   # Block 2: Next 174 predictors
  block3 = 180:324  # Block 3: Next 145 predictors
  
)

## -----------------------------------------------------------------------------
set.seed(123)

fit_bin <- priorityelasticnet(
  X = as.matrix(Pen_Data[, 1:324]), 
  Y = Pen_Data[, 325],
  family = "binomial", 
  alpha = 0.5, 
  type.measure = "auc",
  blocks = blocks,
  standardize = FALSE
)

## -----------------------------------------------------------------------------
predictions <- predict(fit_bin, type = "response")
head(predictions)


## -----------------------------------------------------------------------------
predictions <- predict(fit_bin, newdata = as.matrix(Pen_Data[, 1:324]), type = "response")
head(predictions)

## -----------------------------------------------------------------------------
library(pROC)
roc_curve <- roc(Pen_Data[, 325], predictions[,1])
plot(roc_curve, col = "red", main = "ROC Curve for Binomial Model")
text(0.1, 0.1, labels = paste("AUC =", round(roc_curve$auc, 2)), col = "black", cex = 1.2)



## -----------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Number of observations and predictors
n <- 100  # Number of observations
p <- 50   # Number of predictors
k <- 3    # Number of classes

# Simulate a matrix of predictors
x <- matrix(rnorm(n * p), n, p)

# Simulate a response vector with three classes
y <- factor(sample(1:k, n, replace = TRUE))

## -----------------------------------------------------------------------------
blocks <- list(
  block1 = 1:10,   # First block with predictors 1 to 10
  block2 = 11:30,  # Second block with predictors 11 to 30
  block3 = 31:50   # Third block with predictors 31 to 50
)

## -----------------------------------------------------------------------------
fit_multinom <- priorityelasticnet(
  X = x, 
  Y = y, 
  family = "multinomial", 
  alpha = 0.5, 
  type.measure = "class", 
  blocks = blocks,
  block1.penalization = TRUE,
  lambda.type = "lambda.min",
  standardize = TRUE,
  nfolds = 5
)

## -----------------------------------------------------------------------------
fit_multinom$min.cvm

## -----------------------------------------------------------------------------
fit_multinom$coefficients

## -----------------------------------------------------------------------------
fit_multinom$lambda.min

## -----------------------------------------------------------------------------
fit_no_penalty <-
  priorityelasticnet(
    X,
    Y,
    family = "gaussian",
    type.measure = "mse",
    blocks = blocks,
    block1.penalization = FALSE
  )

## -----------------------------------------------------------------------------
fit_no_penalty

## -----------------------------------------------------------------------------
mcontrol <-missing.control(handle.missingdata = "impute.offset", nfolds.imputation = 5)

fit_missing <- priorityelasticnet(
  X,
  Y,
  family = "gaussian",
  type.measure = "mse",
  blocks = blocks,
  mcontrol = mcontrol
)

## -----------------------------------------------------------------------------
fit_missing

## -----------------------------------------------------------------------------
blocks1 <- list(1:10, 11:30, 31:50)
blocks2 <- list(1:5, 6:20, 21:50)

fit_cvm <-
  cvm_priorityelasticnet(
    X,
    Y,
    blocks.list = list(blocks1, blocks2),
    family = "gaussian",
    type.measure = "mse",
    weights = NULL,
    foldid = NULL
  )

## -----------------------------------------------------------------------------
fit_cvm

## ----eval = FALSE-------------------------------------------------------------
# weightedThreshold(object = fit_bin)

## -----------------------------------------------------------------------------
coef(fit_bin)

## -----------------------------------------------------------------------------
set.seed(123)
X_new <- matrix(rnorm(406 * 324), 406, 324)

predictions < predict(fit_bin, newdata = X_new, type = "response")
head(predictions)

## -----------------------------------------------------------------------------
# Set the random seed for reproducibility
set.seed(1234)

# Simulate high-dimensional data
n <- 200  # Number of observations
p <- 100  # Number of predictors
n_strong <- 10  # Number of strong predictors
n_weak <- 20  # Number of weak predictors

# Design matrix (predictors)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Generate coefficients: strong predictors with large effects, weak with small effects
beta <- c(rep(2, n_strong), rep(0.5, n_weak), rep(0, p - n_strong - n_weak))

# Generate response with Gaussian noise
Y <- X %*% beta + rnorm(n)

## -----------------------------------------------------------------------------
# Define blocks of predictors for the model
blocks <- list(
  strong_block = 1:n_strong,               # Strong predictors
  weak_block = (n_strong + 1):(n_strong + n_weak),  # Weak predictors
  noise_block = (n_strong + n_weak + 1):p  # Noise (irrelevant predictors)
)

## -----------------------------------------------------------------------------
# Run priorityelasticnet with Adaptive Elastic Net
result <- priorityelasticnet(X = X, 
                             Y = Y, 
                             family = "gaussian", 
                             alpha = 0.5, 
                             type.measure = "mse", 
                             blocks = blocks, 
                             adaptive = TRUE,
                             initial_global_weight = FALSE, 
                             verbose = TRUE)

## -----------------------------------------------------------------------------
# Examine the coefficients
cat("Final model coefficients:")
result$coefficients

## -----------------------------------------------------------------------------
# Examine the adaptive weights
cat("Adaptive weights for each predictor:")
result$adaptive_weights

## -----------------------------------------------------------------------------
plot(result$glmnet.fit[[1]], xvar = "lambda", label = TRUE, main = "Coefficient Paths for Strong Block")

## -----------------------------------------------------------------------------
plot(result$glmnet.fit[[2]], xvar = "lambda", label = TRUE, main = "Coefficient Paths for Weak Block")

## -----------------------------------------------------------------------------
plot(result$glmnet.fit[[3]], xvar = "lambda", label = TRUE, main = "Coefficient Paths for Noise Block")

## -----------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Number of observations and predictors
n <- 50  # Number of observations
p <- 300  # Number of predictors

# Number of non-zero coefficients
nzc <- trunc(p / 10)

# Simulate predictor matrix
x <- matrix(rnorm(n * p), n, p)

# Simulate regression coefficients for non-zero predictors
beta <- rnorm(nzc)

# Calculate linear predictor
fx <- x[, seq(nzc)] %*% beta / 3

# Calculate hazard function
hx <- exp(fx)

# Simulate survival times using exponential distribution
ty <- rexp(n, hx)

# Generate censoring indicator (30% censoring probability)
tcens <- rbinom(n = n, prob = .3, size = 1)

# Load survival library and create survival object
library(survival)
y <- Surv(ty, 1 - tcens)

## -----------------------------------------------------------------------------
blocks <- list(
  bp1 = 1:20,    # First block with predictors 1 to 20
  bp2 = 21:200,  # Second block with predictors 21 to 200
  bp3 = 201:300  # Third block with predictors 201 to 300
)

## -----------------------------------------------------------------------------
# Fit Cox model using priorityelasticnet
result_cox <- priorityelasticnet(
  x, 
  y, 
  family = "cox", 
  alpha = 1, 
  type.measure = "deviance", 
  blocks = blocks,
  block1.penalization = TRUE,
  lambda.type = "lambda.min",
  standardize = TRUE,
  nfolds = 5,
  adaptive = TRUE,
  initial_global_weight = FALSE
)

## -----------------------------------------------------------------------------
# Examine the coefficients
cat("Final model coefficients:")
result_cox$coefficients


## -----------------------------------------------------------------------------
result_cox$initial_coeff

## -----------------------------------------------------------------------------
# Examine the adaptive weights
cat("Adaptive weights for each predictor:")
result_cox$adaptive_weights

## -----------------------------------------------------------------------------
# Run priorityelasticnet with Adaptive Elastic Net
result_bin <- priorityelasticnet(X = as.matrix(Pen_Data[, 1:324]), Y = Pen_Data[, 325],
                             family = "binomial", alpha = 0.5, type.measure = "auc",
                             blocks = list(bp1 = 1:5, bp2 = 6:179, bp3 = 180:324),
                             standardize = FALSE,
                             adaptive = TRUE,
                             initial_global_weight = FALSE, 
                             verbose = TRUE)

## -----------------------------------------------------------------------------
result_bin$nzero

## -----------------------------------------------------------------------------
result_bin$min.cvm

## -----------------------------------------------------------------------------
result_bin$lambda.min

## -----------------------------------------------------------------------------
result_bin$adaptive_weights

## -----------------------------------------------------------------------------
result_bin$coefficients

## -----------------------------------------------------------------------------
predictions <- predict(result_bin, newdata = as.matrix(Pen_Data[, 1:324]), type = "response")
head(predictions)

## -----------------------------------------------------------------------------
library(pROC)
roc_curve <- roc(Pen_Data[, 325], predictions[,1])
plot(roc_curve, col = "red", main = "ROC Curve for Binomial Model")
text(0.1, 0.1, labels = paste("AUC =", round(roc_curve$auc, 2)), col = "black", cex = 1.2)



## -----------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(123)

# Number of observations and predictors
n <- 100  # Number of observations
p <- 50   # Number of predictors
k <- 3    # Number of classes

# Simulate a matrix of predictors
x <- matrix(rnorm(n * p), n, p)

# Simulate a response vector with three classes
y <- factor(sample(1:k, n, replace = TRUE))

## -----------------------------------------------------------------------------
blocks <- list(
  block1 = 1:10,   # First block with predictors 1 to 10
  block2 = 11:30,  # Second block with predictors 11 to 30
  block3 = 31:50   # Third block with predictors 31 to 50
)

## -----------------------------------------------------------------------------

# Run priorityelasticnet
result_multinom <- priorityelasticnet(
  X = x, 
  Y = y, 
  family = "multinomial", 
  alpha = 0.5, 
  type.measure = "class", 
  blocks = blocks,
  block1.penalization = TRUE,
  lambda.type = "lambda.min",
  standardize = TRUE,
  nfolds = 10,
  adaptive = TRUE,
  initial_global_weight = FALSE
  
)


## -----------------------------------------------------------------------------
result_multinom$coefficients

## -----------------------------------------------------------------------------
result_multinom$adaptive

## -----------------------------------------------------------------------------
result_multinom$adaptive_weights

## -----------------------------------------------------------------------------
result_multinom$glmnet.fit

## -----------------------------------------------------------------------------
result_multinom$min.cvm

## -----------------------------------------------------------------------------
result_multinom$lambda.min

