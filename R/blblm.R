#' @import purrr
#' @import stats
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c(".","plan","fit","multicore","capture.output","select")
                       )

#' @inheritParams split_data
#' @inheritParams lm_each_subsample
#' @export
blblm <- function(formula, data, m , B ) {

  ans<-readline(prompt ="Do you want to use parallelization?(Yes/No) Please enter: ")
  dat<-selectList(data)

  if(ans=="Yes"||ans=="Y"){
    res<-par_blblm(formula, dat, n=4, m, B)
    class(res) <- "blblm"
    invisible(res)
  }else{
    data_list <- split_data(dat, m)
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = dat, n = nrow(dat), B = B))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}

#' @inheritParams split_data
#' @inheritParams lm_each_subsample
#' use parallelization to do the little bag of bootstraps.
par_blblm <- function(formula, data, n, m, B ){

  suppressWarnings(future::plan(future::multicore, workers = n))

  data_list <- split_data(data, m)

  estimates <- furrr::future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  return(res)
}

#' @inheritParams split_data
#' ask the user to choose specify columns of the dataset that could be
#' used in the little bag of bootstraps algorithm.
selectList<-function(data){
  pos<-readline(prompt = "Would you like to choose a list of file of datasets? (Yes/No)")
  dat<-data
  col<-NULL
  while(pos=="Yes" || pos=="Y"){
    name<-readline(prompt = "Insert one column name that would be choosed: ")
    col<-c(col,name)
    dat<-data[,col]
    pos<-readline(prompt = "Would you like to choose another list of file of datasets? (Yes/No) ")
  }
  return(dat)
}

#' @param data a data set
#' @param m integer
#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' @param B repeat B times
#' @inheritParams lm_each_boot
#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' @param n split n small data sets
#' @inheritParams lm1
#' compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' @param formula fit linear model
#' @param data a data set
#' @param freqs frequency
#' @param fit fitted linear model
#' estimate the regression estimates based on given the number of repetitions
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' @param fit variable Yhat from the fitted linear model
#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}

#' @param fit
#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

#' GLM model
#'
#' @inheritParams split_data
#' @inheritParams glm_each_subsample
#' @inheritParams glm1
#' @export
GLM_blblm <- function(formula, data, m , B ) {

  ans<-readline(prompt ="Do you want to use parallelization?(Yes/No) Please enter: ")
  dat<-selectList(data)

  if(ans=="Yes"||ans=="Y"){
    res<-par_blblm(formula, dat, n=4, m, B)
    class(res) <- "blblm"
    invisible(res)
  }else{
    data_list <- split_data(dat, m)
    estimates <- map(
      data_list,
      ~ glm_each_subsample(formula = formula, data = dat, n = nrow(dat), B = B))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}

#' @inheritParams split_data
#' @inheritParams glm_each_subsample
#' use parallelization to do the little bag of bootstraps.
par_glmblblm <- function(formula, data, n, m, B ){

  suppressWarnings(future::plan(future::multicore, workers = n))

  data_list <- split_data(data, m)

  estimates <- furrr::future_map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  return(res)
}
#' @param B repeat B times
#' @inheritParams glm_each_boot
#' compute the estimates
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}

#' @param n split n small data sets
#' @inheritParams glm1
#' compute the GLM estimates for a blb dataset
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs)
}

#' @param formula fit linear model
#' @param data a data set
#' @param freqs frequency
#' @param fit fitted linear model
#' estimate the regression estimates based on given the number of repetitions
glm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}