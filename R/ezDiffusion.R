#' EZ Diffusion
#'
#' Calculate the parameters of the EZ diffusion model from mrt, vrt, pc...
#'
#' @param proportion_correct Character. The name of the variable that contains proportion of correct trials.
#' @param rt_variance Character. The name of the variable that contains the variance of reaction times for correct decisions
#' @param rt_mean Character. The name of the variable that contains the mean of reaction times for correct decisions.
#' @param data A \code{data.frame} that conitains the above-specified variables.
#' @param s Numeric. A scaling constant. Take a look at the literature, because different
#' authors use different scaling constants.
#'
#' @examples
#' data <- data.frame(
#'   prop_correct = runif(n = 100, min = 0.2, max = .8)
#'   , rt_mean = runif(n = 100, min = .2, max = .8)
#'   , rt_var = runif(n = 100, min = .01, max = .1)
#' )
#'
#' ezDiffusion(
#'   data = data
#'   , proportion_correct = "prop_correct"
#'   , rt_variance = "rt_var"
#'   , rt_mean = "rt_mean"
#'  )
#'
#' @export


ezDiffusion <- function(proportion_correct, rt_variance, rt_mean, s = 0.1, data){

  pc <- data[[proportion_correct]]
  vrt <- data[[rt_variance]]
  mrt <- data[[rt_mean]]

  parameters <- data.frame(
    "v" = rep(NA, nrow(data))
    , "a" = rep(NA, nrow(data))
    , "t0" = rep(NA, nrow(data))
  )
  # apply edge corrections if necessary
  if(any(pc==0)){
    pc[pc==0] <- 1e-9
  }

  if(any(pc==.5)){
    pc[pc==.5] <- .5 + 1e-9 * sample(rep(-1, 1), size = sum(pc==.5), replace = TRUE)
  }

  if(any(pc==1)){
    pc[pc==1] <- 1 - 1e-9
  }

  if(any(vrt==0)){
    warning("variances==0 supplied")
  }

  # The function qlogis calculates the logit.
  L <- qlogis(pc)

  # These give drift rate
  x <- L * (L * pc^2 - L * pc + pc -.5)/vrt
  v <- sign(pc - .5) * s * x^(1/4)

  # This gives boundary separation
  a <- s^2 * L / v

  # This gives non-decision time
  y <- -v*a/s^2
  mdt <- (a/(2*v)) * (1-exp(y))/(1+exp(y))
  t0 <- mrt - mdt

  parameters$v <- v
  parameters$a <- a
  parameters$t0 <- t0

  return(parameters)
}

#' The EZ-diffusion model
#'
#' This is a convenience wrapper function for the EZ-diffusion model.
#'
#' @param data A \code{data.frame} that contains the data.
#' @param id Character. The name of the variable that contains the participant identifier.
#' @param rt Character. The name of the variable that contains the reaction times.
#' @param correct Character. The name of the variable that codes correct vs. erroneous decision. Needs to be 0 = erroneous and 1 = correct.
#' @param within Character. Possibly a vector of names of variables that contain within-subjects factors.
#' @param s Numeric. An optional scaling constant.
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:150, each = 100)
#'   , rt = runif(150 * 100, min = .2, max = .8)
#'   , correct = sample(0:1, size = 150 * 100, replace = TRUE)
#' )
#' ez_diffusion_model(data = data, id = "id", rt = "rt", correct = "correct")

#' @export

ez_diffusion_model <- function(data, id, rt, correct, within = NULL, s = .1){
  correct_trials <- data[data[[correct]]==1, ]
  mrt <- papaja:::fast_aggregate(data = correct_trials, dv = rt, factors = c(id, within), fun = mean)
  vrt <- papaja:::fast_aggregate(data = correct_trials, dv = rt, factors = c(id, within), fun = var)
  pc <- papaja:::fast_aggregate(data = data, dv = correct, factors = c(id, within), fun = mean)

  colnames(mrt)[which(colnames(mrt)==rt)] <- "rt_mean"
  colnames(vrt)[which(colnames(vrt)==rt)] <- "rt_variance"
  colnames(pc)[which(colnames(pc)==correct)] <-"proportion_correct"

  dat_ <- merge(mrt, vrt)
  dat_ <- merge(dat_, pc)

  output <- ezDiffusion(
    data = dat_
    , proportion_correct = "proportion_correct"
    , rt_mean = "rt_mean"
    , rt_variance = "rt_variance"
    , s = s
  )
  cbind(dat_, output)
}

# ez2_moments <- function(data, rt, error, factors){
#
#   means <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors, error), fun = mean)
#   variances <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors, error), fun = var)
#   proportion_correct <- papaja:::fast_aggregate(data = data, dv = error, factors = c(factors), fun = mean)
#
#   data <- data.frame(
#     "pc" = proportion_correct[[error]]
#     , "mrt1" = means[[rt]][means[[error]]==1]
#     , "mrt0" = means[[rt]][means[[error]]==0]
#     , "vrt1" = variances[[rt]][variances[[error]]==1]
#     , "vrt0" = variances[[rt]][variances[[error]]==0]
#   )
#
#   data
# }

# ez2_moments <- function(data, rt, error, factors, within = NULL){
#
#   colnames(data) <- gsub(colnames(data), pattern = " ", replacement = "__")
#   factors <- gsub(factors, pattern = " ", replacement = "__")
#   within <- gsub(within, pattern = " ", replacement = "__")
#   for(i in factors){
#     data[[i]] <- droplevels(as.factor(data[[i]]))
#   }
#   # within <- gsub(within, pattern = " ", replacement = "__")
#
#   means <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors), fun = mean)
#   colnames(means) <- gsub(colnames(means), pattern = rt, replacement = "mrt")
#   variances <- papaja:::fast_aggregate(data = data, dv = rt, factors = c(factors), fun = var)
#   colnames(variances) <- gsub(colnames(variances), pattern = rt, replacement = "vrt")
#   proportion_correct <- papaja:::fast_aggregate(data = data, dv = error, factors = c(factors), fun = mean)
#   colnames(proportion_correct) <- gsub(colnames(proportion_correct), pattern = error, replacement = "pc")
#
#   data <- merge(means, variances)
#   data <- merge(data, proportion_correct)
#
#   data <- data[, c(factors, "pc", "mrt", "vrt")]
#
#
#   if(!is.null(within)){
#     tmp1 <- data[data[[within]]==levels(data[[within]])[1], ]
#     tmp2 <- data[data[[within]]==levels(data[[within]])[2], ]
#     tmp1 <- tmp1[, c(setdiff(factors, within), "mrt", "vrt", "pc")]
#     tmp2 <- tmp2[, c(setdiff(factors, within), "mrt", "vrt", "pc")]
#
#     for(i in c("mrt", "vrt", "pc")){
#       colnames(tmp1)[which(colnames(tmp1)==i)] <- paste0(i, "_", levels(data[[within]])[1])
#       colnames(tmp2)[which(colnames(tmp2)==i)] <- paste0(i, "_", levels(data[[within]])[2])
#     }
#     data <- merge(tmp1, tmp2)
#   }
#
#   colnames(data) <- gsub(colnames(data), pattern = "__", replacement = " ")
#
#   return(data)
# }


