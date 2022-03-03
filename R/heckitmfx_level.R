#' Marginal Effects for Two-Step Sample Selection Models (Heckit) with Outcome Variable in Level
#'
#' @param data The dataset
#' @param selection The selection equation
#' @param outcome The outcome equation
#'
#' @return A dataframe with estimated parameters, partial and total marginal effects
#' @export
#'
#' @import stats
#'
#' @description Partial and total marginal effects for two-step sample selection models (Heckit) with outcome variable in level. Partial marginal effects refer to marginal effects associated with changes in the probability of selection and with changes in the conditional outcome. Total marginal effects refer to marginal effects on the unconditional outcome.
#'
#' @examples
#' # Example following Hoffmann & Kassouf (2005) specification (eqs. 1 & 2)
#' \dontrun{heckitmfx_log(data, L ~ z1 + z2, log(g) ~ x1 + x2)}
#'
#' @examples
#' # Example using simulated data of household tourism expenditures (in level, i.e. monetary units)
#' data(tourexp)
#' selection <- "participation ~ income + education + health"
#' outcome <- "expenditure ~ income + education + tripweather"
#' heckitmfx_level(tourexp, selection, outcome)
#'
#' @references Hoffmann, R. & Kassouf, A. L. (2005) Deriving conditional and unconditional marginal effects in log earnings equations estimated by Heckman's procedure, Applied Economics, 37:11, 1303-1311, DOI: \href{https://www.doi.org/10.1080/00036840500118614}{https://www.doi.org/10.1080/00036840500118614}. (Equation A9).
#'
heckitmfx_level <- function(data, selection, outcome){
  probit <- suppressWarnings(glm(selection, family=binomial(link=probit), data=data, na.action=na.exclude))
  data$alpha <- predict(probit)
  data$densalpha <- dnorm(data$alpha)
  data$probalpha <- pnorm(data$alpha)
  data$lambda <- data$densalpha / data$probalpha
  data$delta <- data$lambda * (data$lambda + data$alpha)
  data$probalpha.delta <- (1 - data$probalpha) * data$delta
  gamma <- probit$coefficients
  estimates <- data.frame(gamma)

  regression <- lm(paste0(outcome, " + lambda"), data=data, na.action = na.exclude)
  data$xb <- predict(regression)
  beta <- regression$coefficients
  betalambda <- regression[["coefficients"]][["lambda"]]
  estimates <- merge(estimates, data.frame(beta), by = "row.names", all = TRUE)

  data$xb.densalpha <- data$xb * data$densalpha
  dpseldx.i <- as.matrix(data$xb.densalpha) %*% t(as.matrix(estimates[,2]))
  dpseldx <- apply(dpseldx.i, 2, FUN = mean, na.rm = TRUE)
  estimates <- cbind(estimates, dpseldx)

  ei1 <- as.matrix(1-data$probalpha) %*% t(as.matrix(estimates[,3]))
  ei2 <- as.matrix(betalambda * data$probalpha.delta) %*% t(as.matrix(estimates[,2]))
  ei1[is.na(ei1)] <- 0
  ei2[is.na(ei2)] <- 0
  ei <- ei1 - ei2
  dyconddx <- apply(ei, 2, FUN = mean, na.rm = TRUE)
  estimates <- cbind(estimates, dyconddx)

  estimates[is.na(estimates)] <- 0
  estimates$dydx <- estimates$dyconddx + estimates$dpseldx
  estimates <- subset(estimates, (Row.names!="lambda" & Row.names!="(Intercept)"))
  names(estimates)[names(estimates) == "Row.names"] <- "variable"
  legend <- c("gamma: Selection coefficients", "beta: Outcome coefficients", "dpseldx: Partial marginal effects on selection", "dyconddx: Partial marginal effects on outcome", "dydx: Total marginal effects")
  return <- list(estimates, legend)

  return(return)
}
