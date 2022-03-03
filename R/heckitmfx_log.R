#' Marginal Effects for Heckit (Log-Tranformed Variable)
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
#' @description Partial and total marginal effects for two-step sample selection models (Heckit) with log-transformed outcome variable. Partial marginal effects refer to marginal effects associated with changes in the probability of selection and with changes in the conditional outcome. Total marginal effects refer to marginal effects on the unconditional outcome.
#'
#' @examples
#' # Example following Hoffmann & Kassouf (2005) specification (eqs. 1 & 2)
#' \dontrun{heckitmfx_log(data, L ~ z1 + z2, log(g) ~ x1 + x2)}
#'
#' @examples
#' # Example using simulated data of household tourism expenditures (in logarithms)
#' data(tourexp)
#' selection <- "participation ~ income + education + health"
#' outcome <- "log(expenditure) ~ income + education + tripweather"
#' heckitmfx_log(tourexp, selection, outcome)
#'
#' @references Hoffmann, R. & Kassouf, A. L. (2005) Deriving conditional and unconditional marginal effects in log earnings equations estimated by Heckman's procedure, Applied Economics, 37:11, 1303-1311, DOI: \href{https://www.doi.org/10.1080/00036840500118614}{https://www.doi.org/10.1080/00036840500118614}. (Equation 14).
#'
heckitmfx_log <- function(data, selection, outcome){
  probit <- suppressWarnings(glm(selection, family=binomial(link=probit), data=data, na.action=na.exclude))
  data$alpha <- predict(probit)
  data$densalpha <- dnorm(data$alpha)
  data$probalpha <- pnorm(data$alpha)
  data$lambda <- data$densalpha / data$probalpha
  data$delta <- data$lambda * (data$lambda + data$alpha)
  gamma <- probit$coefficients
  estimates <- data.frame(gamma)

  regression <- lm(paste0(outcome, " + lambda"), data=data, na.action = na.exclude)
  beta <- regression$coefficients
  estimates <- merge(estimates, data.frame(beta), by = "row.names", all = TRUE)

  dpseldx <- gamma * mean(data$lambda)
  betalambda <- regression[["coefficients"]][["lambda"]]

  ei2 <- gamma * betalambda * mean(data$delta)
  estimates <- merge(estimates, data.frame(ei2, Row.names = names(ei2)), by = "Row.names", all = TRUE)
  estimates <- merge(estimates, data.frame(dpseldx, Row.names=names(dpseldx)), by = "Row.names", all = TRUE)

  estimates[is.na(estimates)] <- 0
  estimates$dyconddx <- estimates$beta - estimates$ei2
  estimates$dydx <- estimates$dyconddx + estimates$dpseldx
  estimates <- subset(estimates, (Row.names!="lambda" & Row.names!="(Intercept)"))
  estimates <- subset(estimates, select = -c(ei2))
  names(estimates)[names(estimates) == "Row.names"] <- "variable"
  legend <- c("gamma: Selection coefficients", "beta: Outcome coefficients", "dpseldx: Partial marginal effects on selection", "dyconddx: Partial marginal effects on outcome", "dydx: Total marginal effects")
  return <- list(estimates, legend)

  return(return)
}

utils::globalVariables(c("Row.names"))
