#' anti.trafo
#'
#' @param x x
#' @param adult.age adult.age
#'
#' @return res
anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) *
    x + adult.age)
}


#' .compute_ages
#' @importFrom stats na.omit
#' @param betas matrix of betas
#' @param coeff coeff
#'
#' @return res
#'
.compute_ages <- function(betas, coeff) {
  ages <- apply(betas, 2, function(x) {
    ## .split_intercept_from_coeff
    intercept <- coeff["(Intercept)"]
    if (is.na(intercept)) {
      intercept <- 0
    }
    interceptless_coeff <- coeff[!names(coeff) %in% "(Intercept)"]
    coef2 <- list(intercept = intercept, coeffs = interceptless_coeff)
    ## .handle_missing
    not_missing <- names(coef2$coeffs) %in% names(na.omit(x))
    coef2$missing_probes <- paste0(na.omit(names(coef2$coeffs)[!not_missing]),
      collapse = ";"
    )
    coef2$n_missing <- length(na.omit(names(coef2$coeffs)[!not_missing]))
    coef2$coeffs <- coef2$coeffs[not_missing]
    data <- x[names(coef2$coeffs)]
    the_sum <- data %*% coef2$coeffs + coef2$intercept
    return(data.frame(the_sum, coef2$n_missing, coef2$missing_probes,
      stringsAsFactors = FALSE
    ))
  })
  ages <- do.call("rbind", ages)
  return(ages)
}

#' agep2
#' @importFrom utils data
#' @param betas Matrix of betas
#' @param method Default is "horvath"
#' @param n_missing Logical, additionally output the number of missing CpGs
#'     for each sample using the specified method or coeff list.
#' @param missing_probes Logical, additionally output the names of missing CpGs
#'     for each sample using the specified method or coeff list.
#'
#' @return Returns matrix of predicted ages per sample.
#' @export
#'
agep2 <- function(betas,
                  method = c("horvath", "hannum", "phenoage", "skinblood", "lin", "all"),
                  n_missing = TRUE,
                  missing_probes = FALSE) {
  method <- match.arg(method)
  # data("ageCoefs", package = "DnaMethyXY", envir = environment())
  ages <- switch(method,
    horvath = {
      pre <- .compute_ages(betas = betas, coeff = ageCoefs[["Horvath"]])
      pre[, 1] <- anti.trafo(pre[, 1], adult.age = 20)
      colnames(pre) <- c(
        "horvath.age", "horvath.n_missing",
        "horvath.missing_probes"
      )
      pre
    },
    hannum = {
      pre <- .compute_ages(betas = betas, coeff = ageCoefs[["Hannum"]])
      colnames(pre) <- c(
        "hannum.age", "hannum.n_missing",
        "hannum.missing_probes"
      )
      pre
    },
    lin = {
      pre <- .compute_ages(betas = betas, coeff = ageCoefs[["Lin"]])
      colnames(pre) <- c("lin.age", "lin.n_missing", "lin.missing_probes")
      pre
    },
    skinblood = {
      pre <- .compute_ages(betas = betas, coeff = ageCoefs[["SkinBlood"]])
      pre[, 1] <- anti.trafo(pre[, 1], adult.age = 20)
      colnames(pre) <- c(
        "skinblood.age", "skinblood.n_missing",
        "skinblood.missing_probes"
      )
      pre
    },
    phenoage = {
      pre <- .compute_ages(betas = betas, coeff = ageCoefs[["PhenoAge"]])
      colnames(pre) <- c(
        "phenoage.age", "phenoage.n_missing",
        "phenoage.missing_probes"
      )
      pre
    },
    all = {
      clocks <- c(
        horvath = "horvath", hannum = "hannum",
        phenoage = "phenoage", skinblood = "skinblood",
        lin = "lin"
      )
      out <- lapply(clocks, function(x, betas) {
        agep2(
          betas = betas, method = x,
          n_missing = n_missing, missing_probes = missing_probes
        )
      }, betas = betas)
      out <- do.call("cbind", out)
      return(out)
    }
  )
  return(ages[, c(TRUE, n_missing, missing_probes)])
}
