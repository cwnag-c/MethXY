#' methxy.dmp
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats p.adjust pt
#' @param matrix beta or M matrix.
#' @param Design.matrix Some users need to adjust for covariates,
#'     so here only the design matrix is input.
#' @param use_eBayes  default F; use eBayes
#' @param adjust.method "bonferroni" or "BH". default "bonferroni"
#'
#' @return Data.frame
#' @export
#'
methxy.dmp <- function(matrix, Design.matrix,
                       use_eBayes = F,
                       adjust.method = "bonferroni") {
  if (use_eBayes) {
    fit <- lmFit(matrix, Design.matrix)
    fitE <- eBayes(fit)
    if (adjust.method == "BH") {
      DMPs <- topTable(fit = fitE, adjust.method = "BH", number = nrow(matrix), coef = 2)
    } else {
      if (adjust.method == "bonferroni") {
        DMPs <- topTable(fit = fitE, adjust.method = "bonferroni", number = nrow(matrix), coef = 2)
      }
    }
    DMPs$CpG <- rownames(DMPs)
    DMPs <- DMPs[, c("CpG", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(DMPs) <- c("CpG", "coefficient", "t", "P.Value", "adj.P.Val")
  } else {
    fit <- lmFit(matrix, Design.matrix)
    t <- fit$coefficients[, 2] / (fit$stdev.unscaled[, 2] * fit$sigma)
    P.Value <- 2 * pt(-abs(t), df = fit$df.residual[1])
    Coefficient <- fit$coefficients[, 2]
    DMPs <- data.frame(
      CpG = rownames(matrix),
      coefficient = Coefficient,
      t = t,
      P.Value = P.Value
    )
    if (adjust.method == "BH") {
      DMPs$adj.P.Val <- p.adjust(DMPs$P.Value, method = "BH")
    } else {
      if (adjust.method == "bonferroni") {
        DMPs$adj.P.Val <- p.adjust(DMPs$P.Value, method = "bonferroni")
      }
    }
  }
  return(DMPs)
}
