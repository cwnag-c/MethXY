#' VarCpG
#' @importFrom s20x levene.test
#' @importFrom stats bartlett.test fligner.test
#' @param row row
#' @param Sample_Group Sample_Group
#' @param Probe_name Probe_name
#'
#' @return res
#' @export
#'
VarCpG <- function(row, Sample_Group, Probe_name) {
  lev <- levene.test(row ~ Sample_Group, show.table = F, digit = 200)
  bar <- bartlett.test(row ~ Sample_Group)
  fli <- fligner.test(row ~ Sample_Group)
  res <- data.frame(
    CpG = Probe_name,
    Lev_F = lev$f.value,
    Lev_P = lev$p.value,
    Bar_K2 = unname(bar$statistic),
    Bar_P = bar$p.value,
    FK_Chi2 = unname(fli$statistic),
    FK_P = fli$p.value
  )
  return(res)
}

#' methxy.vmp
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats residuals var
#' @importFrom dplyr bind_rows inner_join
#' @param matrix matrix
#' @param design design
#' @param Sample_Group Sample_Group
#' @param cores cores
#' @param chunks chunks
#'
#' @return res
#' @export
#'
methxy.vmp <- function(matrix, design, Sample_Group, cores = 10, chunks = 10) {
  tmp <- tempdir()
  if (!dir.exists(tmp)) dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  fit <- lmFit(matrix, design)
  Residuals <- residuals(fit, matrix)
  ## chunk
  rows_per_chunks <- ceiling(nrow(Residuals) / chunks)
  dat_list <- lapply(1:chunks, function(i) {
    start <- (i - 1) * rows_per_chunks + 1
    end <- min(i * rows_per_chunks, nrow(Residuals))
    if (start <= nrow(Residuals)) {
      Residuals[start:end, , drop = FALSE]
    } else {
      NULL
    }
  })
  ## case/ctrl means and variance
  resid_case <- Residuals[, which(Sample_Group == "1")]
  resid_control <- Residuals[, which(Sample_Group == "0")]
  rm(Residuals)
  Mean_cases <- apply(resid_case, 1, function(x) {
    mean(x, na.rm = T)
  })
  Variance_cases <- apply(resid_case, 1, function(x) {
    var(x, na.rm = T)
  })
  Mean_controls <- apply(resid_control, 1, function(x) {
    mean(x, na.rm = T)
  })
  Variance_controls <- apply(resid_control, 1, function(x) {
    var(x, na.rm = T)
  })
  rm(resid_case, resid_control)
  df <- data.frame(
    Mean_cases = Mean_cases,
    Mean_controls = Mean_controls,
    Variance_cases = Variance_cases,
    Variance_controls = Variance_controls
  )
  df$Mean_difference <- df$Mean_cases - df$Mean_controls
  df$Variance_difference <- df$Variance_case - df$Variance_controls
  df$CpG <- rownames(df)
  ## run
  res <- lapply(dat_list, function(l) {
    res_tmp <- foreach(
      cg = rownames(l),
      .combine = rbind,
      .packages = c("s20x", "MethXY"),
      .export = c("Sample_Group")
    ) %dopar% {
      VarCpG(l[cg, ], Sample_Group, cg)
    }
    return(res_tmp)
  })
  stopCluster(cl)
  res <- bind_rows(res)
  ## combind
  res <- inner_join(res, df, by = "CpG")
  return(res)
}
