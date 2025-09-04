#' preprocess.estimate
#' @importFrom EpiSmokEr epismoker
#' @importFrom CETYGO projectCellTypeWithError
#' @importFrom minfi sampleNames getAnnotation
#'     featureNames preprocessQuantile getBeta
#' @importFrom utils write.csv read.csv
#' @param DNAm_env returned by preprocess.filter
#' @param save_dir Generates and saves phenotype data using UniqueID
#'     as the filename.
#' @param cell.proportion "Whether to perform cell composition estimation.
#'     The default value is 'F'
#' @param CellType "Blood" or "Prefrontal_cortex".
#'     Specify which tissue type to use for estimation.
#'     The default value is 'Blood'.
#' @param Referenceset To predict cell composition in blood tissue,
#'     the path to the reference file must be provided.
#' @param predict.smok The default value is 'F'.Whether to predict smoking score
#' @param predict.age The default value is 'F'.Whether to predict biological age
#'
#' @return The predicted information has been
#'     added to the target object in the environment.
#' @export
#'
preprocess.estimate <- function(DNAm_env,
                                save_dir = "./",
                                cell.proportion = FALSE, CellType = "Blood",
                                Referenceset = NULL,
                                predict.smok = FALSE,
                                predict.age = FALSE) {
  UniqueID <- DNAm_env$UniqueID
  targets_file <- paste0(save_dir, "/", UniqueID, "_predicted.csv")
  if (file.exists(targets_file)) {
    message(save_dir, " already contains the file ", "\n",
            paste0(UniqueID, "_predicted.csv"), "\n", "Reading it directly!")
    targets <- as.data.frame(read.csv(targets_file, row.names = 1))
    targets <- targets[match(sampleNames(DNAm_env$GenomicMethylSet),
                             rownames(targets)), ]
    colData(DNAm_env$GenomicMethylSet) <- DataFrame(targets)
  } else {
    ##### need rgSet/beta
    if (cell.proportion) {
      if (CellType == "Blood") {
        referenceset <- load(Referenceset, envir = environment())
        referenceset <- get(referenceset)
        DNAm_env$GenomicMethylSet$Sex <- DNAm_env$GenomicMethylSet$Gender
        cellp <- estimateCellCounts2_GMset(DNAm_env$GenomicMethylSet,
                                           referenceset = referenceset)
        cellp <- cellp$prop
        cellp <- as.data.frame(cellp)
        colData(DNAm_env$GenomicMethylSet) <- cbind(
          colData(DNAm_env$GenomicMethylSet),
          cellp[rownames(colData(DNAm_env$GenomicMethylSet)), , drop = FALSE]
        )
        GenomicMethylSet <- fixMethOutliers(DNAm_env$GenomicMethylSet)
        GRSet <- preprocessQuantile(GenomicMethylSet,
                                    sex = colData(GenomicMethylSet)$Gender,
                                    fixOutliers = FALSE)
        rm(GenomicMethylSet)
        beta <- getBeta(GRSet)
        rm(GRSet)
      }
      if (CellType == "Prefrontal_cortex") {
        GenomicMethylSet <- fixMethOutliers(DNAm_env$GenomicMethylSet)
        GRSet <- preprocessQuantile(GenomicMethylSet,
                                    sex = colData(GenomicMethylSet)$Gender,
                                    fixOutliers = FALSE)
        rm(GenomicMethylSet)
        beta <- getBeta(GRSet)
        rm(GRSet)
        cellp <- projectCellTypeWithError(beta, modelBrainCoef[["IDOL"]][[7]])
        cellp <- cellp[, 1:4]
        cellp <- as.data.frame(cellp)
        colData(DNAm_env$GenomicMethylSet) <- cbind(
          colData(DNAm_env$GenomicMethylSet),
          cellp[rownames(colData(DNAm_env$GenomicMethylSet)), , drop = FALSE]
        )
      }
      if (CellType != "Blood" && CellType != "Prefrontal_cortex") {
        stop("Error: CellType must be either 'Blood' or 'Prefrontal_cortex'.")
      }
    } else {
      GenomicMethylSet <- fixMethOutliers(DNAm_env$GenomicMethylSet)
      GRSet <- preprocessQuantile(GenomicMethylSet,
                                  sex = colData(GenomicMethylSet)$Gender,
                                  fixOutliers = FALSE)
      rm(GenomicMethylSet)
      beta <- getBeta(GRSet)
      rm(GRSet)
    }
    ### need beta
    if (predict.smok) {
      smoking <- epismoker(beta,
        ref.Elliott = Illig_data,
        method = "SSc"
      )
      colData(DNAm_env$GenomicMethylSet) <- cbind(
        colData(DNAm_env$GenomicMethylSet),
        smoking[rownames(colData(DNAm_env$GenomicMethylSet)), , drop = FALSE]
      )
    }
    ### need beta
    if (predict.age) {
      age <- agep2(beta, method = c("horvath"))
      colnames(age) <- c("predictedAge", "horvath.n_missing")
      colData(DNAm_env$GenomicMethylSet) <- cbind(
        colData(DNAm_env$GenomicMethylSet),
        age[rownames(colData(DNAm_env$GenomicMethylSet)), , drop = FALSE]
      )
    }
    targets <- as.data.frame(colData(DNAm_env$GenomicMethylSet))
    write.csv(targets,
              file = paste0(save_dir, "/", UniqueID, "_predicted.csv"),
              row.names = TRUE)
  }
}
