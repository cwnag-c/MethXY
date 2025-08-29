#' read_from_intensity
#' @importFrom data.table fread fwrite
#' @importFrom minfi readGEORawFile
#' @importFrom utils write.csv
#' @param targets_dir argets_dir
#' @param DNAmType DNAmType
#' @param UniqueID UniqueID
#' @param pattern 1 or 3
#' @param detp file_dir
#' @param methyl file_dir
#' @param unmethyl file_dir
#' @param merge file_dir
#' @param nThread 12
#'
#' @return env
#' @export
#'
read_from_intensity <- function(targets_dir, DNAmType, UniqueID,
                                pattern = NULL, detp = NULL, methyl = NULL, unmethyl = NULL, merge = NULL,
                                nThread = 12) {
  ## read file####
  targets <- read.csv(targets_dir)
  GSM_id_pos <- paste0(targets$Sample_Name, "_", targets$Sentrix_ID, "_", targets$Sentrix_Position)
  id_pos <- paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)
  if (pattern == 3) {
    detp <- as.data.frame(fread(detp, nThread = nThread))
    methyl <- as.data.frame(fread(methyl, nThread = nThread))
    unmethyl <- as.data.frame(fread(unmethyl, nThread = nThread))
    rownames(detp) <- detp$V1
    rownames(methyl) <- methyl$V1
    rownames(unmethyl) <- unmethyl$V1
    detp$V1 <- NULL
    methyl$V1 <- NULL
    unmethyl$V1 <- NULL
  }
  if (pattern == 1) {
    merge <- as.data.frame(fread(merge, nThread = nThread))
    rownames(merge) <- merge$V1
    merge$V1 <- NULL
    #### split
    methyl <- merge[, grepl("Methylated", colnames(merge))]
    unmethyl <- merge[, grepl("Unmethylated", colnames(merge))]
    detp <- merge[, !grepl("_Methylated|_Unmethylated", colnames(merge))]
    rm(merge)
  }
  ## process #####
  if (!all(ncol(detp) == c(ncol(methyl), ncol(unmethyl)))) {
    stop("detp, methyl, and unmethyl have unequal numbers of columns.")
  }
  index_tf_detp <- sapply(id_pos, function(x) any(grepl(x, colnames(detp))))
  index_tf_methyl <- sapply(id_pos, function(x) any(grepl(x, colnames(methyl))))
  index_tf_unmethyl <- sapply(id_pos, function(x) any(grepl(x, colnames(unmethyl))))
  prop_true_detp <- mean(index_tf_detp)
  prop_true_methyl <- mean(index_tf_methyl)
  prop_true_unmethyl <- mean(index_tf_unmethyl)
  if (!all(prop_true_detp == c(prop_true_methyl, prop_true_unmethyl))) {
    stop("Possible ID mismatch: detp, methyl, and unmethyl may come from different sources.")
  }
  ## Ensure that each element in id_pos_t matches at least one item, and that no unmatched cases exist.
  id_pos_t <- id_pos[index_tf_detp]
  GSM_id_pos_t <- GSM_id_pos[index_tf_detp]
  # Extract index positions
  index_detp <- sapply(id_pos_t, function(x) grep(x, colnames(detp)))
  index_methyl <- sapply(id_pos_t, function(x) grep(x, colnames(methyl)))
  index_unmethyl <- sapply(id_pos_t, function(x) grep(x, colnames(unmethyl)))
  # Ensure that each element in id_pos_t matches exactly one item.
  if (any(lengths(index_detp) > 1)) {
    stop("Duplicate IDs detected in the detp that match target IDs. Please ensure each target ID uniquely corresponds to one methylation column name.")
  }
  if (any(lengths(index_methyl) > 1)) {
    stop("Duplicate IDs detected in the methyl that match target IDs. Please ensure each target ID uniquely corresponds to one methylation column name.")
  }
  if (any(lengths(index_unmethyl) > 1)) {
    stop("Duplicate IDs detected in the unmethyl that match target IDs. Please ensure each target ID uniquely corresponds to one methylation column name.")
  }

  if (prop_true_detp == 0) {
    stop("The target failed to match the column names of detp.This may be because the column names are not in the 'SentrixID_SentrixPosition' format.")
  } else {
    if (prop_true_detp != 1) {
      message("Warning: IDs in the target and methylation files do not fully match. Only the intersection of both will be used in subsequent analyses. Please verify your input data.")
    }
  }
  if (length(id_pos_t) != ncol(detp)) {
    message("Warning: IDs in the target and methylation files do not fully match. Only the intersection of both will be used in subsequent analyses. Please verify your input data.")
  }
  ## pData
  pData <- targets[index_tf_detp, ]
  detp <- detp[, index_detp]
  methyl <- methyl[, index_methyl]
  unmethyl <- unmethyl[, index_unmethyl]
  # names
  rownames(pData) <- GSM_id_pos_t
  colnames(detp) <- GSM_id_pos_t
  colnames(methyl) <- GSM_id_pos_t
  colnames(unmethyl) <- GSM_id_pos_t
  ###### write and read####
  colnames(methyl) <- paste0(colnames(methyl), " Methylated_signal")
  colnames(unmethyl) <- paste0(colnames(unmethyl), " Unmethylated_signal")
  merged_signals <- cbind(methyl, unmethyl)
  rm(methyl, unmethyl)
  temp_path <- tempdir()
  fwrite(merged_signals,
    paste0(temp_path, "/", "merged_signals.csv"),
    row.names = TRUE, sep = ",", nThread = nThread
  )
  rm(merged_signals)
  # Making RAW GenomicMethylSet
  if (DNAmType == "450k") {
    message("The platform is 450k, proceeding with 450k-specific reading.")
    GenomicMethylSet <- readGEORawFile(
      filename = paste0(temp_path, "/", "merged_signals.csv"),
      pData = pData,
      Uname = "Unmethylated_signal",
      Mname = "Methylated_signal",
      array = "IlluminaHumanMethylation450k",
      sep = ","
    )
  } else {
    if (DNAmType == "EPIC") {
      message("The platform is EPIC, proceeding with EPIC-specific reading.")
      GenomicMethylSet <- readGEORawFile(
        filename = paste0(temp_path, "/", "merged_signals.csv"),
        pData = pData,
        Uname = "Unmethylated_signal",
        Mname = "Methylated_signal",
        array = "IlluminaHumanMethylationEPIC",
        annotation = "ilm10b2.hg19",
        sep = ","
      )
    } else {
      stop("Chip type not recognized. Please confirm whether it is 450k or EPIC.")
    }
  }
  unlink(temp_path, recursive = TRUE)
  preprocessMethod <- "Raw (no normalization or bg correction)"
  names(preprocessMethod) <- "rg.norm"
  GenomicMethylSet@preprocessMethod <- preprocessMethod
  DNAm_env <- new.env()
  DNAm_env$detp <- detp
  DNAm_env$GenomicMethylSet <- GenomicMethylSet
  DNAm_env$UniqueID <- UniqueID
  DNAm_env$Type <- "Intensity"
  return(DNAm_env)
}
