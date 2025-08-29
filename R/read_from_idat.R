#' read_from_idat
#' @importFrom utils write.csv
#' @importFrom minfi read.metharray.exp sampleNames
#' @param targets_dir Location of the targets file
#' @param idat_dir Location of the idat
#' @param UniqueID UniqueID
#' @param pattern "Name_ID_Position" or "ID_Position"
#'
#' @return rgset in env
#' @export
#'
read_from_idat<-function(targets_dir,idat_dir,UniqueID,pattern="GM_ID_Position"
){
  if(pattern=="GM_ID_Position"){
    targets<-read.csv(targets_dir)
    required_cols <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position")
    if(all(required_cols %in% colnames(targets))){
      targets$Basename<-paste0(idat_dir,"/",targets$Sample_Name,"_",targets$Sentrix_ID,"_",targets$Sentrix_Position)
    }else{
      missing_cols <- required_cols[!required_cols %in% colnames(targets)]
      stop(paste("'targets' is missing the following required column(s):", paste(missing_cols, collapse = ", ")))
    }
  }
  if(pattern=="ID_Position"){
    targets<-read.csv(targets_dir)
    required_cols <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position")
    if(all(required_cols %in% colnames(targets))){
      targets$Basename<-paste0(idat_dir,"/",targets$Sentrix_ID,"_",targets$Sentrix_Position)
    }else{
      missing_cols <- required_cols[!required_cols %in% colnames(targets)]
      stop(paste("'targets' is missing the following required column(s):", paste(missing_cols, collapse = ", ")))
    }
  }
  rgSet <- read.metharray.exp(targets = targets,extended = T);gc()
  colData(rgSet)$Basename<-NULL;colData(rgSet)$filenames<-NULL
  #######创建环境###
  DNAm_env <- new.env()
  DNAm_env$rgSet<-rgSet
  DNAm_env$Type<-"IDAT"
  DNAm_env$UniqueID<-UniqueID
  return(DNAm_env)
}
