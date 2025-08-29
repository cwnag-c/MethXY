#' preprocess.norm
#' @importFrom minfi sampleNames getAnnotation featureNames ratioConvert GenomicMethylSet preprocessQuantile
#' @param DNAm_env Output returned by preprocess.filter
#' @param Norm  "preprocessQuantile"
#' @param save_result_locally Save results to local directory; default is FALSE
#' @param save_dir Directory to save results;default: "./"
#' @param return_res Return results; default is TRUE
#'
#' @return The environment contains:
#'    GRSet without sex chromosomes,GRSet with X ,GRSet with Y ,
#'    and the targets table
#' @export
#'
preprocess.norm<-function(DNAm_env,
                          Norm="preprocessQuantile",
                          save_result_locally=F,
                          save_dir="./",
                          return_res=T){
  UniqueID<-DNAm_env$uniqueID
  DNAm_GRSet<-list()
  ##reconstruct
  Annotation<-as.data.frame(getAnnotation(DNAm_env$GenomicMethylSet));Annotation<-Annotation[,c(1,4)]
  lines_X = which(Annotation[,'chr'] == 'chrX')
  probeX = rownames(Annotation[lines_X,])
  DNAm_env$GenomicMethylSet_X<-DNAm_env$GenomicMethylSet[probeX,]
  probeAuto<-setdiff(featureNames(DNAm_env$GenomicMethylSet),probeX)
  DNAm_env$GenomicMethylSet<-DNAm_env$GenomicMethylSet[probeAuto,]# Contains only autosomes

  if(is.null(Norm)){
    DNAm_GRSet$GRSet<-ratioConvert(DNAm_env$GenomicMethylSet,what="M", keepCN = F)
    DNAm_GRSet$GRSet_Y<-ratioConvert(DNAm_env$GenomicMethylSet_Y,what="M", keepCN = F)
    DNAm_GRSet$GRSet_X<-ratioConvert(DNAm_env$GenomicMethylSet_X,what="M", keepCN = F)
    rm(list = c("GenomicMethylSet", "GenomicMethylSet_Y", "GenomicMethylSet_X"), envir = DNAm_env)
  }else{
    if(Norm=="preprocessQuantile"){
      ## Note: For Quantile, only M and CN are returned; beta values are not
      DNAm_GRSet$GRSet<-preprocessQuantile(DNAm_env$GenomicMethylSet,sex = colData(DNAm_env$GenomicMethylSet)$Gender)
      DNAm_GRSet$GRSet_Y<-preprocessQuantile(DNAm_env$GenomicMethylSet_Y,sex = colData(DNAm_env$GenomicMethylSet_Y)$Gender)
      DNAm_GRSet$GRSet_X<-preprocessQuantile(DNAm_env$GenomicMethylSet_X,sex = colData(DNAm_env$GenomicMethylSet_X)$Gender)
      DNAm_GRSet$GRSet@assays@data@listData[["CN"]]<-NULL
      DNAm_GRSet$GRSet_Y@assays@data@listData[["CN"]]<-NULL
      DNAm_GRSet$GRSet_X@assays@data@listData[["CN"]]<-NULL
      rm(list = c("GenomicMethylSet", "GenomicMethylSet_Y", "GenomicMethylSet_X"), envir = DNAm_env)
    }else{
      message("Warning: Invalid input")
    }
  }
  if(save_result_locally){
    DNAm_GRSet$Message<-DNAm_env$Message
    DNAm_GRSet$UniqueID<-DNAm_env$UniqueID
    save(file=paste0(save_dir,"/",UniqueID,".Rdata"),
         list = c("DNAm_GRSet"))
  }

  if(return_res){
    DNAm_env$GRSet<-DNAm_GRSet$GRSet
    DNAm_env$GRSet_Y<-DNAm_GRSet$GRSet_Y
    DNAm_env$GRSet_X<-DNAm_GRSet$GRSet_X
  }
}
