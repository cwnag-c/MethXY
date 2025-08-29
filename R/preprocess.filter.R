
#' preprocess.filter
#' @importFrom minfi getSex addSex sampleNames getAnnotation featureNames detectionP mapToGenome preprocessRaw
#' @importFrom dplyr bind_rows
#' @param DNAm_env The object returned by read_from_idat or read_from_intensity
#'
#' @return An environment containing:
#'    a GenomicMethylSet without Y chromosome ,
#'    a GenomicMethylSet with only Y chromosome ,
#'    the targets file, a UniqueID,
#'    and a message object containing filtering information.
#' @export
#'
preprocess.filter<-function(DNAm_env){
  Message<-list()
  UniqueID<-DNAm_env$UniqueID
  if(DNAm_env$Type=="IDAT"){
    MethylSet<-preprocessRaw(DNAm_env$rgSet)
    GenomicMethylSet<-mapToGenome(MethylSet);rm(MethylSet)
    detP <- detectionP(DNAm_env$rgSet)# run much time
    beadc<-beadcount2(DNAm_env$rgSet)
    rm(rgSet,envir = DNAm_env)
  }
  if(DNAm_env$Type=="Intensity"){
    GenomicMethylSet<-DNAm_env$GenomicMethylSet
    detP <- DNAm_env$detp
    rm(list = c("GenomicMethylSet", "detp"), envir = DNAm_env)
  }
  if(GenomicMethylSet@annotation["array"]=="IlluminaHumanMethylation450k"){
    DNAmType="450k"
  }else{
    if(GenomicMethylSet@annotation["array"]=="IlluminaHumanMethylationEPIC"){
      DNAmType="EPIC"
    }
  }

  #--1.filter samples:discordance between reported and estimated sex----
  #--1.1 predict sex----#
  predictedSex <- getSex(GenomicMethylSet)
  GenomicMethylSet<-addSex(GenomicMethylSet, sex = predictedSex)
  colData(GenomicMethylSet)$xMed<-NULL;colData(GenomicMethylSet)$yMed<-NULL
  #--1.2.First check for samples with missing Sex; if Gender is missing, fill it using the predicted sex.#
  targets<-as.data.frame(colData(GenomicMethylSet))
  sex_merge<-apply(targets,1,function(x){
    tmp <- ifelse(is.na(x["Gender"]), x["predictedSex"], x["Gender"])
    return(tmp)
  })
  if(identical(match(names(sex_merge), rownames(colData(GenomicMethylSet))),
               1:ncol(GenomicMethylSet))){
    colData(GenomicMethylSet)$Gender<-sex_merge
  }else{
    stop("The row names of colData(GenomicMethylSet) do not match the names of sex_merge. Cannot safely assign gender values.")
  }
  targets<-as.data.frame(colData(GenomicMethylSet))
  discordance<-targets[targets$Gender!=targets$predictedSex,]
  targets<-targets[targets$Gender==targets$predictedSex,]
  if(nrow(discordance)==0){
    ##message#
    message("All samples are cordance between reported and predicted sex")
    remove_sample_text<-"All samples are cordance between reported and predicted sex"
    Message$remove_sample_text<-remove_sample_text
  }else{
    #Check whether pd$Sample_Name is a subset of sampleNames(DNAm_env$GenomicMethylSet).
    GenomicMethylSet<-GenomicMethylSet[,match(rownames(targets), sampleNames(GenomicMethylSet))]
    ##message#
    discordance<-discordance$Sample_Name
    number_discordance<-length(discordance)
    message(paste0(number_discordance," samples !! are discordance and removed :"),
            "\n",paste0(discordance,"\n"))
    remove_sample_text<-paste0(paste0(number_discordance," samples !! are discordance and removed :"),
                               "\n",paste0(discordance,"\n"))
    Message$remove_sample_text<-remove_sample_text
  }
  #--2.separate Y chromosomes-----
  Annotation<-as.data.frame(getAnnotation(GenomicMethylSet));Annotation<-Annotation[,c(1,4)]
  lines_Y = which(Annotation[,'chr'] == 'chrY')
  lines_X = which(Annotation[,'chr'] == 'chrX')
  line_atuo_x =  which(Annotation[,'chr'] != 'chrY')
  line_atuo =  which(Annotation[,'chr'] != 'chrY' & Annotation[,'chr'] != 'chrX')
  probeY = rownames(Annotation[lines_Y,])
  probeX = rownames(Annotation[lines_X,])
  probeAtuo = rownames(Annotation[line_atuo,])
  probeAtuo_x = rownames(Annotation[line_atuo_x,])

  GenomicMethylSet_Y<-GenomicMethylSet[probeY,]
  detP_Y<-detP[probeY,]
  if(DNAm_env$Type=="IDAT"){
    beadc_Y<-beadc[probeY,]
  }
  GenomicMethylSet_Y<-GenomicMethylSet_Y[,which(colData(GenomicMethylSet_Y)$Gender=="M")]
  detP_Y<-detP_Y[,match(sampleNames(GenomicMethylSet_Y),colnames(detP_Y))]
  if(DNAm_env$Type=="IDAT"){
    beadc_Y<-beadc_Y[,match(sampleNames(GenomicMethylSet_Y),colnames(beadc_Y))]
  }
  #auto and X
  detP<-detP[probeAtuo_x ,]
  if(DNAm_env$Type=="IDAT"){
    beadc<-beadc[probeAtuo_x,]
  }
  GenomicMethylSet<-GenomicMethylSet[probeAtuo_x,]

  #--3.filter samples and probes using detP and beadcount-----
  #######3.1.filter samples:poor quality----------------#
  #A sample is considered bad if more than 10% of its CpG probes have detection p-values (detP) ≥ 0.01.

  ####Auto+X
  remove_samples<-colnames(detP)[colSums(detP >= 0.01) >= 0.10 * nrow(detP)]
  if(length(remove_samples)==0){#If the number of rows is equal, it means no samples were excluded due to quality issues.
    Message$remove_sample_dueto_detP<-"All samples are retained"
  }else{
    retain_samples<-setdiff(sampleNames(GenomicMethylSet),remove_samples)
    GenomicMethylSet<-GenomicMethylSet[,retain_samples]
    Message$remove_sample_dueto_detP<-remove_samples
  }
  ###Y
  remove_samples<-colnames(detP_Y)[colSums(detP_Y >= 0.01) >= 0.10 * nrow(detP_Y)]
  if(length(remove_samples)==0){#If the number of rows is equal, it means no samples were excluded due to quality issues.
    Message$remove_sample_dueto_detP_Y<-"All samples are retained"
  }else{
    retain_samples<-setdiff(sampleNames(GenomicMethylSet_Y),remove_samples)
    GenomicMethylSet_Y<-GenomicMethylSet_Y[,retain_samples]
    Message$remove_sample_dueto_detP_Y<-remove_samples
  }
  #######3.2.1 filter probes: poor quality(detectionP)--------#
  ##auto+x##
  remove_probes<-rownames(detP)[rowSums(detP >= 0.01) >= 0.20 * ncol(detP)]
  if(length(remove_probes)== 0){
    Message$remove_cpg_dueto_detP<-"None"
    rm(detP)
  }else{
    retain_probes<-setdiff(featureNames(GenomicMethylSet),remove_probes)
    GenomicMethylSet<-GenomicMethylSet[retain_probes,]
    Message$remove_cpg_dueto_detP<-remove_probes
    rm(detP)
  }
  ###Y###
  remove_probes_Y<-rownames(detP_Y)[rowSums(detP_Y >= 0.01) >= 0.20 * ncol(detP_Y)]
  if(length(remove_probes_Y)== 0){
    Message$remove_cpg_dueto_detP_Y<-"None"
    rm(detP_Y)
  }else{
    retain_probes<-setdiff(featureNames(GenomicMethylSet_Y),remove_probes)
    GenomicMethylSet_Y<-GenomicMethylSet_Y[retain_probes,]
    Message$remove_cpg_dueto_detP_Y<-remove_probes_Y
    rm(detP_Y)
  }
  #######3.2.2 filter probes: poor quality(beadcount)--------#
  if(DNAm_env$Type=="IDAT"){
    ##auto+x##
    remove_probes <- rownames(beadc)[rowSums(is.na(beadc)) >= 0.05*(ncol(beadc))]
    if(length(remove_probes)==0){
      Message$remove_cpg_dueto_beadc<-"None"
      rm(beadc)
    }else{
      retain_probes<-setdiff(featureNames(GenomicMethylSet),remove_probes)
      GenomicMethylSet<-GenomicMethylSet[retain_probes,]
      Message$remove_cpg_dueto_beadc<-remove_probes
      rm(beadc)
    }
    ###Y###
    remove_probes_Y <- rownames(beadc_Y)[rowSums(is.na(beadc_Y)) >= 0.05*(ncol(beadc_Y))]
    if(length(remove_probes_Y)==0){
      Message$remove_cpg_dueto_beadc_Y<-"None"
      rm(beadc_Y)
    }else{
      retain_probes<-setdiff(featureNames(GenomicMethylSet_Y),remove_probes)
      GenomicMethylSet_Y<-GenomicMethylSet_Y[retain_probes,]
      Message$remove_cpg_dueto_beadc_Y<-remove_probes_Y
      rm(beadc_Y)
    }
  }
  #######3.3.filter probes: non-CpG probes-------#
  CpG_probes <- rownames(Annotation)[which(substr(rownames(Annotation), 1, 2) == 'cg')]
  GenomicMethylSet<-GenomicMethylSet[intersect(featureNames(GenomicMethylSet),CpG_probes),]
  GenomicMethylSet_Y<-GenomicMethylSet_Y[intersect(featureNames(GenomicMethylSet_Y),CpG_probes),]
  #######3.4.filter probes: non-specific probes, and polymorphic probes(Zhou et al. 2017)----#
  #data("mask_cpg", package = "MethXY", envir = environment())
  if(DNAmType=="450k"){
    GenomicMethylSet<-GenomicMethylSet[setdiff(featureNames(GenomicMethylSet),maskname_450k),]
    GenomicMethylSet_Y<-GenomicMethylSet_Y[setdiff(featureNames(GenomicMethylSet_Y),maskname_450k),]
  }else{
    if(DNAmType=="EPIC"){
      GenomicMethylSet<-GenomicMethylSet[setdiff(featureNames(GenomicMethylSet),maskname_850k),]
      GenomicMethylSet_Y<-GenomicMethylSet_Y[setdiff(featureNames(GenomicMethylSet_Y),maskname_850k),]
    }
  }
  #######3.5.XY:filter probes: X-transposed region/repetitive elements/cancer-testis gene----#
  if(DNAmType=="450k"){
    #data("filter_xy_450k", package = "MethXY", envir = environment())
    filter_xy<-filter_xy_450k
  }else{
    if(DNAmType=="EPIC"){
      #data("filter_xy_850k", package = "MethXY", envir = environment())
      filter_xy<-filter_xy_850k
    }
  }
  remove_probe_Y<-filter_xy[(!(filter_xy$MASK_X_transposed_region==F &
                                 filter_xy$MASK_Cancer_testis_gene==F &
                                 filter_xy$MASK_SINE_LINE==F &
                                 filter_xy$MASK_AnyRepetitive==F))&
                              filter_xy$chr=="chrY",]$probeID
  remove_probe_X<-filter_xy[(!(filter_xy$MASK_X_transposed_region==F &
                                 filter_xy$MASK_Cancer_testis_gene==F &
                                 filter_xy$MASK_SINE_LINE==F &
                                 filter_xy$MASK_AnyRepetitive==F))&
                              filter_xy$chr=="chrX",]$probeID

  GenomicMethylSet<-GenomicMethylSet[setdiff(featureNames(GenomicMethylSet),remove_probe_X),]
  GenomicMethylSet_Y<-GenomicMethylSet_Y[setdiff(featureNames(GenomicMethylSet_Y),remove_probe_Y),]
  Message$remove_probeX<-remove_probe_X
  Message$remove_probeY<-remove_probe_Y
  ######4.Return#########
  DNAm_env$GenomicMethylSet<-GenomicMethylSet
  DNAm_env$GenomicMethylSet_Y<-GenomicMethylSet_Y
  DNAm_env$Message<-Message
  DNAm_env$UniqueID<-UniqueID
}
