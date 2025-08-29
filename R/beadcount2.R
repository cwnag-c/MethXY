#' beadcount2
#'
#' @param x rgSet_extend
#' @importFrom minfi getNBeads getProbeInfo getManifestInfo
#' @return beadcount
#' @export
beadcount2<-function (x)
{
  nb <- getNBeads(x)
  typeIadd <- getProbeInfo(x, type = "I")
  typeImatchA <- match(typeIadd$AddressA, rownames(nb))
  typeImatchB <- match(typeIadd$AddressB, rownames(nb))
  typeIIadd <- getProbeInfo(x, type = "II")
  typeIImatch <- match(typeIIadd$Address, rownames(nb))
  nbcg <- nb
  locusNames <- getManifestInfo(x, "locusNames")
  bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                    dimnames = list(locusNames, sampleNames(x)))
  TypeII.Name <- getProbeInfo(x, type = "II")$Name
  bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,
  ]
  TypeI <- getProbeInfo(x, type = "I")
  bcB <- bc_temp
  bcA <- bc_temp
  bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB, ]
  bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
  bcB3 <- which(bcB < 3)
  bcA3 <- which(bcA < 3)
  rm(bcB,bcA);gc()
  bc_temp[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
  bc_temp[bcA3] <- NA
  bc_temp[bcB3] <- NA
  return(bc_temp)
}
