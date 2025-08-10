
#' Get QTL Exceed 95\% and 99\% Confidence Interval
#'
#' This function will export <outPrefix>.95CI.csv, <outPrefix>.99CI.csv, <outPrefix>.ED4QTL.csv, <outPrefix>.GprimeQTL.csv and figures of target chromosomes.
#'
#' @param x The QTLseq S3 object.
#' @param outPrefix The prefix of outfiles.
#' @param minN The minimum number of SNP in a window.
#'
#' @return
#' NULL
#' @export
#'
#' @examples
#' #' library(easyQTLseq)
#' # Example with sample data from a GATK table.
#' file_path <- system.file("extdata", "subset.table.gz", package = "easyQTLseq")
#' # readr::read_tsv() has a faster speed than read.table() when reading a file.
#' data <- readr::read_tsv(file = file_path)
#' x <- select_sample_and_SNP(data = data, highP = "qY", lowP = "R3", highB = "Y", lowB = "R", popType = "F2", bulkSize = c(30, 30))
#' x_filter <- filterDP(x = x)
#' x_filter <- calc_index_etc(x = x_filter, outPrefix = "outprefix", winSize = 2000000, winStep = 200000)
#' getQTL_and_exportFigure(x = x_filter, outPrefix = "outprefix", minN = 20)
getQTL_and_exportFigure <- function(x, ...) {
  UseMethod("getQTL_and_exportFigure")
}

#' @rdname getQTL_and_exportFigure
#' @export
getQTL_and_exportFigure.WithParent <- function(x, outPrefix, minN) {
  data <- x$slidwin
  CI95 <- getQTL_for_index(data = data, CI = 95, n = minN, export = TRUE, filename = paste(outPrefix, "95CI.csv", sep = "."))
  CI99 <- getQTL_for_index(data = data, CI = 99, n = minN, export = TRUE, filename = paste(outPrefix, "99CI.csv", sep = "."))
  ED4qtl <- getQTL_for_ed4(data = data, n = minN, export = TRUE, filename = paste(outPrefix, "ED4QTL.csv", sep = "."))
  Gprimeqtl <- getQTL_for_Gprime(data = data, n = minN, export = TRUE, filename = paste(outPrefix, "GprimeQTL.csv", sep = "."))

  ## plot target chrom
  options(scipen=200)
  plotTargetChrom_for_index(x = x, CI = 95, minN = minN, outPrefix = outPrefix)
  plotTargetChrom_for_index(x = x, CI = 99, minN = minN, outPrefix = outPrefix)
  plotTargetChrom_for_ED4(x = x, minN = minN, outPrefix = outPrefix)
  plotTargetChrom_for_Gprime(x = x, minN = minN, outPrefix = outPrefix)
}

#' @rdname getQTL_and_exportFigure
#' @export
getQTL_and_exportFigure.WithoutParent <- function(x, outPrefix, minN) {
  data <- x$slidwin
  ED4qtl <- getQTL_for_ed4(data = data, n = minN, export = TRUE, filename = paste(outPrefix, "ED4QTL.csv", sep = "."))
  Gprimeqtl <- getQTL_for_Gprime(data = data, n = minN, export = TRUE, filename = paste(outPrefix, "GprimeQTL.csv", sep = "."))

  ## plot target chrom
  options(scipen=200)
  plotTargetChrom_for_ED4(x = x, minN = minN, outPrefix = outPrefix)
  plotTargetChrom_for_Gprime(x = x, minN = minN, outPrefix = outPrefix)
}
