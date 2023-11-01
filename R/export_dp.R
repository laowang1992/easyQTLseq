
#' Export Allele Depth
#'
#' export_dp writes allele depth information to <outPrefix>.Depth_information.txt and <outPrefix>.Depth_information.csv
#'
#' @param x A QTLseq S3 object.
#' @param outPrefix The prefix of output file.
#'
#' @return NULL
#' @export
#'
#' @examples
#' library(easyQTLseq)
#' # Example with sample data from a GATK table.
#' file_path <- system.file("extdata", "subset.table.gz", package = "easyQTLseq")
#' # readr::read_tsv() has a faster speed than read.table() when reading a file.
#' data <- readr::read_tsv(file = file_path)
#' x <- importData(data = data, highP = "qY", lowP = "R3", highB = "Y", lowB = "R", popType = "F2", bulkSize = c(30, 30))
#' x_filter <- filterDP(x = x)
#' export_dp(x = x, outPrefix = "outfile")
export_dp <- function(x, ...) {
  UseMethod("export_dp")
}

#' @rdname export_dp
#' @export
export_dp.WithParent <- function(x, outPrefix){
  # export allele Depth data
  outtb <- x$data %>%
    select(CHROM, POS, REF, ALT, HB.HP.AD, HB.LP.AD, LB.HP.AD, LB.LP.AD)
  colnames(outtb) <- c("CHROM", "POS", "REF", "ALT",
                       paste(x$highB, x$highP, "depth", sep = "."),
                       paste(x$highB, x$lowP, "depth", sep = "."),
                       paste(x$lowB, x$highP, "depth", sep = "."),
                       paste(x$lowB, x$lowP, "depth", sep = "."))
  write_tsv(x = outtb, file = paste(outPrefix, "Depth_information.txt", sep = "."))
  write_csv(x = outtb, file = paste(outPrefix, "Depth_information.csv", sep = "."))
}
#' @rdname export_dp
#' @export
export_dp.WithoutParent <- function(x, outPrefix){
  # export allele Depth data
  outtb <- x$data %>%
    select(CHROM, POS, REF, ALT, HB.REF.AD, HB.ALT.AD, LB.REF.AD, LB.ALT.AD)
  colnames(outtb) <- c("CHROM", "POS", "REF", "ALT",
                       paste(highB, "ref", "depth", sep = "."),
                       paste(highB, "alt", "depth", sep = "."),
                       paste(lowB, "ref", "depth", sep = "."),
                       paste(lowB, "alt", "depth", sep = "."))
  write_tsv(x = outtb, file = paste(outPrefix, "Depth_information.txt", sep = "."))
  write_csv(x = outtb, file = paste(outPrefix, "Depth_information.csv", sep = "."))
}
