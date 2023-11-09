#' Convert vcfR Object to A Table
#'
#' vcf2table is used to convert a vcfR object to a table, then this table can be used for select_sample_and_SNP() to construct QTLseq object.
#'
#' @param x A vcfR object
#'
#' @return A table
#' @export
#'
#' @examples
#' library(vcfR)
#' library(easyQTLseq)
#' file_path <- system.file("extdata", "A07.SNPs.vcf.gz", package = "easyQTLseq")
#' x <- read.vcfR(file = file_path)
#' data <- vcf2table(x = x)
vcf2table <- function(x) {
  data <- as_tibble(x@fix) %>% select(CHROM, POS, REF, ALT)

  gt <- extract.gt(x = x, element = "GT", return.alleles = TRUE)
  gt[gt=="."] <- "./."

  sample <- colnames(gt)

  ad <- extract.gt(x = x, element = "AD")

  gq <- extract.gt(x = x, element = "GQ")
  gq <- type.convert(gq, as.is = TRUE)

  colnames(gt) <- paste(sample, "GT", sep = ".")
  colnames(ad) <- paste(sample, "AD", sep = ".")
  colnames(gq) <- paste(sample, "GQ", sep = ".")

  df <- cbind(data, gt, ad, gq) %>% as_tibble()
  return(df)
}
