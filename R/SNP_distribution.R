#' Analyzing the Distribution of SNPs on Chromosomes
#'
#' The SNP number in every 1 Mb bins will be shown in histogram along chromosome.
#'
#' @param x The QTLseq S3 object.
#' @param outPrefix The prefix of output file. <outPrefix>.SNP_number_per_chr.txt|csv and <outPrefix>.SNP_distribution_histogram.pdf|png
#' @param targetChr Only chromosomes in targetChr will be analyzed, default is all chromosomes.
#' @param chrLabel The labels of chromosomes, default is the same of targetChr.
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
#' x <- select_sample_and_SNP(data = data, highP = "qY", lowP = "R3", highB = "Y", lowB = "R", popType = "F2", bulkSize = c(30, 30))
#' x_filter <- filterDP(x = x)
#' SNP_distribution(x = x_filter, outPrefix = "outprefix",
#'                  targetChr = c("scaffoldA01", "scaffoldA07", "scaffoldA09"),
#'                  chrLabel = c("A01", "A07", "A09"))
SNP_distribution <- function(x, ...) {
  UseMethod("SNP_distribution")
}

#' @rdname SNP_distribution
#' @export
SNP_distribution.QTLseq <- function(x, outPrefix, targetChr, chrLabel){
  # 计算每条染色体位点个数
  SNPnumber <- x$data %>% group_by(CHROM) %>% count()
  write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
  write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))

  # SNP分布图
  options(scipen = 200)
  colourCount = length(targetChr)
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  chr <- tibble(CHROM = targetChr, LABEL = chrLabel)
  chr$LABEL <- factor(chr$LABEL, levels = chr$LABEL)
  Phist <- chr %>%
    left_join(x$data, by = "CHROM") %>%
    ggplot(aes(x = POS)) +
    geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +
    labs(x = NULL, y = "SNP Count / 1Mb", fill = "Chrom") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 2.1) +
    scale_fill_manual(values = getPalette(colourCount)) +
    theme_half_open() +
    theme(strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color = NA, fill = NA),
          legend.position = "NULL") +
    facet_grid(LABEL ~ .)
  ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.pdf", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5)
  ggsave(Phist, filename = paste(outPrefix, "SNP_distribution_histogram.png", sep = "."), width = 8, height = dim(chr)[[1]] * 0.45 + 0.5, dpi = 500)
}
