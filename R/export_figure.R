
#' Export Figures of Sliding Window SNP Index and Euclidean Distance
#'
#' @param x A QTLseq S3 object with sliding window of SNP index and Euclidean Distance.
#' @param outPrefix The prefix of figures.
#' @param targetChr Target chromosome to be drawn in figures, default is all chromosomes in the data.
#' @param chrLabel The label for chromosome shown in figures, default is chromosome names in the data.
#' @param minN The windows containing SNPs less than minN will be omitted in figures.
#' @param width The width of figures.
#' @param height The height of figures.
#' @param color A vector of colors for chromosomes.
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
#' x_filter <- calc_index_etc(x = x_filter, outPrefix = "outprefix", winSize = 2000000, winStep = 200000)
#' export_figure(x = x_filter, outPrefix = "outprefix", targetChr = c("scaffoldA01", "scaffoldA07", "scaffoldA09"), chrLabel = c("A01", "A07", "A09"), minN = 20, width = 6, height = 3)

export_figure <- function(x, ...) {
  UseMethod("export_figure")
}
#' @rdname export_figure
#' @export
export_figure.WithParent <- function(x, outPrefix, targetChr, chrLabel, minN, width, height,
                          color = c("#4197d8","#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d")) {
  if (missing(targetChr)) {
    targetChr <- unique(x$slidwin$CHROM)
  }
  if (missing(chrLabel)) {
    chrLabel <- targetChr
  }
  highB <- x$highB
  lowB <- x$lowB
  chr <- tibble(CHROM = targetChr, LABEL = chrLabel)
  #COLOR <- tibble(CHROM = targetChr, LABEL = chrLabel, color = rep(color, len = length(targetChr)))
  COLOR <- chr %>% mutate(color = rep(color, len = nrow(chr)))
  newLen <- chr %>% left_join(x$chrLen, by = "CHROM") %>% select(LABEL, Len)

  dataforPlot <- x$slidwin %>% filter(nSNPs >= minN) %>% right_join(chr, by = "CHROM") %>%
    addUp(len = newLen, group = "LABEL", pos = "POS", band = 0.005)


  p <- ggplot(dataforPlot$df, aes(x = POS_addUp, group = CHROM, color = CHROM)) +
    geom_vline(xintercept = dataforPlot$gaps, linetype = "dashed", color = "gray") +
    scale_x_continuous(breaks = dataforPlot$breaks, labels = dataforPlot$labels, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(breaks = COLOR$CHROM, values = COLOR$color) +
    theme_half_open() +
    theme(legend.position = "NULL",
          #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank())

  #
  p_dlt99 <- p +
    geom_line(aes(y = CI99upper), color = "#666666") +
    geom_line(aes(y = CI99lower), color = "#666666") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(y = "Delta SNP index")
  p_dlt99_L <- p_dlt99 + geom_line(aes(y = delta.index), linewidth = 1)
  ggsave(p_dlt99_L, filename = paste(outPrefix, "delta_SNP_index.99CI.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_dlt99_L, filename = paste(outPrefix, "delta_SNP_index.99CI.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_dlt99_P <- p_dlt99 + geom_point(aes(y = delta.index), size = 1.5)
  ggsave(p_dlt99_P, filename = paste(outPrefix, "delta_SNP_index.99CI.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_dlt99_P, filename = paste(outPrefix, "delta_SNP_index.99CI.point.png", sep = "."), width = width, height = height, dpi = 500)

  p_dlt95 <- p +
    geom_line(aes(y = CI95upper), color = "#666666") +
    geom_line(aes(y = CI95lower), color = "#666666") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(y = "Delta SNP index")
  p_dlt95_L <- p_dlt95 + geom_line(aes(y = delta.index), linewidth = 1)
  ggsave(p_dlt95_L, filename = paste(outPrefix, "delta_SNP_index.95CI.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_dlt95_L, filename = paste(outPrefix, "delta_SNP_index.95CI.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_dlt95_P <- p_dlt95 + geom_point(aes(y = delta.index), size = 1.5)
  ggsave(p_dlt95_P, filename = paste(outPrefix, "delta_SNP_index.95CI.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_dlt95_P, filename = paste(outPrefix, "delta_SNP_index.95CI.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_HB <- p +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = paste(highB, "SNP index", sep = " "))
  p_HB_L <- p_HB + geom_line(aes(y = HB.index), linewidth = 1)
  ggsave(p_HB_L, filename = paste(outPrefix, highB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_HB_L, filename = paste(outPrefix, highB, "SNP_index.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_HB_P <- p_HB + geom_point(aes(y = HB.index), size = 1.5)
  ggsave(p_HB_P, filename = paste(outPrefix, highB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_HB_P, filename = paste(outPrefix, highB, "SNP_index.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_LB <- p +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = paste(lowB, "SNP index", sep = " "))
  p_LB_L <- p_LB + geom_line(aes(y = LB.index), linewidth = 1)
  ggsave(p_LB_L, filename = paste(outPrefix, lowB, "SNP_index.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_LB_L, filename = paste(outPrefix, lowB, "SNP_index.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_LB_P <- p_LB + geom_point(aes(y = LB.index), size = 1.5)
  ggsave(p_LB_P, filename = paste(outPrefix, lowB, "SNP_index.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_LB_P, filename = paste(outPrefix, lowB, "SNP_index.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_ED <- p +
    coord_cartesian(ylim = c(0, max(dataforPlot$df$ED, na.rm = TRUE)*1.05)) +
    labs(y = "ED")
  p_ED_L <- p_ED + geom_line(aes(y = ED), linewidth = 1)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED_P <- p_ED + geom_point(aes(y = ED), size = 1.5)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_ED4 <- p +
    geom_line(aes(y = ED4threshold), color = "#666666", linetype = "dashed") +
    coord_cartesian(ylim = c(0-max(dataforPlot$df$ED4, na.rm = TRUE)*0.015, max(dataforPlot$df$ED4, na.rm = TRUE)*1.05)) +
    labs(y = bquote(ED^4))
  p_ED4_L <- p_ED4 +
    geom_line(aes(y = ED4), linewidth = 1)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED4_P <- p_ED4 +
    geom_point(aes(y = ED4), size = 1.5)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_G <- p +
    geom_hline(yintercept = 6.634897, color = "#666666", linetype = "dashed") +
    coord_cartesian(ylim = c(0-max(dataforPlot$df$Gprime, na.rm = TRUE)*0.015, max(dataforPlot$df$Gprime, na.rm = TRUE)*1.05)) +
    labs(y = "G' value")
  p_G_L <- p_G +
    geom_line(aes(y = Gprime), linewidth = 1)
  ggsave(p_G_L, filename = paste(outPrefix, "Gprime.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_G_L, filename = paste(outPrefix, "Gprime.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_G_P <- p_G +
    geom_point(aes(y = Gprime), size = 1.5)
  ggsave(p_G_P, filename = paste(outPrefix, "Gprime.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_G_P, filename = paste(outPrefix, "Gprime.point.png", sep = "."), width = width, height = height, dpi = 500)
}
#' @rdname export_figure
#' @export
export_figure.WithoutParent <- function(x, outPrefix, targetChr, chrLabel, minN, width, height,
                                        color = c("#4197d8","#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d")) {
  if (missing(targetChr)) {
    targetChr <- unique(x$slidwin$CHROM)
    chrLabel <- targetChr
  }
  if (missing(chrLabel)) {
    chrLabel <- targetChr
  }
  highB <- x$highB
  lowB <- x$lowB
  chr <- tibble(CHROM = targetChr, LABEL = chrLabel)
  #COLOR <- tibble(CHROM = targetChr, LABEL = chrLabel, color = rep(color, len = length(targetChr)))
  COLOR <- chr %>% mutate(color = rep(color, len = nrow(chr)))
  newLen <- chr %>% left_join(x$chrLen, by = "CHROM") %>% select(LABEL, Len)

  dataforPlot <- x$slidwin %>% filter(nSNPs >= minN) %>% right_join(chr, by = "CHROM") %>%
    addUp(len = newLen, group = "LABEL", pos = "POS", band = 0.005)

  p <- ggplot(dataforPlot$df, aes(x = POS_addUp, group = CHROM, color = CHROM)) +
    geom_vline(xintercept = dataforPlot$gaps, linetype = "dashed", color = "gray") +
    scale_x_continuous(breaks = dataforPlot$breaks, labels = dataforPlot$labels, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(breaks = COLOR$CHROM, values = COLOR$color) +
    theme_half_open() +
    theme(legend.position = "NULL",
          #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank())

  #
  p_ED <- p +
    coord_cartesian(ylim = c(0, max(dataforPlot$df$ED, na.rm = TRUE)*1.05)) +
    labs(y = "ED")
  p_ED_L <- p_ED + geom_line(aes(y = ED), linewidth = 1)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_L, filename = paste(outPrefix, "ED.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED_P <- p_ED + geom_point(aes(y = ED), size = 1.5)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED_P, filename = paste(outPrefix, "ED.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_ED4 <- p +
    geom_line(aes(y = ED4threshold), color = "#666666", linetype = "dashed") +
    coord_cartesian(ylim = c(0-max(dataforPlot$df$ED4, na.rm = TRUE)*0.015, max(dataforPlot$df$ED4, na.rm = TRUE)*1.05)) +
    labs(y = bquote(ED^4))
  p_ED4_L <- p_ED4 +
    geom_line(aes(y = ED4), linewidth = 1)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_L, filename = paste(outPrefix, "ED4.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_ED4_P <- p_ED4 +
    geom_point(aes(y = ED4), size = 1.5)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_ED4_P, filename = paste(outPrefix, "ED4.point.png", sep = "."), width = width, height = height, dpi = 500)

  #
  p_G <- p +
    geom_hline(yintercept = 6.634897, color = "#666666", linetype = "dashed") +
    coord_cartesian(ylim = c(0-max(dataforPlot$df$Gprime, na.rm = TRUE)*0.015, max(dataforPlot$df$Gprime, na.rm = TRUE)*1.05)) +
    labs(y = "G' value")
  p_G_L <- p_G +
    geom_line(aes(y = Gprime), linewidth = 1)
  ggsave(p_G_L, filename = paste(outPrefix, "Gprime.line.pdf", sep = "."), width = width, height = height)
  ggsave(p_G_L, filename = paste(outPrefix, "Gprime.line.png", sep = "."), width = width, height = height, dpi = 500)
  p_G_P <- p_G +
    geom_point(aes(y = Gprime), size = 1.5)
  ggsave(p_G_P, filename = paste(outPrefix, "Gprime.point.pdf", sep = "."), width = width, height = height)
  ggsave(p_G_P, filename = paste(outPrefix, "Gprime.point.png", sep = "."), width = width, height = height, dpi = 500)
}
