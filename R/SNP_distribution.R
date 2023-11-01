#' @export
SNP_distribution <- function(x, ...) {
  UseMethod("SNP_distribution")
}
#' @export
SNP_distribution.QTLseq <- function(x, outPrefix, targetChr, chrLabel){
  # 计算每条染色体位点个数
  SNPnumber <- x$data %>% group_by(CHROM) %>% count()
  write_tsv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.txt", sep = "."))
  write_csv(x = SNPnumber, file = paste(outPrefix, "SNP_number_per_chr.csv", sep = "."))

  # SNP分布图
  options(scipen = 200)
  colourCount = nrow(chr)
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
