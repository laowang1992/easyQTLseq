#' @export
plotTargetChrom_for_index <- function(x, CI = 95, minN = 0, outPrefix){
  interval <- read_csv(file = paste(outPrefix, paste(CI, "CI", sep = ""), "csv", sep = "."), show_col_types = FALSE)
  #interval <- get(paste("CI", CI, sep = ""))
  if (nrow(interval) > 0) {
    for (i in unique(interval$CHROM)) {
      d <- x$slidwin %>% filter(CHROM == i, nSNPs >= minN) %>%
        select(CHROM, POS = POS,
               CIupper = all_of(paste("CI", CI, "upper", sep = "")),
               CIlower = all_of(paste("CI", CI, "lower", sep = "")),
               delta.index)
      inter <- interval %>% filter(CHROM == i)
      chrLen <- x$chrLen$Len[x$chrLen$CHROM == i]

      p_index <- ggplot() +
        geom_rect(data = inter, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300") +
        geom_line(data = d, aes(x = POS, y = CIupper), color = "gray60") +
        geom_line(data = d, aes(x = POS, y = CIlower), color = "gray60") +
        geom_line(data = d, aes(x = POS, y = delta.index), color = "blue") +
        scale_x_continuous(limits = c(0, chrLen), labels = scales::comma_format(scale = 1e-6, suffix = "Mb")) +
        scale_y_continuous(limits = c(-1, 1)) +
        geom_hline(aes(yintercept=0)) +
        labs(x = NULL, y = "Delta SNP index") +
        theme_half_open()
      ggsave(p_index, filename = paste(outPrefix, i, paste(CI, "CI", sep = ""), "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p_index, filename = paste(outPrefix, i, paste(CI, "CI", sep = ""), "pdf", sep = "."), height = 3.5, width = 4.5)

      print(paste(i, "has been done...", sep = " "))
    }
  }
}

#' @export
plotTargetChrom_for_ED4 <- function(x, minN = 0, outPrefix){
  interval <- read_csv(file = paste(outPrefix, "ED4QTL.csv", sep = "."), show_col_types = FALSE)
  if (nrow(interval) > 0) {
    for (i in unique(interval$CHROM)) {
      d <- x$slidwin %>% filter(CHROM == i, nSNPs >= minN) %>%
        select(CHROM, POS = POS,
               ED, ED4, ED4threshold)
      inter <- interval %>% filter(CHROM == i)
      chrLen <- x$chrLen$Len[x$chrLen$CHROM == i]

      p_ed4 <- ggplot() +
        geom_rect(data = inter, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),fill = '#FF3300', color = "#FF3300") +
        geom_line(data = d, aes(x = POS, y = ED4threshold), color = "gray60", linetype = "dashed") +
        geom_line(data = d, mapping = aes(x = POS, y = ED4), color = "blue") +
        scale_x_continuous(limits = c(0, chrLen), labels = scales::comma_format(scale = 1e-6, suffix = "Mb")) +
        scale_y_continuous(limits = c(0, max(d$ED4, na.rm = TRUE))) +
        labs(x = NULL, y = bquote(ED^4)) +
        theme_half_open()
      ggsave(p_ed4, filename = paste(outPrefix, i, "ED4", "png", sep = "."), height = 3.5, width = 4.5, dpi = 500)
      ggsave(p_ed4, filename = paste(outPrefix, i, "ED4", "pdf", sep = "."), height = 3.5, width = 4.5)

      print(paste(i, "has been done...", sep = " "))
    }
  }
}
