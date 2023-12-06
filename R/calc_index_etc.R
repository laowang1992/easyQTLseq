
#' Calculate Smoothed SNP Index and Euclidean Distance Using Sliding Window
#'
#' @param x The QTLseq S3 object.
#' @param outPrefix The prefix of <outPrefix>.SlidingWindow.txt and <outPrefix>.SlidingWindow.csv
#' @param winSize Window size for sliding window.
#' @param winStep Window step for sliding window.
#'
#' @return A QTLseq S3 object with sliding window of SNP index and Euclidean Distance.
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
#' x_filter <- calc_index_etc(x = x_filter, outPrefix = "outprefix", winSize = 2000000, winStep = 20000)
calc_index_etc <- function(x, ...) {
  UseMethod("calc_index_etc")
}
#' @rdname calc_index_etc
#' @export
calc_index_etc.WithParent <- function(x, outPrefix, winSize, winStep){
  dp <- x$data %>% select(HB.DP, LB.DP) %>% distinct()
  dltIndex_CI <- getQTLseqCI(df = dp, popType = x$popType, bulkSizeH = x$bulkSize[1], bulkSizeL = x$bulkSize[2], repN = 10000)
  df <- x$data %>%
    mutate(HB.index = HB.HP.AD / HB.DP,
                  LB.index = LB.HP.AD / LB.DP,
                  delta.index = HB.index- LB.index,
                  ED = sqrt((HB.index - LB.index)^2 + ((1-HB.index) - (1-LB.index))^2)) %>%
    left_join(dltIndex_CI, by = c("HB.DP", "LB.DP"))

  # 滑窗统计，windowscanr这个包好久不更新，说不定那天就不弄用了，还是自己写一个滑窗统计函数吧
  x$slidwin <- slidingWindow(df = df,
                           winSize = winSize,
                           winStep = winStep,
                           groups = "CHROM",
                           position = "POS",
                           values = c("HB.index", "LB.index", "delta.index", "ED",
                                      "CI95upper", "CI95lower", "CI99upper", "CI99lower"),
                           fun = "mean") %>% as_tibble() %>%
    select(CHROM, win_start, win_end,
                  HB.index = HB.index_mean, LB.index = LB.index_mean, delta.index = delta.index_mean,
                  CI95upper = CI95upper_mean, CI95lower = CI95lower_mean,
                  CI99upper = CI99upper_mean, CI99lower = CI99lower_mean,
                  ED= ED_mean, nSNPs = N) %>%
    mutate(ED4 = ED^4, POS = win_start/2 + win_end/2)
  write_tsv(x = x$slidwin %>% select(-POS) %>%
              rename(!!paste(x$highB, "index", sep = ".") := HB.index, !!paste(x$lowB, "index", sep = ".") := LB.index),
            file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
  write_csv(x = x$slidwin %>% select(-POS) %>%
              rename(!!paste(x$highB, "index", sep = ".") := HB.index, !!paste(x$lowB, "index", sep = ".") := LB.index),
            file = paste(outPrefix, "SlidingWindow.csv", sep = "."))
  x
}
#' @rdname calc_index_etc
#' @export
calc_index_etc.WithoutParent <- function(x, outPrefix, winSize, winStep){
  df <- x$data %>%
    mutate(HB.index = HB.REF.AD / HB.DP,
                  LB.index = LB.REF.AD / LB.DP,
                  ED = sqrt((HB.index - LB.index)^2 + ((1-HB.index) - (1-LB.index))^2))
  # 滑窗统计，windowscanr这个包好久不更新，说不定那天就不弄用了，还是自己写一个滑窗统计函数吧
  slidwin <- slidingWindow(df = df,
                           winSize = winSize,
                           winStep = winStep,
                           groups = "CHROM",
                           position = "POS",
                           values = c("ED"),
                           fun = "mean") %>% as_tibble() %>%
    select(CHROM, win_start, win_end,
                  ED= ED_mean, nSNPs = N) %>%
    mutate(ED4 = ED^4, POS = win_start/2 + win_end/2)
  write_tsv(x = slidwin %>% select(-POS), file = paste(outPrefix, "SlidingWindow.txt", sep = "."))
  write_csv(x = slidwin %>% select(-POS), file = paste(outPrefix, "SlidingWindow.csv", sep = "."))
  slidwin
}
