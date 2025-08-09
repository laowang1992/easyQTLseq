#' @export
getQTL_for_index <- function(data, chr = "CHROM", pos = "POS",
                   index = "delta.index", nSNPs = "nSNPs", CI = 95, n = 10,
                   export = TRUE, filename = "SignificantQTL_for_index.csv"){
  ## 整理格式
  df <- data %>%
    select(
      chr = all_of(chr),
      pos = all_of(pos),
      index = all_of(index),
      nSNPs = all_of(nSNPs),
      CIupper = all_of(paste("CI", CI, "upper", sep = "")),
      CIlower = all_of(paste("CI", CI, "lower", sep = ""))
    )

  ## 初始化结果变量
  interval <- tibble(CHROM = character(0), Start = numeric(0), End = numeric(0), Length = numeric(0))
  ##
  for (chromosome in unique(df$chr)) {
    # 提取子集
    subdf <- df %>% filter(chr == chromosome & nSNPs >= n)
    #subdf <- slidwin %>% filter(CHROM == "scaffoldA03")
    # 初始化变量
    start <- NULL
    end <- NULL
    sig <- FALSE
    for (i in seq_along(subdf$chr)) {
      if ((subdf$index[i]>subdf$CIupper[i] | subdf$index[i]<subdf$CIlower[i]) & (!sig)) {
        start <- subdf$pos[i]
        sig <- TRUE
      }
      if (((subdf$index[i]<subdf$CIupper[i] & subdf$index[i]>subdf$CIlower[i]) | i == nrow(subdf)) & sig) {
        if (i == nrow(subdf)) {
          end = subdf$pos[i]
        } else {
          end = subdf$pos[i-1]
        }
        sig = FALSE
        new_interval <- tibble(CHROM = chromosome, Start = start, End = end, Length = end - start)
        interval <- rbind(interval, new_interval)
      }
    }
  }

  # Peak information
  interval <- interval %>% mutate(PeakPos = numeric(nrow(interval)), PeakDeltaIndex = numeric(nrow(interval)))
  for (i in seq_along(interval$CHROM)) {
    peak <- df %>% filter(chr == interval$CHROM[i], pos >= interval$Start[i], pos <= interval$End[i]) %>% arrange(desc(index))
    # 如果QTL位点超过upperCI，也就是峰向上
    if (peak$index[1] > peak$CIupper[1]) {
      interval$PeakPos[i] <- peak %>% head(1) %>% pull(pos)
      interval$PeakDeltaIndex[i] <- peak %>% head(1) %>% pull(index)
    # 如果QTL位点超过lowerCI，也就是峰向下
    }else{
      interval$PeakPos[i] <- peak %>% tail(1) %>% pull(pos)
      interval$PeakDeltaIndex[i] <- peak %>% tail(1) %>% pull(index)
    }
  }
  if (export == TRUE) {
    write_csv(x = interval, file = filename)
  }
  return(interval)
}

#' @export
getQTL_for_ed4 <- function(data, chr = "CHROM", pos = "POS", nSNPs = "nSNPs", n = 10,
                          export = TRUE, filename = "SignificantQTL_for_ED4.csv"){
  ## 整理格式
  df <- data %>%
    select(
      chr = all_of(chr),
      pos = all_of(pos),
      ED4, ED4threshold,
      nSNPs = all_of(nSNPs)) %>%
    filter(nSNPs >= n)

  ## 初始化结果变量
  interval <- tibble(CHROM = character(0), Start = numeric(0), End = numeric(0), Length = numeric(0))
  ##
  for (chromosome in unique(df$chr)) {
    # 提取子集
    subdf <- df %>% filter(chr == chromosome)
    #subdf <- slidwin %>% filter(CHROM == "scaffoldA03")
    # 初始化变量
    start <- NULL
    end <- NULL
    sig <- FALSE
    for (i in seq_along(subdf$chr)) {
      if (subdf$ED4[i] > subdf$ED4threshold[i] & (!sig)) {
        start <- subdf$pos[i]
        sig <- TRUE
      }
      if ((subdf$ED4[i] < subdf$ED4threshold[i] | i == nrow(subdf)) & sig) {
        if (i == nrow(subdf)) {
          end = subdf$pos[i]
        } else {
          end = subdf$pos[i-1]
        }
        sig = FALSE
        new_interval <- tibble(CHROM = chromosome, Start = start, End = end, Length = end - start)
        interval <- rbind(interval, new_interval)
      }
    }
  }

  # Peak information
  interval <- interval %>% mutate(PeakPos = numeric(nrow(interval)), PeakED4 = numeric(nrow(interval)))
  for (i in seq_along(interval$CHROM)) {
    peak <- df %>% filter(chr == interval$CHROM[i], pos >= interval$Start[i], pos <= interval$End[i]) %>% arrange(desc(ED4))
    # 只有ED4大于阈值的情况
    interval$PeakPos[i] <- peak %>% head(1) %>% pull(pos)
    interval$PeakED4[i] <- peak %>% head(1) %>% pull(ED4)
  }
  if (export == TRUE) {
    write_csv(x = interval, file = filename)
  }
  return(interval)
}

#' @export
getQTL_for_Gprime <- function(data, chr = "CHROM", pos = "POS", nSNPs = "nSNPs", n = 10,
                           export = TRUE, filename = "SignificantQTL_for_Gprime.csv"){
  ## 整理格式
  df <- data %>%
    select(
      chr = all_of(chr),
      pos = all_of(pos),
      Gprime, p_for_Gprime, adj_p_for_Gprime,
      nSNPs = all_of(nSNPs)) %>%
    filter(nSNPs >= n) %>%
    mutate(Gprimethreshold = 6.634897)

  ## 初始化结果变量
  interval <- tibble(CHROM = character(0), Start = numeric(0), End = numeric(0), Length = numeric(0))
  ##
  for (chromosome in unique(df$chr)) {
    # 提取子集
    subdf <- df %>% filter(chr == chromosome)
    #subdf <- slidwin %>% filter(CHROM == "scaffoldA03")
    # 初始化变量
    start <- NULL
    end <- NULL
    sig <- FALSE
    for (i in seq_along(subdf$chr)) {
      if (subdf$Gprime[i] > subdf$Gprimethreshold[i] & (!sig)) {
        start <- subdf$pos[i]
        sig <- TRUE
      }
      if ((subdf$Gprime[i] < subdf$Gprimethreshold[i] | i == nrow(subdf)) & sig) {
        if (i == nrow(subdf)) {
          end = subdf$pos[i]
        } else {
          end = subdf$pos[i-1]
        }
        sig = FALSE
        new_interval <- tibble(CHROM = chromosome, Start = start, End = end, Length = end - start)
        interval <- rbind(interval, new_interval)
      }
    }
  }

  # Peak information
  interval <- interval %>% mutate(PeakPos = numeric(nrow(interval)), PeakGprime = numeric(nrow(interval)))
  for (i in seq_along(interval$CHROM)) {
    peak <- df %>% filter(chr == interval$CHROM[i], pos >= interval$Start[i], pos <= interval$End[i]) %>% arrange(desc(Gprime))
    # 只有Gprime大于阈值的情况
    interval$PeakPos[i] <- peak %>% head(1) %>% pull(pos)
    interval$PeakGprime[i] <- peak %>% head(1) %>% pull(Gprime)
  }
  if (export == TRUE) {
    write_csv(x = interval, file = filename)
  }
  return(interval)
}
