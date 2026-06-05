#' @export
getQTLseqCI <- function(df, popType, bulkSizeH, bulkSizeL, repN = 10000){
  dltIndex_CI <- df %>% mutate(CI95upper = NA, CI95lower = NA, CI99upper = NA, CI99lower = NA)
  cat(date(), ", start calculating confidence interval ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(dltIndex_CI), clear = FALSE, width = width
  )
  for (i in 1:nrow(dltIndex_CI)) {
    depthH <- dltIndex_CI[i, 1][[1]]
    depthL <- dltIndex_CI[i, 2][[1]]

    if (popType == "RIL") {
      PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
      PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
    } else if (popType == "F2") {
      PH <- apply(rmultinom(repN, bulkSizeH, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeH, 2, sum)
      PL <- apply(rmultinom(repN, bulkSizeL, c(1, 2, 1)) * c(1, 0.5, 0) / bulkSizeL, 2, sum)
    }

    indexH <- rbinom(repN, depthH, PH) / depthH
    indexL <- rbinom(repN, depthL, PL) / depthL
    dltIndex <- indexH - indexL

    dltIndex_CI[i, c("CI95upper", "CI95lower", "CI99upper", "CI99lower")] <- t(
      c(quantile(dltIndex, 0.975), quantile(dltIndex, 0.025),
        quantile(dltIndex, 0.995), quantile(dltIndex, 0.005))
    )

    # 打印进度条
    pb$tick()
  }
  cat(date(), ", finish calculating 95% and 99% confidence interval.\n", sep = "")
  return(dltIndex_CI)
}

#' @export
#getQTLseqCI_fast <- function(popType, bulkSizeH, bulkSizeL, minDepth = 5, maxDepth = 150, repN = 10000){
getQTLseqCI_fast <- function(df, popType, bulkSizeH, bulkSizeL, repN = 10000){
  #df <- tibble(HB.DP = ceiling(runif(n = 199, 5, 150)), LB.DP = ceiling(runif(n = 199, 10, 180))) %>% distinct() %>% arrange(HB.DP, LB.DP)
  # 1. 计算所有独立的深度值
  #depths <- minDepth:maxDepth
  #n_depths <- length(depths)
  #total_comb <- n_depths * n_depths
  depths_H <- df %>% pull(HB.DP) %>% unique()
  depths_L <- df %>% pull(LB.DP) %>% unique()
  n_depths_H <- length(depths_H)
  n_depths_L <- length(depths_L)
  total_comb <- nrow(df)

  cat(date(), " | Step 1: calculating allele frequency (PH & PL)...\n", sep = "")

  # 【优化一】将 PH 和 PL 移出循环，只计算一次！
  # 同时用高速的矩阵算术代替 F2 中的 apply(..., 2, sum)
  if (popType == "RIL") {
    PH <- rbinom(repN, bulkSizeH, 0.5) / bulkSizeH
    PL <- rbinom(repN, bulkSizeL, 0.5) / bulkSizeL
  } else if (popType == "F2") {
    mH <- rmultinom(repN, bulkSizeH, c(1, 2, 1))
    PH <- (mH[1, ] + 0.5 * mH[2, ]) / bulkSizeH

    mL <- rmultinom(repN, bulkSizeL, c(1, 2, 1))
    PL <- (mL[1, ] + 0.5 * mL[2, ]) / bulkSizeL
  } else {
    stop("Unknown population type，select 'RIL' or 'F2'")
  }

  cat(date(), " | Step 2: Precalculating SNP-index matrices across different depths...\n", sep = "")

  # 【优化二】预先生成所有独立深度的测序抽样矩阵 (行=repN, 列=独立深度数)
  # 这样 rbinom 的调用次数从 38416 次骤降到 196 次
  #indexH_mat <- matrix(NA, nrow = repN, ncol = n_depths)
  #indexL_mat <- matrix(NA, nrow = repN, ncol = n_depths)
  indexH_mat <- matrix(NA, nrow = repN, ncol = n_depths_H)
  indexL_mat <- matrix(NA, nrow = repN, ncol = n_depths_L)

  #for(j in 1:n_depths) {
  #  indexH_mat[, j] <- rbinom(repN, depths[j], PH) / depths[j]
  #  indexL_mat[, j] <- rbinom(repN, depths[j], PL) / depths[j]
  #}
  for(j in 1:n_depths_H) {
    indexH_mat[, j] <- rbinom(repN, depths_H[j], PH) / depths_H[j]
  }
  for(j in 1:n_depths_L) {
    indexL_mat[, j] <- rbinom(repN, depths_L[j], PL) / depths_L[j]
  }

  cat(date(), " | Step 3: calculating delta index and confidence interval...\n", sep = "")

  # 【优化三】利用矩阵索引和预分配内存，彻底告别 Tibble 按行写入
  #h_indices <- rep(1:n_depths, times = n_depths)
  #l_indices <- rep(1:n_depths, each = n_depths)
  h_indices <- df %>% mutate(pos = match(HB.DP, depths_H)) %>% pull(pos)
  l_indices <- df %>% mutate(pos = match(LB.DP, depths_L)) %>% pull(pos)


  # 创建原代码对应的基础网格
  #dltIndex_CI <- tibble(
  #  HB.DP = rep(depths, times = n_depths),
  #  LB.DP = rep(depths, each = n_depths)
  #)
  dltIndex_CI <- df

  # 预分配存储置信区间的矩阵
  res_cl <- matrix(NA, nrow = total_comb, ncol = 4)

  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = total_comb, clear = FALSE, width = options()$width
  )

  # 核心循环：现在内部只有极快的矩阵列减法和分位数计算
  for (i in 1:total_comb) {
    dltIndex <- indexH_mat[, h_indices[i]] - indexL_mat[, l_indices[i]]

    # names = FALSE 可以略微加速 quantile 函数
    res_cl[i, ] <- quantile(dltIndex, probs = c(0.975, 0.025, 0.995, 0.005), names = FALSE)
    pb$tick()
  }

  # 将计算结果快速合并回 tibble
  dltIndex_CI <- dltIndex_CI %>%
    mutate(
      CI95upper = res_cl[, 1],
      CI95lower = res_cl[, 2],
      CI99upper = res_cl[, 3],
      CI99lower = res_cl[, 4]
    )

  cat(date(), " | finish！\n", sep = "")
  return(dltIndex_CI)
}
