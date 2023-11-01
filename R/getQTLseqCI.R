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
