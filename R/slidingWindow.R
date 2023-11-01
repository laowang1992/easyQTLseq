#' @export
slidingWindow <- function(df, winSize, winStep, groups, position, values, fun){
  # winStep应小于等于winSize
  if (winStep>winSize) {
    stop("winStep should not be bigger than winSize!")
  }
  df <- df %>% dplyr::select(groups = all_of(groups), position = all_of(position), all_of(values))
  chrLen <- df %>% group_by(groups) %>% summarise(Len = max(position))
  # 生成windows
  wid <- tibble()
  for (i in 1:nrow(chrLen)) {
    if (chrLen$Len[i]<=winSize) {
      start <- 1
      end <- chrLen$Len[i]
      #end <- winSize
    }else{
      end <- seq(from = winSize, to = chrLen$Len[i], by = winStep)
      start <- end - winSize + 1
      if (end[length(end)] < chrLen$Len[i]) {
        start <- c(start, start[length(start)]+winStep)
        end <- c(end, chrLen$Len[i])
      }
      #start = seq(from = 1, to = chrLen$Len[i]-winSize+winStep+1, by = winStep)
      #end = seq(from = winSize, to = chrLen$Len[i]+winStep, by = winStep)
    }
    wid_tmp <- tibble(groups = chrLen$groups[i], Start = start, End = end)
    wid <- rbind(wid, wid_tmp)
  }
  # 添加统计列（N，sum）
  # 分窗口统计
  cat(date(), ", Sliding window statistics start ...\n", sep = "")
  # 定义一个进度条
  width <- options()$width
  pb <- progress::progress_bar$new(
    format = 'Progress [:bar] :percent eta: :eta',
    total = nrow(wid), clear = FALSE, width = width
  )
  for (i in 1:nrow(wid)) {
    df_tmp <- df %>% filter(groups == wid$groups[i], position >= wid$Start[i], position <= wid$End[i]) %>% na.omit()
    wid[i, "N"] <- nrow(df_tmp)
    x <- apply(df_tmp[, values], 2, fun)
    wid[i, paste(values, fun, sep = "_")] <- t(x[values])

    # 打印一个进度条
    ## 使用progress扩展包打印进度条
    pb$tick()
    ## 以下是一个手动打印进度条的代码
    #########################################################
    #cat('[', paste0(rep('#', i/nrow(wid)*width), collapse=''),
    #    paste0(rep('-', width - i/nrow(wid)*width), collapse=''),
    #    ']',
    #    round(i/nrow(wid)*100),'%')
    #if(i==nrow(wid))cat('\n', date(), ' DONE!\n', sep = "")
    #else cat('\r')
    ##########################################################
  }
  colnames(wid)[1:3] <- c(groups, "win_start", "win_end")
  return(wid)
}
