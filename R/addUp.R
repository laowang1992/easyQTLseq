#' @export
addUp <- function(df, len = NULL, group, pos, band = 0.01){
  df_tmp <- df
  df <- df %>% ungroup() %>% dplyr::select(group = all_of(group), all_of(pos))
  # 取最大值
  if (is.null(len)) {
    len <- df %>% mutate(max = apply(df[,-1], 1, max)) %>% group_by(group) %>% summarise(Len = max(max))
  }else{
    len <- len %>% dplyr::select(group = all_of(group), Len)
  }
  accu <- rep(0, rep(nrow(len)))
  for (i in seq_along(accu)[-1]) {
    accu[i] <- sum(len$Len[1:i-1]) + sum(len$Len)*band*(i-1)
  }
  names(accu) <- len$group
  breaks <- accu + len$Len/2
  gaps <- accu[-1] - sum(len$Len)*band/2
  labels <- names(accu)

  for (p in pos) {
    dd <- df %>% select(group, p = all_of(p)) %>% mutate(p_addUp = 0)
    for (c in names(accu)) {
      dd <- dd %>% mutate(p_addUp = if_else(group == c, p+accu[c], p_addUp))
    }
    colnames(dd) <- c("group", p, paste(p, "addUp", sep = "_"))
    df <- df %>% left_join(dd, by = c("group", p))
  }

  colnames(df)[1] <- group
  df <- df %>% left_join(df_tmp, by = c(group, pos))
  outList <- list(df = df, breaks = breaks, labels = labels, gaps = gaps)
  return(outList)
}
