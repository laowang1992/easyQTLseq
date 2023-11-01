#' Filter SNP according depth
#'
#' For low coverage depth SNP may have low reliability and accuracy, and extremely high coverage depth may be derived from repetitive sequence.
#' This function filter low coverage depth and extremely high coverage depth SNP.
#'
#' @param x A QTLseq S3 object.
#' @param minHPdp Minimum coverage depth of high phenotype parent.
#' @param minLPdp Minimum coverage depth of low phenotype parent.
#' @param minHBdp Minimum coverage depth of high phenotype bulk.
#' @param minLBdp Minimum coverage depth of low phenotype bulk.
#' @param maxHPdp Maximum coverage depth of high phenotype parent, default if average+3*sd.
#' @param maxLPdp Maximum coverage depth of low phenotype parent, default if average+3*sd.
#' @param maxHBdp Maximum coverage depth of high phenotype bulk, default if average+3*sd.
#' @param maxLBdp Maximum coverage depth of low phenotype bulk, default if average+3*sd.
#'
#' @return
#' filtered QTLseq S3 object.
#' @export
#'
#' @examples
#' library(easyQTLseq)
#' # Example with sample data from a GATK table.
#' file_path <- system.file("extdata", "subset.table.gz", package = "easyQTLseq")
#' # readr::read_tsv() has a faster speed than read.table() when reading a file.
#' data <- readr::read_tsv(file = file_path)
#' x <- importData(data = data, highP = "qY", lowP = "R3", highB = "Y", lowB = "R", popType = "F2", bulkSize = c(30, 30))
#' x_filter <- filterDP(x = x)
filterDP <- function(x, ...) {
  UseMethod("filterDP")
}

#' @rdname filterDP
#' @export
filterDP.BothParent <- function(x, minHPdp = 6, minLPdp = 6, minHBdp = 6, minLBdp = 6, maxHPdp, maxLPdp, maxHBdp, maxLBdp){
  # 根据深度过滤，上限由以前的自定义改为设置成ave+3*sd
  dp <- x$data %>% select(HP = any_of("HP.DP"), LP = any_of("LP.DP"), HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  stat <- dp %>% group_by(sample) %>%
    summarise(ave = mean(depth), sd = sd(depth), upper = ave+3*sd, outer = sum(depth>upper), `outerRate` = paste(round(x = outer/n()*100, digits = 3), "%"))
  print(stat)
  if(missing(maxHPdp)) maxHPdp <- ceiling(stat$upper[stat$sample %in% "HP"])
  if(missing(maxLPdp)) maxLPdp <- ceiling(stat$upper[stat$sample %in% "LP"])
  if(missing(maxHBdp)) maxHBdp <- ceiling(stat$upper[stat$sample %in% "HB"])
  if(missing(maxLBdp)) maxLBdp <- ceiling(stat$upper[stat$sample %in% "LB"])
  # 根据深度过滤
  if (x$highP == "REF") {
    x$data <- x$data %>%
      filter(LP.DP >= minLPdp, LP.DP <= maxLPdp,
             HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  } else if (x$lowP == "REF") {
    x$data <- x$data %>%
      filter(HP.DP >= minHPdp, HP.DP <= maxHPdp,
             HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  } else {
    x$data <- x$data %>%
      filter(HP.DP >= minHPdp, HP.DP <= maxHPdp,
             LP.DP >= minLPdp, LP.DP <= maxLPdp,
             HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  }
  return(x)
}

#' @rdname filterDP
#' @export
filterDP.HighParent <- function(x, minHPdp = 6, minHBdp = 6, minLBdp = 6, maxHPdp, maxHBdp, maxLBdp){
  # 根据深度过滤，上限由以前的自定义改为设置成ave+3*sd
  dp <- x$data %>% select(HP = any_of("HP.DP"), HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  stat <- dp %>% group_by(sample) %>%
    summarise(ave = mean(depth), sd = sd(depth), upper = ave+3*sd, outer = sum(depth>upper), `outerRate` = paste(round(x = outer/n()*100, digits = 3), "%"))
  print(stat)
  if(missing(maxHPdp)) maxHPdp <- ceiling(stat$upper[stat$sample %in% "HP"])
  if(missing(maxHBdp)) maxHBdp <- ceiling(stat$upper[stat$sample %in% "HB"])
  if(missing(maxLBdp)) maxLBdp <- ceiling(stat$upper[stat$sample %in% "LB"])
  # 根据深度过滤
  if (x$highP == "REF") {
    x$data <- x$data %>%
      filter(HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  } else {
    x$data <- x$data %>%
      filter(HP.DP >= minHPdp, HP.DP <= maxHPdp,
             HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  }
  return(x)
}

#' @rdname filterDP
#' @export
filterDP.LowParent <- function(x, minLPdp = 6, minHBdp = 6, minLBdp = 6, maxLPdp, maxHBdp, maxLBdp){
  # 根据深度过滤，上限由以前的自定义改为设置成ave+3*sd
  dp <- x$data %>% select(LP = any_of("LP.DP"), HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  stat <- dp %>% group_by(sample) %>%
    summarise(ave = mean(depth), sd = sd(depth), upper = ave+3*sd, outer = sum(depth>upper), `outerRate` = paste(round(x = outer/n()*100, digits = 3), "%"))
  print(stat)
  if(missing(maxLPdp)) maxLPdp <- ceiling(stat$upper[stat$sample %in% "LP"])
  if(missing(maxHBdp)) maxHBdp <- ceiling(stat$upper[stat$sample %in% "HB"])
  if(missing(maxLBdp)) maxLBdp <- ceiling(stat$upper[stat$sample %in% "LB"])
  # 根据深度过滤
  if (x$lowP == "REF") {
    x$data <- x$data %>%
      filter(HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  } else {
    x$data <- x$data %>%
      filter(LP.DP >= minHPdp, LP.DP <= maxHPdp,
             HB.DP >= minHBdp, HB.DP <= maxHBdp,
             LB.DP >= minLBdp, LB.DP <= maxLBdp)
  }
  return(x)
}

#' @rdname filterDP
#' @export
filterDP.WithoutParent <- function(x, minHBdp = 6, minLBdp = 6, maxHBdp, maxLBdp){
  # 根据深度过滤，上限由以前的自定义改为设置成ave+3*sd
  dp <- x$data %>% select(HB = HB.DP, LB = LB.DP) %>%
    gather(key = "sample", value = "depth")
  stat <- dp %>% group_by(sample) %>%
    summarise(ave = mean(depth), sd = sd(depth), upper = ave+3*sd, outer = sum(depth>upper), `outerRate` = paste(round(x = outer/n()*100, digits = 3), "%"))
  print(stat)
  if(missing(maxHBdp)) maxHBdp <- ceiling(stat$upper[stat$sample %in% "HB"])
  if(missing(maxLBdp)) maxLBdp <- ceiling(stat$upper[stat$sample %in% "LB"])
  # 根据深度过滤
  x$data <- x$data %>%
    filter(HB.DP >= minHBdp, HB.DP <= maxHBdp,
           LB.DP >= minLBdp, LB.DP <= maxLBdp)
  return(x)
}
