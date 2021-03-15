
#' Draws a sunset plot showing the completeness of a dataset.
#'
#' @param mat              A \code{matrix}.
#' @param title            A \code{character} to specify the title of your plot.
#' @param col.pal          A \code{character} vector of length 3 matching R
#'                         colors to used as a palette for the plot.
#'                         Alternatively, col.pal can use the palettes provides
#'                         within the script of the function:
#'                         pal_sunset, pal_westworld, pal_startrek,
#'                         pal_margesimpson, pal_eighties, pal_insta\cr
#'                         (Default: col.pal = pal_sunset).
#' @param horizontal       A \code{logical} to specify whether the plot should
#'                         be drawn vertically or horizontally\cr
#'                         (Default: horizontal = FALSE).
#' @param reverse          A \code{logical} to specify whether the order
#'                         displayed of the bins should reversed\cr
#'                         (Default: reverse = FALSE).
#' @param keep_2nd_ticks   A \code{logical} to specify whether the ticks of the
#'                         secondary axis should be displayed or not\cr
#'                         (Default: keep_2nd_ticks = FALSE).
#' @param n.grad           An \code{integer} specifying the number of graduation
#'                         to use on the Y-Axis\cr(Default: n.grad = 15).
#' @param display.cutoff   A \code{double} specifying the minimum height of
#'                         a bin to label the number of samples on it. Below
#'                         this cutoff the number of samples is displayed in the
#'                         side margin\cr
#'                         (Default: display.cutoff = 0.03).
#' @param display.num.smpl A \code{double} specifying the minimum height of a
#'                         bin to display the number of samples in the side
#'                         margin\cr(Default: display.num.smpl = 0.01).
#'
#' @return A \code{gg} sunset plot.
#' @author Yoann Pageaud.
#' @export

sunset <- function(
  mat, title = "Number of rows without missing data",
  col.pal = c("red","gold","blue4"), horizontal = FALSE, reverse = FALSE,
  keep_2nd_ticks = FALSE, n.grad = 15, display.cutoff = 0.03,
  display.num.smpl = 0.01, lgd.pos = "bottom"){

  if(!is.matrix(mat)){ stop("'mat' must be a matrix.") }
  #Calculate number values not being NAs for each row in the matrix
  pos.cov <- as.integer(rowSums(!is.na(mat)))
  #Get total number of samples
  N <- dim(mat)[2]
  #Get total number of rows
  N.rows <- dim(mat)[1]
  #Tabulate data
  tab.data <- tabulate(pos.cov + 1L, N + 1L)

  #Compute size of bins
  dt.data <- data.table(
    "sample.amount" = c(0, seq(N)), "data.covered" = tab.data,
    "sample.string" = paste0(
      c("NONE", seq(N)), " (", round((tab.data/sum(tab.data))*100), "%", ")"))
  #Add index to dt.data
  dt.data[, index := .I]

  # x <- tabulate(pos.cov + 1L, N + 1L)
  # smpl_string <- paste0(
  #   c("NONE", seq(N)), " (", round((x/sum(x))*100), "%", ")")
  # Sample.Amt <- c(0, seq(N))

  #Reverse option
  if(reverse){
    dt.data[, c("sample.amount", "data.covered", "sample.string") := .(
      rev(sample.amount), rev(data.covered), rev(sample.string)) ]
    # x <- rev(x)
    # smpl_string <- rev(smpl_string)
    # Sample.Amt <- rev(Sample.Amt)
  }
  #Add label position, percentage, white_line position, difference bins &
  # cumulated sum
  dt.data[, c(
    "label.pos", "diff.bins", "cumulated", "right.y.data.covered",
    "right.y.sample.string", "right.y.label.pos","right.y.diff.bins") := .(
      cumsum(data.covered)-0.5*data.covered, data.covered/sum(data.covered),
      cumsum(data.covered), data.covered, sample.string,
      cumsum(data.covered)-0.5*data.covered, data.covered/sum(data.covered))]
  # #Create dataframes
  # dframe <- data.frame(Sample.Amt = Sample.Amt,
  #                      CpG.Covered = x,
  #                      label.pos = cumsum(x)-0.5*x,
  #                      percent = smpl_string,
  #                      white_lines = cumsum(x),
  #                      diff_bins = x/sum(x),
  #                      cumulated = cumsum(x),
  #                      stringsAsFactors = F)
  # df_right_y <- data.frame(CpG.Covered = x,
  #                          cumulated_smpl = cumsum(x)-0.5*x,
  #                          right_Y = smpl_string,
  #                          diff_bins = x/sum(x),
  #                          stringsAsFactors = F)

  #Apply display cut-off on main labels
  dt.data[data.covered < display.cutoff * N.rows, hiden := TRUE]
  #Apply display cut-off on right y labels
  dt.data[right.y.data.covered > display.cutoff * N.rows, right.y.hiden := TRUE]
  dt.data[right.y.diff.bins < display.num.smpl, right.y.hiden := TRUE]
  #Set white lines display
  dt.data[diff.bins < display.num.smpl, white.line.hide := TRUE]
  #Set to FALSE all NAs in hiden, right.y.hiden & white.line.hide
  dt.data[is.na(hiden), hiden := FALSE]
  dt.data[is.na(right.y.hiden), right.y.hiden := FALSE]
  dt.data[is.na(white.line.hide), white.line.hide := FALSE]
  #Hide last white line if one
  if(tail(dt.data$white.line.hide, n = 1) == FALSE){
    dt.data[max(index), white.line.hide := TRUE]
  }
  #Keep white line of level before the one having a white line
  invisible(lapply(X = dt.data[, .I], FUN = function(i){
    if(!all(is.na(dt.data[i+1,]))){
      if((dt.data[i+1, ]$hiden == FALSE | dt.data[i+1, ]$right.y.hiden == FALSE) &
         dt.data[i, ]$white.line.hide != FALSE){
        dt.data[i, white.line.hide := FALSE]
      }
    }
  }))

  # dframe[dframe$CpG.Covered < display.cutoff * max(dframe$CpG.Covered),
  # ]$percent <- " "
  # dframe[dframe$CpG.Covered/sum(x) < display.sep,]$white_lines <- NA
  # df_right_y[df_right_y$CpG.Covered > display.cutoff * max(dframe$CpG.Covered),
  # ]$right_Y <- " "
  # df_right_y[df_right_y$diff_bins < display.num.smpl,]$right_Y <- " "

  #Get positions of strings that will not be displayed
  pos_str <- data.table(dt.data[hiden == TRUE & right.y.hiden == TRUE,
                                which = TRUE])
  # pos_str <- match(smpl_string[-(sort(match(c(
  #   df_right_y$right_Y, dframe$percent), smpl_string)))], smpl_string)

  #Get groups of following samples
  smpl_intervals <- do.call(
    rbind, by(pos_str, cumsum(c(0, diff(pos_str$V1) != 1)), function(g){
      data.table(start = min(g$V1), end = max(g$V1),
                 width = diff(range(g$V1)) + 1)}))
  # smpl_intervals <- do.call(c, lapply(pos_str, function(i){
  #   if(i+1 - i == 1){ #If first sample
  #     if(!is.na(match(i+1, pos_str))){
  #       IRanges::IRanges(start = i, end = i+1)
  #     }
  #   }
  # }))
  # smpl_intervals <- IRanges::reduce(smpl_intervals)

  #Get groups average label.pos cumulative diff_bins
  if(nrow(smpl_intervals) != 0){
    invisible(lapply(smpl_intervals[, .I], function(i){
      dt_grp <- dt.data[c(smpl_intervals[i, ]$start:smpl_intervals[i, ]$end), ]
      # df_grp <- dframe[
      #   c(smpl_intervals[i, ]$start:smpl_intervals[i, ]$end), c(1:3, 6, 7)]
      data_cov_grp <- sum(dt_grp$data.covered)
      # cpg_cov_grp <- sum(df_grp$CpG.Covered)
      per_grp <- sum(dt_grp$diff.bins)
      # per_grp <- sum(df_grp$diff_bins)
      str_grp <- paste0(
        dt_grp[index == smpl_intervals[i, ]$start]$sample.amount, " to ",
        dt_grp[index == smpl_intervals[i, ]$end]$sample.amount, " (",
        round(per_grp*100), "%)")
      # str_grp <- paste0(
      #   df_grp[as.character(smpl_intervals[i, ]$start),]$Sample.Amt, " to ",
      #   df_grp[as.character(smpl_intervals[i, ]$end),]$Sample.Amt, " (",
      #   round(per_grp*100), "%)")
      if(dt_grp[1, ]$index == 1){
        pos_lab_grp <- tail(dt_grp, n = 1)$cumulated/2
        amnt_grp <- tail(dt_grp$sample.amount, n = 1)
      } else {
        pos_lab_grp <- (dt.data[dt_grp[1, ]$index - 1,]$cumulated +
                          tail(dt_grp, n = 1)$cumulated)/2
        amnt_grp <- tail(dt_grp$sample.amount, n = 1) - head(
          dt_grp$sample.amount, n = 1) + 1
      }
      if(data_cov_grp > display.cutoff * N.rows){
        hiden_grp <- FALSE
        right.y.hiden_grp <- TRUE
      } else {
        hiden_grp <- TRUE
        if(per_grp < display.num.smpl){ right.y.hiden_grp <- TRUE
        } else { right.y.hiden_grp <- TRUE }
      }
      new.dt <- data.table(
        "sample.amount" = amnt_grp, "data.covered" = data_cov_grp,
        "sample.string" = str_grp, "index" = tail(dt_grp$index, n = 1),
        "label.pos" = pos_lab_grp, "diff.bins" = per_grp,
        "cumulated" = tail(dt_grp$cumulated, n = 1),
        "right.y.data.covered" = data_cov_grp,
        "right.y.sample.string" = str_grp, "right.y.label.pos" = pos_lab_grp,
        "right.y.diff.bins" = per_grp, "hiden" = hiden_grp,
        "right.y.hiden" = right.y.hiden_grp)
      # if(per_grp > display.num.smpl){
      #   new.df <- data.frame(
      #     "CpG.Covered" = cpg_cov_grp, "cumulated_smpl" = pos_lab_grp,
      #     "right_Y" = str_grp, "diff_bins" = per_grp)
      # }
      # return(new.df)
      dt.data[index == new.dt$index, c(
        "sample.string", "label.pos", "right.y.sample.string",
        "right.y.label.pos", "hiden", "right.y.hiden") := .(
          new.dt$sample.string, new.dt$label.pos, new.dt$sample.string,
          new.dt$label.pos, new.dt$hiden, new.dt$right.y.hiden)]
    }))
  }
  #Rbind all new data.tables
  # new.dts <- rbindlist(l = ls.new.df)
  #
  #   #Rbind dt.data with new.dts
  #   dt.data <- rbind(dt.data, new.dts)
  #   #Order dt.data following index
  #   dt.data <- dt.data[order(index)]

  # if(length(smpl_intervals) != 0){
  #   invisible(lapply(seq_along(smpl_intervals), function(i){
  #     df_grp<-dframe[c(IRanges::start(smpl_intervals[i]):IRanges::end(smpl_intervals[i])),
  #                    c(1:3,6,7)]
  #     cpg_cov_grp <- sum(df_grp$CpG.Covered)
  #     per_grp<-sum(df_grp$diff_bins)
  #     str_grp<-paste0(df_grp[as.character(start(smpl_intervals[i])),1]," to ",
  #                     df_grp[as.character(end(smpl_intervals[i])),1]," (",
  #                     round(per_grp*100),"%)")
  #     if(rownames(head(df_grp,n=1L))=="1"){
  #       pos_lab_grp<-tail(df_grp$cumulated,n=1L)/2
  #     } else {
  #       pos_lab_grp<-(dframe[as.integer(rownames(head(df_grp,n=1L)))-1,
  #       ]$cumulated + tail(df_grp$cumulated,n=1L))/2
  #     }
  #     if(per_grp > display.num.smpl){
  #       df_right_y<<-rbind(df_right_y,data.frame(
  #         "CpG.Covered" = cpg_cov_grp, "cumulated_smpl" = pos_lab_grp,
  #         "right_Y" = str_grp, "diff_bins" = per_grp))
  #     }
  #   }))
  # }

  # if(keep_2nd_ticks == FALSE) {
  #   df_right_y <- df_right_y[df_right_y$right_Y != " ", ]
  # }

  #Vector for white lines
  # white_lines <- dframe[!is.na(dframe$white_lines), ]$white_lines

  #"Sunset" Plot of the Amount of CpGs Covered by Number of Samples
  sun.plt <- ggplot() + theme_gray() +
    theme(legend.title.align = 0.5,
          legend.text = element_text(size = 12),
          legend.position = lgd.pos) +
    geom_bar(data = dt.data,
             mapping = aes(x = 0, y = data.covered, fill = sample.amount),
             stat = "identity") +
    geom_text(data = dt.data[hiden == FALSE], mapping = aes(
      x = 0, y = label.pos, label = sample.string), vjust = 0.35, hjust = 0.5,
      color = "white", size = 5) +
    scale_fill_gradient2(low = col.pal[1], mid = col.pal[2], high = col.pal[3],
                         midpoint = round(N/2)) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq(0, N.rows, length.out = n.grad),
      labels = function(x) format(x, digits = 2, scientific = TRUE),
      sec.axis = sec_axis(
        trans = ~.,
        breaks = dt.data[right.y.hiden == FALSE]$right.y.label.pos,
        labels = dt.data[right.y.hiden == FALSE]$sample.string)) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_hline(data = dt.data[white.line.hide == FALSE],
               mapping = aes(yintercept = cumulated), color = "white")

  # Sunset <- ggplot(data = dframe,
  #                  aes(x = 0, y = CpG.Covered, fill = Sample.Amt)) +
  #   theme_gray() +
  #   theme(legend.title.align=0.5,
  #         legend.text=element_text(size=12),
  #         legend.position=lgd.pos) +
  #   geom_bar(stat = "identity") +
  #   geom_text(aes(y = label.pos, label = percent), vjust = 0.25, hjust = 0.5,
  #             colour = "white", size = 5) +
  #   scale_fill_gradient2(low = col.pal[1], mid = col.pal[2], high = col.pal[3],
  #                        midpoint = round(N/2)) +
  #   scale_y_continuous(
  #     expand = c(0,0),
  #     breaks = seq(0, sum(dframe$CpG.Covered), length.out = n.grad),
  #     labels = function(x) format(x, digits = 2, scientific = TRUE),
  #     sec.axis = sec_axis(trans = ~., breaks = df_right_y$cumulated_smpl,
  #                         labels = df_right_y$right_Y)) +
  #   scale_x_continuous(expand = c(0, 0)) +
  #   geom_hline(yintercept = white_lines, color = "white")

  if(horizontal) {
    sun.plt <- sun.plt +
      # Sunset <- Sunset +
      ggtitle(title) +
      theme(plot.title = element_text(title, size = 15, hjust = 0.5),
            axis.title = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 14, vjust = 0.8),
            plot.margin = margin(0, 0, 0, 20)) +
      labs(fill = "Number of Samples") +
      coord_flip()
  } else {
    sun.plt <- sun.plt +
      # Sunset <- Sunset +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 15),
            legend.title = element_text(size = 14, vjust = 1),
            legend.title.align = 0) +
      ylab(title) +
      labs(fill = "Number\nof samples")
  }
  sun.plt
}
