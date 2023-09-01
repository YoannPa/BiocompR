
#' Draws a craviola plot (half-splitted and percentile-binned violin plot).
#'
#' @param data           A \code{data.frame}. 2 formats of data.frame are
#'                       supported:
#' \itemize{
#'  \item{A "complete" data.frame of 4 columns:
#'   \itemize{
#'    \item{column 1 must be the "group" column. There must not be
#'          more than 2 conditions in each group.}
#'    \item{column 2 must be the "condition" column.}
#'    \item{column 3 must be the "values" column.}
#'    \item{column 4 is optional. It should be another "values" column to be
#'          used if bins = TRUE. When 'bins' is TRUE, the function specified in
#'          'bin.fun' is applied on the column 4 for each bin of each
#'          distribution. The result obtain for each bin is then associated with
#'          a specific shade of the filling colors.}
#'   }
#'  }
#'  \item{A "minimal" data.frame of 2 columns:
#'   \itemize{
#'    \item{column 1 must be the "samples" column. This column can take a
#'    maximum of 2 possible levels. Based on this column ggcraviola() will
#'    extrapolate the "grouping" and the "filling".}
#'    \item{column 2 must be the "values" column.}
#'   }
#'   The "minimal" data.frame is useful to easily plot a single craviola to
#'   compare 2 distributions. You cannot plot more than 1 craviola with the
#'   "minimal" data.frame format. the "minimal" format also do not support the
#'   binning.
#'  }
#' }
#' @param craviola.width A \code{double} value to specify the width of
#'                       craviolas\cr(Default: craviola.width = 1).
#' @param boxplots       A \code{logical} specifying if boxplots should be
#'                       displayed or not\cr(Default: boxplots = TRUE).
#' @param boxplot.width  A \code{double} specifying the width of boxplots\cr
#'                       (Default: boxplot.width = 0.04).
#' @param mean.value     A \code{logical} specifying if the red dot of the mean
#'                       value should be displayed or not\cr
#'                       (Default: mean.value = TRUE).
#' @param bins           A \code{logical} specifying if the distributions should
#'                       be binned following specific quantiles to be displayed
#'                       following different opacity thanks to the values in the
#'                       5th column of the data data.frame\cr
#'                       (Default: bins = FALSE).
#' @param bins.quantiles A \code{double} vector to define the limits between the
#'                       bins as percentiles of distributions\cr
#'                       (Default: bins.quantiles = seq(0.1,0.9,0.1)).
#' @param bin.fun        A \code{character} to specify the function to apply on
#'                       values stored in the column 4 of the complete
#'                       data.frame, for each bin of each distribution.\cr
#'                       Examples of functions supported: bin.fun = c(
#'                       "stats::sd", "stats::mad", "base::mean",
#'                       "stats::median"); Default: bin.fun = "stats::sd".
#' @param lines.col      A \code{character} matching a color to use for the
#'                       border lines of craviolas' bins\cr
#'                       (Default: lines.col = NA).
#' @param verbose        A \code{logical} to display information about the
#'                       step-by-step processing of the data if TRUE
#'                       (Default: verbose = FALSE).
#' @return A \code{gg} craviola plot.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=`
#' @export
#' @examples
#' # Using a 'minimal' data.frame:
#' df.minimal = data.frame(
#'   Samples = rep(paste0("Sample", c(1:2)), each = 1000),
#'   Values = c(rnorm(1000, 0), rnorm(1000, 0.5)))
#' ggcraviola(data = df.minimal, lines.col = "black")
#'
#' # Using a 'complete' data.frame:
#' df.complete = data.frame(
#'     Groups = rep(c('A', 'B', 'C'), each = 2000),
#'     Conditions = rep(c('I', 'J'), each = 1000,3),
#'     Values = c(rnorm(1000, 0), rnorm(1000, 0.5),
#'                rnorm(1000, 3), rnorm(1000, 3.5),
#'                rnorm(1000, -3), rnorm(1000, -3.5)))
#' ggcraviola(data = df.complete, lines.col = "black")
#'
#' # Using a 'complete' data.frame with support of an additional variable:
#' df.complete2 = data.frame(
#'     Groups = rep(c('A', 'B', 'C'), each = 2000),
#'     Conditions = rep(c('I', 'J'), each = 1000, 3),
#'     Values = c(rnorm(1000, 0), rnorm(1000, 0.5),
#'                rnorm(1000, 3), rnorm(1000, 3.5),
#'                rnorm(1000, -3), rnorm(1000, -3.5)),
#'     Scnd.Var = rep(rep(x = c(60, 50, 40, 30, 20, 10, 30, 40, 50, 60),
#'                        each = 100), 6))
#' ggcraviola(data = df.complete2, lines.col = "black", bins = TRUE)
#'
#' #Use ggplot2 to add components and customize the Craviola plot
#' ggcraviola(data = df.complete2, lines.col = "black", bins = TRUE) +
#' labs(alpha = "SD(2nd Variable)") + # Rename alpha legend
#' ggtitle("This is a Craviola plot!") + # Add title
#' theme(plot.title = element_text(hjust = 0.5),
#'       axis.text = element_text(size = 14, color = "black"),# Custom axis text
#'       axis.title = element_text(size = 15),
#'       legend.title = element_text(size = 13), # Change legend font size
#'       legend.text = element_text(size = 12),
#'       panel.background = element_blank(), # Change panel appearance
#'       panel.grid.major.x = element_blank(),
#'       panel.grid.minor.x = element_blank(),
#'       panel.grid.major.y = element_line(color = "grey"),
#'       panel.grid.minor.y = element_line(color = "grey")) +
#' scale_y_continuous(expand = c(0, 0)) + #Expand fully plot panel on Y-axis
#' scale_fill_manual(
#'   labels = c("Control", "Case"), # Rename conditions
#'   values = biopalette(name = "BiocompR_cond3", mute = TRUE)) # Change colors

ggcraviola <- function(
  data, craviola.width = 1, boxplots = TRUE, boxplot.width = 0.04,
  mean.value = TRUE, bins = FALSE, bins.quantiles = seq(0.1, 0.9, 0.1),
  bin.fun = "stats::sd", lines.col = NA, linewidth = 0.5, verbose = FALSE){
  # Fix BiocCheck() complaining about these objects initialization

  # Samples <- NULL
  # dens.curv <- NULL
  # y.pos <- NULL
  # Var.col <- NULL
  # Var1 <- NULL
  # bin <- NULL
  # bin.av.val2 <- NULL
  # x <- NULL
  # low <- NULL
  # mid <- NULL
  # top <- NULL

  if(verbose){ cat("Formatting data and annotation table...") }
  #Check if data is a data.table and convert if not
  if(!data.table::is.data.table(data)){
    data <- data.table::as.data.table(data)
  } else {
    dt.new <- data.table::copy(data) # Make to avoid modifying original object
    data <- dt.new
  }
  # Rename columns
  if(ncol(data) > 3){
    copydt <- data.table::copy(data[NA])
    old_cols <- colnames(copydt)
    data.table::setnames(x = data, old = 1, new = "groups")
    data.table::setnames(x = data, old = 2, new = "conds")
    data.table::setnames(x = data, old = 3, new = "vals")
    data.table::setnames(x = data, old = 4, new = "opac")
    data[, samples := interaction(groups, conds)]
    data[, c("groups", "conds") := .(
      as.factor(groups), as.factor(conds))]
  } else if(ncol(data) == 3){
    copydt <- data.table::copy(data[NA])
    old_cols <- colnames(copydt)
    data.table::setnames(x = data, old = 1, new = "groups")
    data.table::setnames(x = data, old = 2, new = "conds")
    data.table::setnames(x = data, old = 3, new = "vals")
    data[, samples := interaction(groups, conds)]
    data[, c("groups", "conds") := .(as.factor(groups), as.factor(conds))]
  } else if(ncol(data) == 2){
    copydt <- data.table::copy(data[NA])
    old_cols <- c(colnames(copydt), "Groups", "Conditions")
    data.table::setnames(x = data, old = 1, new = "samples")
    data.table::setnames(x = data, old = 2, new = "vals")
    data[, samples := as.factor(samples)]
    if(length(levels(data$samples)) != 2){
      stop("there should be only 2 possible values in column 1.")
    } else {
      data[, groups := as.factor(1)]
      data[, conds := samples]
    }
  } else { stop("not enough columns in data.") }
  # print(old_cols)
  if(length(levels(data$conds)) > 2){
    stop("more than 2 conditions detected. Only 2 conditions tolerated.")
  }
  # Make annotation table
  annot_table <- unique(x = data, by = "samples")
  annot_table <- annot_table[, -c("vals"), ]
  if("opac" %in% colnames(annot_table)){
    annot_table <- annot_table[, -c("opac"), ]
  }
  if(any(annot_table[, data.table::.N, by = "groups"]$N > 2)){
    stop(paste(
      "At least 1 group contains data for more than 2 samples.",
      "Cannot represent more than 2 distributions per group."))
  }
  # unique.dt <- unique(x = data, by = 1)
  if(verbose){ cat("Done.\n") }

  # if(ncol(data) == 2){
  #   if(nrow(annot_table) == 2) {
  #     # colnames(data)[1] <- "Samples"
  #     # colnames(unique.dt)[1] <- "Samples"
  #     # if(!is.factor(unique.dt[[1]])){ #Convert as factor column 1
  #     #   unique.dt[, Samples := as.factor(Samples)]
  #     # }
  #     annot_table[, groups := ]
  #     Annot.table <- data.table::data.table(
  #       unique.dt$samples, "Groups" = as.factor(c(1, 1)),
  #       "Conditions" = unique.dt[[1]], unique.dt[, 2])
  #     data <- merge(
  #       x = Annot.table[, -4], y = data, by = "Samples", all.y = TRUE)
  #   } else{ stop("Missing columns in the data provided.") }
  # } else {
  #   Annot.table <- unique.dt[, seq(3), with = FALSE]
  #   data.table::setnames(
  #     x = Annot.table, old = colnames(Annot.table)[2], new = "Groups")
  #   data.table::setnames(
  #     x = Annot.table, old = colnames(Annot.table)[1], new = "Samples")
  # }
  # if(nrow(unique(Annot.table, by = 3)) > 2){
  #   stop("More than 2 conditions inputed. Only 2 conditions tolerated.")
  # }
  # if(any(Annot.table[, data.table::.N, by = "Groups"]$N > 2)){
  #   stop(paste(
  #     "At least 1 group contains data for more than 2 samples.",
  #     "A craviola cannot represent more than 2 distributions per group."))
  # }
  # #Check if all conditional variables are factors
  # if(!all(Annot.table[, lapply(X = .SD, FUN = is.factor)] == TRUE)){
  #   #Convert annotations and data conditional variable into factors
  #   if(verbose){ cat("Converting conditional variables as factors...") }
  #   Annot.table <- Annot.table[, lapply(X = .SD, FUN = as.factor)]
  #   cols <- colnames(data)[seq(3)]
  #   data[, (cols) := lapply(
  #     X = .SD, FUN = as.factor), .SDcols = cols]
  #   if(verbose){ cat("Done.\n") }
  # }
  # original.var.col <- levels(Annot.table[[3]])
  # if(length(levels(Annot.table[[3]])) > 2) {
  #   stop("More levels than possible values. Only 2 conditions tolerated. Remove the excess levels.")
  # } else {
  #   levels(Annot.table[[3]])[1] <- "1"
  #   levels(Annot.table[[3]])[2] <- "2"
  # }
  # if (length(Annot.table[[1]]) < length(levels(Annot.table[[1]]))){
  #   stop("More levels than matching values found in column 1 of data.")
  # }
  # amount.grp <- length(unique(Annot.table[[2]]))
  # if(amount.grp > 1){
  # data.table::setattr(
  #   x = Annot.table[[2]], name = "levels", value = as.character(
  #     seq_along(levels(Annot.table[[2]]))))
  # }
  if(verbose){ cat("Splitting dataset on samples...") }
  # mylist_data <- split(x = data, f = data[[1]])
  mylist_data <- split(x = data, f = data$samples)
  # list_val1 <- lapply(X = mylist_data, FUN = subset, select = 4)
  list_val1 <- lapply(X = mylist_data, FUN = subset, select = "vals")
  rm(mylist_data)
  list_vect.val1 <- lapply(X = list_val1, FUN = unlist, use.names = FALSE)
  rm(list_val1)
  if(verbose){ cat("Done.\n") }
  # Create stats plots
  if(verbose){ cat("Compute boxplot stats...") }
  list.bp.stat <- lapply(X = seq_along(list_vect.val1), FUN = function(i){
    qiles <- stats::quantile(list_vect.val1[[i]])
    means <- mean(list_vect.val1[[i]])
    # if(Annot.table[Annot.table[[1]] == names(list_vect.val1)[i], 3] == 1){
    if(annot_table[
      samples == names(list_vect.val1)[i]
    ]$conds == levels(annot_table$conds)[1]){
      x.pos <- -boxplot.width
    } else { x.pos <- boxplot.width }
    # if(Annot.table[Annot.table[[1]] == names(list_vect.val1)[i], 2] != 1){
    if(annot_table[
      samples == names(list_vect.val1)[i]
    ]$groups != levels(annot_table$groups)[1]){
      # x.pos <- x.pos + (as.integer(Annot.table[
      #   Annot.table[[1]] == names(list_vect.val1)[i], 2])-1)
      x.pos <- x.pos + (as.integer(annot_table[
        samples == names(list_vect.val1)[i]]$groups) - 1)
    }
    # data.table::data.table(
    #   Annot.table[Annot.table[[1]] == names(list_vect.val1)[i], 3], "x" = x.pos,
    #   "min" = qiles[1], "low" = qiles[2], "mid" = qiles[3], "top" = qiles[4],
    #   "max" = qiles[5], "mean" = means, "pos.crav" = round(x.pos))
    data.table::data.table(
      annot_table[samples == names(list_vect.val1)[i], conds], "x" = x.pos,
      "min" = qiles[1], "low" = qiles[2], "mid" = qiles[3], "top" = qiles[4],
      "max" = qiles[5], "mean" = means, "pos.crav" = round(x.pos))
  })
  box.dframe <- data.table::rbindlist(list.bp.stat)
  if(verbose){ cat("Done.\n") }
  #Create Craviola plot
  if(verbose){ cat("Compute density on each distribution...") }
  # density.scaler <- mean(abs(data[[4]]), na.rm = TRUE)/2 # Deflt craviola scale
  density.scaler <- mean(abs(data$vals), na.rm = TRUE)/2 # Deflt craviola scale
  list_dens.res <- lapply(X = list_vect.val1, FUN = stats::density)
  list_dens.df <- lapply(X = list_dens.res, FUN = function(i){
    # data.frame("y.pos" = i$x, "dens.curv" = i$y*density.scaler*craviola.width)
    data.table::data.table(
      "y.pos" = i$x, "dens.curv" = i$y*density.scaler*craviola.width)
  })
  if(verbose){ cat("Done.\n") }

  list_oriented_dens <- lapply(X = names(list_dens.df), FUN = function(i){
    dens.df <- list_dens.df[[i]]
    if(annot_table[samples == i]$conds == levels(annot_table$conds)[1]){
      dens.df[, dens.curv := dens.curv * -1]
    }
    if(annot_table[samples == i]$groups != levels(annot_table$groups)[1]){
      dens.df[, dens.curv := dens.curv + as.integer(
        annot_table[samples == i]$groups) - 1]
    }
    dens.df
  })
  # list_oriented_dens <- lapply(names(list_dens.df), function(i){
  #   if(as.integer(Annot.table[Annot.table[[1]] == i, 3]) == 1){
  #     base::`<<-` (
  #       list_dens.df[[i]]$dens.curv, list_dens.df[[i]]$dens.curv * -1)
  #     list_dens.df[[i]]
  #   } else { list_dens.df[[i]] }
  #   if(as.integer(Annot.table[Annot.table[[1]] == i, 2]) > 1){
  #     base::`<<-` (list_dens.df[[i]]$dens.curv, list_dens.df[[i]]$dens.curv +
  #                    as.integer(Annot.table[Annot.table[[1]] == i, 2]) - 1)
  #     list_dens.df[[i]]
  #   } else { list_dens.df[[i]] }
  # })

  # Remove density values outside the extrema
  if(verbose){ cat("Cutting density at the extrema...") }
  xtrems <- BiocompR:::ls.quantile(ls = list_vect.val1, qtiles = c(0, 1))
  bined.xtrm.dens <- BiocompR:::bin.polygons(
    list_oriented_dens = list_oriented_dens, list.quant.lim = xtrems,
    annot_table = annot_table)
  # bined.xtrm.dens <- BiocompR:::bin.polygons(
  #   list_oriented_dens = list_oriented_dens, list.quant.lim = xtrems,
  #   Annot.table = Annot.table)
  #Keep only bin 1 for each sample
  list_oriented_dens <- lapply(
    X = seq_along(bined.xtrm.dens), FUN = function(i){
    # df <- bined.xtrm.dens[[i]][bined.xtrm.dens[[i]]$bin == 1, c(2, 3)]
    df <- bined.xtrm.dens[[i]][bin == 1, c("y.pos", "dens.curv")]
    # rownames(df) <- NULL
    df
  })
  names(list_oriented_dens) <- names(list_dens.df)
  if(verbose){ cat("Done.\n") }

  # # Reorder Annotation table following order of dataframes
  # Annot.table <- Annot.table[order(match(Samples, names(list_oriented_dens)))]

  # Create Bins based on a third variable
  if(bins){ #Create bin polygons
    if(verbose){ cat("Cutting density by quantiles...") }
    #Bin polygons
    list.quant.lim <- BiocompR:::ls.quantile(
      ls = list_vect.val1, qtiles = bins.quantiles)
    list.dfs <- BiocompR:::bin.polygons(
      list_oriented_dens = list_oriented_dens, list.quant.lim = list.quant.lim,
      annot_table = annot_table)
    names(list.dfs) <- names(list.quant.lim)

    # Compute bin.fun on the 3rd variable for each bin
    list.fun.val2 <- lapply(X = seq_along(list.quant.lim), FUN = function(i){
      smpl.data <- data[samples == names(list.quant.lim)[i], ]
      # Check if external quantiles are min and max
      new_quant.lim <- list.quant.lim[[i]]
      if(min(list_vect.val1[[i]]) != new_quant.lim[1]){
        # Add minimum value at the beginning of the vector
        new_quant.lim <- c(min(list_vect.val1[[i]]), new_quant.lim)
      }
      if(max(list_vect.val1[[i]]) != rev(new_quant.lim)[1]){
        # Add maximum value at the end of the vector
        new_quant.lim <- c(new_quant.lim, max(list_vect.val1[[i]]))
      }
      smpl.data[, bin.groups := findInterval(
        vals, new_quant.lim, all.inside = TRUE) - 1]
      unlist(lapply(X = sort(unique(smpl.data$bin.groups)), FUN = function(j){
        eval(parse(text = paste0(
          bin.fun,"(smpl.data[smpl.data$bin.groups == j]$opac, na.rm = TRUE)")))
      }))
    })
    # list.fun.val2 <- lapply(X = seq_along(list.quant.lim), FUN = function(i){
    #   smpl.data <- data[data[[1]] == names(list.quant.lim)[i], ]
    #   # Check if external quantiles are min and max
    #   if(min(list_vect.val1[[i]]) != list.quant.lim[[i]][1]){
    #     # Add minimum value at the beginning of the vector
    #     base::`<<-` (
    #       list.quant.lim[[i]], c(min(list_vect.val1[[i]]), list.quant.lim[[i]]))
    #   }
    #   if(max(list_vect.val1[[i]]) != rev(list.quant.lim[[i]])[1]){
    #     # Add maximum value at the end of the vector
    #     base::`<<-` (
    #       list.quant.lim[[i]], c(list.quant.lim[[i]], max(list_vect.val1[[i]])))
    #   }
    #   smpl.data[["bin.groups"]] <- findInterval(
    #     smpl.data[[4]], list.quant.lim[[i]], all.inside = TRUE)-1
    #   if(bin.fun == "mean"){
    #     unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
    #       mean(smpl.data[smpl.data$bin.groups == j][[5]], na.rm = TRUE)
    #     }))
    #   } else if(bin.fun == "sd"){
    #     unlist(lapply(X = sort(unique(smpl.data$bin.groups)), FUN = function(j){
    #       stats::sd(smpl.data[smpl.data$bin.groups == j][[5]], na.rm = TRUE)
    #     }))
    #   } else if(bin.fun == "mad"){
    #     unlist(lapply(sort(unique(smpl.data$bin.groups)), function(j){
    #       stats::mad(smpl.data[smpl.data$bin.groups == j][[5]], na.rm = TRUE)
    #     }))
    #   } else {
    #     stop("Unsupported function. Supported functions: bin.fun = c('mean','sd' and 'mad').")
    #   }
    # })

    # vec_av.val2 <- unlist(list.fun.val2) # Make vector average val2
    # Map Bins average value on 3rd variable to the dataframe list
    list.dfs <- lapply(X = seq_along(list.dfs), FUN = function(i){
      list.bins <- split(x = list.dfs[[i]], f = list.dfs[[i]]$bin)
      list.bins <- Map(cbind, list.bins, bin.av.val2 = list.fun.val2[[i]])
      do.call(rbind, list.bins)
    })
    names(list.dfs) <- names(list.quant.lim)
    dframe <- data.table::rbindlist(l = list.dfs, idcol = "samples")
    # # Add Sample IDs, Var.grp and Var.col
    # list.dframes <- Map(
    #   cbind, Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
    #   Var.col = Annot.table[[3]], list.dfs)
    if(verbose){ cat("Done.\n") }
  } else { # No bin polygons
    dframe <- data.table::rbindlist(l = list_oriented_dens, idcol = "samples")
    # Add Sample IDs, Var.grp and Var.col
    # list.dframes <- Map(
    #   cbind, Var1 = Annot.table[[1]], Var.grp = Annot.table[[2]],
    #   Var.col = Annot.table[[3]], list_oriented_dens[Annot.table[[1]]])
  }
  dframe <- merge(x = annot_table, y = dframe, by = "samples")
  rm(list_vect.val1)
  # # Make data.frame
  # dframe <- do.call(rbind, list.dframes)

  # Plot
  craviola.plot <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(dframe$groups)) - 1,
      labels = levels(data$groups),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(
      x = old_cols[which(colnames(data) == "groups")],
      y = old_cols[which(colnames(data) == "vals")],
      fill = old_cols[which(colnames(data) == "conds")],
      alpha = paste0(
          bin.fun, "(", old_cols[which(colnames(data) == "opac")], ")")) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_fill_manual(
      values = BiocompR::biopalette(name = "BiocompR_cond"),
      labels = levels(annot_table$conds))
  # craviola.plot <- ggplot2::ggplot() +
  #   ggplot2::scale_x_continuous(
  #     breaks = as.integer(levels(dframe$Var.grp))-1,
  #     labels = levels(data[[2]]),
  #     expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  #   ggplot2::xlab(colnames(Annot.table)[2]) + ggplot2::ylab(colnames(data)[4]) +
  #   ggplot2::labs(fill = colnames(Annot.table)[3], color = "Extrema",
  #                 alpha = colnames(data)[5]) +
  #   ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
  #   ggplot2::scale_fill_manual(
  #     values = BiocompR::biopalette(name = "BiocompR_cond", mute = TRUE),
  #     labels = original.var.col)

  # Plot Options
  if(bins){ # bins  = TRUE
    alpha_lim <- c(
      floor(min(dframe$bin.av.val2)), ceiling(max(dframe$bin.av.val2)))
    craviola.plot <- craviola.plot +
      ggplot2::geom_polygon(
        data = dframe, mapping = ggplot2::aes(
          x = dens.curv, y = y.pos, fill = conds,
          group = interaction(samples, bin), alpha = bin.av.val2),
        colour = lines.col) +
      ggplot2::guides(alpha = ggplot2::guide_legend(order = 2)) +
      ggplot2::scale_alpha_continuous(
        limits = alpha_lim,
        breaks = round(seq(alpha_lim[1], alpha_lim[2], length.out = 5)))
    # craviola.plot <- craviola.plot +
    #   ggplot2::geom_polygon(data = dframe, mapping = ggplot2::aes(
    #     dens.curv, y.pos, fill = Var.col, group = interaction(Var1, bin),
    #     alpha = bin.av.val2), colour = lines.col) +
    #   ggplot2::guides(alpha = ggplot2::guide_legend(order = 2)) +
    #   ggplot2::scale_alpha_continuous(
    #     limits = c(floor(min(vec_av.val2)), ceiling(max(vec_av.val2))),
    #     breaks = round(seq(floor(min(vec_av.val2)), ceiling(max(vec_av.val2)),
    #                        length.out = 5)))
  } else { # bins = FALSE
    craviola.plot <- craviola.plot +
      ggplot2::geom_polygon(
        data = dframe, mapping = ggplot2::aes(
          x = dens.curv, y = y.pos, fill = conds, group = samples),
        colour = lines.col, linewidth = linewidth)
    # craviola.plot <- craviola.plot +
    #   ggplot2::geom_polygon(
    #     data = dframe, mapping = ggplot2::aes(
    #       dens.curv, y.pos, fill = Var.col, group = interaction(Var.col, Var1)),
    #     colour = lines.col, linewidth = linewidth)
  }
  rm(data)
  if(boxplots){ #boxplots TRUE
    craviola.plot <- craviola.plot +
      ggplot2::geom_boxplot(
        data = box.dframe,
        mapping = ggplot2::aes(
          x = x, ymin = low, lower = low, middle = mid, upper = top, ymax = top,
          group = x), stat = "identity")
  }
  if(mean.value){ #mean.value TRUE
    craviola.plot <- craviola.plot +
      ggplot2::geom_point(data = box.dframe, mapping = ggplot2::aes(
        x = x, y = mean), size = 2, color = "red")
  }
  return(craviola.plot)
}
