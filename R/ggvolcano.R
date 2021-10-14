
#' Checks the data.table provided to ggpanel.corr() and build.ggvolcano().
#'
#' @param data      A \code{data.table} with 3 to 5 columns:
#'                  \itemize{
#'                   \item{column 1 - The labels.}
#'                   \item{column 2 - Y-axis values.}
#'                   \item{column 3 - P-values.}
#'                   \item{column 4 (optional) - groups to be used to color
#'                   points. If none provided, you can use 'x.col.sign' to color
#'                   dots based on the X-axis category they belong to.}
#'                   \item{column 5 (optional) - values to be used to define
#'                   dots sizes based on another variable.}
#'                  }
#' @param data.type A \code{character} to specify the type of test used to
#'                  generate 'data'.
#'                  (Supported: data.type = c("t.test", "corr")).
#' @return A valid \code{data.table} with new column names.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

chk.dt <- function(data, data.type){
  #Check col2
  if(is.numeric(data[[2]])){
    if(data.type == "corr"){
      if(any(data[, 2] < -1) | any(data[, 2] > 1)){
        stop(
          "Correlation values in column 2 must be comprised between -1 and 1.")
      }
    }
  } else { stop("Column 2 type must be numeric.") }
  #Check col3
  if(is.numeric(data[[3]])){
    if(any(data[, 3] < 0) | any(data[, 3] > 1)){
      stop("P-values in column 3 must be comprised between 0 and 1.")}
  } else { stop("Column 3 type must be numeric.") }
  #Check data.type and select appropriate colnames
  if(data.type == "corr"){
    vec.colnames <- c("labels", "corr", "pval", "grp", "size")
  } else if(data.type == "t.test"){
    vec.colnames <- c("labels", "fold", "pval", "grp", "size")
  } else if(data.type == "free"){
    vec.colnames <- c("labels", "xval", "pval", "grp", "size")
  }
  #Check ncol(data) and change column names
  if(ncol(data) < 3){ stop("Data should contain at least 3 columns.")
  } else if(ncol(data) == 3){ colnames(data) <- vec.colnames[1:3]
  } else if(ncol(data) == 4){ colnames(data) <- vec.colnames[1:4]
  } else if(ncol(data) == 5){
    if(!is.numeric(data[[5]])){ stop("Column 5 type must be numeric.") }
    colnames(data) <- vec.colnames[1:5]
  } else if (ncol(data) > 5){
    stop("Too many columns. ncol(data) must be <= 5.")
  }
  return(data)
}


#' Checks cut-off value(s) given and computes the negative and positive
#' cut-offs.
#'
#' @param cutoff    A \code{numeric} vector of 1 or 2 values, to be used as a
#'                  cut-off on Y-axis values (Default: cutoff = 0).
#'                  \itemize{
#'                   \item{If 1 value is given, it will be defined as the
#'                   minimum cut-off on absolute Y-axis values
#'                   (positive and negative ones).}
#'                   \item{If 2 values are given: The smallest one will be used
#'                   as a maximum cut-off on negative Y-axis values. The biggest
#'                   one will be used as a minimum cut-off on positive Y-axis
#'                   values. The smallest value must be inferior or equal to 0.
#'                   The biggest value must be superior or equal to 0.}
#'                  }
#' @param data.type A \code{character} to specify the type of test used to
#'                  generate 'data'
#'                  (Supported: data.type = c("t.test", "corr")).
#' @return A \code{numeric} vector of length 2 containing the negative cut-off
#'         and the positive cut-off for Y-axis values.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

chk.cutoff <- function(cutoff = 0, data.type){
  #Assign cutoff values
  if(length(cutoff) == 2){
    neg.cutoff <- min(cutoff)
    pos.cutoff <- max(cutoff)
  } else if(length(cutoff) == 1){
    neg.cutoff<- -cutoff
    pos.cutoff<- cutoff
  } else { stop("'cutoff' must be either of length 1 or 2.") }
  #Check cutoff values
  if(data.type == "corr"){
    if(!(-1 <= neg.cutoff & neg.cutoff <= 0)){
      stop("min(cutoff) must be comprised between -1 and 0.")
    }
    if(!(0 <= pos.cutoff & pos.cutoff <= 1)){
      stop("max(cutoff) must be comprised between 0 and 1.")
    }
  } else if(data.type %in% c("t.test", "free")){
    if(!(neg.cutoff <= 0)){ stop("min(cutoff) must be inferior to 0.") }
    if(!(pos.cutoff >= 0)){ stop("max(cutoff) must be superior to 0.") }
  }
  return(c("negative.cutoff" = neg.cutoff, "positive.cutoff" = pos.cutoff))
}


#' Checks parameters for tests functions
#'
#' @param data         A \code{data.table} with 3 to 5 columns:
#'                     \itemize{
#'                      \item{column 1 - The labels.}
#'                      \item{column 2 - Y-axis values.}
#'                      \item{column 3 - P-values.}
#'                      \item{column 4 (optional) - groups to be used to color
#'                      points. If none provided, you can use 'x.col.sign' to
#'                      color dots based on the X-axis category they belong
#'                      to.}
#'                      \item{column 5 (optional) - values to be used to define
#'                      dots sizes based on another variable.}
#'                     }
#' @param data.type    A \code{character} to specify the type of test used to
#'                     generate 'data'
#'                     (Supported: data.type = c("t.test", "corr")).
#' @param label.cutoff A \code{numeric} vector of 1 or 2 values, to be used as a
#'                     cut-off on Y-axis values (Default: label.cutoff = 0).
#'                     \itemize{
#'                      \item{If 1 value is given, it will be defined as the
#'                      minimum cut-off on absolute Y-axis values
#'                      (positive and negative ones).}
#'                      \item{If 2 values are given: The smallest one will be
#'                      used as a maximum cut-off on negative Y-axis values. The
#'                      biggest one will be used as a minimum cut-off on
#'                      positive Y-axis values. The smallest value must be
#'                      inferior or equal to 0. The biggest value must be
#'                      superior or equal to 0.}
#'                     }
#'                     Dots outside these limits will be labeled. Dots within
#'                     the range of these limits will remain unlabeled.
#' @return A \code{list} containing:
#'         \itemize{
#'          \item{The original columns names from 'data'.}
#'          \item{The modified data.}
#'          \item{The label cut-off positive and negative values.}
#'         }
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

chk.param <- function(data, data.type, label.cutoff){
  #Get colnames
  orig.cnames <- colnames(data)
  #Convert as a data.table
  if(!data.table::is.data.table(data)){
    data <- data.table::as.data.table(data)
  }
  #Check & format data.table
  data <- BiocompR::chk.dt(data = data, data.type = data.type)
  #Check label.cutoff values
  cutoff.values <- BiocompR::chk.cutoff(
    cutoff = label.cutoff, data.type = data.type)
  return(list(orig.cnames, data, cutoff.values))
}


#' Plots results of correlation test between a single variable and multiple
#' others as jittered scatter plot divided into 4 different panels.
#'
#' @param data              A \code{data.table} with 3 to 5 columns:
#'                          \itemize{
#'                           \item{column 1 - strings to be used as labels for
#'                           individual dots.}
#'                           \item{column 2 - correlation values.}
#'                           \item{column 3 - correlation p-values.}
#'                           \item{column 4 (optional) - groups to be used for
#'                           coloring the dots.}
#'                           \item{column 5 (optional) - values to be used to
#'                           define dots sizes. It can be the sample size used
#'                           for the calculation of the correlation test.}
#'                          }
#' @param p.cutoff          A \code{numeric} between 0 and 1 to be used as a
#'                          maximum cut-off on p-values
#'                          (Default: p.cutoff = 0.01).
#' @param label.cutoff      A \code{numeric} vector of 1 or 2 values, between -1
#'                          and 1 to be used as the minimum cut-off on positive
#'                          and negative correlation values.
#'                          \itemize{
#'                           \item{If 1 value is given, it will be define as the
#'                           minimum cut-off on absolute correlation values
#'                           (positive and negative ones).}
#'                           \item{If 2 values are given, the smallest one will
#'                           be used as a maximum cut-off on negative
#'                           correlation values. the biggest one will be used as
#'                           a minimum cut-off on positive correlation values.
#'                           The smallest value must be comprised between -1 and
#'                           0. The biggest value must be comprised between 0
#'                           and 1.}
#'                          }
#' @param jitter.height A \code{numeric} to specify the amount of vertical
#'                      jitter (Default: jitter.height = 0.4).
#' @return A \code{gg} plot object with 4 panels:
#'         \itemize{
#'          \item{1 panel with significant positive correlation values.}
#'          \item{1 panel with significant negative correlation values.}
#'          \item{1 panel with non-significant positive correlation values.}
#'          \item{1 panel with non-significant positive correlation values.}
#'         }
#'         Each panel displays results as jittered scatter plots.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=`
#' @export

ggpanel.corr <- function(
  data, p.cutoff = 0.01, label.cutoff = 0, jitter.height = 0.4){
  #Check test parameters
  res.param <- BiocompR::chk.param(
    data = data, data.type = "corr", label.cutoff = label.cutoff)
  orig.cnames <- res.param[[1]]
  data <- res.param[[2]]
  cutoff.values <- res.param[[3]]

  #Define P-value intervals
  data$P.value <- cut(
    x = data$pval, breaks = c(min(data$pval), p.cutoff, max(data$pval)),
    labels = c(paste0("<= ", p.cutoff,":\n[ ", formatC(
      x = min(data$pval), format = "e", digits = 2), ", ", p.cutoff, " ]"),
      paste0("> ", p.cutoff, ":\n] ",p.cutoff,", ",
             formatC(x = max(data$pval), format = "e", digits = 2), " ]")),
    include.lowest = TRUE)
  #Define correlation value intervals
  data$cor.cat <- cut(x = data$corr, breaks = c(-Inf, 0, +Inf),
                      labels = c("Negative Correlation","Positive Correlation"))
  #Remove labels out of cut-offs
  data[!(corr >= cutoff.values[2] | corr <= cutoff.values[1]), labels := ""]
  #Seed position for using geom_label_repel() with geom_jitter().
  pos <- ggplot2::position_jitter(seed = 1, height = jitter.height)
  #Plot Spearman correlation results
  if(ncol(data) == 5){
    ggpan <- ggplot2::ggplot(
      data = data, mapping = ggplot2::aes(
        x = corr, y = P.value, label = labels))
  } else if(ncol(data) == 6){
    ggpan <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(
      x = corr, y = P.value, label = labels, color = grp))
  } else{
    ggpan <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(
      x = corr, y = P.value, label = labels, color = grp, size = size))
  }
  ggpan <- ggpan +
    ggplot2::geom_jitter(position = pos) +
    ggplot2::facet_grid(
      P.value ~ cor.cat, scales = "free", space = "free", switch = "y") +
    ggrepel::geom_label_repel(position = pos, size = 4) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey"),
      panel.border = ggplot2::element_rect(
        color = "black", fill = NA, size = 1),
      axis.title = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(size = 12),
      strip.background = ggplot2::element_rect(
        fill = NA, colour = 'black', size = 1),
      legend.key = ggplot2::element_blank()) +
    ggplot2::labs(
      x = orig.cnames[2], color = orig.cnames[4], size = orig.cnames[5])
  return(ggpan)
}


#' Core function that builds a volcano plot using ggplot2.
#'
#' @param data           A \code{data.table} with 3 to 5 columns:
#'                       \itemize{
#'                        \item{column 1 - The labels.}
#'                        \item{column 2 - X-axis values.}
#'                        \item{column 3 - P-values.}
#'                        \item{column 4 (optional) - groups to be used to color
#'                        points. If none provided, you can use 'x.col.sign' to
#'                        color dots based on the X-axis category they belong
#'                        to.}
#'                        \item{column 5 (optional) - values to be used to
#'                        define dots sizes based on another variable.}
#'                       }
#' @param data.type      A \code{character} to specify the type of test used to
#'                       generate 'data'
#'                       (Supported: data.type = c("t.test", "corr")).
#' @param label.cutoff   A \code{numeric} vector of 1 or 2 values, to be used as
#'                       a cut-off on Y-axis values (Default: label.cutoff = 0).
#'                       \itemize{
#'                        \item{If 1 value is given, it will be defined as the
#'                        minimum cut-off on absolute Y-axis values
#'                        (positive and negative ones).}
#'                        \item{If 2 values are given: The smallest one will be
#'                        used as a maximum cut-off on negative Y-axis values.
#'                        The biggest one will be used as a minimum cut-off on
#'                        positive Y-axis values. The smallest value must be
#'                        inferior or equal to 0. The biggest value must be
#'                        superior or equal to 0.}
#'                       }
#'                       Dots outside these limits will be labeled. Dots within
#'                       the range of these limits will remain unlabeled.
#' @param p.cutoff       A \code{numeric} between 0 and 1 to be used as a
#'                       maximum cut-off on p-values (Default: p.cutoff = 0.01).
#' @param x.cutoff       A \code{numeric} vector of 1 or 2 values, to be used as
#'                       a cut-off on Y-axis values (Default: x.cutoff = NULL).
#'                       \itemize{
#'                        \item{If 1 value is given, it will be defined as the
#'                        minimum cut-off on absolute Y-axis values
#'                        (positive and negative ones).}
#'                        \item{If 2 values are given: The smallest one will be
#'                        used as a maximum cut-off on negative Y-axis values.
#'                        The biggest one will be used as a minimum cut-off on
#'                        positive Y-axis values. The smallest value must be
#'                        inferior or equal to 0. The biggest value must be
#'                        superior or equal to 0.}
#'                       }
#'                       x.cutoff value(s) will be used to draw vertical lines
#'                       on the volcano plot.
#' @param title.x.cutoff A \code{character} to specify how you wish to label the
#'                       Y-axis vertical lines on the volcano plot
#'                       (Default: title.cutoff = "X-Axis cutoff").
#' @param x.col.sign     A \code{logical} to activate automatic coloring of data
#'                       based on their correlation category. If TRUE, elements
#'                       will be divided into 3 categories:
#'                       \itemize{
#'                        \item{Blue - elements below the negative 'x.cutoff'.}
#'                        \item{Red - elements above the positive 'x.cutoff'.}
#'                        \item{Grey - any other element not meeting the
#'                        requirements in the first nor the second categories.}
#'                       }
#'                       If FALSE, dots will be colored using groups from data
#'                       (if any provided). If no groups are passed within data,
#'                       dots will remain black.
#' @param force.label    A \code{character} vector matching elements in the
#'                       first column of 'data' that you wish to label,
#'                       independently from their actual significance
#'                       (Default: force.label = NULL). force.label overrides
#'                       other labeling parameters.
#' @return A \code{gg} volcano plot of your statistical test results.
#' @author Yoann Pageaud.
#' @export
#' @importFrom data.table `:=`
#' @keywords internal

build.ggvolcano <- function(
  data, data.type, label.cutoff = 0, p.cutoff = 0.01, x.cutoff = NULL,
  title.x.cutoff = "X-Axis cutoff", x.col.sign = FALSE, force.label = NULL){
  #Fix BiocCheck() complaining about these objects initialization
  corr <- NULL
  grp <- NULL
  fold <- NULL
  xval <- NULL
  pval <- NULL
  size <- NULL
  `P-value` <- NULL
  #Make a copy of the data.table to work on
  dt.data <- data.table::as.data.table(data)
  #Check test parameters
  res.param <- BiocompR::chk.param(
    data = dt.data, data.type = data.type, label.cutoff = label.cutoff)
  orig.cnames <- res.param[[1]]
  dt.data <- res.param[[2]]
  cutoff.values <- res.param[[3]]
  #Check x.cutoff values
  x.cutoff.values <- BiocompR::chk.cutoff(
    cutoff = x.cutoff, data.type = data.type)
  #If x.col.sign TRUE add grp color to data
  if(x.col.sign){
    if(data.type == "corr"){
      groups <- c("Negatively correlated", "Insufficiently correlated",
                  "Positively correlated")
      dt.data[corr <= x.cutoff.values[1], grp := groups[1]]
      dt.data[corr >= x.cutoff.values[2], grp := groups[3]]
      dt.data[corr > x.cutoff.values[1] & corr < x.cutoff.values[2],
           grp := groups[2]]
    } else if(data.type == "t.test"){
      groups <- c("Negative log2(Fold change)", "Insufficient log2(Fold change)",
                  "Positive log2(Fold change)")
      dt.data[fold <= x.cutoff.values[1], grp := groups[1]]
      dt.data[fold >= x.cutoff.values[2], grp := groups[3]]
      dt.data[fold > x.cutoff.values[1] & fold < x.cutoff.values[2],
           grp := groups[2]]
    } else if(data.type == "free"){
      groups <- c("Negatively associated", "Insufficiently associated",
                  "Positively associated")
      dt.data[xval <= x.cutoff.values[1], grp := groups[1]]
      dt.data[xval >= x.cutoff.values[2], grp := groups[3]]
      dt.data[xval > x.cutoff.values[1] & xval < x.cutoff.values[2],
           grp := groups[2]]
    }
    dt.data[, grp := as.factor(grp)]
    dt.data[, grp := factor(grp, levels = levels(grp)[
      order(match(levels(grp), groups))])]
  }
  #Create shading conditions
  dt.data[pval > p.cutoff, "P-value" := as.factor(paste0("> ", p.cutoff))]
  dt.data[pval <= p.cutoff, "P-value" := paste0("<= ", p.cutoff)]
  #Build scatter plot
  ggvol <- ggplot2::ggplot()
  if(ncol(dt.data) == 6){
    ggvol <- ggvol + ggplot2::geom_point(
      data = dt.data, mapping = ggplot2::aes(
        x = dt.data[[2]], y = -log10(pval), color = grp, size = size,
        alpha = `P-value`))
  } else if(ncol(dt.data) == 5){
    ggvol <- ggvol + ggplot2::geom_point(
      data = dt.data, mapping = ggplot2::aes(
        x = dt.data[[2]], y = -log10(pval), color = grp, alpha = `P-value`))
  } else {
    ggvol <- ggvol + ggplot2::geom_point(
      data = dt.data, mapping = ggplot2::aes(
        x = dt.data[[2]], y = -log10(pval), alpha = `P-value`))
  }
  if(x.col.sign){
    ggvol <- ggvol + ggplot2::scale_color_manual(
      values = c("darkblue", "grey", "darkred"))
  }
  #Make p-value cut-off
  ggvol <- ggvol + ggplot2::geom_hline(
    yintercept = -log10(p.cutoff), color = 'black')
  #Make Y-Axis cut-off
  if(!is.null(x.cutoff)){
    #Make negative and positive cut-off
    ggvol <- ggvol +
      ggplot2::geom_vline(xintercept = x.cutoff.values[1], color = 'darkblue') +
      ggplot2::geom_vline(xintercept = x.cutoff.values[2], color = 'darkred')
  }
  #Create labels table
  if(is.null(force.label)){
    if(data.type == "corr"){
      dt.label <- dt.data[(corr >= cutoff.values[2] | corr <= cutoff.values[1]) &
                         pval <= p.cutoff]
    } else if(data.type == "t.test"){
      dt.label <- dt.data[(fold >= cutoff.values[2] | fold <= cutoff.values[1]) &
                         pval <= p.cutoff]
    } else if(data.type == "free"){
      dt.label <- dt.data[(xval >= cutoff.values[2] | xval <= cutoff.values[1]) &
                         pval <= p.cutoff]
    } else { stop("Unsupported 'data.type'.") }
  } else {
    #Keep only some labels of interest
    if(all(force.label %in% dt.data[["labels"]])){
      #Override label display to force labeling of specific data
      dt.label <- dt.data[labels %in% force.label]
    } else {
      warning(
        "Some labels specified in 'force.label' are not part of the dataset.")
      force.label <- force.label[force.label %in% dt.data[["labels"]]]
      dt.label <- dt.data[labels %in% force.label]
    }
  }
  #Create labels
  if(ncol(dt.data) > 4){
    ggvol <- ggvol + ggrepel::geom_label_repel(
      data = dt.label, mapping = ggplot2::aes(
        x = dt.label[[2]], y = -log10(pval), label = labels, color = grp),
      size = 4.5, max.overlaps = Inf)
  } else if(ncol(dt.data) == 4){
    ggvol <- ggvol + ggrepel::geom_label_repel(
      data = dt.label, mapping = ggplot2::aes(
        x = dt.label[[2]], y = -log10(pval), label = labels),
      size = 4.5, max.overlaps = Inf)
  }
  #Create P-value label
  ggvol <- ggvol + ggrepel::geom_label_repel(
    data = data.frame(), mapping = ggplot2::aes(
      x = -Inf, y = -log10(p.cutoff), fontface = 1, label = "P-value cut-off"),
    color = "black", direction = "x", size = 4)
  #Create Y-axis label
  if(!is.null(x.cutoff)){
    #Make negative and positive cut-off label
    ggvol <- ggvol +
      ggrepel::geom_label_repel(data = data.frame(), mapping = ggplot2::aes(
        x = x.cutoff.values[1], y = Inf, fontface = 1, label = title.x.cutoff),
        color = "darkblue", direction = "y", size = 4) +
      ggrepel::geom_label_repel(data = data.frame(), mapping = ggplot2::aes(
        x = x.cutoff.values[2], y = Inf, fontface = 1, label = title.x.cutoff),
        color = "darkred", direction = "y", size = 4)
  }
  #Add default volcano theme
  ggvol <- ggvol +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey"),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.key = ggplot2::element_blank()) +
    ggplot2::labs(x = orig.cnames[2], y = paste0(
      "-log10(", orig.cnames[3], ")"), color = orig.cnames[4],
      size = orig.cnames[5])
  #Return plot removing warning about discrete variable given to alpha
  BiocompR::warn.handle(
    pattern = "Using alpha for a discrete variable is not advised.",
    print(ggvol))
}


#' Plots results of correlation test between a single variable and multiple
#' others as volcano plot.
#'
#' @param data         A \code{data.table} with 3 to 5 columns:
#'                     \itemize{
#'                      \item{column 1 - The labels.}
#'                      \item{column 2 - correlation values.}
#'                      \item{column 3 - P-values.}
#'                      \item{column 4 (optional) - groups to be used to color
#'                      points. If none provided, you can use 'x.col.sign' to
#'                      color dots based on the correlation category they belong
#'                      to.}
#'                      \item{column 5 (optional) - values to be used to define
#'                      dots sizes based on another variable.}
#'                     }
#' @param p.cutoff     A \code{numeric} between 0 and 1 to be used as a maximum
#'                     cut-off on p-values (Default: p.cutoff = 0.01).
#' @param corr.cutoff  A \code{numeric} vector of 1 or 2 values, comprised
#'                     between -1 and 1, to be used as a cut-off on correlation
#'                     values (Default: corr.cutoff = NULL).
#'                     \itemize{
#'                      \item{If 1 value is given, it will be defined as the
#'                      minimum cut-off on absolute correlation values
#'                      (positive and negative ones).}
#'                      \item{If 2 values are given: The smallest one will be
#'                      used as a maximum cut-off on negative correlation
#'                      values. The biggest one will be used as a minimum
#'                      cut-off on positive correlation values.
#'                      The smallest value must be inferior or equal to 0. The
#'                      biggest value must be superior or equal to 0.}
#'                     }
#'                     corr.cutoff value(s) will be used to draw vertical lines
#'                     on the volcano plot.
#' @param title.cutoff A \code{character} to specify how you wish to label the
#'                     correlation limits on the volcano plot
#'                     (Default: title.cutoff = "Correlation cut-off").
#' @param label.cutoff A \code{numeric} vector of 1 or 2 values, comprised
#'                     between -1 and 1, to be used as a cut-off on correlation
#'                     values (Default: corr.cutoff = NULL).
#'                     \itemize{
#'                      \item{If 1 value is given, it will be defined as the
#'                      minimum cut-off on absolute correlation values
#'                      (positive and negative ones).}
#'                      \item{If 2 values are given: The smallest one will be
#'                      used as a maximum cut-off on negative correlation
#'                      values. The biggest one will be used as a minimum
#'                      cut-off on positive correlation values.
#'                      The smallest value must be inferior or equal to 0. The
#'                      biggest value must be superior or equal to 0.}
#'                     }
#'                     Dots outside these limits will be labeled. Dots within
#'                     the range of these limits will remain unlabeled.
#' @param x.col.sign   A \code{logical} to activate automatic coloring of data
#'                     based on their correlation category. If TRUE, elements
#'                     will be divided into 3 categories:
#'                     \itemize{
#'                      \item{Blue - elements below the negative 'corr.cutoff'.}
#'                      \item{Red - elements above the positive 'corr.cutoff'.}
#'                      \item{Grey - any other element not meeting the
#'                      requirements in the first nor the second categories.}
#'                     }
#'                     If FALSE, dots will be colored using groups from data
#'                     (if any provided). If no groups are passed within data,
#'                     dots will remain black.
#' @param force.label  A \code{character} vector matching elements in the
#'                     first column of 'data' that you wish to label,
#'                     independently from their actual significance
#'                     (Default: force.label = NULL). force.label overrides
#'                     other labeling parameters.
#' @return A \code{gg} volcano plot of your correlation test results.
#' @author Yoann Pageaud.
#' @export

ggvolcano.corr <- function(
  data, p.cutoff = 0.01, corr.cutoff = NULL,
  title.cutoff = "Correlation cut-off", label.cutoff = 0, x.col.sign = FALSE,
  force.label = NULL){
  #Build volcano plot for correlation data
  BiocompR::build.ggvolcano(
    data = data, data.type = "corr", label.cutoff = label.cutoff,
    p.cutoff = p.cutoff, x.cutoff = corr.cutoff, title.x.cutoff = title.cutoff,
    x.col.sign = x.col.sign, force.label = force.label)
}


#' Plots results of statistical tests as a volcano plot.
#'
#' @param data         A \code{data.table} with 3 to 5 columns:
#'                     \itemize{
#'                      \item{column 1 - The labels.}
#'                      \item{column 2 - log2(fold changes) values.}
#'                      \item{column 3 - P-values.}
#'                      \item{column 4 (optional) - groups to be used to color
#'                      points. If none provided, you can use 'x.col.sign' to
#'                      color dots based on the fold change category they belong
#'                      to.}
#'                      \item{column 5 (optional) - values to be used to define
#'                      dots sizes based on another variable.}
#'                     }
#' @param p.cutoff     A \code{numeric} between 0 and 1 to be used as a maximum
#'                     cut-off on p-values (Default: p.cutoff = 0.01).
#' @param l2fc.cutoff  A \code{numeric} vector of 1 or 2 values, to be used as
#'                     a cut-off on log2(fold change) values
#'                     (Default: l2fc.cutoff = NULL).
#'                     \itemize{
#'                      \item{If 1 value is given, it will be defined as the
#'                      minimum cut-off on absolute log2(fold change) values
#'                      (positive and negative ones).}
#'                      \item{If 2 values are given: The smallest one will be
#'                      used as a maximum cut-off on negative log2(fold change)
#'                      values. The biggest one will be used as a minimum
#'                      cut-off on positive log2(fold change) values.
#'                      The smallest value must be inferior or equal to 0. The
#'                      biggest value must be superior or equal to 0.}
#'                     }
#'                     l2fc.cutoff value(s) will be used to draw vertical lines
#'                     on the volcano plot.
#' @param title.cutoff A \code{character} to specify how you wish to label the
#'                     log2(fold change) limits on the volcano plot
#'                     (Default: title.cutoff = "L2FC cut-off").
#' @param label.cutoff A \code{numeric} vector of 1 or 2 values, to be used as
#'                     a cut-off on log2(fold change) values
#'                     (Default: label.cutoff = 0).
#'                     \itemize{
#'                      \item{If 1 value is given, it will be defined as the
#'                      minimum cut-off on absolute log2(fold change) values
#'                      (positive and negative ones).}
#'                      \item{If 2 values are given: The smallest one will be
#'                      used as a maximum cut-off on negative log2(fold change)
#'                      values. The biggest one will be used as a minimum
#'                      cut-off on positive log2(fold change) values.
#'                      The smallest value must be inferior or equal to 0. The
#'                      biggest value must be superior or equal to 0.}
#'                     }
#'                     Dots outside these limits will be labeled. Dots within
#'                     the range of these limits will remain unlabeled.
#' @param x.col.sign   A \code{logical} to activate automatic coloring of data
#'                     based on their fold change category. If TRUE, elements
#'                     will be divided into 3 categories:
#'                     \itemize{
#'                      \item{Blue - elements below the negative 'l2fc.cutoff'.}
#'                      \item{Red - elements above the positive 'l2fc.cutoff'.}
#'                      \item{Grey - any other element not meeting the
#'                      requirements in the first nor the second categories.}
#'                     }
#'                     If FALSE, dots will be colored using groups from data
#'                     (if any provided). If no groups are passed within data,
#'                     dots will remain black.
#' @param l2.transform A \code{logical} to specify whether the data in the
#'                     second column should be log2 transformed
#'                     (l2.transform = TRUE) or if the log2 transformation has
#'                     already been applied and no further transformation is
#'                     needed (Default: l2.transform = TRUE).
#' @param force.label  A \code{character} vector matching elements in the
#'                     first column of 'data' that you wish to label,
#'                     independently from their actual significance
#'                     (Default: force.label = NULL). force.label overrides
#'                     other labeling parameters.
#' @return A \code{gg} volcano plot of your statistical test results.
#' @author Yoann Pageaud, Verena Bitto.
#' @importFrom data.table `:=`
#' @export

ggvolcano.test <- function(
  data, p.cutoff = 0.01, l2fc.cutoff = NULL, title.cutoff = "L2FC cut-off",
  label.cutoff = 0, x.col.sign = FALSE, l2.transform = FALSE,
  force.label = NULL){
  #Apply log2 transformation if needed
  if(l2.transform){ data[, (2) := log2(x = data[[2]])] }
  #Build volcano plot for t-test data
  BiocompR::build.ggvolcano(
    data = data, data.type = "t.test", label.cutoff = label.cutoff,
    p.cutoff = p.cutoff, x.cutoff = l2fc.cutoff, title.x.cutoff = title.cutoff,
    x.col.sign = x.col.sign, force.label = force.label)
}


#' Plots any king of results with P-values that can be displayed as a volcano
#' plot.
#'
#' @param data           A \code{data.table} with 3 to 5 columns:
#'                       \itemize{
#'                        \item{column 1 - The labels.}
#'                        \item{column 2 - X-axis values.}
#'                        \item{column 3 - P-values.}
#'                        \item{column 4 (optional) - groups to be used to color
#'                        points. If none provided, you can use 'x.col.sign' to
#'                        color dots based on the Y-axis category they belong
#'                        to.}
#'                        \item{column 5 (optional) - values to be used to
#'                        define dots sizes based on another variable.}
#'                       }
#' @param label.cutoff   A \code{numeric} vector of 1 or 2 values, to be used as
#'                       a cut-off on Y-axis values (Default: label.cutoff = 0).
#'                       \itemize{
#'                        \item{If 1 value is given, it will be defined as the
#'                        minimum cut-off on absolute Y-axis values
#'                        (positive and negative ones).}
#'                        \item{If 2 values are given: The smallest one will be
#'                        used as a maximum cut-off on negative Y-axis values.
#'                        The biggest one will be used as a minimum cut-off on
#'                        positive Y-axis values. The smallest value must be
#'                        inferior or equal to 0. The biggest value must be
#'                        superior or equal to 0.}
#'                       }
#'                       Dots outside these limits will be labeled. Dots within
#'                       the range of these limits will remain unlabeled.
#' @param p.cutoff       A \code{numeric} between 0 and 1 to be used as a
#'                       maximum cut-off on p-values (Default: p.cutoff = 0.01).
#' @param x.cutoff       A \code{numeric} vector of 1 or 2 values, to be used as
#'                       a cut-off on Y-axis values (Default: x.cutoff = NULL).
#'                       \itemize{
#'                        \item{If 1 value is given, it will be defined as the
#'                        minimum cut-off on absolute Y-axis values
#'                        (positive and negative ones).}
#'                        \item{If 2 values are given: The smallest one will be
#'                        used as a maximum cut-off on negative Y-axis values.
#'                        The biggest one will be used as a minimum cut-off on
#'                        positive Y-axis values. The smallest value must be
#'                        inferior or equal to 0. The biggest value must be
#'                        superior or equal to 0.}
#'                       }
#'                       x.cutoff value(s) will be used to draw vertical lines
#'                       on the volcano plot.
#' @param title.x.cutoff A \code{character} to specify how you wish to label the
#'                       Y-axis vertical lines on the volcano plot
#'                       (Default: title.cutoff = "X-Axis cutoff").
#' @param x.col.sign     A \code{logical} to activate automatic coloring of data
#'                       based on their correlation category. If TRUE, elements
#'                       will be divided into 3 categories:
#'                       \itemize{
#'                        \item{Blue - elements below the negative 'x.cutoff'.}
#'                        \item{Red - elements above the positive 'x.cutoff'.}
#'                        \item{Grey - any other element not meeting the
#'                        requirements in the first nor the second categories.}
#'                       }
#'                       If FALSE, dots will be colored using groups from data
#'                       (if any provided). If no groups are passed within data,
#'                       dots will remain black.
#' @param force.label    A \code{character} vector matching elements in the
#'                       first column of 'data' that you wish to label,
#'                       independently from their actual significance
#'                       (Default: force.label = NULL). force.label overrides
#'                       other labeling parameters.
#' @return A \code{gg} volcano plot of your results.
#' @author Yoann Pageaud.
#' @importFrom data.table `:=`
#' @export

ggvolcano.free <- function(
  data, label.cutoff = 0, p.cutoff = 0.01, x.cutoff = NULL,
  title.x.cutoff = "X-Axis cutoff", x.col.sign = FALSE, force.label = NULL){
  #Build volcano plot for free data
  BiocompR::build.ggvolcano(
    data = data, data.type = "free", label.cutoff = label.cutoff,
    p.cutoff = p.cutoff, x.cutoff = x.cutoff, title.x.cutoff = title.x.cutoff,
    x.col.sign = x.col.sign, force.label = force.label)
}
