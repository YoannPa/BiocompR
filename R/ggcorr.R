
#' Checks the data.table provided to ggpanel.corr() and ggvolcano.corr().
#'
#' @param data A \code{data.table} which contains correlation data.
#' @return A \code{data.table} with new column names.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

chk.dt<-function(data){
  #Check col2
  if(is.numeric(data[[2]])){
    if(any(data[,2] < -1) | any(data[,2] > 1)){
      stop("Correlation values in column 2 must be comprised between -1 and 1.")
    }
  } else {stop("Column 2 type must be numeric.")}
  #Check col3
  if(is.numeric(data[[3]])){
    if(any(data[,3] < 0) | any(data[,3] > 1)){
      stop("P-values in column 3 must be comprised between 0 and 1.")}
  } else { stop("Column 3 type must be numeric.") }
  #Check ncol(data)
  if(ncol(data) < 3){ stop("Data should contain at least 3 columns.")
  } else if(ncol(data) == 3){
    colnames(data) <- c("labels","corr","pval")
  } else if(ncol(data) == 4){
    colnames(data) <- c("labels","corr","pval","grp")
  } else if(ncol(data) == 5){
    colnames(data) <- c("labels","corr","pval","grp","size")
    if(!is.numeric(data[[5]])){ stop("Column 5 type must be numeric.") }
  } else if (ncol(data) > 5){
    stop("Too many columns. ncol(data) must be <= 5.")
  }
  return(data)
}

#' Checks the correlation cut-off value(s) given and computes the negative and
#' positive correlation cut-offs.
#'
#' @param corr.label.cutoff A \code{numeric} vector of 1 or 2 values, between -1
#'                          and 1 to be used as the minimum cut-off on
#'                          positive and negative correlation values.
#'                          \itemize{
#'                           \item{If 1 value is given, it will be define as the
#'                           minimum cut-off on absolute correlation values
#'                           (positive and negative ones).}
#'                           \item{If 2 values are given, the smalest one will
#'                           be used as a maximum cut-off on negative
#'                           correlation values. the biggest one will be used as
#'                           a minimum cut-off on positive correlation values.
#'                           The smalest value must be comprised between -1 and
#'                           0. The biggest value must be comprised between 0
#'                           and 1.}
#'                          }
#' @return A \code{numeric} vector of length 2 containing the negative cut-off
#'         and the positive cut-off for correlation values.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

chk.corr.lbl.cutoff<-function(corr.label.cutoff){
  #Check corr.label.cutoff values
  if(length(corr.label.cutoff) == 2){
    neg.corr.label.cutoff<-min(corr.label.cutoff)
    pos.corr.label.cutoff<-max(corr.label.cutoff)
    if(!(-1 <= neg.corr.label.cutoff & neg.corr.label.cutoff <= 0)){
      stop("min(corr.label.cutoff) must be comprised between -1 and 0.")
    }
    if(!(0 <= pos.corr.label.cutoff & pos.corr.label.cutoff <= 1)){
      stop("max(corr.label.cutoff) must be comprised between 0 and 1.")
    }
  } else if(length(corr.label.cutoff) == 1){
    neg.corr.label.cutoff<- -corr.label.cutoff
    pos.corr.label.cutoff<- corr.label.cutoff
  } else {
    stop("'corr.label.cutoff' must be either of length 1 or 2.")
  }
  return(c("negative.cutoff" = neg.corr.label.cutoff,
           "positive.cutoff" = pos.corr.label.cutoff))
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
#' @param corr.label.cutoff A \code{numeric} vector of 1 or 2 values, between -1
#'                          and 1 to be used as the minimum cut-off on
#'                          positive and negative correlation values.
#'                          \itemize{
#'                           \item{If 1 value is given, it will be define as the
#'                           minimum cut-off on absolute correlation values
#'                           (positive and negative ones).}
#'                           \item{If 2 values are given, the smalest one will
#'                           be used as a maximum cut-off on negative
#'                           correlation values. the biggest one will be used as
#'                           a minimum cut-off on positive correlation values.
#'                           The smalest value must be comprised between -1 and
#'                           0. The biggest value must be comprised between 0
#'                           and 1.}
#'                          }
#' @param jitter.height     A \code{numeric} to specify the amount of vertical
#'                          jitter (Default: jitter.height = 0.4).
#' @return A \code{gg} plot object with 4 panels:
#'         \itemize{
#'          \item{1 panel with significant positive correlation values.}
#'          \item{1 panel with significant negative correlation values.}
#'          \item{1 panel with non-significant positive correlation values.}
#'          \item{1 panel with non-significant positive correlation values.}
#'         }
#'         Each panel displays results as jittered scatter plots.
#' @author Yoann Pageaud.
#' @export

ggpanel.corr<-function(
  data, p.cutoff = 0.01, corr.label.cutoff, jitter.height = 0.4){

  #Get colnames
  orig.cnames <- colnames(data)
  #Convert as a data.table
  if(!is.data.table(data)){ data<-as.data.table(data) }
  #Check & format data.table
  data <- chk.dt(data = data)
  #Check corr.label.cutoff values
  cutoff.values <- chk.corr.lbl.cutoff(corr.label.cutoff = corr.label.cutoff)
  #Define P-value intervals
  data$P.value <- cut(
    x = data$pval, breaks=c(min(data$pval), p.cutoff, max(data$pval)),
    labels = c(paste0("<= ", p.cutoff,":\n[ ", formatC(x=min(data$pval),
                                                      format="e",digits=2),
                      ", ", p.cutoff, " ]"),
               paste0("> ", p.cutoff, ":\n] ",p.cutoff,", ",
                      formatC(x=max(data$pval), format="e",digits=2)," ]")),
    include.lowest = TRUE)
  #Define correlation value intervals
  data$cor.cat <- cut(x = data$corr, breaks = c(-Inf, 0, +Inf),
                      labels = c("Negative Correlation","Positive Correlation"))
  #Remove labels out of cut-offs
  data[!(corr >= cutoff.values[2] | corr <= cutoff.values[1]), labels := ""]
  #Seed position for using geom_label_repel() with geom_jitter().
  pos<-position_jitter(seed = 1, height = jitter.height)
  #Plot Spearman correlation results
  if(ncol(data) == 5){
    ggpan <- ggplot(data = data,
                    mapping = aes(x = corr, y = P.value, label = labels))
  } else if(ncol(data) == 6){
    ggpan <- ggplot(data = data, mapping = aes(
      x = corr, y = P.value, label = labels, color = grp))
  } else{
    ggpan <- ggplot(data = data, mapping = aes(
      x = corr, y = P.value, label = labels, color = grp, size = size))
  }
  ggpan <- ggpan +
    geom_jitter(position = pos) +
    facet_grid(P.value ~ cor.cat, scales = "free", space = "free", switch="y") +
    ggrepel::geom_label_repel(position = pos, size = 4) +
    theme(axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.x = element_line(colour = "grey"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.title = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = NA, colour = 'black', size=1),
          legend.key = element_blank()) +
    labs(x = orig.cnames[2], color = orig.cnames[4], size = orig.cnames[5])
  return(ggpan)
}

#' Plots results of correlation test between a single variable and multiple
#' others as volcano plot.
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
#' @param corr.cutoff       A \code{numeric} value between -1 and 1 to be used
#'                          as a cut-off on correlation values.
#'                          \itemize{
#'                           \item{if corr.cutoff < 0: it will be a maximum
#'                           cut-off on negative correlation values.}
#'                           \item{if corr.cutoff > 0: it will be a minimum
#'                           cut-off on positive correlation values.}
#'                          }
#'                          It is not used as an absolute cut-off on positive
#'                          and negative correlation values as if the sign of
#'                          the correlation value is different from the one of
#'                          the cut-off, it might not necessary make sense apply
#'                          the cut-off to it.
#' @param title.corr.cutoff A \code{character} to be used as the label of the
#'                          correlation cut-off vertical blue line
#'                          (Default: title.corr.cutoff="Correlation cut-off").
#' @param corr.label.cutoff A \code{numeric} vector of 1 or 2 values, between -1
#'                          and 1 to be used as the minimum cut-off on
#'                          positive and negative correlation values.
#'                          \itemize{
#'                           \item{If 1 value is given, it will be define as the
#'                           minimum cut-off on absolute correlation values
#'                           (positive and negative ones).}
#'                           \item{If 2 values are given, the smalest one will
#'                           be used as a maximum cut-off on negative
#'                           correlation values. the biggest one will be used as
#'                           a minimum cut-off on positive correlation values.
#'                           The smalest value must be comprised between -1 and
#'                           0. The biggest value must be comprised between 0
#'                           and 1.}
#'                          }
#' @return A \code{gg} volcano plot object.
#' @author Yoann Pageaud.
#' @export

ggvolcano.corr<-function(
  data, p.cutoff = 0.01, corr.cutoff = NULL,
  title.corr.cutoff = "Correlation cut-off", corr.label.cutoff){

  #Get colnames
  orig.cnames <- colnames(data)
  #Convert as a data.table
  if(!is.data.table(data)){ data<-as.data.table(data) }
  #Check & format data.table
  data <- chk.dt(data = data)
  #Check corr.label.cutoff values
  cutoff.values<-chk.corr.lbl.cutoff(corr.label.cutoff = corr.label.cutoff)

  #Create shading conditions
  data[pval > p.cutoff, "P-value" := as.factor(paste0("> ", p.cutoff))]
  data[pval <= p.cutoff, "P-value" := paste0("<= ", p.cutoff)]
  #Scatter plot Correlation of L1 Methylation with L1 transposition
  ggvol <- ggplot()
  if(ncol(data) == 6){
    ggvol <- ggvol +
      geom_point(data = data, aes(x = corr, y = -log10(pval), color = grp,
                                  size = size, alpha = `P-value`))
  } else if(ncol(data) == 5){
    ggvol <- ggvol +
      geom_point(data = data, aes(x = corr, y = -log10(pval), color = grp,
                                  alpha = `P-value`))
  } else {
    ggvol <- ggvol +
      geom_point(data = data, aes(x = corr, y = -log10(pval), alpha=`P-value`))
  }
  ggvol <- ggvol + geom_hline(yintercept = -log10(p.cutoff), color = 'red')

  if(!is.null(corr.cutoff)){
    ggvol <- ggvol +
      geom_vline(xintercept = corr.cutoff, color = 'blue')
  }
  if(ncol(data) > 4){
    ggvol <- ggvol + ggrepel::geom_label_repel(
      data = data[(corr >= cutoff.values[2] | corr <= cutoff.values[1]) &
                    pval <= p.cutoff],
      mapping = aes(x = corr, y = -log10(pval), label = labels, color = grp),
      size = 4.5)
  } else if(ncol(data) == 4){
    ggvol <- ggvol + ggrepel::geom_label_repel(
      data = data[(corr >= cutoff.values[2] | corr <= cutoff.values[1]) &
                    pval <= p.cutoff],
      mapping = aes(x = corr, y = -log10(pval), label = labels), size = 4.5)
  }
  ggvol <- ggvol + ggrepel::geom_label_repel(
    data = data.frame(), aes(x = -Inf, y = -log10(p.cutoff), fontface = 1,
                             label = "P-value cut-off"),
    color = "red", direction = "x", size = 4)
  if(!is.null(corr.cutoff)){
    ggvol <- ggvol + ggrepel::geom_label_repel(
      data = data.frame(), aes(x = corr.cutoff, y = Inf, fontface = 1,
                               label = title.corr.cutoff),
      color = "blue", direction = "y", size = 4)
  }
  ggvol <- ggvol +
    theme(axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey"),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.key = element_blank()) +
    labs(x = orig.cnames[2], y = paste0("-log10(",orig.cnames[3],")"),
         color = orig.cnames[4], size = orig.cnames[5])
  #Return plot removing warning about discrete variable given to alpha
  warn.handle(pattern = "Using alpha for a discrete variable is not advised.",
              print(ggvol))
}
