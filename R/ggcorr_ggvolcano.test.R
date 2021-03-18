#' Plots results of statistical as volcano plot.
#'
#' @param data              A \code{data.table} with 3 to 5 columns:
#'                          \itemize{
#'                           \item{column 1 - strings to be used as labels for
#'                           individual dots.}
#'                           \item{column 2 - log2 fold changes.}
#'                           \item{column 3 - p-values.}
#'                           \item{column 4 (optional) - groups to be used for
#'                           coloring the dots. If not provided the dots will be
#'                           split into two groups 1) those how meet the log2 fc.cutoff
#'                           and p.cutoff, and 2) those who not.
#'                           }
#'                          }
#' @param p.cutoff          A \code{numeric} between 0 and 1 to be used as a
#'                          maximum cut-off on p-values
#'                          (Default: p.cutoff = 0.01).
#' @param fc.cutoff        A \code{numeric} to be used as a cut-off on fold changes
#'                          (Default: fc.cutoff = 0.0).
#' @param label            A \code{character} vector to be used to define which strings of
#'                        \code{data.table}, column 1 are displayed to label the dots.
#'                          (Default: label = NULL).
#' @param group.label       A \code{boolean} to be used to define whether labeling is done on
#'                         individual (label.level = FALSE) or group level (label.level = TRUE).
#'                         Parameter is only used if \code{data table} column 4 is set.
#'                         (Default: group.label = FALSE).
#' @return A \code{gg} volcano plot object.
#' @author Verena Bitto.
#' @export

ggvolcano.test <- function(data, p.cutoff = 0.01, fc.cutoff = 0.0, label = NULL, group.label = FALSE){

   if(!is.numeric(data[[2]])){
    stop("Log2 fold changes in column 2 must be numeric.")
  }
  ##duplicate
  #Check col3
  if(is.numeric(data[[3]])){
    if(any(data[,3] < 0) | any(data[,3] > 1)){
      stop("P-values in column 3 must be comprised between 0 and 1.")}
  } else { stop("Column 3 type must be numeric.") }
  #
  ##

  ggvol <- ggplot()

  if(ncol(data) == 4){
    colnames(data) <- c("gene", "logFC", "pval", "grp")
    # duplicate - I like the shading idea, but implemented a little bit different
    data$shading <- ifelse(data$pval < p.cutoff, paste0("P-value < ", p.cutoff), paste0("P-value >= ", p.cutoff))

    ggvol <- ggplot(data = data, aes(x = logFC, y = -log10(pval)))

    if(!group.label){
      ggvol <- ggvol +
        ggrepel::geom_text_repel(min.segment.length = 0,
                                 aes(label = ifelse(gene %in% label, gene,'')))
    } else if(group.label) {
      ggvol <-  ggvol +
        ggrepel::geom_text_repel(min.segment.length = 0,
                                 aes(label = ifelse(grp %in% label, gene,'')))
    }
  } else {
    colnames(data) <- c("gene", "logFC", "pval")
    cutoff <- data$pval < p.cutoff & abs(data$logFC) >= fc.cutoff
    data$grp <- ifelse(cutoff,
                       paste0("P-value < ", p.cutoff, " AND \n log2FC >= |", fc.cutoff, "|"),
                       paste0("P-value >= ", p.cutoff, " OR \n log2FC < |", fc.cutoff, "|"))
    # duplicate - I like the shading idea, but implemented a little bit different
    data$shading <- ifelse(data$pval < p.cutoff, paste0("P-value < ", p.cutoff), paste0("P-value >= ", p.cutoff))

    ggvol <- ggplot(data = data, aes(x = logFC, y = -log10(pval)))

    if(missing(label)){
      ggvol <-  ggvol +
        ggrepel::geom_text_repel(min.segment.length = 0,
                                 aes(label = ifelse(cutoff, gene,'')))
    } else {
      ggvol <-  ggvol +
        ggrepel::geom_text_repel(min.segment.length = 0,
                                 aes(label = ifelse(gene %in% label, gene,'')))
    }

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
    labs(x = "log2 foldchange", col = "Significance groups", alpha = "P-value") +
    geom_point(aes(color = grp, fill = grp, alpha=shading)) +
      scale_alpha_discrete(range = c(0.9, 0.2)) +
      scale_fill_discrete(guide = FALSE)

  return(ggvol)
}
