
#' Calculate the variance accounted of an eigenvalue.
#'
#' @param colname          A \code{character} matching a column name.
#' @param eigen.values     A \code{numeric} vector of all eigenvalues.
#' @return A \code{character} containg the value of the eigenvalue variance
#'         accounted for its associated eigen vector.
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

var.accounted <- function(colname, eigen.values){
    vect.num <- gsub(pattern = "[^0-9.]", replacement = "", x = colname)
    var.accounted <- round(
        eigen.values[as.integer(vect.num)]/length(eigen.values)*100, 2)
    paste0("Eigenvector ", vect.num,
           " (Variance accounted = ", var.accounted, "%)")
}

#' Creates an eigenvector plot using ggplot2.
#'
#' @param data        A \code{data.frame} containg:
#'                    N eigen vectors, one vector by column, with column names
#'                    numbered from "V1" to "VN",
#'                    a column "Groups" for the variable groups,
#'                    and a column "Labels" for the variable names.
#' @param xcol        A \code{character} matching the column name of an
#'                    eigenvector, for which the values will be used as X axis
#'                    coordinates.
#' @param ycol        A \code{character} matching the column name of an
#'                    eigenvector, for which the values will be used as Y axis
#'                    coordinates.
#' @param eigenvalues A \code{double} vector containing the eigenvalues of the
#'                    matrix.
#' @param colors      A \code{character} vector of colors for the eigenvectors.
#'                    The length of this vector has to match the number of
#'                    different groups existing.
#' @param title       A \code{character} that will be used as a title for the
#'                    plot.
#' @return A \code{gg} plot of the variables following the 2 eigenvectors
#'         selected.
#' @author Yoann Pageaud.
#' @export

ggeigenvector <- function(
    #TODO: Replace aes_string() by aes()
    data, xcol, ycol, eigenvalues, colors, label = TRUE, title){
    ev.plot <- ggplot2::ggplot(data = data) +
        ggplot2::ggtitle(title) +
        ggplot2::geom_point(
            ggplot2::aes_string(x = xcol, y = ycol, color = "Groups")) +
        ggplot2::geom_segment(ggplot2::aes_string(
            xend = xcol, yend = ycol, color = "Groups"), x = 0, y = 0, size = 1,
            arrow = grid::arrow(length = grid::unit(0.3, "cm"))) +
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::labs(
            x = BiocompR::var.accounted(
                colname = xcol, eigen.values = eigenvalues),
            y = BiocompR::var.accounted(
                colname = ycol, eigen.values = eigenvalues)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 13),
                       legend.text = ggplot2::element_text(size = 12)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
    if(label){
        ev.plot <- ev.plot + ggrepel::geom_label_repel(
            ggplot2::aes_string(x = xcol, y = ycol, label = "Labels"))
    }
    ev.plot
}

#' Computes eigenvectors, PC scores and correlations from a correlation test.
#'
#' @param data        A \code{matrix} or a \code{data.frame} containing
#'                    variables by columns and values to be used for the
#'                    correlation test by rows.
#' @param use         A \code{character} to specify how to handle missing values
#'                    when calculating a correlation. Possible values are
#'                    'pairwise' and 'complete'. 'pairwise' is the default value
#'                    and will do pairwise deletion of cases. 'complete' will
#'                    select just complete cases.
#' @param method      The correlation method to use as a \code{character}
#'                    matching one of these: 'pearson','spearman','kendall'.
#' @param adjust      A \code{character} specifying what adjustment for multiple
#'                    tests should be used. Possible values are: "holm",
#'                    "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and
#'                    "none".
#' @param var.min     A \code{double} setting the minimum variance accountable
#'                    for an eigenvector to be considered in the plots
#'                    generated.
#' @param groups      A \code{character} vector of groups to which variables
#'                    belong for eigenvector annotations. The length of this
#'                    vector has to match the number/variables of columns in the
#'                    data.
#' @param colors      A \code{character} vector of colors for the eigenvectors.
#'                    The length of this vector has to match the number of
#'                    different groups existing.
#' @return A \code{list} containing the principal components correlations,
#'         a list of the Eigenvector Plots generated and principal components
#'         scores of the scaled data.
#' @author Yoann Pageaud.
#' @export

EVA <- function(data, use = "pairwise", method = "pearson", adjust = "none",
                var.min = 0.01, groups = as.character(seq(ncol(data))),
                colors = grDevices::rainbow(n = ncol(data)), label = TRUE){
    #Compute a correlation test on the data
    M <- psych::corr.test(
        x = data, use = use, method = method, adjust = adjust)$r
    #Get eigenvalues from the correlation matrix
    eigvals <- eigen(M)$values
    #Get percentage of variance accounted by each eigen values
    var.acc <- eigvals/length(eigvals)
    #How many eigenvalues are above the minimum threshold for variance accounted
    true.eigvals <- length(var.acc[var.acc >= var.min])
    if(true.eigvals == 0){
        stop(paste("No eigenvalues accounting for more than the var.min",
                   "minimum of the variance. Please set a lower var.min."))
    }
    #Get eigenvectors
    eigvects <- as.data.frame(eigen(M)$vectors[, seq(true.eigvals)])
    dframe <- as.data.frame(cbind(eigvects, "Groups" = as.factor(groups),
                                  "Labels" = colnames(data)))
    #Create all possible combinations
    combs <- expand.grid(seq(true.eigvals), seq(true.eigvals))
    #Remove duplicate and order
    combs <- combs[combs$Var1 != combs$Var2, ]
    combs <- combs[order(combs[, 1]), ]
    my_cols <- colnames(dframe)[!colnames(dframe) %in% c("Groups", "Labels")]
    #Generate all Eigenvector Plots
    list_EVplots <- lapply(my_cols, function(col1){
        lapply(my_cols, function(col2){
            if(col1 == col2){ return(NULL) } else {
                BiocompR::ggeigenvector(
                    data = dframe, eigenvalues = eigvals, xcol = col1,
                    ycol = col2,colors = colors, label = label, title = paste(
                        "Eigenvector Plot - Pairwise", method,
                        "correlation with", adjust, "adjustment"))
            }
        })
    })
    #Flatten list
    list_EVplots <- unlist(list_EVplots, recursive = FALSE)
    #Remove NULL elements
    list_EVplots <- Filter(Negate(is.null), list_EVplots)
    names(list_EVplots) <- paste(combs[,1], "&", combs[, 2])
    #Scaling the matrix values.
    data <- scale(data[stats::complete.cases(data), ])
    #Matricial product of scaled values and eigenvectors.
    pca.scores <- data %*% eigen(M)$vectors
    colnames(pca.scores) <- paste("EV", seq(ncol(data)), sep = "")
    #Get correlation between samples and Principal components.
    PC.cor <- psych::corr.test(data, pca.scores, use = use, method = method)$r
    return(list("PC.cor" = PC.cor, "EV.plots" = list_EVplots,
                "PC.scores" = as.data.frame(pca.scores)))
}
