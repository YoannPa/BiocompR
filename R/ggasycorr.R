
#' Draws an asymmetrical pairwise correlation plot.
#'
#' @param psych_obj          A \code{psych} object containing correlation
#'                           results, generated by psych::corr.test().
#' @param pval_cutoff        A \code{numeric} comprised within [0;1] specifying
#'                           the maximum cut-off for a P-value to be considered
#'                           significant (Default: pval_cutoff = 0.05).
#' @param show_pval_below    A \code{numeric} comprised within [0;1] specifying
#'                           the maximum cut-off on P-values for results to be
#'                           displayed (Default: show_pval_below = 1).
#' @param show_abs_cor_above A \code{numeric} comprised within [0;1] specifying
#'                           the minimum cut-off on absolute correlation for
#'                           results to be displayed
#'                           (Default: show_abs_cor_above = 0).
#' @return A \code{gg} asymmetrical pairwise correlation plot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Most basic asymetrical pairwise correlation plot using New York air quality
#' # measurements (Pearson's correlation with Holm's correction applied)
#' cor_res <- psych::corr.test(airquality[1:2], airquality[3:6])
#' ggasycorr(psych_obj = cor_res)
#' # There is a significant positive correlation between solar radiation
#' # intensities and temperatures measured.
#' # There is a significant negative correlation between ozone levels and
#' # wind speeds.
#' # There is a significant positive correlation between ozone levels and
#' # temperatures measured.
#' #
#' # Do Spearman's correlation instead
#' cor_res <- psych::corr.test(
#'     airquality[1:2], airquality[3:6], method = "spearman")
#' ggasycorr(psych_obj = cor_res)
#' # The positive correlation between solar radiation intensities and
#' # temperatures measured is not significant anymore.
#' #
#' # Do Spearman's correlation and apply Benjamini–Hochberg's correction
#' cor_res <- psych::corr.test(
#'     airquality[1:2], airquality[3:6], method = "spearman", adjust = "BH")
#' ggasycorr(psych_obj = cor_res)
#' # The positive correlation between solar radiation intensities and
#' # temperatures measured is significant again with Benjamini–Hochberg's
#' # procedure
#' #
#' # Apply no P-values correction procedure
#' cor_res <- psych::corr.test(
#'     airquality[1:2], airquality[3:6], method = "spearman", adjust = "none")
#' ggasycorr(psych_obj = cor_res)
#' # Set a P-value cut-off at 0.01 instead of 0.05 (default)
#' cor_res <- psych::corr.test(
#'     airquality[1:2], airquality[3:6], method = "spearman", adjust = "BH")
#' ggasycorr(psych_obj = cor_res, pval_cutoff = 0.01)
#' # The positive correlation between solar radiation intensities and
#' # temperatures measured is no longer significant at P-value <= 0.01
#' #
#' # You can also choose to display only results with a P-value below a maximum
#' # value:
#' ggasycorr(psych_obj = cor_res, show_pval_below = 0.1)
#' # And/or results above a minimum absolute correlation:
#' ggasycorr(psych_obj = cor_res, show_abs_cor_above = 0.5)
#' # Customize X & Y labels and add title using ggplot2 grammar
#' ggasycorr(psych_obj = cor_res, pval_cutoff = 0.01) +
#'     labs(
#'         x = "Weather & Time of the year", y = "Conditions",
#'         title = "Spearman's correlation of ozone vs. wind & temperature")

ggasycorr <- function(
    psych_obj, pval_cutoff = 0.05, show_pval_below = 1, show_abs_cor_above = 0){
    # Check display filters
    if(show_pval_below > 1 | show_pval_below < 0){
        stop("Invalid value for 'show_pval_below'. Must be comprised in [0;1].")
    }
    if(show_abs_cor_above > 1 | show_abs_cor_above < 0){
        stop("Invalid value for 'show_abs_cor_above'. Must be comprised in [0;1].")
    }
    # Get corr parameters
    if(any(as.character(psych_obj$Call) == "pearson")){
        corr_name <-"Pearson's correlation"
    } else if(any(as.character(psych_obj$Call) == "spearman")){
        corr_name <-"Spearman's correlation"
    } else if(any(as.character(psych_obj$Call) == "kendall")){
        corr_name <-"Kendall's correlation"
    } else { corr_name <-"Pearson's correlation" }
    # Get adj method
    if(any(as.character(psych_obj$Call) %in% c(
        "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))){
        adj_name <- "-log10(adj. P-values)"
    } else if(any(as.character(psych_obj$Call) == "none")){
        adj_name <- "-log10(P-values)"
    } else { adj_name <- "-log10(adj. P-values)" }
    # Get corr results
    dt_adj <- as.data.table(
        x = -log10(psych_obj$p.adj), keep.rownames = "Rows")
    dt_cor <- as.data.table(x = psych_obj$r, keep.rownames = "Rows")
    dt_t <- as.data.table(x = psych_obj$t, keep.rownames = "Rows")
    dt_se <- as.data.table(x = psych_obj$se, keep.rownames = "Rows")
    # Melt dts
    dt_adj <- melt.data.table(
        data = dt_adj, id.vars = "Rows", variable.name = "Cols",
        value.name = "-log10(adj. p-values)")
    dt_cor <- melt.data.table(
        data = dt_cor, id.vars = "Rows", variable.name = "Cols",
        value.name = "Correlation")
    dt_t <- melt.data.table(
        data = dt_t, id.vars = "Rows", variable.name = "Cols",
        value.name = "t-test")
    dt_se <- melt.data.table(
        data = dt_se, id.vars = "Rows", variable.name = "Cols",
        value.name = "Std. Err.")
    # Cbind all
    dt_res <- cbind(
        dt_cor, dt_adj[, "-log10(adj. p-values)"], dt_t[, "t-test"],
        dt_se[, "Std. Err."])
    # Add significance column
    dt_res[
        `-log10(adj. p-values)` >= -log10(pval_cutoff), c(
            "P-value", "border_col") := .("Significant", "black")]
    dt_res[
        `-log10(adj. p-values)` < -log10(pval_cutoff), c(
            "P-value", "border_col") := .("Not significant", "grey75")]
    dt_res[, border_col := as.factor(border_col)]
    dt_res[, border_col := factor(
        border_col, levels = c("grey75", "black"))]
    # Set max -log10(pval)
    dt_res[`-log10(adj. p-values)` >= 3, `-log10(adj. p-values)` := 3]
    # Set display cutoff on p-values
    dt_res <- dt_res[`-log10(adj. p-values)` >= -log10(show_pval_below)]
    # Set display cutoff on absolute correlation values
    dt_res <- dt_res[abs(Correlation) >= abs(show_abs_cor_above)]
    # Remove unused levels
    dt_res[, border_col := droplevels(border_col)]
    # Plot asymetrical pairwise correlation results
    cor_asy_plt <- ggplot() +
        geom_point(
            data = dt_res,
            mapping = aes(
                x = Cols, y = Rows, fill = Correlation, color = `P-value`,
                size = `-log10(adj. p-values)`), shape = 21, stroke = 1) +
        scale_radius(
            range = c(3, 7), limits = c(0, 3), breaks = seq(1, 3, by = 1),
            labels = c("1", "2", ">= 3")) +
        scale_fill_gradientn(
            limits = c(-1, 1),
            colors = BiocompR::biopalette(name = "RColorBrewer_RdBu8")) +
        scale_color_manual(values = levels(dt_res$border_col)) +
        guides(
            color = guide_legend(
                override.aes = list(geom = "point", shape = 21, size = 4),
                order = 3),
            size = guide_legend(
                title = adj_name, override.aes = list(
                    geom = "point", shape = 21), order = 2),
            fill = guide_colorbar(
                title = corr_name, ticks.linewidth = 0.3,
                ticks.colour = "black", frame.colour = "black",
                frame.linewidth = 0.3, order = 1)) +
        theme(
            axis.text.y = element_text(size = 11, color = "black"),
            axis.text.x = element_text(
                angle = -30, hjust = 0, vjust = 0.5, size = 11,
                color = "black"),
            axis.title = element_text(size = 12),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0.1, "cm"),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(
                fill = "transparent", color = "black"),
            panel.grid = element_line(color = "black", linewidth = 0.3),
            panel.spacing = unit(0, "cm"),
            legend.text = element_text(size = 10),
            legend.key = element_rect(fill = "transparent"))
    return(cor_asy_plt)
}
