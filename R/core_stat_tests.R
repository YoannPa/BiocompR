
#' Converts to uppercase the first letter in a character string.
#'
#' @param x A \code{character}.
#' @return A \code{character} with the first letter converted to uppercase.
#' @author Yoann Pageaud.
#' @examples simplecap(x = "test")
#' @references \href{https://www.R-project.org/}{R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.}
#' @keywords internal

simplecap <- function(x){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
}

#' Returns the alias of a correlation test metric.
#'
#' @param metric A \code{character} matching a specific metric return by the
#'               \link[psych]{corr.test} function.
#' @return A \code{character} being the alias of the metric given.
#' @author Yoann Pageaud.
#' @examples
#' # Get the alias of the metric 'r'
#' metric_alias(metric = "r")
#' @references \href{https://www.scholars.northwestern.edu/en/publications/psych-procedures-for-personality-and-psychological-research}{William R Revelle, psych: Procedures for Personality and Psychological Research. Northwestern University, Evanston, Illinois, USA (2017).}
#' @keywords internal

metric_alias <- function(metric = "r"){
    if(metric == "r"){
        alias <- "correlations"
    } else if(metric == "n"){
        alias <- "n. cases"
    } else if(metric == "t"){
        alias <- "t-test values"
    } else if(metric == "p"){
        alias <- "-log10(P-values)"
    } else if(metric == "se"){
        alias <- "standard error"
    }
    return(alias)
}

#' Extracts Kolmogorov-Smirnov test results and return them as a matrix.
#'
#' @param df.ks.tests A \code{data.frame} containing raw results from multiple
#'                    Kolmogorov-Smirnov tests.
#' @param statistic   A \code{character} specifying the type of statistic to
#'                    retrieve from the test\cr
#'                    (Supported: statistic = c('n','stat','p')).
#' @return A \code{matrix} containing the KS test's statistic values of
#'         interest.
#' @author Yoann Pageaud.
#' @keywords internal

get.ks.stat <- function(table_combinations, df.ks.tests, statistic){
    ks.res <- df.ks.tests[[statistic]]
    tbl.ks <- cbind(table_combinations, ks.res)
    #Cast a molten tables into a matrix
    mat.ks <- data.table::dcast(
        data = tbl.ks, formula = Var1 ~ Var2, value.var = "ks.res")
    rownames(mat.ks) <- mat.ks$Var1
    mat <- as.matrix(mat.ks[,-1])
    return(mat)
}

#' Computes a Kolmogrov-Smirnov test between all columns of a data.frame.
#'
#' @param data      A \code{data.frame} of numerical values.
#' @param statistic A \code{character} specifying the type of statistic to
#'                  retrieve from the test\cr
#'                  (Supported: statistic = c('n','stat','p')).
#' @param ncores    An \code{integer} to specify the number of cores/threads to
#'                  be used to parallel-run tests.
#' @return A \code{list} containing a matrix of the statistic values, a
#'         data.frame of the pairwise KS raw test results, and a table of
#'         pairwise combinations.
#' @author Yoann Pageaud.
#' @keywords internal

pairwise.ks <- function(data, statistic, ncores){
    table_combinations <- expand.grid(colnames(data), colnames(data))
    List_ks.tests <- parallel::mclapply(
        seq(nrow(table_combinations)), mc.cores = ncores, function(i){
            #Compute KS test
            ks.res <- stats::ks.test(x = data[,table_combinations[i,1]],
                                     y = data[,table_combinations[i,2]])
            #Create table with all statistics of the KS test
            data.frame('n' = nrow(data[stats::complete.cases(
                data[, c(table_combinations[i, 1], table_combinations[i, 2])]),
                ]), 'stat' = ks.res$statistic, 'p' = ks.res$p.value,
                row.names = NULL)
        })
    #Create Table
    df.ks.tests <- do.call("rbind", List_ks.tests)
    #Get matrix of the stat of interest
    mat <- BiocompR::get.ks.stat(
        table_combinations = table_combinations, df.ks.tests = df.ks.tests,
        statistic = statistic)
    return(list("res.statistic" = mat, "res.test" = df.ks.tests,
                "table_combinations" = table_combinations))
}

#' Converts a list of type 'psych' objects into a data.table.
#'
#' @param psych.list A \code{psych} list of correlation results.
#' @return A \code{data.table} of the correlation results with columns:
#'         \itemize{
#'          \item{'var.name' - the name of the psych in the list.}
#'          \item{'cor' - the value of the correlation test.}
#'          \item{'pvalue' - the P-value associated to the correlation test.}
#'          \item{'nsample' - the number of samples used for the pairwise
#'          correlation.}
#'          \item{'stderr' - the standard error associated to the correlation
#'          test.}
#'         }
#' @author Yoann Pageaud.
#' @references \href{https://www.scholars.northwestern.edu/en/publications/psych-procedures-for-personality-and-psychological-research}{William R Revelle, psych: Procedures for Personality and Psychological Research. Northwestern University, Evanston, Illinois, USA (2017).}
#' @keywords internal

ls.psych.as.dt <- function(psych.list){
    dt <- data.table::rbindlist(l = lapply(X = psych.list, FUN = function(i){
        data.table::data.table(cor = i[["r"]], pvalue = i[["p"]],
                               nsample = i[["n"]], stderr = i[["se"]])
    }), use.names = TRUE, fill = TRUE)
    dt <- cbind(var.name = names(psych.list), dt)
    return(dt)
}

#' Checks if a function exists and package of origin.
#'
#' @param fun        A \code{character} string matching a function you are
#'                   looking for. The syntax must either be 'function' or
#'                   'package::function'.
#' @param param.name A \code{character} specifying the name of the parameter to
#'                   which the function is related. This name is used in error
#'                   messages.
#' @param ncores     An \code{integer} to specify the number of cores/threads to
#'                   be used to parallel-compute distances for dendrograms.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Looking for the function "abbreviate"
#' check_fun(fun = "abbreviate") #Result should be: "base::abbreviate"

check_fun <- function(fun, param.name = "fun", ncores = 1){
    fun.str <- unlist(strsplit(x = fun, "::"))
    if(length(fun.str) == 1){
        # Look first through functions of loaded packages
        ls.packfun <- parallel::mclapply(
            X = search(), mc.cores = ncores, FUN = function(p){
                ls.fun <- as.character(utils::lsf.str(p))
                if(any(ls.fun == fun.str)){
                    paste(unlist(strsplit(x = p, split = ":"))[2], fun.str,
                          sep = "::")
                } else { NULL }
            })
        packfun <- unlist(ls.packfun[!sapply(ls.packfun, is.null)])
        if(length(ls.packfun) != 1){
            # Look through all packages installed
            ls.packfun2 <- parallel::mclapply(
                X = .packages(all.available = TRUE), mc.cores = ncores,
                FUN = function(p){
                    if(any(getNamespaceExports(p) == fun.str)){
                        paste(p, fun.str, sep = "::")
                    } else { NULL }
                })
            ls.packfun2 <- unlist(ls.packfun2[!sapply(ls.packfun2, is.null)])
            if(length(ls.packfun2) == 1){
                packfun <- ls.packfun2
            } else {
                stop(paste0(
                    "Too many matches for the function in '", param.name,
                    "'. Please precise the package of the function with the",
                    " syntax 'package::function'"))
            }
        }
    } else if(length(fun.str) == 2){
        if(any(.packages(all.available = TRUE) == fun.str[1])){
            if(any(getNamespaceExports(fun.str[1]) == fun.str[2])){
                packfun <- paste(fun.str, collapse = "::")
            } else {
                stop(paste(
                    "Function", fun.str[2], "from package", fun.str[1],
                    "not found."))
            }
        } else { stop(paste("Package", fun.str[1], "not found.")) }
    } else {
        stop(paste0(
            param.name, " must be a 'function', or a 'package::function'."))
    }
    return(packfun)
}

#' Tests association of an annotation with another one or with a PC.
#'
#' @param x           A \code{vector} of numerics or factors containing data
#'                    from one annotation.
#' @param y           A \code{vector} of numerics or factors containing data
#'                    from a PCA principal component or another annotation.
#' @param perm.matrix A \code{matrix} containing permutations of the order of
#'                    rows, as integers, based on annotations length. The
#'                    permutation matrix can optionally be used to test the
#'                    significance of a correlation value between 2 annotations
#'                    or an annotation and a PCA principal component.
#' @return A \code{data.table} containing all results from an association test
#'         between x and y.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Most basic statistic test between 2 vectors of integers
#' test.annots(x = 1:15, y = 24:10)

test.annots <- function (x, y, perm.matrix = NULL){
    # Set nominal parameter
    nominal <- TRUE
    #If annotation x or y is a date then convert it into integers
    if(class(x) == "Date"){ x <- as.integer(x) }
    if(class(y) == "Date"){ y <- as.integer(y) }
    inds <- which(!(is.na(x) | is.na(y)))
    if (length(inds) < 2) {
        nominal <- FALSE
        warning("Not enough common values.")
    }
    x <- x[inds]
    if (is.factor(x)) {
        x <- as.factor(as.character(x))
        if (nlevels(x) < 2) {
            nominal <- FALSE
            warning("Not enough categories.")
        }
    }
    y <- y[inds]
    if (is.factor(y)) {
        y <- as.factor(as.character(y))
        if (nlevels(y) < 2) {
            nominal <- FALSE
            warning("Not enough categories.")
        }
    }
    if(nominal){ # If a stat. test can be run
        # Retrieve test p-value
        get.p <- function(expr) {
            tryCatch(suppressWarnings(expr$p.value), error = function(er){
                as.double(NA)
            })
        }
        if(is.factor(x)){
            if(is.factor(y)) {
                # If both vectors are factors do Non-parametric Fisher test with
                # 50 000 replicates for Monte Carlo test.
                simulate <- (nlevels(x) > 2 || nlevels(y) > 2)
                dt.res <- data.table::data.table(
                    test = "Fisher", correlation = NA, pvalue = get.p(
                        fisher.test(
                            x, y, conf.int = FALSE, simulate.p.value = simulate,
                            B = 50000)))
            } else if (nlevels(x) == 2) {
                # If Y is a numeric vector and X has only 2 possible values do
                # Non-parametric Wilcoxon–Mann–Whitney test
                # aka Wilcoxon rank-sum test
                values <- tapply(y, x, identity)
                dt.res <- data.table::data.table(
                    test = "Wilcoxon-Mann-Whitney", correlation = NA,
                    pvalue = get.p(wilcox.test(
                        values[[1]], values[[2]], alternative = "two.sided"))
                )
            } else {
                # If X has multiple levels & Y is a numeric vector do
                # Non-parametric Kruskal-Wallis test
                dt.res <- data.table::data.table(
                    test = "Kruskal-Wallis", correlation = NA, pvalue = get.p(
                        kruskal.test(y, x)))
            }
        } else if(is.factor(y)) {
            if (nlevels(y) == 2) {
                # If X is a numeric vector and Y has only 2 possible values do
                # Non-parametric Wilcoxon–Mann–Whitney test
                # aka Wilcoxon rank-sum test
                values <- tapply(x, y, identity)
                dt.res <- data.table::data.table(
                    test = "Wilcoxon-Mann-Whitney", correlation = NA,
                    pvalue = get.p(wilcox.test(
                        values[[1]], values[[2]], alternative = "two.sided")))
            } else {
                # If X is a numeric vector and Y has multiple levels do
                # Non-parametric Kruskal-Wallis test
                dt.res <- data.table::data.table(
                    test = "Kruskal-Wallis", correlation = NA, pvalue = get.p(
                        kruskal.test(x, y)))
            }
        } else {
            # If both X & Y are numeric vectors do a correlation test
            N <- length(inds)
            if(is.null(perm.matrix)){
                cor_res <- cor(x, y) # Calculate Pearson correlation
                # Get T statistic
                t_stat <- (cor_res*sqrt(N-2))/(sqrt(1-cor_res^2))
                dt.res <- data.table::data.table(
                    test = "Pearson Corr.", correlation = cor_res,
                    pvalue = 2*pt(-abs(t_stat), N-2))
            } else {
                values <- apply(X = perm.matrix, MARGIN = 2, FUN = function(i){
                    cor(x[i[i <= N]], y)
                })
                abs_val <- abs(values)
                # Correlation p-value here is the proportion of times where the
                # 1st correlation value is inferior or equal to calculated
                # correlations in all permutations. 10 000 permutations to make
                # sure that 1 result out of 10 000 will be significant enough if
                # correlated.
                dt.res <- data.table::data.table(
                    test = "Pearson Corr.", correlation = values[1],
                    pvalue = mean(abs_val[1] <= abs_val))
            }
        }
    } else { # Empty data.table result when no test possible
        dt.res <- data.table::data.table(
            test = NA, correlation = NA, pvalue = NA)
    }
    return(dt.res)
}

#' Prepares annotations to be tested for associations.
#'
#' @param annot.table A \code{data.frame} containing all annotations, 1
#'                    annotation per column.
#' @param verbose     A \code{logical} to display information about the
#'                    step-by-step processing of the data if TRUE
#'                    (Default: verbose = FALSE).
#' @return A \code{list} containing updated annotations, the number of
#'         annotations available, and the annotation table itself.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create mtcars annotation table
#' mtcars_annot <- data.table::as.data.table(
#'     mtcars[, c("cyl", "vs", "am", "gear", "carb")], keep.rownames = "cars")
#' # Prepare the annotations for testing
#' ls_annot <- prepare_annot_asso(annot.table = mtcars_annot, verbose = TRUE)

prepare_annot_asso <- function(annot.table, verbose = FALSE){
    annots <- lapply(X = annot.table, FUN = function(x) {
        # Exclude columns containing the same value for all rows
        # (except when missing)
        if(length(unique(na.omit(x))) < 2){ return(NULL) }
        # Convert as factor character variables
        if(is.character(x)){ x <- as.factor(x) }
        # Convert as factor TRUE/FALSE variables
        if(is.logical(x)){ x <- as.factor(x) }
        if(is.factor(x)){ if(anyDuplicated(na.omit(x)) == 0) { return(NULL) } }
        x
    })
    annots <- annots[!vapply(
        X = annots, FUN = is.null, FUN.VALUE = logical(length = 1L))]
    n.annot <- length(annots)
    if(n.annot == 0){
        stop("No suitable annotation found for association tests.")}
    if(verbose){
        cat(c("Testing the following annotations for associations:\n",
              paste(names(annots), collapse = ", "), ".\n"))
    }
    res <- list(
        "annotations" = annots, "n.annot" = n.annot,
        "annot.table" = annot.table)
    return(res)
}

#' Tests associations between a set of annotations and PCs from a prcomp object.
#'
#' @param annot.table A \code{data.frame} containing all annotations, 1
#'                    annotation per column.
#' @param prcomp.res A PCA result of classes \code{prcomp} or
#'                   \code{irlba_prcomp} resulting from stats::prcomp() or
#'                   irlba::prcomp_irlba().
#' @param perm.count An \code{integer} specifying the number of permutations to
#'                   realize on a vector, for the permutations matrix
#'                   initialization, to be used for calculating the significance
#'                   of a correlation test (Default: perm.count = 10000).
#' @param max.PCs    An \code{integer} specifying the maximum number of
#'                   principal components to consider for association tests with
#'                   annotations (Default: max.PCs = 8).
#' @param verbose    A \code{logical} to display information about the
#'                   step-by-step processing of the data if TRUE
#'                   (Default: verbose = FALSE).
#' @return A \code{data.table} containing all association test results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create mtcars annotation table
#' mtcars_annot <- data.table::as.data.table(
#'     mtcars[, c("cyl", "vs", "am", "gear", "carb")], keep.rownames = "cars")
#' # Create mtcars matrix with scaled and transposed data
#' mat_mtcars <- scale(as.matrix(
#'     mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")]))
#' # Compute PCA on mtcars matrix
#' pcs_mtcars <- prcomp(x = mat_mtcars)
#' # Test association between annotations and principal components
#' asso_res <- test_asso_annot_pc(
#'     annot.table = mtcars_annot, prcomp.res = pcs_mtcars, verbose = TRUE)
#' # Table summarizing tests' results
#' asso_res

test_asso_annot_pc <- function(
    annot.table, prcomp.res, perm.count = 10000, max.PCs = 8, verbose = FALSE){
    prep.annot.table <- prepare_annot_asso(
        annot.table = annot.table, verbose = verbose)
    annots <- prep.annot.table$annotations
    n.annot <- prep.annot.table$n.annot
    annot.table <- prep.annot.table$annot.table

    if(perm.count != 0 && (
        (!is.null(prcomp.res)) || sum(!vapply(
            X = annots, FUN = is.factor,
            FUN.VALUE = logical(length = 1L))) >= 2)) {
        # Create the random permutation matrix
        perm.matrix <- mapply(
            FUN = sample, rep(nrow(annot.table), times = perm.count))
        perm.matrix[, 1] <- 1:nrow(perm.matrix)
    } else {
        warning("Cannot initialize the permutations matrix.")
        perm.matrix <- NULL
    }

    pc.association.count <- max.PCs
    # Get PCs % variance explained
    PCA_metrics <- BiocompR::prepare_pca_data(
        prcomp.res = prcomp.res, dt.annot = annot.table, PCs = 1:max.PCs,
        scale = 1)
    dpoints <- prcomp.res$x
    if(!is.null(dpoints)){
        if (ncol(dpoints) > pc.association.count) {
            # Reduce principal components coordinates to the maximum number of
            # dimensions wanted
            dpoints <- dpoints[, 1:pc.association.count]
        }
        # Test all annotations against all PCs
        ls_allres <- lapply(X = seq(n.annot), FUN = function(i){
            if(verbose){ cat("Testing association of", names(annots)[i], "&\n") }
            ls_tres <- lapply(X = seq(ncol(dpoints)), FUN = function(j){
                if(verbose){ cat("\t", colnames(dpoints)[j], "\n") }
                t.result <- BiocompR::test.annots(
                    x = annots[[i]], y = dpoints[, j],
                    perm.matrix = perm.matrix)
                t.result[, c("annotation", "PC", "var.explained") := .(
                    names(annots)[i], colnames(dpoints)[j],
                    (PCA_metrics$var.explained*100)[j])]
                t.result
            })
            data.table::rbindlist(l = ls_tres)
        })
        # Rbind all results
        dt_allres <- data.table::rbindlist(l = ls_allres)
        dt_allres[, log_trans_pval := -log10(pvalue)]
        rm(ls_allres)
        # Convert PC as factor to keep the right order
        dt_allres[, PC := as.factor(x = PC)]
        dt_allres[, PC := factor(
            x = PC, levels = paste0("PC", seq(pc.association.count)))]
    }
    rm(dpoints)
    # Return association test results
    return(dt_allres)
}

#' Tests associations between all annotations stored in a data.frame.
#'
#' @param annot.table A \code{data.frame} containing all annotations, 1
#'                    annotation per column.
#' @param perm.count An \code{integer} specifying the number of permutations to
#'                   realize on a vector, for the permutations matrix
#'                   initialization, to be used for calculating the significance
#'                   of a correlation test (Default: perm.count = 10000).
#' @param verbose    A \code{logical} to display information about the
#'                   step-by-step processing of the data if TRUE
#'                   (Default: verbose = FALSE).
#' @return A \code{data.table} containing all association test results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create mtcars annotation table
#' mtcars_annot <- data.table::as.data.table(
#'     mtcars[, c("cyl", "vs", "am", "gear", "carb")], keep.rownames = "cars")
#' # Test association between all annotations
#' asso_res <- test_asso_all_annot(
#'     annot.table = mtcars_annot, perm.count = 10000, verbose = TRUE)
#' # Table summarizing tests' results
#' asso_res

test_asso_all_annot <- function(
    annot.table, perm.count = 10000, verbose = FALSE){
    prep_res <- BiocompR::prepare_annot_asso(
        annot.table = annot.table, verbose = verbose)
    annots <- prep_res$annotations
    n.annot <- prep_res$n.annot
    annot.table <- prep_res$annot.table
    if(perm.count != 0 && sum(!vapply(
        X = annots, FUN = is.factor, FUN.VALUE = logical(length = 1L))) >= 2) {
        # Create the random permutation matrix
        perm.matrix <- mapply(
            FUN = sample, rep(nrow(annot.table), times = perm.count))
        perm.matrix[, 1] <- 1:nrow(perm.matrix)
    } else {
        warning("Cannot initialize the permutations matrix.")
        perm.matrix <- NULL
    }
    if (n.annot > 1) {
        # Create matrix of tests combinations
        test_matrix <- utils::combn(x = names(annots), m = 2)
        # Test association between all annotations available
        ls_annotres <- apply(X = test_matrix, MARGIN = 2, FUN = function(i){
            if(verbose){ cat("Testing association of", i[1], "&", i[2], "\n") }
            t.result <- BiocompR::test.annots(
                x = annots[[i[1]]], y = annots[[i[2]]],
                perm.matrix = perm.matrix)
            t.result[, c("annotation1", "annotation2") := .(i[1], i[2])]
        })
        dt_annotres <- data.table::rbindlist(l = ls_annotres)
        #Duplicate results for the full table
        dt_annotres_bis <- data.table::copy(dt_annotres)
        data.table::setnames(
            x = dt_annotres_bis, old = "annotation1", new = "annotation2_new")
        data.table::setnames(
            x = dt_annotres_bis, old = "annotation2", new = "annotation1")
        data.table::setnames(
            x = dt_annotres_bis, old = "annotation2_new", new = "annotation2")
        # Rbind all results
        dt_annotres <- rbind(dt_annotres, dt_annotres_bis, use.names = TRUE)
        rm(dt_annotres_bis, perm.matrix)
        dt_annotres[, log_trans_pval := -log10(pvalue)]
    } else {
        stop(paste(
            "only 1 annotation usable. Cannot compute association of an",
            "annotation against itself alone."))
    }
    return(dt_annotres)
}
