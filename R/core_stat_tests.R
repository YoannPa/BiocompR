
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
#' @export
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
#' @export
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
#' @export
#' @keywords internal

ls.psych.as.dt <- function(psych.list){
    dt <- data.table::rbindlist(l = lapply(X = psych.list, FUN = function(i){
        data.table::data.table(cor = i[["r"]], pvalue = i[["p"]],
                               nsample = i[["n"]], stderr = i[["se"]])
    }), use.names = TRUE, fill = TRUE)
    dt <- cbind(var.name = names(psych.list), dt)
    return(dt)
}

#' Checks if a function exists.
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
                ls.fun <- as.character(lsf.str(p))
                if(any(as.character(lsf.str(p)) == fun.str)){
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
