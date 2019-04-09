
##IMPORTS
Imports = c('psych','ggplot2','Hmisc','ggrepel','gtools')
lapply(Imports, library, character.only = T)

##FUNCTIONS

# ggEigenvector ################################################################

#' Create an Eigenvector Plot using GGplot.
#' 
#' @param dataframe   A \code{data.frame} containg:
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

ggEigenvector<- function(dataframe, xcol, ycol, eigenvalues, colors, title){
  ggplot(data=dataframe) +
    ggtitle(title) +
    geom_point(aes_string(x = xcol, y = ycol, color = "Groups")) +
    geom_segment(aes_string(xend = xcol, yend = ycol, color = "Groups"),
                 x = 0, y = 0, size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) +
    geom_label_repel(aes_string(x = xcol, y = ycol, label = "Labels")) +
    scale_color_manual(values=colors) +
    xlab(paste0("Eigenvector ",gsub("[^0-9.]","",xcol),
                " (Variance accounted = ",
                round(eigenvalues[as.integer(gsub("[^0-9.]","",xcol))]
                      /nrow(dataframe)*100,2),
                "%)")) +
    ylab(paste0("Eigenvector ",gsub("[^0-9.]","",ycol),
                " (Variance accounted = ",
                round(eigenvalues[as.integer(gsub("[^0-9.]","",ycol))]
                      /nrow(dataframe)*100,2),
                "%)")) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 13),
          legend.text = element_text(size=12)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_vline(xintercept = 0, linetype="dashed")
}

# EVA - EigenVector Analysis ###################################################

#' From a Correlation test return Eigenvectors, Principal components scores and
#' principal components correlations with the data.  
#' 
#' @param data        A \code{matrix} or a \code{data.frame} containing
#'                    variables by columns and values to be used for the
#'                    correlation test.
#' @param use         A \code{character} to specify how to handle missing values
#'                    when calculating a correlation. Possible types are
#'                    'pairwise' and 'complete'. 'pairwise' is the default value
#'                    and will do pairwise deletion of cases. 'complete' will
#'                    select just complete cases.
#' @param method      The correlation method to use as a \code{character}
#'                    matching one of these: 'pearson','spearman','kendall'.
#' @param adjust      A \code{character} specifying what adjustment for multiple
#'                    tests should be used. Possible values are: "holm",
#'                    "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#'                    "none".
#' @param var.min     A \code{double} setting the minimum variance accountable
#'                    for an eigenvector to be considered in the plots
#'                    generated.
#' @param groups      A \code{character} vector of groups to which variable
#'                    belongs for eigenvector annotations. The length of this
#'                    vector has to match the number/variables of columns in the
#'                    data.
#' @param colors      A \code{character} vector of colors for the eigenvectors.
#'                    The length of this vector has to match the number of
#'                    different groups existing.
#' @return A \code{list} object containing the principal components
#'         correlations, a list of the Eigenvector Plots generated and principal
#'         components scores of the data.
#' @author Yoann Pageaud.

EVA<-function(data, use = "pairwise", method = "pearson", adjust = "none",
              var.min = 0.01, groups = as.character(seq(ncol(data))),
              colors = rainbow(n = ncol(data))){
  #Compute a correlation test on the data
  M<-corr.test(x = data, use = use, method = method, adjust = adjust)$r
  #Get eigenvalues from the correlation matrix
  eigvals <- eigen(M)$values
  
  #Get percentage of variance accounted by each eigen values
  var.acc<-eigvals/length(eigvals)
  #How many eigenvalues are above the minimum threshold for variance accounted
  true.eigvals<-length(var.acc[var.acc >= var.min])
  
  #Get eigenvectors
  eigvects <- as.data.frame(eigen(M)$vectors[,c(1:true.eigvals)])
  dframe<-as.data.frame(cbind(eigvects,"Groups" = groups,
                              "Labels" = colnames(data)))
  dframe$Groups<-factor(dframe$Groups,levels = unique(dframe$Groups))
  
  #Create all possible combinations
  combs<-expand.grid(seq(true.eigvals),seq(true.eigvals))
  #Remove duplicate and order
  combs<-combs[combs$Var1 != combs$Var2,]
  combs<-combs[order(combs[,1]),]

  my_cols<-rev(rev(colnames(dframe))[-c(1,2)])
  
  #Generate all Eigenvector Plots
  list_EVplots<-lapply(my_cols, function(col1){
    lapply(my_cols, function(col2){
      
      if(col1 == col2){
        return(NULL)
      }
      
      ggEigenvector(dataframe = dframe, eigenvalues = eigvals, xcol = col1,
                    ycol = col2, colors = colors,
                    title = paste("Eigenvector Plot - Pairwise",
                                  capitalize(method),"correlation with", adjust,
                                  "adjustment"))
    })
  })
  
  #Flatten list
  list_EVplots<-unlist(list_EVplots,recursive = F)
  #Remove NULL elements
  list_EVplots<-Filter(Negate(is.null), list_EVplots)
  names(list_EVplots)<-paste(combs[,1],"&",combs[,2])

  #Scaling the matrix values.
  data<-scale(data[complete.cases(data),])
  #Matricial product of scaled values and eigenvectors.
  pca.scores<- data %*% eigen(M)$vectors
  colnames(pca.scores)<-paste("EV",seq(ncol(data)), sep="")
  #Get correlation between samples and Principal components.
  PC.cor<-corr.test(data,pca.scores,use = use, method = method)$r
  return(list("PC.cor" = PC.cor, "EV.plots" = list_EVplots,
              "PC.scores" = as.data.frame(pca.scores)))
}