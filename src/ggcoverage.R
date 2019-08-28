
##IMPORTS
Imports = c("data.table","gtools","ggplot2")
lapply(Imports, library, character.only = T)

##FUNCTIONS
# ggcoverage ###################################################################

#' @description Plot an annotated barplot.
#' 
#' @param data       A \code{data.frame} containing labels in column 1, total
#'                   amount in column 2, subset amount in column 3.
#' @param round.unit An \code{integer} to specify how much decimals should be
#'                   kept when rounding percentages.
#' @param horizontal A \code{logical} to specify whether the plot should be
#'                   plotted horizontally or vertically.
#' @value a \code{type} object returned description.
#' @author Yoann Pageaud.

ggcoverage<-function(data, round.unit=2, rev.stack=FALSE, invert.percent=FALSE,
                     horizontal=FALSE,log.scaled=FALSE,decreasing.order=FALSE){
  colnames(data)[1:3]<-c("IDs","Total","Subset")
  suppressWarnings(data<-data[, lapply(.SD, na.replace, replace = 0)])
  data[, remainings:=.(Total-Subset)]
  if(invert.percent){
    data[, percents:=.(round((Subset/Total)*100,2))]  
  } else {
    data[, percents:=.(round((remainings/Total)*100,2))]
  }
  data[, percents:=.(paste0(percents,"%"))]
  if(log.scaled){
    data[Subset == 0, Subset:=1]
    data[,c("logTotal","logSubset"):=.(log10(Total),log10(Subset))]
    data[, logremainings:=.(logTotal-logSubset)]
    data<-melt(
      data,
      id.vars = c("IDs","Total","Subset","logTotal","remainings","percents"),
      measure.vars = c("logSubset","logremainings"))
  } else {
    data<-melt(data, id.vars = c("IDs","Total","percents"),
               measure.vars = c("Subset","remainings"))  
  }
  if(!rev.stack){
  data<-data[order(variable,decreasing = TRUE)] #Change order for cumsum
  }
  data[, label_ypos:=.(cumsum(value) - 0.5*value), by=IDs]
  if(log.scaled){
    data[["value.char"]]<-as.character(melt(unique(data[
      ,c("IDs","remainings","Subset"),]), id.vars = "IDs")[
        order(variable,decreasing = TRUE)]$value)
  } else {
    data[, value.char:=.(as.character(value))]
  }
  #Set the orientation of interest
  if (isTRUE(horizontal)) {
    display.count.cutoff<-0.04
    coeff.max.margin<-0.1
  } else {
    display.count.cutoff<-0.02
    coeff.max.margin<-0.05
  }
  data$IDs<-
    factor(data$IDs,levels = unique(data[
      order(Total,decreasing = decreasing.order)]$IDs))
  if(rev.stack){
    data$variable<-factor(data$variable, levels = rev(levels(data$variable)))  
  }
  #Removing duplicated strings to not display it
  data[, filter.val:= .(value - display.count.cutoff*max(data$value))]
  data[filter.val<0, value.char := " "]
  data[variable == "Subset", percents := " "]
  if(log.scaled){
    data[, "Total":=.(logTotal)]  
  }
  #Barplot
  ggcov<-ggplot(data=data, aes(x=IDs, y=value, fill=variable)) +
    geom_bar(stat = "identity") +
    geom_text(aes(y=label_ypos, label=value.char), vjust=0.5, color="white",
              size=4,fontface = "bold")
  if (isTRUE(horizontal)) {
    ggcov<-ggcov + geom_text(aes(y=Total, label=percents), hjust=-0.1) + coord_flip()
  } else {
    ggcov<-ggcov + geom_text(aes(y=Total, label=percents), vjust=-1, hjust=0.38)
  }
  return(ggcov)
}
