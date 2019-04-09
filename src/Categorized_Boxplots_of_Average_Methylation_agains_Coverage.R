##Categorized Boxplots of Average Methylation against Coverage 

setwd("/media/yoann/3CAD87DD271F7BEC/PhD_2018-2021/Hematopoiesis_DNMT1_Hypomorph")


## L I B R A R I E S ###########################################################

library(ggplot2)
library(RColorBrewer)


## F U N C T I O N S ###########################################################

#' Boxplot of categories
#' 
#' Creates a boxplot per category.
#' 
#' @param dframe   Two-column \code{data.frame} containing values to categorize
#' and values to show the distributions of.
#' @param cat.step Step size for the categorization. The first category starts
#' with value \code{0} (if )
#' @param cat.max  Maximum bin to consider for categories. Any values lying
#' outside this bin are put into an extra bin
#'                 named \code{"over ..."}.
#' @return The newly created plot as an instance of \code{\link{ggplot}}.
#' 
#' @author Yassen Assenov, Yoann Pageaud
#' @export

plot.box.categories <- function(dframe, cat.step = 10L, cat.max = 10L) {
  #Create Categories
  dframe$`category` <- as.integer(floor(dframe[, 1] / cat.step))
  #Assign CpGs to last category if there category is above the maximum category
  dframe[dframe[, "category"] > cat.max, "category"] <- cat.max
  #Get Number of CpGs by categories
  cat.sizes <- tabulate(dframe$`category` + 1L, cat.max + 1L)
  cat.breaks <- cumsum(rep(cat.step, cat.max))
  cat.first <- ifelse(is.integer(dframe[, 1]), "1", "0")
  #Create X Axis Labels
  cat.labels <- paste0("[", c(cat.first, cat.breaks[-length(cat.breaks)]),
                       ", ", cat.breaks, "]")
  cat.labels <- c(cat.labels, paste("over", tail(cat.breaks, 1)))
  #Convert categories as factors
  dframe$`category` <- factor(dframe$`category`,
                              levels = c(0:(length(cat.breaks))))
  #Add levels
  levels(dframe$`category`) <- cat.labels
  #Get boxplot statistics
  dfr <- tapply(dframe[, 2], dframe[, 3],
                function(x) { boxplot.stats(x)$stats })
  #Create new dataframe of boxplot stats with categories and amount of CpGs
  dfr <- data.frame(
    "x" = factor(rep(names(dfr), each = length(dfr[[1]])), levels = names(dfr)),
    "y" = unlist(dfr, use.names = FALSE),
    "amount" = factor(rep(cat.sizes, each = length(dfr[[1]])),
                      levels = cat.sizes))
  #Categories as Levels
  levels(dfr$x) <- paste0(levels(dfr$x), "\n",
                          as.integer(round(cat.sizes / 1000)))
  #Return dataframe 
  dfr
}


## G L O B A L S ###############################################################

## File storing the summary objects
FILE.SUMMARY <- paste0("data/","CpG_Table_DNMT1.RDS")

## Load the statistics table
tbl <- readRDS(FILE.SUMMARY)
tbl$pos.cov <- as.integer(tbl$pos.cov)
tbl$pos.cov.ten <- as.integer(tbl$pos.cov.ten)
# av.cov: average coverage
# pos.cov: number of samples with positive coverage
# av.pos.cov: average coverage amongth samples with pos. coverage
# min.cov: minimum coverage
# fith.perc.cov: fifth percentile coverage
# fiftyth.perc.cov: median coverage
# ninetyfith.perc.cov: 95th percentile coverage
# pos.cov.ten: number of samples with coverage >= 10
# av.meth: average methylation (among samples with positive coverage)
# sd.meth: SD of methylation (among samples with positive coverage)
rm(FILE.SUMMARY)

#Subset table
dframe <- tbl[, c("av.cov", "min.cov", "av.meth")]
#Remove CpGs for which there is no coverage and so no methylation across samples
dframe <- dframe[dframe[, 1] != 0L, ]
#Give colnames
colnames(dframe) <- c("[Average Coverage Ranges]\nAmount of CpGs in thousands",
                      "[Minimum Coverage Ranges]\nAmount of CpGs in thousands",
                      "Average Methylation")


## Plot Categorized Average Methylation against Average Coverage ###############
#Generate dataframe of Average Methylation and Average Coverage
dfr <- plot.box.categories(dframe[, c(1, 3)])

#Create Coverage Palette
cov_palette<-colorRampPalette(c("aliceblue","deepskyblue3"))

#Get Category Sizes 
cat.sizes<-as.numeric(as.character(unique(dfr$amount)))

#Plot Boxplots
ggplot(dfr, aes_string('x', 'y')) +
  geom_boxplot(outlier.size = NA, aes(fill = dfr$x)) +
  labs(x = colnames(dframe)[1], y = colnames(dframe)[3]) +
  scale_y_continuous(breaks = seq(0,1,0.2), labels = seq(0,1,0.2))+
  scale_fill_manual(values=cov_palette(max(cat.sizes))[cat.sizes]) + 
  theme_minimal() +
  theme(legend.position="none",
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13))
  

## Plot Categorized Average Methylation against Minimum Coverage ###############
#Generate dataframe of Average Methylation and Minimum Coverage
dfr <- plot.box.categories(dframe[, c(2, 3)],cat.step = 5L)

#Create Coverage Palette
cov_palette<-colorRampPalette(c("aliceblue","deepskyblue3"))

#Get Category Sizes 
cat.sizes<-as.numeric(as.character(unique(dfr$amount)))

#Plot Boxplots
ggplot(dfr, aes_string('x', 'y')) +
  geom_boxplot(outlier.size = NA, aes(fill = dfr$x)) +
  labs(x = colnames(dframe)[2], y = colnames(dframe)[3]) +
  scale_y_continuous(breaks = seq(0,1,0.2), labels = seq(0,1,0.2))+
  scale_fill_manual(values=cov_palette(max(cat.sizes))[cat.sizes]) + 
  theme_minimal() +
  theme(legend.position="none",
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13))
