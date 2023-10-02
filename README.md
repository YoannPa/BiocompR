# BiocompR - Advanced visualizations for data comparison

![GitHub repo size](https://img.shields.io/github/repo-size/YoannPa/BiocompR)
![GitHub issues](https://img.shields.io/github/issues-raw/YoannPa/BiocompR)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/YoannPa/BiocompR)  

_**BiocompR** is an R package built upon ggplot2, and using data.table. It improves some visualisations commonly used in biology and genomics for data comparison and dataset exploration, introduces new kind of plots, provides a toolbox of functions to work with ggplot2 and grid objects, and ultimately, allows users to customize plots produced into publication ready figures._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**How to cite:** _Pageaud Y. et al., BiocompR - Advanced visualizations for data comparison._  

![GitHub R package version](https://img.shields.io/github/r-package/v/YoannPa/BiocompR?label=Package%20version&logo=RStudio&logoColor=white&style=for-the-badge)  
<img src="https://img.shields.io/static/v1?label=compatibility&message=4.3.1&color=blue&logo=R&logoColor=white&style=for-the-badge" />  
![GitHub last commit](https://img.shields.io/github/last-commit/YoannPa/BiocompR?logo=git&style=for-the-badge)  
![GitHub](https://img.shields.io/github/license/YoannPa/BiocompR?color=brightgreen&style=for-the-badge)  

### Ackowledgment
_I would like to thank every people who contributed to the development of this package, with their code, their test datasets, their advices and feedbacks. - Yoann_.  
**Contributors:** Dr. Schefzik R.<sup>2</sup>; Mr. Hruska D.<sup>1</sup>; Mrs. Bitto V.<sup>1</sup>; Dr. Kurilov R.<sup>1</sup>; Mr. Beumer N.<sup>1</sup>; Mrs. Wursthorn A.<sup>3</sup>; Dr. Feuerbach L.<sup>1</sup>.  
**1.** [**DKFZ - Division of Applied Bioinformatics, Germany.**](https://www.dkfz.de/en/applied-bioinformatics/index.php)  
**2.** [**Klinik für Anästhesiologie und Operative Intensivmedizin, Medizinische Fakultät Mannheim, Universität Heidelberg, Germany.**](https://www.umm.de/klinik-fuer-anaesthesiologie-und-operative-intensivmedizin/)  
**3.** [**DKFZ - Clinical Cooperation Unit Translational Radiation Oncology, Germany.**](https://www.dkfz.de/en/molekulare-radioonkologie/index.php)  

## Installing BiocompR
In R execute the following command:
```R
devtools::install_github("YoannPa/BiocompR")
```

## Content
Currently the package BiocompR contains **39 exported functions**:

* `basic.sidebar()` - Draws a ggplot2 of a basic sidebar.  
* `biopalette()` - A color palette advisor for biology plots.  
* `bivar.plot()` - Draws boxplots or violins from a variable values against ranges of a 2nd one.  
* `build_legends_layout()` - Builds legends layout.  
* `check_fun()` - Checks if a function exists and package of origin.  
* `cross.biplot()` - Computes and draws biplots for multiple principal components at once.  
* `EVA()` - Computes eigenvectors, principal component scores and correlations from a correlation test.  
* `fused.plot()` - Creates a plot summarizing results from 2 different pairwise comparisons.  
* `fused.view()` - Displays 2 matrices of results as a fused plot.  
* `gg2heatmap()` - Creates a custom heatmap with dendrograms and annotations.  
* `ggbipca()` - Computes and draws a custom PCA biplot.  
* `ggcirclart()` - Circlizes ggplot2 objects.  
* `ggcoverage()` - Plots an annotated barplot.  
* `ggcraviola()` - Draws a craviola plot (half-splitted and percentile-binned violin plot).  
* `ggdend()` - Creates a dendogram in ggplot2.  
* `ggdensity_map()` - Plots a density color map from a matrix or a molten data.frame.  
* `ggeigenvector()` - Creates an eigenvector plot using ggplot2.  
* `ggfusion.free()` - Displays 2 triangle matrices fused together in a single plot.  
* `gghist()` - Plots an histogram using ggplot2 from a numeric or character vector.  
* `ggpanel.corr()` - Plots results of correlation test between a single variable and multiple others as jittered scatter plot divided into 4 different panels.  
* `ggstackbar()` - Draws stacked barplots from an annotation table.  
* `ggsunset()` - Draws a sunset plot showing the completeness of a dataset.  
* `ggtriangle()` - Draws a triangle ggplot for a basic molten triangle matrix.  
* `ggvolcano.corr()` - Plots results of correlation test between a single variable and multiple others as volcano plot.  
* `ggvolcano.free()` - Plots any kind of results with P-values that can be displayed as a volcano plot.  
* `ggvolcano.test()` - Plots results of a Plots results of statistical tests as volcano plot.  
* `ks.plot()` - Computes pairwise Kolmogorov-Smirnov tests on a matrix and display results in a fused plot.  
* `manage.na()` - Keeps, removes or imputes missing values in a matrix or a data.frame based on sample groups.  
* `plot.col.sidebar()` - Creates a colored side annotation bars in ggplot2.  
* `plot_asso_all_annot()`	- Draws association test results between all columns from a data.frame.  
* `plot_asso_annot_PC()` - Plots association tests' results between some annotations and some PCs.  
* `raster.gg2grob()` - Rasterize a gg plot into a raster grob.  
* `raster.ggplot.to.grob()` - Rasterize a gg plot into a raster grob.  
* `resize.grob.oneway()` - Resizes heights or widths of a grob based on the dimensions of another grob.  
* `resize.grobs()` - Resizes heights or widths of multiple grobs based on a given grob dimensions.  
* `test.annots()` - Tests association of an annotation with another one or with a PC.  
* `test_asso_all_annot()` - Write function description here.  
* `test_asso_annot_pc()` - Tests associations between a set of annotations and PCs from a prcomp object.  
* `warn.handle()` - Filters irrelevant warnings matching a regular expression.  

## Problems ? / I need help !
For any questions **Not related to bugs or development** please check the section "**Known Issues**" available below. If the issue you experience is not adressed in the known issues you can write me at [y.pageaud@dkfz.de](y.pageaud@dkfz.de).

### Known Issues 
**❎ Error in UseMethod("depth")**  
```R
Error in UseMethod("depth") : 
  no applicable method for 'depth' applied to an object of class "NULL"
```
This error seems to happen randomly when executing code using the ggplot2 and/or grid packages. Usually executing one more time the chunck of code solve the error. The current statues of this issue can be tracked [**here**](https://github.com/tidyverse/ggplot2/issues/2514).


**❎ Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto), : Viewport has zero dimension(s)**  
```R
Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto), :
  Viewport has zero dimension(s)
```
This error can arise when using the `ggbipca()` function: if you define a legend with too many values, the plotting area becomes too small to print the plot in plotting panel of RStudio.  
When it happens, you can try to manually increase the size of the plotting panel in your RStudio interface. If doing this doesn't solve the error, then it is advised to define a legend with fewer values for colors and/or shapes.  

**⚠️ Reached elapsed time limit.**
```
Warning message:
In grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  :
  reached elapsed time limit.
```
This warning seems to happen randomly when executing code using the ggplot2 and/or grid packages. Usually after executing one more time the chunck of code the warning does not display anymore. The current statues of this issue can be tracked [**here**](https://stackoverflow.com/questions/51247102/reached-elapsed-time-limit-errors-in-r).

**⚠️ Using alpha for a discrete variable is not advised.**
```R
Warning message:
Using alpha for a discrete variable is not advised. 
```
This warning can arise when using the function `ggvolcano.corr()` with additionnal ggplot2 components. It doesn't compromise the printing of the plot, however you might feel annoyed by it.  
A quick fix to suppress specifically this warning is to use the function `warn.handle()`, which filters out annoying warnings using pattern matching, as following:  
```R
#Create your correlation volcano plot (this will also print the 'default' volcano plot)
my_volcano <- ggvolcano.corr(
  data = dfrm_my_correlation_res, p.cutoff = 0.01, corr.cutoff = 0.1,
  title.corr.cutoff = "Samples default correlation",
  corr.label.cutoff = c(-0.35,0.40)) +
  scale_color_manual(values = ggsci::pal_npg("nrc", alpha = 1)(10)) +
  xlab("Spearman correlation") + ylab("Spearman P-value") +
  ggtitle("Spearman correlation between multiple variables and my variable of interest")

#Print your volcano plot without displaying the annoying warning
warn.handle(
  pattern = "Using alpha for a discrete variable is not advised.",
  print(my_volcano)) 
```
Nevertheless, using `ggvolcano.corr()` without additionnal ggplot2 components
should not raise this warning.  

**⚠️ ggrepel: ## unlabeled data points (too many overlaps). Consider increasing max.overlaps.**
```R
Warning message:
ggrepel: ## unlabeled data points (too many overlaps). Consider increasing max.overlaps
```
This warning can arise when using the function `ggbipca()`. If the scale is too small, and you want to display too many loadings labels, then those overlapping will not be displayed, and this warning will be printed. Using this function has shown that this specific warning can persist, and be printed randomly afterward when running other commands. It is unclear why this is happening. But it can be fixed by executing the following command, once you ran ggbipca():  
```R
assign("last.warning", NULL, envir = baseenv())
```
The current statues of this issue can be tracked [**here**](https://github.com/slowkow/ggrepel/issues/187).  

**⚠️ In min(x) : no non-missing arguments to min; returning Inf / In max(x) : no non-missing arguments to max; returning -Inf**

```R
Warning messages:
1: In min(x) : no non-missing arguments to min; returning Inf
2: In max(x) : no non-missing arguments to max; returning -Inf
```
This warning can arise when using the function `sunset()`when there is only 1 label displayed on the right Y axis. This warning does not compromise the result and should be ignored.  
The current statues of this issue can be tracked [**here**](https://github.com/tidyverse/ggplot2/issues/4368).  



## Technical questions / Development / Feature request
If you encounters issues or if a feature you would expect is not available in a BiocompR function, please check if an existing issue adresses your point [here](https://github.com/YoannPa/BiocompR/issues/). If not, create a [new issue here](https://github.com/YoannPa/BiocompR/issues/new).  

## References
⚠️ **Work in progress !**  

1. [_Share a legend between two ggplot2 graphs - Mara Averick_](https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs)  
2. [_Align two plots on a page - Mara Averick_](https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page)  
3. [_ggfortify: Data Visualization Tools for Statistical Analysis Results_](https://cran.r-project.org/web/packages/ggfortify/index.html)  
4. [_ggfortify: Plotting PCA (Principal Component Analysis_](https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html)  
5. [_Loadings vs eigenvectors in PCA: when to use one or another?_](https://stats.stackexchange.com/questions/143905/loadings-vs-eigenvectors-in-pca-when-to-use-one-or-another)  
6. [_What is the proper association measure of a variable with a PCA component?_](https://stats.stackexchange.com/questions/119746/what-is-the-proper-association-measure-of-a-variable-with-a-pca-component-on-a/)  

## Licence
BiocompR is currently under the GPL-3.0 licence.  

