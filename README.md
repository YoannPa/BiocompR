# BiocompR - Advanced visualizations for data comparison
_**BiocompR** is an R package built upon ggplot2 to improve commonly used plots dedicated to data comparison and dataset exploration and ultimately provides users with versatile and customizable graphics._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**Contributors: SCHEFZIK R.<sup>2</sup>; HRUSKA D.<sup>1</sup>; BITTO V.<sup>1</sup>; KURILOV R.<sup>1</sup>; BEUMER N.<sup>1</sup>; WURSTHORN A.<sup>4</sup>; MAYAKONDA A.<sup>3</sup>; FEUERBACH L.<sup>1</sup>; TOTH R.<sup>3</sup>.**  
**1-** [**DKFZ - Division of Applied Bioinformatics, Germany.**](https://www.dkfz.de/en/applied-bioinformatics/index.php)  
**2-** [**Klinik für Anästhesiologie und Operative Intensivmedizin, Medizinische Fakultät Mannheim, Universität Heidelberg, Germany.**](https://www.umm.de/klinik-fuer-anaesthesiologie-und-operative-intensivmedizin/)  
**3-** [**DKFZ - Computational Cancer Epigenomics, Germany.**](https://www.dkfz.de/en/CanEpi/CompEpigen/index.html)  
**4-** [**DKFZ - Clinical Cooperation Unit Translational Radiation Oncology, Germany.**](https://www.dkfz.de/en/molekulare-radioonkologie/index.php)  

**Version: 0.0.152 (Beta)**  
**R Compatibility: Version 4.0.5**  
**Last Update: 14/10/2021**  
**How to cite:** _Pageaud Y. et al., BiocompR - Advanced visualizations for data comparison._  

## Content
Currently the package BiocompR contains **25 functions**:

* `basic.sidebar()` - Draws a ggplot2 of a basic sidebar.  
* `bivar.plot()` - Computes boxplots or violins from 1 variable values against ranges of a 2nd one.  
* `build.layout()` - Builds legends layout.  
* `cross.biplot()` - Computes and draws biplots for multiple principal components at once.  
* `EVA()` - Computes eigenvectors, principal component scores and correlations from a correlation test.  
* `fancy.hist()` - Computes in parallel and plot an histogram using ggplot2 from a given vector of values.  
* `fused.plot()` - Creates a plot summarizing results from 2 different pairwise comparisons.  
* `fused.view()` - Displays 2 matrices of results as a fused plot.  
* `gg2heatmap()` - Creates a custom heatmap with dendrograms and annotations.  
* `ggbipca()` - Computes and draws a custom PCA biplot.  
* `ggcoverage()` - Plots an annotated barplot.  
* `ggcraviola()` - Draws a craviola plot (half-splitted and percentile-binned violin plot).  
* `ggdend()` - Creates a dendogram in ggplot2.  
* `ggeigenvector()` - Creates an eigenvector plot using ggplot2.  
* `ggpanel.corr()` - Plots results of correlation test between a single variable and multiple others as jittered scatter plot divided into 4 different panels.  
* `ggvolcano.corr()` - Plots results of correlation test between a single variable and multiple others as volcano plot.  
* `ggvolcano.test()` - Plots results of a Plots results of statistical tests as volcano plot.  
* `ks.plot()` - Computes pairwise Kolmogorov-Smirnov tests on a matrix and display results in a fused plot.  
* `load.palettes()` - Loads pre-defined palettes.  
* `plot.col.sidebar()` - Creates a colored side annotation bars in ggplot2.  
* `raster.ggplot.to.grob()` - Rasterize a gg plot into a raster grob.  
* `resize.grob.oneway()` - Resizes heights or widths of a grob based on the dimensions of another grob.  
* `resize.grobs()` - Resizes heights or widths of multiple grobs based on a given grob dimensions.  
* `sunset()` - Draws a sunset plot showing the completeness of a dataset.  
* `warn.handle()` - Filters unrelevant warnings matching a regular expression.  

## Prerequisites
### Install Linux dependencies
In the terminal do:  
```bash
sudo apt-get install libcurl4-openssl-dev libxml2-dev libmagick++-dev
```
### Install Bioconductor dependencies
In R do:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IRanges")
```
### Install CRAN dependencies
```R
inst.pkgs = c('corrplot', 'data.table', 'devtools', 'fastcluster', 'ggdendro',
	      'ggplot2', 'ggrepel', 'grid', 'gridExtra', 'psych', 'parallel',
	      'parallelDist', 'quantmod', 'magick')
install.packages(inst.pkgs)
```

## Installing
1. In the Git repository click on "Clone or Download".
2. Copy the HTTPS link.
3. Open a terminal and type:
```bash
git clone https://github.com/YoannPa/BiocompR.git
```
4. Open the folder BiocompR and open the "BiocompR.Rproj" file in RStudio.
5. In the RStudio console, type:
```R
devtools::install()
```

## Problems ? / I need help !
For any questions **Not related to bugs or development** please check the section "**Known Issues**" available below. If the issue you experience is not adressed in the known issues you can write me at [y.pageaud@dkfz.de](y.pageaud@dkfz.de).

### Known Issues 
**❎ Error in UseMethod("depth")**  
```R
Error in UseMethod("depth") : 
  no applicable method for 'depth' applied to an object of class "NULL"
```
This error seems to happen randomly when executing code using the ggplot2 and/or grid packages. Usually executing one more time the chuck of code solve the error. The current statues of this issue can be tracked [**here**](https://github.com/tidyverse/ggplot2/issues/2514).


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
This warning seems to happen randomly when executing code using the ggplot2 and/or grid packages. Usually after executing one more time the chuck of code the warning does not display anymore. The current statues of this issue can be tracked [**here**](https://stackoverflow.com/questions/51247102/reached-elapsed-time-limit-errors-in-r).

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

## Licence
The repository BiocompR is currently under the GPL-3.0 licence.  

