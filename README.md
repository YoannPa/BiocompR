# BiomcompR - Advanced visualizations for data comparison
_**BiocompR** is an R package built upon ggplot2 to improve commonly used plots dedicated to data comparison and dataset exploration and ultimately provides users with versatile and customizable graphics._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**Contributors: SCHEFZIK R.<sup>2</sup>; HRUSKA D.<sup>1</sup>; KURILOV R.<sup>1</sup>; WURSTHORN A.<sup>4</sup>; MAYAKONDA A.<sup>3</sup>; FEUERBACH L.<sup>1</sup>; TOTH R.<sup>3</sup>.**  
**1-** [**DKFZ - Division of Applied Bioinformatics, Germany.**](https://www.dkfz.de/en/applied-bioinformatics/index.php)  
**2-** [**Klinik für Anästhesiologie und Operative Intensivmedizin, Medizinische Fakultät Mannheim, Universität Heidelberg, Germany.**](https://www.umm.de/klinik-fuer-anaesthesiologie-und-operative-intensivmedizin/)  
**3-** [**DKFZ - Computational Cancer Epigenomics, Germany.**](https://www.dkfz.de/en/CanEpi/CompEpigen/index.html)  
**4-** [**DKFZ - Clinical Cooperation Unit Translational Radiation Oncology, Germany.**](https://www.dkfz.de/en/molekulare-radioonkologie/index.php)  

**Version: 0.0.76 (Beta)**  
**R Compatibility: Version 3.6.3**  
**Last Update: 18/08/2020**  
**How to cite:** _Pageaud Y. et al., BiomcompR - Advanced visualizations for data comparison._  

## Content
Currently the package BiocompR contains **18 functions**:

* `basic.sidebar()` - Draws a ggplot2 of a basic sidebar.  
* `bivar.plot()` - Computes boxplots or violins from 1 variable values against ranges of a 2nd one.  
* `EVA()` - Computes eigenvectors, principal component scores and correlations from a correlation test.  
* `fancy.hist()` - Computes in parallel and plot an histogram using ggplot2 from a given vector of values.  
* `fused.plot()` - Creates a plot summarizing results from 2 different pairwise comparisons.  
* `fused.view()` - Displays 2 matrices of results as a fused plot.  
* `gg2heatmap()` - Creates a custom heatmap with dendrograms and annotations.
* `ggcoverage()` - Plots an annotated barplot.  
* `ggcraviola()` - Draws a craviola plot (half-splitted and percentile-binned violin plot).  
* `ggdend()` - Creates a dendogram in ggplot2.  
* `ggeigenvector()` - Creates an eigenvector plot using ggplot2.  
* `ggpanel.corr()` - Plots results of correlation test between a single variable and multiple others as jittered scatter plot divided into 4 different panels.  
* `ggvolcano.corr()` - Plots results of correlation test between a single variable and multiple others as volcano plot.  
* `ks.plot()` - Computes pairwise Kolmogorov-Smirnov tests on a matrix and display results in a fused plot.  
* `plot.col.sidebar()` - Creates a colored side annotation bars in ggplot2.  
* `resize.grobs()` - Resizes heights or widths of multiple grobs based on a given grob dimensions.  
* `sunset()` - Draws a sunset plot showing the completeness of a dataset.  
* `warn.handle()` - Filters unrelevant warnings matching a regular expression.  

## Prerequisites
### Install Linux dependencies
In the terminal do:  
```bash
sudo apt-get install libcurl4-openssl-dev libxml2-dev
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
	      'parallelDist', 'quantmod')
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

## Technical questions / Development / Feature request
If you encounters issues or if a feature you would expect is not available in of BiocompR functions, please check if an existing issue adresses your point [here](https://github.com/YoannPa/BiocompR/issues/). If not, create a [new issue here](https://github.com/YoannPa/BiocompR/issues/new).  

## References
⚠️ **Work in progress !**  

## Licence
The repository BiocompR is currently under the GPL-3.0 licence.  

