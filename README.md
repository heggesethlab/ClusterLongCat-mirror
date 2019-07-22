# Clustering Longitudinal Categorical Data

This shiny app is designed as a tool for comparing methods for clustering longitudinal categorical data. It supports doing cluster analysis on such data using a variety of appropriate clustering methods. For each method, it produces visualizations and statistics to help interpret the clustering assignments. We also allow for comparing different clusterings and producing clustering comparison statistics such as the adjusted Rand index. We include two data sets, but users could upload and analyze datasets of their choosing.

## Getting Started

A pared down version of our project can be found on our [shinyapps.io page](https://heggesethlab.shinyapps.io/ShinyApp/). The full version should be run locally due to the computational complexity of performing clusterings.

### Prerequisites  

You'll need a copy of [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/) in order to run the app locally.  

Additionally, you'll need several R packages, which can be installed by inputting the following into the RStudio console:
```
install.packages(c("abind",
                   "cluster",
                   "data.table",
                   "dendextend",
                   "dplyr",
                   "ggplot2",
                   "htmlTable",
                   "igraph",
                   "markdown",
                   "partitionComparison",
                   "purrr",
                   "rmarkdown",
                   "scales",
                   "seqHMM",
                   "seriation",
                   "shiny",
                   "shinybusy",
                   "shinydashboard",
                   "shinyMatrix",
                   "stringr",
                   "tidyr",
                   "TraMineR",
                   "viridis"))

install.packages("devtools")
devtools::install_github("cran/bayesMCClust")
```

### Running the app

Open `ShinyApp/shiny_prototype` in RStudio, then press *Run App* in the upper left hand corner. A window with the app running in it should pop up.


## Authors
* Dr. Brianna Heggeseth
* Ellen Graham
* Zuofu Huang
* Kieu-Giang Nguyen
