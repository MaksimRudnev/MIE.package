# MIE -- Measurement invariance explorer <a href="https://maksimrudnev.github.io/MIE.package/"><img src="man/figures/logo.png" alt="MIE website" align="right" height="120"/></a>

To install, type in the RStudio console: `devtools::install_github("maksimrudnev/MIE.package", dependencies = TRUE)`

Invariance explorer helps to find groups that might demonstrate measurement invariance of latent factors. Instead of repeatedly running the models to find a subset of groups that do show invariance, MIE makes it easy to plot the groups based on a set of invariance metrics. First, you choose a metric of invariance in `getMetrics()`, second, you plot them with `plotDistances()`. There are many options available, a set of convenience functions to, e.g. put together and extract results of Mplus alignment models, as well as a friendly Shiny app to do all the steps interactively.


### Covariance-based approach

Usually, assessment of measurement invariance is based on a measurement model (which is assumed to be a true model). However, one may employ completely exploratory approach, using covariances between observed variables to find groups with the most similar structures. This is the default (and the quickest) option in `MIE` package. Of course, at some point one has to specify a measurement model, but it may be done *after* finding a set of groups with similar covariance structures.

If, or when, a measurement model is given, one can use it in `MIE` to find sets of groups that are closer to each other using either model parameters or model fit indices.

### Parameter-based approach

You need to specify a measurement model using [lavaan syntax](http://lavaan.ugent.be/tutorial/syntax1.html). `MIE` computes configural model to compare factor loadings and metric invariance model (with fixed factor loadings) to compare intercepts. It is recommended first be sure that you found a set of groups that have similar loadings before comparing intercepts.


### Fit indices-based

When using fit indices as a measure of group proximity, `MIE` computes for each pair of groups either configural and metric models, or metric and scalar models, and plots the groups based on a difference between less and more restricted models, as it was suggested by [Chen, 2007](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.459.8501&rep=rep1&type=pdf). Chen also suggested a threshold of .01. For convenience, a circle .01x.01 is added to the plot. If all the groups are inside this circle, it is very likely that measurement invariance holds for these groups. When using this approach, please check the actual fit indices (e.g., CFI should be above .9), as the implied model is assumed to be true.

### Main purpose

This approach is useful when there are either clusters of groups or outlier groups (in terms of implied model/correlations). This approach is not very suitable for finding model misspecifications.


Please read the vignettes for more details: <https://maksimrudnev.github.io/MIE.package/articles/>


### Cautions

-   As follows from its name, `MIE` is exploratory. The plots are based basically on averaging a bunch of information, so it doesn't guarantee that closer groups have measurement invariance.
-   When using implied models, pay attention to the fit indices themselves, not only to their differences. Model in general should fit the data well enough.
-   The use of some fit indices in plots based on multidimensional scaling and in cluster analysis might violate assumptions, as they have not been proven to have a known distribution.
-   With big datasets and complex models calculation of metrics could take *a considerable amount of time*.



