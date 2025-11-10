# MMRR

## Multiple matrix regression with randomization (MMRR)

``` r

library(algatr)
```

``` r

# Install required packages
mmrr_packages()
```

``` r

library(here)
library(raster)
```

**If using MMRR, please cite the following: Wang I.J. (2013) Examining
the full effects of landscape heterogeneity on spatial genetic
variation: a multiple matrix regression approach for quantifying
geographic and ecological isolation. Evolution 67(1):3403-3411.
<DOI:10.1111/evo.12134>**

Multiple matrix regression with randomization ([Wang
2013](https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12134))
examines how heterogeneity in the landscape (both geographic and
environmental) affects spatial genetic variation, allowing a user to
determine the relative contributions of isolation by distance (i.e., an
association between genetic and geographic distances) and isolation by
environment (i.e., an association between genetic and environmental
distances). MMRR can provide information concerning how dependent
variables (in our case, genetic data) change with respect to multiple
independent variables (environmental and geographic distances).

The randomization aspect of MMRR is that significance testing is
performed through random permutations of the rows and columns of the
dependent matrix (the genetic distances), which is necessary because of
the non-independence of values in pairwise distance matrices. Once
significance testing is completed, MMRR provides individual regression
coefficients and p-values for each dependent variable and a “coefficient
ratio,” which is the ratio between regression coefficients, which thus
provides an idea of the relative contributions of IBD and IBE in
explaining variation in the genetic distances in your data.

There are a few assumptions built into this function that users should
be aware of: (1) the coordinates and genetic distance files MUST have
the same ordering of individuals; there isn’t a check for this, and (2)
this function assumes that each individual has its own sampling
coordinates (even if population-based sampling was performed).

#### Read in and process genetic data

Let’s load the example dataset and extract environmental variables based
on our sampling coordinates. To perform MMRR, we need a genetic distance
matrix; please refer to the genetic distances vignette for information
on how to generate this in algatr using a vcf file.

``` r

load_algatr_example()
#> 
#> ---------------- example dataset ----------------
#>  
#> Objects loaded: 
#> *liz_vcf* vcfR object (1000 loci x 53 samples) 
#> *liz_gendist* genetic distance matrix (Plink Distance) 
#> *liz_coords* dataframe with x and y coordinates 
#> *CA_env* RasterStack with example environmental layers 
#> 
#> -------------------------------------------------
#> 
#> 
Y <- as.matrix(liz_gendist)
```

#### Process environmental data

Let’s calculate environmental distances with our extracted environmental
values using the
[`env_dist()`](https://thewanglab.github.io/algatr/reference/env_dist.md)
function, and add on geographic distances to this matrix using the
[`geo_dist()`](https://thewanglab.github.io/algatr/reference/geo_dist.md)
function. For more information on how these two functions work, please
refer to algatr’s environmental data processing vignette.

``` r

# Extract values from our environmental raster (CA_env)
env <- raster::extract(CA_env, liz_coords)
# Calculate environmental distances
X <- env_dist(env)
# Add geographic distance to X
X[["geodist"]] <- geo_dist(liz_coords)
```

### Run MMRR

------------------------------------------------------------------------

We can run MMRR using the
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md)
function. We can either run the full model (with no variable selection)
by setting `model = "full"`, or by running a variable selection
procedure using backwards elimination by setting `model = "best"` within
this function. The arguments to
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md)
are as follows:

- `Y`: matrix of genetic distances

- `X`: matrix of environmental distances (including geographic
  distances) at each sampling coordinate

- `nperm`: number of permutations to perform

- `stdz`: whether matrices are standardized

#### MMRR with all variables

First, let’s run the full MMRR model with all variables included (i.e.,
no variable selection).

``` r

set.seed(10)
results_full <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "full")
```

The results from running the “full” MMRR model contains four elements:

- `coeff_df`: a dataframe with statistics relating to each variable’s
  distance related to genetic distance, including coefficient values for
  each environmental variable and geographic distance

- `mod`: a dataframe containing statistics for the results of the model,
  including an R^2 value, and F statistics

- `X` and `Y`: the input data

#### MMRR with model selection

Now, let’s run MMRR with variable selection.

``` r

# Run MMRR with all variables
set.seed(01)
results_best <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "best")
```

Once again, the results from running
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md)
with model selection contains the same elements as “full”, except that
the “best” model result now also has an `X_best` element, which contains
values for significant variables; in this case, the only significant
variable found was geographic distance. You can also see that the
`coeff_df` element is reduced from all four variables to only those
variables that were significantly associated with genetic distance (in
the test example case, only geographic distance).

### Visualizing MMRR results

------------------------------------------------------------------------

#### Plotting MMRR results with `mmrr_plot()`

Let’s take a look at the results from our MMRR analyses above. We can
produce several plots to visualize results with the
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md)
function, including plotting single variable relationships
(`plot_type = "vars"`), plotting the fitted relationship that compares
predicted against the observed genetic distances
(`plot_type = "fitted"`), or plotting covariances between the predictor
variables (`plot_type = "cov"`). Lastly, you can set `plot_type = "all"`
to produce all three of these plots for your MMRR run.

``` r

# Single variable plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "vars", stdz = TRUE)
```

![](MMRR_vignette_files/figure-html/mmrr%20plots-1.png)

``` r

# Fitted variable plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "fitted", stdz = TRUE)
```

![](MMRR_vignette_files/figure-html/fitted%20plot-1.png)

``` r

# Covariance plot
mmrr_plot(Y, X, mod = results_full$mod, plot_type = "cov", stdz = TRUE)
```

![](MMRR_vignette_files/figure-html/cov%20plot-1.png)

How do the above results compare to those from the “best” model?

``` r

mmrr_plot(Y, X, mod = results_best$mod, plot_type = "all", stdz = TRUE)
```

![](MMRR_vignette_files/figure-html/mmrr%20best%20plots-1.png)![](MMRR_vignette_files/figure-html/mmrr%20best%20plots-2.png)![](MMRR_vignette_files/figure-html/mmrr%20best%20plots-3.png)

#### Look at MMRR output statistics with `mmrr_table()`

We can also take a closer look at the output statistics from our MMRR
runs using the
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md)
function. This function will also provide summary statistics if we set
`summary_stats` to TRUE (the default). The estimate column contains the
coefficients for each of the variables. *N.B.: significance values may
change slightly due to the permutation procedure.*

``` r

mmrr_table(results_full, digits = 2, summary_stats = TRUE)
```

[TABLE]

### Running MMRR with `mmrr_do_everything()`

------------------------------------------------------------------------

The algatr package also has an option to run all of the above
functionality in a single function,
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md).
This function runs MMRR while specifying whether a user would like to
run variable selection or not (using the `model` argument), producing a
table and all plots automatically. It will automatically extract
environmental variables based on coordinates and calculate environmental
distances from these data. Users can also choose to include geographic
distances within the independent variable matrix by setting `geo = TRUE`
(the default setting). **Please be aware that the `do_everything()`
functions are meant to be exploratory. We do not recommend their use for
final analyses unless certain they are properly parameterized.**

#### MMRR with all variables

``` r

set.seed(01)
mmrr_full_everything <- mmrr_do_everything(liz_gendist, liz_coords, env = CA_env, geo = TRUE, model = "full")
```

![](MMRR_vignette_files/figure-html/all%20vars-1.png)![](MMRR_vignette_files/figure-html/all%20vars-2.png)![](MMRR_vignette_files/figure-html/all%20vars-3.png)

[TABLE]

#### Extending MMRR’s functionality with resistance distances

The relationships between geographic distance and/or environmental
variables on spatial genetic variation result in patterns of isolation
by distance (IBD) or isolation by environment (IBE), respectively.
However, one other key phenomenon that can occur - and one in which
landscape genomicists are often interested in examining - is one in
which landscape features act to restrict gene flow, resulting in a
pattern known as isolation by resistance (IBR). IBR is typically
investigated by first generating a resistance surface which provides
information about how the landscape “resists” gene flow given the study
organism. Then, given the resistance surface, resistance distances can
be calculated using (once again) a variety of methods, but most
typically use [circuit
theory](https://conbio.onlinelibrary.wiley.com/doi/abs/10.1111/cobi.13230)
to model connectivity across the landscape.

The resistance surface can be generated using species distribution
models, habitat suitability models, or other types of approaches that
simulate forward evolution given sampling localities and genotypes
(e.g., resistanceGA; [Peterman
2018](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12984)).
Building a reliable resistance surface can be quite complex, requiring
many decisions on the part of the researcher and, in some cases,
extensive parameterization and model testing to generate this type of
surface. These sets of decisions are precisely why algatr does not
explicitly implement any IBR methods, but below we describe one way in
which the package can be extended to do so.

As was briefly mentioned in the environmental data vignette, the
[`geo_dist()`](https://thewanglab.github.io/algatr/reference/geo_dist.md)
function is able to calculate resistance-based distances, largely using
the
[gdistance](https://cran.r-project.org/web/packages/gdistance/index.html)
package ([van Etten
2017](https://www.jstatsoft.org/article/view/v076i13)). Importantly, a
user must provide a resistance surface using whichever their preferred
method is by specifying the `lyr` argument. The function then creates a
transition surface and calculates circuit distances (random-walk commute
distances) using the
[`gdistance::commuteDistance()`](https://AgrDataSci.github.io/gdistance/reference/commuteDistance-methods.html)
function.

Although algatr does not have the explicit functionality to quantify and
visualize IBR, we can investigate IBR using MMRR given that this method
conducts a multiple Mantel test given some distance metric. Thus, it is
possible to run
[`geo_dist()`](https://thewanglab.github.io/algatr/reference/geo_dist.md)
specifying `lyr` with your resistance raster and use the resulting
distances as input into MMRR. The following is a *toy* example of how
this could be done. For this example we’ll pretend our resistance
surface is one of our environmental layers (where greater values
indicate more resistance) and use it to calculate our resistance
distances.

``` r

# here we aggregate the layer for computational speed
lyr <- aggregate(CA_env$CA_rPCA1, 50)
plot(lyr)
points(liz_coords)
```

![](MMRR_vignette_files/figure-html/toy%20resist%20surface-1.png)

``` r


# Recreate MMRR input with resistance distances
# Calculate environmental distances
X <- env_dist(env)
# Add geographic distance to X
X[["resistdist"]] <- geo_dist(liz_coords, type = "resistance", lyr = lyr)
#> Calculating resistance distances... This can be time consuming with many points and large rasters.
```

Now, we can run MMRR as we did above but including resistance-based
distances in our input data.

``` r

# Run MMRR with resistance distances
results_resist <- mmrr_run(Y, X, nperm = 99, stdz = TRUE, model = "full")
mmrr_plot(Y, X, mod = results_resist$mod, plot_type = "all", stdz = TRUE)
```

![](MMRR_vignette_files/figure-html/mmrr%20resist-1.png)![](MMRR_vignette_files/figure-html/mmrr%20resist-2.png)![](MMRR_vignette_files/figure-html/mmrr%20resist-3.png)

``` r

mmrr_table(results_resist)
```

[TABLE]

### Additional documentation and citations

------------------------------------------------------------------------

|  | Citation/URL | Details |
|----|----|----|
| Manuscript with method | [Wang 2013](https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12134) | Paper describing MMRR method |
| Associated code | [MMRR tutorial](https://nature.berkeley.edu/wanglab/wp-content/uploads/2019/09/MMRRtutorial.zip) | Walkthrough of MMRR |
