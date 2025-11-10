# Tracy–Widom test

This function is adapted from the `tw` function in the archived CRAN
package *AssocTests*. The implementation and original documentation were
authored by Lin Wang, Wei Zhang, and Qizhai Li. The AssocTests package
is archived on CRAN; see the CRAN archive for details.

## Usage

``` r
tw(eigenvalues, eigenL, criticalpoint = 2.0234)
```

## Arguments

- eigenvalues:

  A numeric vector of eigenvalues (sorted decreasing).

- eigenL:

  The number of eigenvalues.

- criticalpoint:

  A numeric value corresponding to the significance level. If the
  significance level is 0.05, 0.01, 0.005, or 0.001, set `criticalpoint`
  to `0.9793`, `2.0234`, `2.4224`, or `3.2724`, respectively. Default:
  `2.0234`.

## Value

A list with class `"htest"`:

|               |     |     |                                            |
|---------------|-----|-----|--------------------------------------------|
| `statistic`   |     |     | a vector of Tracy–Widom statistics.        |
| `alternative` |     |     | description of the alternative hypothesis. |
| `method`      |     |     | test name.                                 |
| `data.name`   |     |     | name of the input data.                    |
| `SigntEigenL` |     |     | number of significant eigenvalues.         |

## Details

Find the significant eigenvalues of a matrix (from AssocTests).

This function implements the Tracy–Widom test to determine how many
leading eigenvalues of a matrix are statistically significant. This
function and its documentation are from the `tw` function in the
archived CRAN package *AssocTests* (version 1.0-1).

The input `eigenvalues` should be sorted in descending order. The test
statistic is computed for each leading eigenvalue and compared to a
supplied Tracy–Widom critical value.

## References

Lin Wang, Wei Zhang, and Qizhai Li. AssocTests: An R Package for Genetic
Association Studies. *Journal of Statistical Software*. 2020; 94(5):
1–26.

N. Patterson, A. L. Price, and D. Reich. Population Structure and
Eigenanalysis. *PLoS Genetics*. 2006; 2(12): 2074–2093.

C. A. Tracy and H. Widom. Level-Spacing Distributions and the Airy
Kernel. *Communications in Mathematical Physics*. 1994; 159(1): 151–174.

A. Bejan. Tracy–Widom and Painlevé II: Computational Aspects and
Realisation in S-Plus. In *First Workshop of the ERCIM Working Group on
Computing and Statistics*. 2008, Neuchâtel, Switzerland.

A. Bejan. Largest eigenvalues and sample covariance matrices. *MSc
Dissertation, University of Warwick*. 2005. (This function was written
by A. Bejan and made publicly available.)

## Author

Original authors: Lin Wang, Wei Zhang, and Qizhai Li (AssocTests).
Documentation adaptation: package maintainers of this project.

## Examples

``` r
tw(eigenvalues = c(5, 3, 1, 0), eigenL = 4, criticalpoint = 2.0234)
#> 
#>  Tracy-Widom test
#> 
#> data:  c(5, 3, 1, 0)
#> TW1 = -0.82427, TW2 = -0.60186, TW3 = -0.55525, TW4 = NaN
#> alternative hypothesis: the eigenvalue is significant
#> 
```
