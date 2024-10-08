# mediation.test

**More Powerful  Tests of the  Composite Null Hypothesis Arising  in Mediation
Analysis**

The `mediation.test` `R` package implements more powerful tests of the composite null hypothesis arising in mediation analysis. It is developed by A.  Chambaz (MAP5, UMR CNRS 8145, Université Paris Cité) \& Caleb Miles (Department of Biostatistics, Columbia University Mailman School of Public Health).

## Introduction

The `mediation_test()` function of the `mediation.test` package tests the composite null hypothesis <em>"&delta;<sub>x</sub> * &delta;<sub>y</sub> = 0"</em> against its alternative <em>"&delta;<sub>x</sub> * &delta;<sub>y</sub> &ne; 0"</em>.

```r
> library("mediation.test")
> example(mediation_test)
```

## Citation

To cite the package, see 

```r
> citation("mediation.test")
> toBibtex(citation("mediation.test"))
```

## Installation 

Both a stable version and a development version are available via [GitHub](https://github.com/achambaz/mediation.test) and can be installed in R as:

```r 
devtools::install_github("achambaz/mediation.test", ref = "main")
```

or 

```r 
devtools::install_github("achambaz/mediation.test", ref = "develop")
```
