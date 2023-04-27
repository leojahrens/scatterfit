# SCATTERFIT v1.6

Scatterfit includes two commands for Stata that produce a wide range of scatter plots with overlaid fit lines: scatterfit and slopefit. To install the package, execute the following command in Stata

```
net install scatterfit, from(https://raw.githubusercontent.com/leojahrens/scatterfit/master) replace
```

Once installed, please see `help scatterfit` and `help slopefit` for syntax and the whole range of options.
There are two tutorials below that showcase the possibilities of the two commands. 

# Creating plots with scatterfit
## Basics

Scatterfit uses `graph twoway` to plot (1) scatter points, (2) fit lines, and optionally (3) confidence intervals. The standard syntax is `scatterfit y x [,options]`. In its simplest form, it is essentially a `graph twoway (scatter y x) (lfit y x)` command, although much better looking by default. Let's see it in action, using data from the German European Social Survey sample from 2018-19.

```
scatterfit pinc_net age
```
<img src="./examples/gr1.png" height="300">

While this graph is easy to rebuild with `graph twoway`, scatterfit has a lot of added functionalities that not only make it better looking than standard Stata graphs but also much more flexible, with very simple syntax.

One of the big upsides of scatterfit is its ability to create binned scatter plots. The data are categorized into a number of bins along the x-dimension (in our example age), after which mean values of y and x are plotted as scatter points.

```
scatterfit pinc_net age, binned
```
<img src="./examples/gr2.png" height="300">

The advantages of binned scatter plots should become apparent. You can get a much better feel of how the variables are related, especially when the original data have a lot of datapoints. Another good use case for binned scatter plots are discrete y variables. By default, these are bundled on a limited number of y-values, such as a 5 point Likert scale. Scatter plots are usually useless with these variables, but binned scatter plots work greatly. Compare the two figures below. 

```
scatterfit redistr pinc_net
```
```
scatterfit redistr pinc_net, binned
```
<img src="./examples/gr3.png" height="200"><img src="./examples/gr4.png" height="200">

By default, the x variable is binned according to quantile cutoff-points so that the data are categorized into (nearly) 20 equally sized groups. Using the `nquantiles()` option changes the number of equally sized groups. 

```
scatterfit redistr pinc_net, binned nquantiles(50)
```
<img src="./examples/gr5.png" height="300">

It is also possible to let scatterfit categorize the x variable into equally spaced groups with the option `unibin()`, so that the distance in cutoff points remains constant across the x variable (for example categorizing a variable from 0 to 100 into two equally spaced groups would result in the two categories 0-50, >50-100). As a result, the number of observations within each bin may differ.
 
```
scatterfit pinc_net age, binned unibin(30)
```
<img src="./examples/gr6.png" height="300">

Scatterfit can also treat the x variable as discrete, using each distinct value of x as a bin.

```
scatterfit pinc_net age, binned discrete
```
<img src="./examples/gr7.png" height="300">

Use the `mweighted` option to visualize the differing amount of observations within the bins.

```
scatterfit pinc_net age, binned discrete
```
<img src="./examples/gr8.png" height="300">

The last option is to use `binvar()` to define a variable that already contains all the bins.

## Binary dependent variables

The command also works with binary dependent variables. It will automatically detect if the y variable only has two distinct values (after accounting for sample reductions via listwise deletion or if/in). No matter the underlying scale, scatterfit will transform the variable into a dummy where the original value with the higher scale point gets the "1" coding. Of course, scatter plots make little sense for binary dependent variables because all points are bundled on 0 and 1. But, you guessed it, binned scatter plots work just fine. The points and fit line show the proportion of respondents in the "1" rather than the "0" category.

```
scatterfit right pinc_net, binned
```
<img src="./examples/gr9.png" height="300">

The fit line in the binary DV case is derived from linear predictions from a logit model that is calculated internally.

## Fit types

The command supports several fit lines, manipulated via the `fit()` option. The standard setting is `fit(lfit)` for a linear fit with no confidence intervals. Other possible settings are `fit(lfitci)` to include confidence intervals, `fit(qfit)` or `fit(qfitci)` for quadratic fits, `fit(poly)` or `fit(polyci)` for local polynomial fits, and `fit(lowess)` for a lowess smoother. To control how fine-grained the smoothing is for local polynomials and lowess, use the `bwidth()` option. Here are some examples.

<img src="./examples/gr10.png" height="200"><img src="./examples/gr11.png" height="200">
<img src="./examples/gr12.png" height="200"><img src="./examples/gr13.png" height="200">

## Multiple plots along a by-dimension

The values of y can be plotted separately for values across a by-dimension. Here is the development of income by age separately for men and women.
```
scatterfit pinc_net age, binned fit(lfitci) by(gender)
```
<img src="./examples/gr14.png" height="300">

This is also a good opportunity to showcase the `mlabel()` option, which replaces the usual scatter markers with a string contained in a separate string variable or a numeric variable with labeled values.
```
scatterfit pinc_net age, binned by(gender) mlabel(genderlab)
```
<img src="./examples/gr15.png" height="300">

## Control variables

Work in progress...





