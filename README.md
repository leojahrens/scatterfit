# SCATTERFIT v1.6

Scatterfit includes two commands for Stata that produce a wide range of scatter plots with overlaid fit lines. `scatterfit` visualizes the relationship between two variables. `slopefit` is for interaction models and visualizes the relationship between x and y conditional on another variable z.

To install the package, execute the following command in Stata

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

Scatterfit makes it possible to plot scatterplots with fit after accounting for control variables. If at all, researchers mostly use fitted scatterplots as a step preliminary to the "real" data analysis. The results are treated as descriptive, mostly because scatterplots appear insufficient for multivariate analyses. However, this is not the case. The data can be pre-treated before plotting them so that they only show the residual covariation between x and y after accounting for control variables. This is achieved by first regressing the x and y variables on the control variables and then using the residuals for the plots. If binned scatterplots are involved, this becomes slightly more complicated; see "Cattaneo et al. 2022 - On Binscatter". Scatterfit uses the simple residualization method for all plots without bins and the Cattaneo et al. method for all binned scatterplots.

The `controls()` option is used to control linearly for continuous variables (no factor-notation such as i.x allowed!), and the `fcontrols()` option is used for categorical / factorized variables. So, let's see it in action. 

```
scatterfit pinc_net age, binned controls(eduyrs) fcontrols(educ emp)
```
<img src="./examples/gr16.png" height="300">

You will see that, even though the variables are residualized, both x and y appear to be on their original scale. This is because scatterfit adds the mean of the variables back onto the variables are residualization so that the variables stay on a familiar scale. If you are unhappy with this behavior, use the `standardize` option to z-standardize both x and y in the plot.

## Regression parameters

The `regparameters()` option allows the user to print parameters from a regression between y and x (and, possibly, covariates) into the plot. Allowed are 'coef' for the slope coefficient, 'se' for its standard error, 'pval' for the p-value of a null-hypothesis test, 'sig' for significance stars on the coefficient (*<.1, **<.05, ***<.01), 'r2' for r-squared, 'adjr2' for adjusted r-squared, and 'nobs' for the number of observations. 
```
scatterfit pinc_net age, binned regparameters(coef se pval sig adjr2 nobs) parpos(.9 85)
```
<img src="./examples/gr17.png" height="300">

The program tries to find a good position for the parameters, but this is not always successful. In this case, I used the `parpos()` option to specify the position of the text box.

## Interaction models

If `by()` is specified, the usual behavior of scatterfit is to stratify the sample and plot the relationship between x and y separately for each by-level. Using the `bymethod(interact)` option rather specifies that the results are obtained from an interaction regression model where x is interacted with the by-factor. When no control variables are specified, this leads to identical results, while the results may vary somewhat in the presence of controls.

The interesting part is that it now becomes possible to plot regression parameters for the interaction(s). Just add the argument 'int' to the `regparameters()` option to receive the specified parameters for the interaction terms as well. 

```
scatterfit redistr pinc_net, binned by(labf) bymethod(interact) regparameters(coef pval sig int) parpos(2.9 2)
```
<img src="./examples/gr17.png" height="300">

In this case, for example, it becomes evident that those who are part of the labor force have a significantly stronger slope coefficient of income on support for redistribution.

## Changing the overall look

While scatterfit comes with a predefined look that differs considerably from the gory Stata standard graphs, the user may prefer a different look. The `plotscheme()` option changes the overall look. Put a Stata native or user-written scheme in the brackets to use it. And `colorscheme()` changes the colors. Either put a defined color palette here, such as s2color or tableau, or define your own colors using a list of hex colors (for example, `colorscheme(#E6007E #009FE3 #312783)`). 


# Creating plots with slopefit

Work in progress









