# SCATTERFIT v1.6

Scatterfit includes two commands for Stata that produce a wide range of scatter plots with overlaid fit lines: scatterfit and slopefit. To install the package, execute the following command in Stata

```
net install scatterfit, from(https://raw.githubusercontent.com/leojahrens/scatterfit/master) replace
```

Once installed, please see `help scatterfit` and `help slopefit` for syntax and the whole range of options.
There are two tutorials below that showcase the possibilities of the two commands. 

# Creating plots with scatterfit

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