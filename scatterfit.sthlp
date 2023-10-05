{smcl}
{title:help scatterfit}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:scatterfit}}Scatter plots with fit lines{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:scatterfit}
yvar xvar {ifin}
{weight} {cmd:,} [options] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Fit}
{synopt :{opt fit(str)}}Functional form of the fit line. May be {it:linear}, {it:quadratic} for a second polynomial, {it:cubic} for a third polynomial, {it:lpoly} for a local polynomial fit, or {it:lowess} for a lowess smoother.{p_end}
{synopt :{opt fitm:odel(model)}}Changes the regression model used to estimate the fit line. See the table below the options for a full list of supported models.{p_end}
{synopt :{opt bw:idth(num)}} Bandwidth for the local polynomial / lowess smoother fit line.{p_end}

{syntab:Confidence intervals and standard errors}
{synopt :{opt ci}}Includes 95% confidence intervals.{p_end}
{synopt :{opt vce(string)}}Specifies the estimated standard errors, which is relevant for confidence intervals and {it:p}-values. Supports all possible vce() options of the respective regression model as well as {it:vce(wbcluster clustervar)} for wild cluster bootstrap.{p_end}
{synopt :{opt l:evel(num)}}Changes the confidence level of the CIs.{p_end}

{syntab:Bins}
{synopt :{opt bin:ned}}Divides {it:xvar} into equally sized bins based on quantile cutoff points and plots mean values of {it:yvar} and {it:xvar} within these bins. {p_end}
{synopt :{opt nq:uantiles(num)}}Choose the number of equally sized bins / quantiles.{p_end}
{synopt :{opt disc:rete}}Treats {it:xvar} as discrete and plots the means of {it:yvar} within each distinct value of {it:xvar}.{p_end}
{synopt :{opt unib:in(num)}}Divides {it:xvar} into {it:num} equally spaced bins and plots the means of {it:yvar} and {it:xvar} within them.{p_end}
{synopt :{opt binv:ar(varlist)}}Plots the means of {it:yvar} and {it:xvar} within each distinct value of {it:varlist}. Used when the bins are already defined in the dataset.{p_end}
{synopt :{opt binm:ethod(str)}}Specifies how data points are aggregated within bins. May be {it:mean} (the standard) or {it:median}.{p_end}

{syntab:Subsample analysis}
{synopt :{opt by(var)}}Draws separate scatterplots and fit lines for each value of {it:var}.{p_end}
{synopt :{opt onef:it}}Only plots a single fit line when {it:by()} is specified.{p_end}
{synopt :{opt bym:ethod(str)}}Defines whether the results are derived from separate analyses of stratified samples ({opt bymethod(stratify)}) or a unified regression model containing interaction terms ({opt bymethod(interact)}).{p_end}

{syntab:Conditioning on covariates}
{synopt :{opt c:ontrols(varlist)}}Adjusts the relationship between {it:xvar} and {it:yvar} linearly for continuous covariates. Achieved by residualization or, when {opt binned} is used, the method by Cattaneo et al. (2022).{p_end}
{synopt :{opt fc:ontrols(varlist)}}Controls for factor variables, i.e. categorical variables such as gender or region.{p_end}

{syntab:Regression parameters}
{synopt :{opt regp:arameters(str)}}Prints regression parameters into the plot. May contain {it:coef}, {it:se}, {it:pval}, {it:sig} (*<.1, **<.05, ***<.001), {it:int} (for {opt bymethod(interact)}), {it:r2} / {it:adjr2}, {it:nobs} (N).{p_end}
{synopt :{opt parp:os(numlist)}}Overrides the position of the parameters within the plot. For example, {opt parpos(0 5)} choses 0 as the y-coordinate and 5 as the x-coordinate.{p_end}

{syntab:Scatter points}
{synopt :{opt mw:eighted}}Adjusts the marker size depending on the number of observations at distinct values of {it:yvar} and {it:xvar}. Useful with {opt discrete} and {opt unibin()}.{p_end}
{synopt :{opt ml:abel(var)}}Uses the strings / value labels stored in {it:var} as scatter markers instead of the usual circles, diamonds, etc.{p_end}
{synopt :{opt jit:ter(num)}}Randomly varies the location of scatter points.{p_end}

{syntab:Distribution of x variable}
{synopt :{opt xdistr:ibution(string)}}Adds a kernel density plot of {it:xvar}. May be set to {it:auto} or {it:a b}, where {it:a} contains the position and {it:b} the height of the density.{p_end}
{synopt :{opt xdistrbw(num)}}Specifies the bandwidth of the kernel density estimator.{p_end}

{syntab:Scheme and colors}
{synopt :{opt plots:cheme(str)}}Defines an alternative graph scheme, such as {opt plotscheme(white_tableau)}. See the collection in https://github.com/asjadnaqvi/stata-schemepack.{p_end}
{synopt :{opt col:orscheme(str)}}Defines a custom color palette (e.g. {opt colorscheme(tableau)}). To define your own colors, use a list of hex colors (e.g. {opt colorscheme(#E6007E #009FE3 #312783)}).{p_end}
{synopt :{opt cint:ensity(num)}}Changes the color intensity. Higher values make colors darker.{p_end}

{syntab:Plot and element size}
{synopt :{opt scale(num)}}Resizes all elements in the plot. For example, {opt scale(1.1)} increases text, marker, and lines sizes by 10%.{p_end}
{synopt :{opt ms:ize(num)}}Resizes the scatter markers.{p_end}
{synopt :{opt pars:ize(num)}}Resizes the printed regression parameters.{p_end}
{synopt :{opt legs:ize(num)}}Resizes the legend.{p_end}
{synopt :{opt xys:ize(num)}}Specifies the relative width to height. For example, {opt xysize(1.2)} draws a plot with 20% more width than height.{p_end}

{syntab:Other}
{synopt :{opt stand:ardize}}Standardizes {it:xvar} and {it:yvar} to a standard deviation of one and a mean of zero.{p_end}
{synopt :{opt legin:side}}Places the legend inside the plot region.{p_end}
{synopt :{opt slow}}Uses (slower) native Stata commands instead of packages such as Gtools.{p_end}
{synopt :{opt opts(str)}}Passes on options to the twoway plot. For example, {opt opts(xtitle("Example"))} changes the x axis title. See {opt help twoway}.{p_end}


{marker rmodel}{...}
{title:Supported regression models}

{synopt:{it:model}}Description{p_end}
{synoptline}
{synopt:{it:ols}}Linear regression estimated via ordinary least squares. The standard for continuous DVs.{p_end}
{synopt:{it:poisson}}Poisson regression for count data{p_end}
{synopt:{it:quantile}}Quantile regression. Standard setting is a regression at the 50th percentile. For other quantiles, specify {opt fitmodel(quantile, p(x))}, where {it:x} is the respective quantile between 1 and 99.{p_end}
{synopt:{it:randomint}}Multilevel model with random intercepts. The higher-level cluster(s) must be specified with the following syntax: {opt fitmodel(randomint, cluster(x))}, where  {it:x} is the list of cluster variables.{p_end}
{synopt:{it:randomint, reml}}Estimates the model with REML and DF adjustion for unbiased estimation of SEs with few clusters (see Elff et al. in BJPS). Use the following syntax: {opt fitmodel(randomint, reml cluster(x))}{p_end}
{synopt:{it:flogit} / {it:fprobit}}Fractional logit/probit model for 0-1 DVs{p_end}
{synopt:{it:logit} / {it:probit}}Logistical regression with the logit or probit link function. {it:logit} it the standard for binary DVs.{p_end}
{synopt:{it:lpm}}Linear probability model for binary DVs, estimated with OLS{p_end}
{synoptline}


{title:Examples}

A detailed tutorial with examples is available at https://github.com/leojahrens/scatterfit


