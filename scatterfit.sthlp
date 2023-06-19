{smcl}
{title:Title}

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

{syntab:Fit line and CIs}
{synopt :{opt fit(str)}}Fit line. May be {it:lfit} for a linear fit, {it:lfitci} to include CIs, {it:qfit} / {it:qfitci} for a quadratic fit, {it:poly} / {it:polyci} for a local polynomial fit, or {it:lowess} for a lowess smoother.{p_end}
{synopt :{opt bw:idth(num)}} Bandwidth for the local polynomial / lowess smoother fit line. {p_end}

{syntab:Multiple plots}
{synopt :{opt by(var)}}Draws separate scatterplots and fit lines for each value of {it:var}.{p_end}
{synopt :{opt bym:ethod(str)}}Defines whether the results are derived from separate analyses of stratified samples ({opt bymethod(stratify)}) or a unified regression model containing interaction terms ({opt bymethod(interact)}).{p_end}

{syntab:Bins}
{synopt :{opt bin:ned}}Divides {it:xvar} into equally sized bins based on quantile cutoff points and plots mean values of {it:yvar} and {it:xvar} within these bins. {p_end}
{synopt :{opt nq:uantiles(num)}}Choose the number of equally sized bins / quantiles.{p_end}
{synopt :{opt disc:rete}}Treats {it:xvar} as discrete and plots the means of {it:yvar} within each distinct value of {it:xvar}.{p_end}
{synopt :{opt unib:in(num)}}Divides {it:xvar} into {it:num} equally spaced bins and plots the means of {it:yvar} and {it:xvar} within them.{p_end}
{synopt :{opt binv:ar(varlist)}}Plots the means of {it:yvar} and {it:xvar} within each distinct value of {it:varlist}. Used when the bins are already defined in the dataset.{p_end}

{syntab:Conditioning on covariates}
{synopt :{opt c:ontrols(varlist)}}Adjusts the relationship between {it:xvar} and {it:yvar} linearly for continuous covariates. Achieved by residualization or, when {opt binned} is used, the method by Cattaneo et al. (2022).{p_end}
{synopt :{opt fc:ontrols(varlist)}}Controls for factor variables, i.e. categorical variables such as gender or region.{p_end}

{syntab:Uncertainty estimates}
{synopt :{opt vce(string)}}Specifies the estimated standard errors, which is relevant for confidence intervals and {it:p}-values. Supports all possible vce() options of {opt reghdfe} (continuous {it:y}) or {opt logit} (binary {it:y}).{p_end}
{synopt :{opt l:evel(num)}}Specifies the confidence level for the CIs. Standard is 95.{p_end}

{syntab:Regression parameters}
{synopt :{opt regp:arameters(str)}}Prints regression parameters into the plot. May contain {it:coef}, {it:se}, {it:pval}, {it:sig} (*<.1, **<.05, ***<.001), {it:int} (for {opt bymethod(interact)}), {it:r2} / {it:adjr2}, {it:nobs} (N).{p_end}
{synopt :{opt parp:os(numlist)}}Overrides the position of the parameters within the plot. For example, {opt parpos(0 5)} choses 0 as the y-coordinate and 5 as the x-coordinate.{p_end}

{syntab:Scatter points}
{synopt :{opt mw:eighted(num)}}Adjusts the marker size depending on the number of observations at distinct values of {it:yvar} and {it:xvar}. Useful with {opt discrete} and {opt unibin()}.{p_end}
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
{synopt :{opt opts(str)}}Passes on options to the twoway plot. For example, {opt opts(xtitle("Example"))} changes the x axis title. See {opt help twoway}.{p_end}
{synopt :{opt binarym:odel(str)}}Specifies the regression model used for binary dependent variables. Must be {it:logit} or {it:probit}.


{title:Examples}

A detailed tutorial with examples is available at https://github.com/leojahrens/scatterfit

{hline}
{pstd}Load data{p_end}
{phang2}{cmd:. sysuse auto}{p_end}

{pstd}Simple scatter plot with a linear fit{p_end}
{phang2}{cmd:. scatterfit weight length}{p_end}
{hline}

{pstd}Including CIs{p_end}
{phang2}{cmd:. scatterfit weight length, fit(lfitci)}{p_end}
{hline}

{pstd}Local polynomial fit{p_end}
{phang2}{cmd:. scatterfit weight length, fit(polyci) bw(10)}{p_end}
{hline}

{pstd}Plots by foreign{p_end}
{phang2}{cmd:. scatterfit weight length, by(foreign)}{p_end}
{hline}

{pstd}With bins{p_end}
{phang2}{cmd:. scatterfit weight length, binned nq(20)}{p_end}
{hline}

{pstd}With control variables and printed coefficient{p_end}
{phang2}{cmd:. scatterfit weight length, controls(trunk) fcontrols(foreign) regparameters(coef sig pval)}{p_end}
{hline}





