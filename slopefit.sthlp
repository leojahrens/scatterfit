{smcl}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:scatterfit}}Scatter plots with fit lines that visualize the effect of {it:xvar} on {it:yvar} at different values of {it:zvar}{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:scatterfit}
yvar xvar zvar {ifin}
{weight} {cmd:,} [options] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Overall slopes}
{synopt :{opt fit(str)}}Specifies the fit line that plots the overall interaction model results. May be {it:lfit} for a linear fit, {it:lfitci} to include CIs, or {it:qfit} / {it:qfitci} for a quadratic fit.{p_end}

{syntab:Individual slopes}
{synopt :{opt m:ethod(str)}}Determines for what values of {it:zvar} individual slopes of {it:xvar} are plotted. May be {it:quantiles} to divide {it:z} into quantiles, {it:unibin} for uniformly spaced bins, or {it:discrete}.{p_end}
{synopt :{opt binv:ar(varlist)}}Plots the effects of {it:xvar} for each within each distinct value of {it:varlist}. Used when the bins are already defined in the dataset.{p_end}
{synopt :{opt indslopesm:ethod(str)}}Determines how the individual slopes are computed. May be {it:interact} for a unified regression models including interactions or {it:stratify} for separate regressions.{p_end}
{synopt :{opt indslopesci}}Includes 95% confidence intervals for the individual slopes.{p_end}

{syntab:Conditioning on covariates}
{synopt :{opt c:ontrols(varlist)}}Adjusts the relationship between {it:xvar} and {it:yvar} linearly for continuous covariates. Achieved by residualization.{p_end}
{synopt :{opt fc:ontrols(varlist)}}Controls for factor variables, i.e. categorical variables such as gender or region.{p_end}

{syntab:Uncertainty estimates}
{synopt :{opt vce(string)}}Specifies the estimated standard errors, which is relevant for confidence intervals and {it:p}-values. Supports all possible vce() options of {opt reghdfe} (continuous {it:y}) or {opt logit} (binary {it:y}).{p_end}
{synopt :{opt l:evel(num)}}Specifies the confidence level for the CIs. Standard is 95.{p_end}

{syntab:Regression parameters}
{synopt :{opt regp:arameters(str)}}Prints parameters into the plot. May contain {it:coef} for information on how the slope of {it:xvar} changes when {it:zvar} increases by one unit, {it:pval} for the associated p-value, {it:sig} for significance (*<.1, **<.05, ***<.001), and {it:nobs} for N.{p_end}
{synopt :{opt parp:os(numlist)}}Overrides the position of the parameters within the plot. For example, {opt parpos(0 5)} choses 0 as the y-coordinate and 5 as the x-coordinate.{p_end}

{syntab:Scatter points}
{synopt :{opt mw:eighted(num)}}Adjusts the marker size depending on the number of observations at distinct values of {it:yvar} and {it:xvar}. Useful with {opt discrete} and {opt unibin()}.{p_end}
{synopt :{opt ml:abel(var)}}Uses the strings / value labels stored in {it:var} as scatter markers instead of the usual circles, diamonds, etc.{p_end}
{synopt :{opt jit:ter(num)}}Randomly varies the location of scatter points.{p_end}

{syntab:Distribution of z variable}
{synopt :{opt zdistr:ibution(string)}}Adds a kernel density plot of {it:zvar}. May be set to {it:auto} or {it:a b}, where {it:a} contains the position and {it:b} the height of the density.{p_end}
{synopt :{opt zdistrbw(num)}}Specifies the bandwidth of the kernel density estimator.{p_end}

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
{synopt :{opt binarym:odel(str)}}Specifies the regression model used for binary dependent variables. Must be {it:logit} or {it:probit}.{p_end}
{synopt :{opt noy:line}}Removes the line at y=0.{p_end}


{title:Examples}

A detailed tutorial with examples is available at https://github.com/leojahrens/scatterfit

{phang2}{cmd:. sysuse auto}{p_end}
{phang2}{cmd:. slopefit weight length turn}{p_end}
{phang2}{cmd:. slopefit weight length turn, fit(lfitci)}{p_end}
{phang2}{cmd:. slopefit weight length turn, fit(lfitci) method(discrete)}{p_end}




