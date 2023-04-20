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
{synopt :{opt fit(str)}}Fit line. May be {opt lfit} for a linear fit, {opt lfitci} to include CIs, {opt qfit} / {opt qfitci} for a quadratic fit, {opt poly} / {opt polyci} for a local polynomial fit, or {opt lowess} for a lowess smoother.{p_end}
{synopt :{opt bw(num)}} Bandwidth for the local polynomial / lowess smoother fit line. {p_end}

{syntab:Multiple plots}
{synopt :{opt by(var)}}Draws separate scatterplots and fit lines for each value of {it:var}.{p_end}

{syntab:Bins}
{synopt :{opt bin:ned}}Divides {it:xvar} into equally sized bins based on quantile cutoff points and plots mean values of {it:yvar} and {it:xvar} within these bins. {p_end}
{synopt :{opt nq:uantiles(num)}}Choses the number of equally sized bins / quantiles. {p_end}
{synopt :{opt disc:rete}}Treats {it:xvar} as discrete and plots the means of {it:yvar} within each distinct value of {it:xvar}.{p_end}
{synopt :{opt unib:in(num)}}Divides {it:xvar} into {it:num} equally spaced bins and plots the means of {it:yvar} and {it:xvar} within them.{p_end}
{synopt :{opt binv:ar(varlist)}}Plots the means of {it:yvar} and {it:xvar} within each distinct value of {it:varlist}. Used when the bins are already defined in the dataset.{p_end}

{syntab:Conditioning on covariates}
{synopt :{opt c:ontrols(varlist)}}Adjusts the relationship between {it:xvar} and {it:yvar} linearly for the covariates specified in {it:varlist}. Achieved by residualization or, when {opt binned} is specified, the method by Cattaneo et al. (2022).{p_end}
{synopt :{opt fc:ontrols(varlist)}}Used for controlling for factor variables (i.e. categorical variables), such as gender.{p_end}

{syntab:Uncertainty estimates}
{synopt :{opt vce(string)}}Specifies the estimated standard errors, which is relevant for confidence intervals and the estimated {it:p}-value. Supports all possible vce() options of {opt reghdfe} / {opt logit}.{p_end}

{syntab:Regression parameters}
{synopt :{opt coef}}Prints the beta coefficient of {it:xvar} and the associated {it:p}-value into the plot. Obtained from OLS or Logit.{p_end}
{synopt :{opt coefp:lace(numlist)}}Overrides the placement of {opt coef} within the plot. For example, {opt coefplace(0 5)} choses 0 as the y-coordinate and 5 as the x-coordinate.{p_end}

{syntab:Scatter points}
{synopt :{opt mw:eighted(num)}}Adjusts the marker size depending on the number of observations at distinct values of {it:yvar} and {it:xvar}.{p_end}
{synopt :{opt ms:ize(num)}}Resizes the scatter markers. For example, {opt msize(.9)} decreases their size by 10%.{p_end}
{synopt :{opt jit:ter(num)}}Randomly varies the location of scatter points.{p_end}

{syntab:Plot and element size}
{synopt :{opt xys:ize(num)}}Specifies the relative width to height. For example, {opt xysize(1.2)} draws a plot with 20% more width than height.{p_end}
{synopt :{opt scale(num)}}Resizes all elements of the plot. For example, {opt scale(1.1)} makes all elements 10% larger.{p_end}

{syntab:Scheme and colors}
{synopt :{opt plots:cheme(str)}}Specifies an alternative graph scheme, e.g. {opt plotscheme(white_tableau)}.{p_end}
{synopt :{opt col:orscheme(str)}}Uses a custom color palette, e.g. {opt plotscheme(tableau)}. If this option is not specified, a high-intensity version of Set1 is used.{p_end}
{synopt :{opt cint:ensity(num)}}Intesity of the colors. Higher values make all colors darker.{p_end}

{syntab:Other}
{synopt :{opt stand:ardize}}Standardizes {it:xvar} and {it:yvar} so that they have a standard deviation of one and a mean of zero.{p_end}
{synopt :{opt legin:side}}Places the legend inside the plot region.{p_end}
{synopt :{opt opts(str)}}Passes on options to the twoway plot. For example, {opt opts(xtitle("Example"))} changes the x axis title. See {opt help twoway}.{p_end}




{title:Examples}

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
{phang2}{cmd:. scatterfit weight length, controls(trunk) fcontrols(foreign) coef}{p_end}
{hline}





