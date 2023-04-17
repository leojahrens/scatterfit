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

{syntab:Binned scatter plot}
{synopt :{opt bin:ned}}Divides the {it:xvar} into equally sized bins based on quantile cutoff points and plots mean values of {it:yvar} and {it:xvar} within these bins. {p_end}
{synopt :{opt nq:uantiles(num)}}Choses the number of equally sized bins / quantiles. {p_end}
{synopt :{opt disc:rete}}Treats {it:xvar} as a discrete variable and plots the means of {it:yvar} within each distinct value of {it:xvar}.{p_end}
{synopt :{opt binv:ar(varlist)}}Plots the means of {it:yvar} and {it:xvar} within each distinct value of {it:varlist}. Used when the variable (or variable combination) that defines bins already exists.{p_end}

{syntab:Conditioning on covariates}
{synopt :{opt c:ontrols(varlist)}}Adjusts the relationship between {it:xvar} and {it:yvar} linearly for the covariates specified in {it:varlist}. Achieved by residualization or, when {opt binned} is specified, a method described in Cattaneo et al.: {it:On Binscatter}.{p_end}
{synopt :{opt fc:ontrols(varlist)}}Used for controlling for factor variables (i.e. categorical variables), such as gender.{p_end}

{syntab:Beta coefficient}
{synopt :{opt coef}}Prints the beta coefficient of {it:xvar} and the associated {it:p}-value into the plot. Obtained from OLS.{p_end}
{synopt :{opt coefp:lace(numlist)}}Overrides the placement of {opt coef} within the plot. Must be a numlist of length 2. For example, {opt coefplace(0 5)} choses 0 as the y-coordinate and 5 as the x-coordinate.{p_end}

{syntab:Other}
{synopt :{opt vce(string)}}Specifies the estimated standard errors, which is relevant for confidence intervals and the estimated {it:p}-value. May be {opt robust} or {opt cluster(clustervar)}.{p_end}
{synopt :{opt jit:ter(num)}}Randomly varies the location of scatter points. Useful for strongly clustered scatter plots, e.g. with a discrete variable such as education years.{p_end}
{synopt :{opt stand:ardize}}Standardizes {it:xvar} and {it:yvar} so that they have a standard deviation of one and a mean of zero.{p_end}
{synopt :{opt opts(str)}}Passes on options to the twoway plot. For example, {opt opts(xtitle("Example"))} changes the x axis title. See {opt help twoway}.{p_end}
{synopt :{opt plots:cheme(str)}}Specifies a graph scheme to be used, e.g. {opt plotscheme(white_tableau)}. If this option is not specified, a custom scheme defined by the scatterfit program is used.{p_end}
{synopt :{opt col:orscheme(str)}}Specifies a custom color scheme, e.g. {opt plotscheme(tableau)}. If this option is not specified, a high-intensity version of Set1 is used.{p_end}



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





