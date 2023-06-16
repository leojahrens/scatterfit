*! version 1.6.4   Leo Ahrens   leo@ahrensmail.de

program define scatterfit
version 15
	
*-------------------------------------------------------------------------------
* syntax and options
*-------------------------------------------------------------------------------
	
#delimit ;

syntax varlist(min=2 max=2)	[if] [in] [aweight fweight] , [

fit(string) BWidth(numlist max=1)  
by(string) BYMethod(string)
BINned DISCrete NQuantiles(numlist max=1) BINVar(varlist max=1) UNIBin(numlist max=1)
Controls(varlist) Fcontrols(varlist)     
BINARYModel(string)
REGParameters(string) PARPos(string) PARSize(string)
vce(string) Level(numlist max=1)
JITter(numlist max=1)
STANDardize 
LEGINside LEGSize(string)
MWeighted MSize(string) MLabel(varlist max=1)
scale(string) XYSize(string)
PLOTScheme(string asis) COLorscheme(string asis) CINTensity(numlist max=1)
opts(string asis)

/* legacy */
polybw(numlist max=1) 
COVariates(varlist) ABSorb(varlist)
coef COEFPos(string) COEFPLace(string)
] ;

#delimit cr

*-------------------------------------------------------------------------------
* install dependencies
*-------------------------------------------------------------------------------

local gtoolsu = 0
local paletteu = 0
foreach package in reghdfe gtools ftools {
	capture which `package'
	if _rc==111 & "`package'"=="gtools" local gtoolsu = 1 
	if _rc==111 ssc install `package', replace
}
if `gtoolsu'==1 gtools, upgrade
capture which colorpalette
if _rc==111 local paletteu = 1
if `paletteu'==1 ssc install colrspace, replace
if `paletteu'==1 ssc install palettes, replace
capture set scheme plotplain
if _rc==111 ssc install blindschemes, replace
capture which labmask
if _rc==111 ssc install labutil, replace


*-------------------------------------------------------------------------------
* prep
*-------------------------------------------------------------------------------

// suppress output
quietly {

// declare x and y variables
tokenize `varlist'
local y `1'
local x `2'

// weight local
if ("`weight'"!="") local w [`weight'`exp']
if ("`weight'"!="") local weightname = subinstr("`exp'","=","",.)


*-------------------------------------------------------------------------------
* check if options are correct & output errors
*-------------------------------------------------------------------------------

// report that legacy options were specified and changed internally
if "`polybw'"!=""  di as error "You specified the legacy option {bf:polybw()}, which is now {bf:bwidth()}. The setting {bf:bwidth(`polybw')} is assumed."
if "`absorb'"!="" di as error "You specified the legacy option {bf:absorb()}, which is now {bf:fcontrols()}. The setting {bf:fcontrols(`absorb')} is assumed"
if "`covariates'"!="" di as error "You specified the legacy option {bf:covariates()}, which is now {bf:controls()}. The setting {bf:controls(`covariates')} is assumed"
if "`coef'"!="" di as error "You specified the legacy option {bf:coef}, which is now governed by {bf:regparameters()}. The setting {bf:regparameters(coef pval)} is assumed."
if "`coefpos'"!="" | "`coefplace'"!=""  di as error "You specified the legacy option {bf:coefpos()} or {bf:coefplace()}, which is now {bf:parpos()}. The setting {bf:parpos(`coefplace'`coefpos')} is assumed"
if strpos("`regparameters'","beta")  di as error "You specified the legacy option {bf:regparameters(beta)}, which is now {bf:regparameters(coef)}. The setting {bf:regparameters(coef)} is assumed."

// harmonize legacy options
if "`polybw'"!="" & "`bwidth'"=="" local bwidth `polybw'
if "`absorb'"!="" & "`fcontrols'"=="" local fcontrols `absorb'
if "`covariates'"!="" & "`controls'"=="" local controls `covariates'
if "`coef'"!="" & "`regparameters'"=="" local regparameters coef pval
if "`coefplace'"!="" & "`parpos'"=="" local parpos `coefplace'
if "`coefpos'"!="" & "`parpos'"=="" local parpos `coefpos'
if strpos("`regparameters'","beta") local regparameters `regparameters' coef

// x and y variables
foreach mm in x y {
	capture confirm numeric variable ``mm''
	if _rc {
		di as error "{it:`mm'var} must be numeric."
		exit 498
	}
}
if "`binarymodel'"!="" & !("`binarymodel'"=="logit" | "`binarymodel'"=="probit") {
	di as error "{bf:binarymodel()} must contain logit or probit."
	exit 498
}
qui levelsof `y' `if', local(yval)
local yvalcount: word count `yval'
if `yvalcount'==2 {
	if "`controls'`fcontrols'"!="" {
		di as error "The use of covariates with a binary dependent variable is not supported."
		exit 498
	}
	if "`fit'"=="lowess" | strpos("`fit'","poly") {
		di as error "Local polynomial and lowess fits are not supported for binary dependent variables."
		exit 498
	}
	if "`binned'"=="" {
		di as error "It is advised to use the binned option with binary dependent variables."
	}
}
if "`binarymodel'"!="" & `yvalcount'!=2 {
	di as error "Logit/probit models require a binary dependent variable."
	exit 498
}

// by method
if "`bymethod'"!="" {
	if "`bymethod'"!="interact" & "`bymethod'"!="stratify" & "`bymethod'"!="int" & "`bymethod'"!="strat" {
		di as error "Only {bf:bymethod({it:stratify})} and {bf:bymethod({it:interact})} are allowed."
		exit 498
	}
	if "`by'"=="" {
		di as error "{bf:bymethod()} requires {bf:by()}."
		exit 498
	}
}

// fit specification and option combinations correct?
if !("`fit'"=="" | "`fit'"=="lfit" | "`fit'"=="lfitci" | "`fit'"=="qfit" | "`fit'"=="qfitci" | "`fit'"=="poly" | "`fit'"=="polyci" | "`fit'"=="lowess") {
	di as error "fit() must be lfit, lfitci, qfit, qfitci, poly, or polyci."
	exit 498
}
if "`bwidth'"!="" & !(strpos("`fit'","poly") | "`fit'"=="lowess") {
	di as error "The bwidth() option requires a local polynomial fit or a lowess smother as the fit line."
	exit 498
}

// binned options
if ("`discrete'"!="" | "`binvar'"!="" | "`unibin'"!="" | "`nquantiles'"!="") & "`binned'"=="" {
	di as error "The discrete, unibin(), binvar(), and nquatiles() options require the binned option."
	exit 498
}
if "`unibin'"!="" & ("`discrete'"!="" | "`binvar'"!="" | "`nquantiles'"!="") {
	di as error "The unibin() option cannot be combined with discrete, nquantiles(), or binvar()."
	exit 498
}

// CIs and standard errors
if "`controls'`fcontrols'"!="" {
	if strpos("`fit'","poly") | "`fit'"=="lowess" {
		di as error "Polynomial and lowess fits cannot be plotted when control variables are specified."
		exit 498
	}
}
if "`vce'"!="" {
	if !strpos("`fit'","ci") & "`regparameters'"=="" {
		di as error "The vce() option requires that confidence intervals are drawn by fit(lfitci) / fit(qfitci) or that regression parameters are plotted via the regparameters() option."
		exit 498
	}
	if strpos("`fit'","poly") |  "`fit'"=="lowess" {
		dis as error "The vce() option is incompatible with a polynomial or lowess fit."
		exit 498
	}
}

// coefficient print
if "`regparameters'"!="" {
	if strpos("`fit'","qfit") | strpos("`fit'","poly") |  "`fit'"=="lowess" {
		dis as error "Regression parameters can only be plotted for a linear fit."
		exit 498
	}
	if (strpos("`regparameters'","sig") | strpos("`regparameters'","pval") | strpos("`regparameters'","se")) & !strpos("`regparameters'","coef") {
		dis as error "{bf:regparameters({it:sig})}, {bf:regparameters({it:se})}, and {bf:regparameters({it:pval})} can only be used together with {bf:regparameters({it:coef})}"
		exit 498
	}
	if strpos("`regparameters'","int") & "`bymethod'"!="interact" {
		dis as error "Regression parameters of interaction terms can only be plotted if by() and bymethod(interact) is specified."
		exit 498
	}
}
if "`parpos'"!="" & "`regparameters'"=="" {
	di as error "The parpos() option requires the regparameters() option."
	exit 498
}


*-------------------------------------------------------------------------------
* drop superfluous observations and variables
*-------------------------------------------------------------------------------

// preserve original data
preserve

// clean dataset
if "`controls'`fcontrols'"!="" | "`binvar'"!="" | "`mlabel'"!="" {
	foreach v of varlist `controls' `fcontrols' `binvar' `mlabel' {
		local covdrop `covdrop' | mi(`v')
	}
}
if "`by'"!="" local bydrop | mi(`by')
drop if mi(`x') | mi(`y') `covdrop' `bydrop'
if "`if'"!="" keep `if'
if "`in'"!="" keep `in'
if strpos("`vce'","cluster") {
	local clustervar: subinstr local vce "cluster" "", all
}
keep `x' `y' `by' `controls' `fcontrols' `binvar' `mlabel' `weightname' `clustervar'


*-------------------------------------------------------------------------------
* prep x, y, and by variables
*-------------------------------------------------------------------------------

// check if suitable by variable is specified, generate one otherwise
local isthereby = 0
if "`by'"!="" {
	levelsof `by', local(byvals)
	local byvalcount: word count `byvals'
	if `byvalcount'!=1 {
		local isthereby = 1
		local byparen by(`by')
	}
}
if `isthereby'==0 & "`by'"!="" dis as error "The by() variable does not vary across {it:xvar} and {it:yvar}."
if `isthereby'==0 {
	gen sfitbyvar = 1
	local by sfitbyvar
}

// make by variable numeric if necessary
capture confirm numeric variable `by'
if _rc {
	rename `by' oldby
	egen `by' = group(oldby)
	labmask `by', values(oldby)
}

// set bymethod
if `isthereby'==1 {
	if "`bymethod'"=="" local bymethod stratify
	if strpos("`bymethod'","int") local bymethod interact
	if strpos("`bymethod'","strat") local bymethod stratify
}

// by locals
levelsof `by', local(bynum)
if `isthereby'==1 & "`bymethod'"=="stratify" {
	local bynumalt `bynum'
}
else {
	local bynumalt xxx
}

// check if by-variable is labeled
if `isthereby'==1 {
	local isbyvarlabeled = 1
	gen bylabcheck = .
	tostring bylabcheck, replace
	foreach bynum2 in `bynum' {
		local bylabrun: label (`by') `bynum2'
		replace bylabcheck = "`bylabrun'" if `by'==`bynum2'
	}
	destring bylabcheck, replace
	capture confirm numeric variable bylabcheck
	if !_rc { 
		local isbyvarlabeled = 0
	}
}

// check if dependent variable is binary & transform into a dummy if required
local binarydv = 0
if "`yval'"=="0 1" local binarydv = 1
capture confirm numeric variable `y'
if `yvalcount'==2 & ("`yval'"!="0 1" | _rc) {
	local binarydv = 1
	local ylab: variable label `y'
	rename `y' old`y'
	egen `y' = group(old`y')
	replace `y' = `y'-1
	lab var `y' "`ylab'"
}

// retrieve names and labels from variables
local xlab: variable label `x'
local xtitle xtitle("`xlab'")
local ylab: variable label `y'
local ytitle ytitle("`ylab'")
if `binarydv'==1 {
	cap confirm variable old`y'
	if _rc {
		local ylabber `y'
	}
	else {
		local ylabber old`y'
	}
	levelsof `ylabber', local(oldyvals)
	foreach kk in `oldyvals' {
		local ylab: label (`ylabber') `kk'
	}
	gen ytitletest = "`ylab'"
	destring ytitletest,replace
	capture confirm numeric variable ytitletest
	if !_rc { 
		local ylab: variable label `y'
	}
	else {
		local ylab Pr(`ylab')
	}
	local ytitle ytitle("`ylab'")
}
if "`ylab'"=="" local ytitle ytitle("`y'")
if "`xlab'"=="" local xtitle xtitle("`x'")

if `isthereby'==1 {
	foreach bynum2 in `bynum' {
		local bylab`bynum2': label (`by') `bynum2'
	}
}

// standardize x and y
if "`standardize'"!="" {
	if `binarydv'==0 gstats transform (standardize) `y' `x' `w', replace
	if `binarydv'==1 gstats transform (standardize) `x' `w', replace
}


*-------------------------------------------------------------------------------
* generate / specify bin variable
*-------------------------------------------------------------------------------

if "`binned'"!="" {
	if "`nquantiles'"=="" local nquantiles 30
	if "`binvar'"=="" {
		if "`discrete'"=="" & "`unibin'"=="" { // quantiles
			gquantiles `x'_q = `x' `w', xtile nq(`nquantiles') `byparen'  
		}
		else if "`unibin'"=="" {  // discrete
			clonevar `x'_q = `x'  
		}
		else {  // uniform bin
			local `unibin' = `unibin'+1
			gen `x'_q = .
			local bincount = 0
			sum `x' `w'
			range binrange r(min) r(max) `unibin'
			foreach bb of numlist 1/`unibin' {
				local bincount = `bincount'+1
				qui sum binrange if _n==`bb'+1
				local binrange1 = r(mean)
				qui sum binrange if _n==`bb'
				replace `x'_q = `bb' if `x'>=r(mean) & `x'<`binrange1'
			}
		}
	}
	if "`binvar'"!="" {
		clonevar `x'_q = `binvar'
	}
}

if "`binned'"=="" & "`mweighted'"!="" { 
	gen `x'_q = `x'
}

if "`binned'"!="" | "`mweighted'"!="" {
	egen xbin = group(`x'_q `by')
}


*-------------------------------------------------------------------------------
* specify if own predictions & CIs should be used in final plot, and which
*-------------------------------------------------------------------------------

if ((`isthereby'==1 & "`bymethod'"=="interact") | "`controls'`fcontrols'"!="" | "`vce'"!="" | `binarydv'==1) & !strpos("`fit'","poly") & "`fit'"!="lowess" {
	local ownpred_plot = 1
}
else {
	local ownpred_plot = 0
}

// confidence level
if strpos("`fit'","ci") {
	if "`level'"=="" local level 95
	local cilvl level(`level')
	local cilvltext `level'%
}

*-------------------------------------------------------------------------------
* estimate and store regression
*-------------------------------------------------------------------------------

if `ownpred_plot'==1 | "`regparameters'"!="" {

// prepare specification
	local `x'marg c.`x'
	if strpos("`fit'","qfit") local `x'marg ``x'marg'##c.`x'
	if "`bymethod'"=="interact" local `x'marg ``x'marg'##i.`by'
	local hdfeabsorb noabsorb
	if "`fcontrols'"!="" local hdfeabsorb absorb(`fcontrols')
	if "`vce'"!="" local vce vce(`vce')
	if "`binarymodel'"=="" local binarymodel logit
	
	// estimate
	if "`bymethod'"=="stratify" {
		foreach bynum2 in `bynum' {
			if `binarydv'==0 reghdfe `y' ``x'marg' `controls' if `by'==`bynum2' `w', `hdfeabsorb' `vce' 
			if `binarydv'==1 `binarymodel' `y' ``x'marg' if `by'==`bynum2' `w', `vce' 
			est sto regmodel`bynum2'
		}
	}
	else {
		if `binarydv'==0 reghdfe `y' ``x'marg' `controls' `w', `hdfeabsorb' `vce' 
		if `binarydv'==1 `binarymodel' `y' ``x'marg' `w', `vce' 
		est sto regmodel
	}
}

	
*-------------------------------------------------------------------------------
* gather regression parameters
*-------------------------------------------------------------------------------

if "`regparameters'"!="" {

	// gather regression parameters
	foreach bynum2 in `bynumalt' {
		if `isthereby'==1 & "`bymethod'"=="stratify" local bynumname `bynum2'

		* r2, nobs
		est res regmodel`bynumname'
		local nobs`bynumname' = e(N)
		if `binarydv'==0 {
			local r2`bynumname' = e(r2)
			local adjr2`bynumname' = e(r2_a)
		}
		if `binarydv'==1 {
			local r2`bynumname' = e(r2_p)
			local adjr2`bynumname' = e(r2_p)
		}
	}

	* coef, se, pval
	foreach bynum2 in `bynum' {
		if `isthereby'==1 & "`bymethod'"=="stratify" local bynumname `bynum2'
		if `isthereby'==1 & "`bymethod'"=="interact" local atbymarg at(`by'=`bynum2')
		est res regmodel`bynumname'
		margins, dydx(`x') `atbymarg' 
		local coef`bynum2' = r(table)[1,1]
		local se`bynum2' = r(table)[2,1]
		local pval`bynum2' = r(table)[4,1]
		
		if strpos("`regparameters'","sig") {
			local siglevel`bynum2' = 99
			if `pval`bynum2''<.1 & `pval`bynum2''>.05 local siglevel`bynum2' = .1
			if `pval`bynum2''<.05 & `pval`bynum2''>.01 local siglevel`bynum2' = .05
			if `pval`bynum2''<.01 local siglevel`bynum2' = .01
		}
	}

	* coefs, pval, and sig of interaction terms 
	if `isthereby'==1 & "`bymethod'"=="interact" & strpos("`regparameters'","int") {
		local count = 0
		local intmargchecker " "
		est res regmodel
		margins, dydx(`x') at(`by'=(`bynum')) post 
		foreach bynum2 in `bynum' {
			local count = `count'+1
			local count2 = 0
			foreach bynum3 in `bynum' {
				local count2 = `count2'+1
				if `bynum2'!=`bynum3' {
					if !strpos("`intmargchecker'","`bynum3'`bynum2'") {
						local intmargchecker "`intmargchecker' `bynum2'`bynum3'"
						local coef`bynum2'_`bynum3' = _b[`count'._at]-_b[`count2'._at]
						test _b[`count'._at] = _b[`count2'._at]
						local pval`bynum2'_`bynum3' = r(p)
						local parbyintcompiler `parbyintcompiler' coef`bynum2'_`bynum3' pval`bynum2'_`bynum3'
						if strpos("`regparameters'","sig") {
							local siglevel`bynum2'_`bynum3' = 99
							if `pval`bynum2'_`bynum3''<.1 & `pval`bynum2'_`bynum3''>.05 local siglevel`bynum2'_`bynum3' = .1
							if `pval`bynum2'_`bynum3''<.05 & `pval`bynum2'_`bynum3''>.01 local siglevel`bynum2'_`bynum3' = .05
							if `pval`bynum2'_`bynum3''<.01 local siglevel`bynum2'_`bynum3' = .01
						}
					}
				}

			}
		}
	}

	// round and store parameters
	foreach bynum2 in `bynum' {
		foreach par in coef pval se {
			local parbycompiler `parbycompiler' `par'`bynum2'
		}
	}	
	foreach bynum2 in `bynumalt' {
		if `isthereby'==1 & "`bymethod'"=="stratify" local bynumname `bynum2'
		foreach par in r2 adjr2 nobs {
			local parbycompiler `parbycompiler' `par'`bynumname'
		}
	}
	
	foreach par in `parbycompiler' `parbyintcompiler' {
		local smallround`par' = 0
		if ``par''>=10 | ``par''<=-10 {
			local `par'round "1"
		}
		else if ``par''>=1 | ``par''<=-1 {
			local `par'round ".1"
		}
		else {
			local roundcount = 0
			local `par'string = subinstr("``par''","-","",.)
			local `par'round ".0"
			foreach rr of numlist 2/6 {
				if substr("``par'string'",`rr',1)!="0" {
					local `par'round "``par'round'1"
					continue, break
				}
				else {
					local `par'round "``par'round'0"
					local roundcount = `roundcount'+1
				}
			}
		}
		cap if strpos("``par'string'","e") & ``par''>0 local smallround`par' = 1
		cap if strpos("``par'string'","e") & ``par''<0 local smallround`par' = -1
		local `par' = round(``par'',``par'round')
		if `smallround`par''==0 {
			local `par' "=``par''"
		}
		else if `smallround`par''==1 {
			local `par' "<.00001"
		}
		else {
			local `par' "{&cong}0"
		}
		if strpos("``par''","000000") & "``par''"!="<.00001" {
			foreach zz of numlist 1/9 {
				if substr("``par''",`zz',1)=="." local dotpos = `zz'
			}
			if "`dotpos'"!="" {
				if "``par'round'"=="1" local `par' = substr("``par''",1,`dotpos'-1)
				if "``par'round'"==".1" local `par' = substr("``par''",1,`dotpos'+1)
				if "``par'round'"==".01" local `par' = substr("``par''",1,`dotpos'+2)
				if "``par'round'"==".001" local `par' = substr("``par''",1,`dotpos'+3)
				if "``par'round'"==".0001" local `par' = substr("``par''",1,`dotpos'+4)
				if "``par'round'"==".00001" local `par' = substr("``par''",1,`dotpos'+5)
			}
		}
	}
}


*-------------------------------------------------------------------------------
* gather prediction & CIs from regression
*-------------------------------------------------------------------------------

if `ownpred_plot'==1 {
	
	// get x var
	if `binarydv'==1 | "`binned'"!="" gen `x'2 = `x'
	if `binarydv'==0 & "`binned'"=="" {
		reghdfe `x' `controls' `w', `hdfeabsorb' res(`x'2)
		sum `x' `w'
		replace `x'2 = `x'2 + r(mean)
	}
	
	// options for linear prediction estimation
	if `binarydv'==0 local atmean atmeans
	if !strpos("`fit'","qfit") 	local margpoints 30
	if strpos("`fit'","qfit") 	local margpoints 50
	
	foreach bynum2 in `bynum' {
		
		// gen empty variables for storage
		foreach vv in cix ciu cil pe {
			gen `vv'`bynum2' = .
		}

		// prepare settings for margins command
		local mcount = 0
		sum `x'2 if `by'==`bynum2' `w'
		range range`bynum2' r(min) r(max) `margpoints'
		foreach p of numlist 1/`margpoints' {
			local mcount = `mcount'+1
			qui sum range`bynum2' if _n==`p'
			local ranger`mcount' = r(mean)
			local margat`bynum2' `margat`bynum2'' `ranger`mcount''
		}
		
		// predict and store values
		if "`bymethod'"=="stratify" local regmodelname `bynum2'
		if `isthereby'==1 & "`bymethod'"=="interact" local margbysetting `by'=(`bynum2')
		est res regmodel`regmodelname'	
		margins, at(`x'=(`margat`bynum2'') `margbysetting') `atmean' `cilvl'  post
		foreach num of numlist 1/`mcount' {
			replace cix`bynum2' = `ranger`num'' if _n==`num'
			replace cil`bynum2' = r(table)[5,`num'] if _n==`num'
			replace ciu`bynum2' = r(table)[6,`num'] if _n==`num'
			replace pe`bynum2' = r(table)[1,`num'] if _n==`num'
		}
	}
}


*-------------------------------------------------------------------------------
* covariate adjustment of scatter points
*-------------------------------------------------------------------------------

if "`controls'`fcontrols'"!="" {

	if "`binned'"=="" {
		foreach v in `y' `x' {
			gen `v'_r = .
			foreach bynum2 in `bynumalt' {
				if "`bymethod'"=="stratify" local ifbysetting if `by'==`bynum2'
				reghdfe `v' `controls' `ifbysetting' `w', `hdfeabsorb' res(`v'_r`bynum2')
				sum `v' `ifbysetting' `w'
				replace `v'_r = `v'_r`bynum2' + r(mean) `ifbysetting'
			}
		}
		local y `y'_r 
		local x `x'_r
	}

	if "`binned'"!="" {
		gen `y'_r = .
		foreach bynum2 in `bynumalt' {
			if "`bymethod'"=="stratify" local ifbysetting if `by'==`bynum2'
			reghdfe `y' `controls' i.xbin `ifbysetting' `w', `hdfeabsorb'
			predict `y'_r`bynum2' if e(sample), xb
			if "`controls'"!="" {
				foreach v of varlist `controls' {
					replace `y'_r`bynum2' = `y'_r`bynum2' - _b[`v']*`v' `ifbysetting'
				}
			}
			sum `y' `w' `ifbysetting'
			local adjm = r(mean)
			sum `y'_r`bynum2' `w' `ifbysetting'
			replace `y'_r`bynum2' = `y'_r`bynum2' + (`adjm'-r(mean)) `ifbysetting'
			replace `y'_r = `y'_r`bynum2' `ifbysetting'
		}
		local y `y'_r
	}
}


*-------------------------------------------------------------------------------
* mean within bins
*-------------------------------------------------------------------------------

if "`binned'"!="" {
	gegen `y'_mean = mean(`y') `w', by(xbin)
	gegen `x'_mean = mean(`x') `w', by(xbin)
	if "`mweighted'"!="" gegen scw = count(`y'), by(xbin)
	egen tag = tag(xbin) if !mi(`y'_mean)
	replace `y'_mean = . if tag!=1
}


*-------------------------------------------------------------------------------
* specify variable to be plotted depending on binned / non-binned
*-------------------------------------------------------------------------------

if "`binned'"=="" {
	local yplot `y'
	local xplot `x'
}
else {
	local yplot `y'_mean
	local xplot `x'_mean
}

count if !mi(`xplot') & !mi(`yplot') & !mi(`by')
local n = r(N)


*-------------------------------------------------------------------------------
* color palette
*-------------------------------------------------------------------------------

if "`colorscheme'"=="" {
	local cpal `" "210 0 0" "49 113 166" "15 137 1" "255 127 14" "169 58 228" "41 217 231" "250 238 22"  "222 115 50" "'
}
else {
	local cpal `colorscheme'
}
if "`cintensity'"=="" & "`colorscheme'"=="" local cpalo int(1.2)
if "`cintensity'"=="" & "`colorscheme'"!="" local cpalo int(1)
if "`cintensity'"!="" local cpalo int(`cintensity')

if "`colorscheme'"!="" | "`plotscheme'"=="" {
	colorpalette `cpal', `cpalo' nograph local(,prefix(c) nonames)
	foreach i of numlist 25 30 50 75 {
		colorpalette `cpal', `cpalo' op(`i') nograph local(,prefix(c) suffix(o`i') nonames) 
	}
}

*-------------------------------------------------------------------------------
* overall plot scheme
*-------------------------------------------------------------------------------

if "`msize'"!="" {
	if strpos("`msize'","*") local msize = subinstr("`msize'","*","",.)
	local mresize *`msize'
}

if "`plotscheme'"=="" {
	local scheme scheme(plotplain) graphregion(lc(white) lw(vthick)) title(,size(medium)) ///
	ysc(lc(gs5) lw(thin)) ylab(#6, labs(*1.05) tlc(gs5) tlw(thin) glc(gs13) glp(solid) glw(thin) gmin gmax) ///
	xsc(lc(gs5) lw(thin)) xlab(#6, labs(*1.05) tlc(gs5) tlw(thin) glc(gs13) glp(solid) glw(thin) gmin gmax)
	
	foreach i of numlist 1/8 {
		local mlines`i' lc(`c`i'') lw(medthick)
		local mlinesci`i' acol(`c`i'o50') alw(none) clc(`c`i'') clw(medthick)
		local ciareas`i' lw(none) fc(`c`i'o30')
		local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') mlabc(`c`i'') mlw(thin)
		local efullscatterm`i' `mfullscatterm`i''
		if `i'<5 local howthin mlw(vthin)
		if `i'>=5 local howthin mlw(medium)
		local mscatter75m`i' mfc(`c`i'o50') mlc(`c`i'o75') mlabc(`c`i'o80') `howthin' mlalign(inside)
		local escatter75m`i' `mscatter75m`i''
		foreach h of numlist 25 50  {
			if `i'<5 {
				local hollow mlw(none) mlalign(outside)
			}
			else {
				local hollow mlw(medium) mlalign(inside) mlc(`c`i'o`h'')
			}
			local mscatter`h'm`i' mfc(`c`i'o`h'') mlabc(`c`i'o`h'') `hollow'
			local escatter`h'm`i' `mscatter`h'm`i''
		}
	}

	local nf = `n'
	local nfmin = 18 
	local nfmax = 500
	if `nf'<`nfmin' local nf = `nfmin'
	if `nf'>`nfmax' local nf = `nfmax'
	local globmsize = ((((`nfmax'+1-`nf')*(1/`nfmax'))^30)+(`nfmax'*.0005))*(2.6^1.5)

	local osize = 1
	local tsize = .98
	local ssize = .85
	local dsize = .8
	
	if "`mweighted'"!="" local mweightedresize *.3
	if `isthereby'==1 local mbyresize *.8
	local hollowresize = .875

	foreach sizeloc in osize tsize ssize dsize {
	    local e`sizeloc' = ``sizeloc''*`globmsize'
		local m`sizeloc' = ``sizeloc''*`globmsize'`mresize'`mweightedresize'`mbyresize'
	    local he`sizeloc' = `e`sizeloc''*`hollowresize'
		local hm`sizeloc' = `m`sizeloc''*`hollowresize'
	}

	foreach en in e m {
		if `isthereby'==0 {
			local m1 m(d) msize(*``en'dsize')
			local m2 m(o) msize(*``en'osize')
		}
		else {
			local m1 m(o) msize(*``en'osize')
			local m2 m(d) msize(*``en'dsize')
		}
		foreach g in `en'fullscatterm `en'scatter25m `en'scatter50m `en'scatter75m { 
			local `g'1 ``g'1' `m1'
			local `g'2 ``g'2' `m2'
			local `g'3 ``g'3' m(t) msize(*``en'tsize')
			local `g'4 ``g'4' m(s) msize(*``en'ssize')
			local `g'5 ``g'5' m(oh) msize(*`h`en'osize')
			local `g'6 ``g'6' m(dh) msize(*`h`en'dsize')
			local `g'7 ``g'7' m(th) msize(*`h`en'tsize')
			local `g'8 ``g'8' m(sh) msize(*`h`en'dsize')
		}
	}

	foreach g in mlines mlinesci {
		local `g'1 ``g'1' lp(solid)
		local `g'2 ``g'2' lp(dash)
		local `g'3 ``g'3' lp(shortdash)
		local `g'4 ``g'4' lp("-.")
		local `g'5 ``g'5' lp("_-")
		local `g'6 ``g'6' lp(longdash)
		local `g'7 ``g'7' lp("_.")
		local `g'8 ``g'8' lp("--.")
	}
}

else {
	local scheme scheme(`plotscheme')
	if "`colorscheme'"!="" {
		foreach i of numlist 1/8 {
			local mlines`i' lc(`c`i'') lw(medthick)
			local mlinesci`i' acol(`c`i'o50') alw(none) clc(`c`i'') clw(medthick)
			local ciareas`i' lw(none) fc(`c`i'o30')
			local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') mlabc(`c`i'') mlw(thin) msize(`mresize')
			local efullscatterm`i' `mfullscatterm`i''
			local mscatter75m`i' mfc(`c`i'o50') mlc(`c`i'o75') mlabc(`c`i'o80') mlw(vthin) mlalign(inside) msize(`mresize')
			local escatter75m`i' `mscatter75m`i''
			foreach h of numlist 25 50  {
				local mscatter`h'm`i' mfc(`c`i'o`h'') mlabc(`c`i'o`h'') mlalign(outside) mlw(none) msize(`mresize')
				local escatter`h'm`i' `mscatter`h'm`i''
			}
		}
	}
}


*-------------------------------------------------------------------------------
* refine scatter markers
*-------------------------------------------------------------------------------

// marker size weight
if "`mweighted'"!="" {
	if "`binned'"=="" gegen scw = count(`y'), by(xbin)
	sum scw 
	replace scw = scw/r(mean)
	local scw [w=scw]
}

// different scatter marker opacity depending on number of data points
foreach en in e m {
	if `n'<=100 local `en'scattermarkers `en'fullscatterm
	if `n'>100 & `n'<=300 local `en'scattermarkers `en'scatter75m
	if `n'>300 & `n'<=2000 local `en'scattermarkers `en'scatter50m
	if `n'>2000 local `en'scattermarkers `en'scatter25m
}

// labeled scatter markers 
if "`mlabel'"!="" {
	if "`plotscheme'"=="" {
		local mlabsizeloc = `globmsize'`mresize'*.7
		local mlabsizeloc mlabsize(*`mlabsizeloc')
	}
	else {
		if "`msize'"!="" local mlabsizeloc mlabsize(`mresize')
	}
	local markerlab mlabel(`mlabel') `mlabsizeloc' mlabp(0) ms(none)
}


*-------------------------------------------------------------------------------
* fit line
*-------------------------------------------------------------------------------

if `ownpred_plot'==0 {
	if "`fit'"=="lfit" | "`fit'"=="" local fittype lfit
	if "`fit'"=="lfitci" local fittype lfitci
	if "`fit'"=="qfitci" local fittype qfitci
	if "`fit'"=="qfit" local fittype qfit
	if "`fit'"=="polyci" local fittype lpolyci
	if "`fit'"=="poly" local fittype lpoly
	if "`fit'"=="lowess" local fittype lowess
}

if "`fit'"=="" | !strpos("`fit'","ci") {
	foreach ff of numlist 1/9 {
		local o`ff' `mlines`ff''
	}
}
else {
	foreach ff of numlist 1/9 {
		local o`ff' `mlinesci`ff''
	}
}
	

*-------------------------------------------------------------------------------
* put regression parameters in plot
*-------------------------------------------------------------------------------

local wherecoef = 0
if "`regparameters'"!="" {
	
	// figure out where to position the box
	sum `yplot',d
	local ymax = r(max)
	local ymin = r(min)
	local y25 = r(p25)
	local y75 = r(p75)
	count if `yplot'>`y75'
	local y75larger = r(N)
	count if `yplot'<`y25'
	local y25smaller = r(N)
	sum `xplot',d
	local xmax = r(max)
	local xmin = r(min)
	local xmean = r(mean)
	local x75 = r(p75)
	count if `xplot'>`x75'
	local x75larger = r(N)
	count if `xplot'<`x75'
	local x75smaller = r(N)

	count if `yplot'>`y75' & `xplot'<`y75'
	if (r(N)/`n')<.03 {
		local textplacey `ymax'
		local textplacex `xmean'
	}
	else {
		count if `yplot'>`y75' & `xplot'>`y75'
		if (r(N)/`n')<.03 {
			local textplacey `ymax'
			local textplacex `xmax'
		}
		else {
			count if `yplot'<`y25' & `xplot'>`y75'
			if (r(N)/`n')<.03 {
				local textplacey `ymin'
				local textplacex `xmax'
				local wherecoef = 5
			}
			else {
				local textplacey `ymax'
				local textplacex `xmean'
			}
		}
	}
	local textplace `textplacey' `textplacex'
	if "`parpos'"!="" local textplace `parpos'
	

	// text box size
	local parresize size(*.8)
	if "`parsize'"!="" {
	    if  strpos("`parsize'","*") local parsize = subinstr("`parsize'","*","",.)
		local parresize size(*`parsize')
	}
	
	// adapt to by-setting 
	if `isthereby'==1 {
		foreach bynum2 in `bynum' {
			if "`bymethod'"=="stratify" {
				if `isbyvarlabeled'==1 local parbylab`bynum2' {sub:`bylab`bynum2''}
				if `isbyvarlabeled'==0 local parbylab`bynum2' {sub:`by'=`bynum2'}
				if `isbyvarlabeled'==1 local parbylabn`bynum2' {sub:`bylab`bynum2''}
				if `isbyvarlabeled'==0 local parbylabn`bynum2' {sub:`by'=`bynum2'}
			}
			else {
				if `isbyvarlabeled'==1 local parbylab`bynum2' {sub:`bylab`bynum2''}
				if `isbyvarlabeled'==0 local parbylab`bynum2' {sub:`by'=`bynum2'}
				if `isbyvarlabeled'==1 local parbylabn`bynum2'
				if `isbyvarlabeled'==0 local parbylabn`bynum2'
			}
			foreach bynum3 in `bynum' {
				if `isbyvarlabeled'==1 local parbylab`bynum2'_`bynum3' {sub:`bylab`bynum2''*`bylab`bynum2''}
				if `isbyvarlabeled'==0 local parbylab`bynum2'_`bynum3' {sub:`by'=`bynum2'*`by'=`bynum3'}
			}
		}
	}

	// compile the parameters 
	foreach bynum2 in `bynumalt' {
		if `isthereby'==1 & "`bymethod'"=="stratify" local bynumname `bynum2'
		
		* r2
		if strpos("`regparameters'","r2") & !strpos("`regparameters'","adjr2") {
			if `binarydv'==0 local r2par `r2par' {it:R2`parbylabn`bynum2''}`r2`bynumname''
			if `binarydv'==1 local r2par `r2par' {it:Pseudo R2`parbylabn`bynum2''}`r2`bynumname''
		}
		* adj r2
		if strpos("`regparameters'","adjr2") {
			if `binarydv'==0 local adjr2par `adjr2par' {it:Adj. R2`parbylabn`bynum2''}`adjr2`bynumname''
			if `binarydv'==1 local adjr2par `adjr2par' {it:Pseudo R2`parbylabn`bynum2''}`r2`bynumname''
		}
		* observations 
		if strpos("`regparameters'","nobs") local nobspar `nobspar' {it:N`parbylabn`bynum2''}`nobs`bynumname''
	}

	foreach bynum2 in `bynum' {
		* coef, p, se
		if strpos("`regparameters'","coef") {
			if strpos("`regparameters'","sig") {
				if `siglevel`bynum2''==.1 local sigstar`bynum2' *
				if `siglevel`bynum2''==.05 local sigstar`bynum2' **
				if `siglevel`bynum2''==.01 local sigstar`bynum2' ***
			}
			if strpos("`regparameters'","se") | strpos("`regparameters'","pval") {
				if strpos("`regparameters'","se") local separ`bynum2' {it:se}`se`bynum2''
				if strpos("`regparameters'","pval") local pvalpar`bynum2' {it:p}`pval`bynum2''
				if strpos("`regparameters'","se") & strpos("`regparameters'","pval") local leerz ", "
				local sepvalpar`bynum2' " (`separ`bynum2''`leerz'`pvalpar`bynum2'')"
			}
			
			if `binarydv'==0 local whatcoef {it:{&beta}`parbylab`bynum2''}
			if `binarydv'==1 local whatcoef {it:{&delta}Pr/{&delta}x`parbylab`bynum2''}
			local coefpar`bynum2' `whatcoef'`coef`bynum2''`sigstar`bynum2''
			local jcoefpar`bynum2' `coefpar`bynum2''`sepvalpar`bynum2''
			local jcoefparcollect `jcoefparcollect' jcoefpar`bynum2'
			local coefsepvalparcollect `coefparcollect' par`bynum2'
		}
	}
	
	* interaction parameters
	if strpos("`regparameters'","int") & "`bymethod'"=="interact" {
		foreach bynum2 in `bynum' {
			local count = `count'+1
			local count2 = 0
			foreach bynum3 in `bynum' {
				local count2 = `count2'+1
				if `bynum2'!=`bynum3' {
					if !strpos("`intmargchecker'","`bynum3'`bynum2'") {
						local intmargchecker "`intmargchecker' `bynum2'`bynum3'"
						
						if strpos("`regparameters'","int") {
							if strpos("`regparameters'","sig") {
								if `siglevel`bynum2'_`bynum3''==.1 local sigstar`bynum2'_`bynum3' *
								if `siglevel`bynum2'_`bynum3''==.05 local sigstar`bynum2'_`bynum3' **
								if `siglevel`bynum2'_`bynum3''==.01 local sigstar`bynum2'_`bynum3' ***
							}
							if strpos("`regparameters'","pval") local intpvalpar`bynum2'_`bynum3' " ({it:p}`pval`bynum2'_`bynum3'')"
							if `binarydv'==0 local whatcoef {it:{&beta}`parbylab`bynum2''}-{it:ß`parbylab`bynum3''}
							if `binarydv'==1 local whatcoef {it:{&delta}Pr/{&delta}x`parbylab`bynum2''}-{it:{&delta}Pr/{&delta}x`parbylab`bynum3''}
							
							local intcoefpar`bynum2'_`bynum3' `whatcoef'`coef`bynum2'_`bynum3''`sigstar`bynum2'_`bynum3''`intpvalpar`bynum2'_`bynum3''
							local intparcollect `intparcollect' intcoefpar`bynum2'_`bynum3'
						}
					}
				}
			}
		}
	}

	if strpos("`regparameters'","r2") & !strpos("`regparameters'","adjr2") local r2par `" "`r2par'" "'
	if strpos("`regparameters'","adjr2") local adjr2par `" "`adjr2par'" "'	
	if strpos("`regparameters'","nobs") local nobspar `" "`nobspar'" "'
	foreach ee in `intparcollect' {
		local intpar `" `intpar' "``ee''" "'
	}
	if `isthereby'==1  {
		foreach ee in `jcoefparcollect' {
			local coefpar `" `coefpar' "``ee''" "'
		}
	}
	else {
		foreach ee in `coefsepvalparcollect' {
			local coefpar `" `coefpar' "`coef`ee''" "'
			if strpos("`regparameters'","se") local coefpar `" `coefpar' "`se`ee''" "'
			if strpos("`regparameters'","pval") local coefpar `" `coefpar' "`pval`ee''" "'
		}
	}

	// final text box
	local printcoef2 text(`textplace' `coefpar' `intpar' `r2par' `adjr2par' `nobspar', ///
	placement(center) `parresize' box fc(white) lc(gs5) lw(thin) la(outside) margin(vsmall) alignment(middle) linegap(.3))
	
}

*-------------------------------------------------------------------------------
* legend
*-------------------------------------------------------------------------------

// overall legend options
if "`legsize'"=="" local legresize size(*1.05)
if "`legsize'"!="" {
	if strpos("`legsize'","*") local legsize = subinstr("`legsize'","*","",.)
	local legresize size(*`legsize')
}

local legtype region(lc(white)) pos(3)
if "`leginside'"!="" {
    local leginsidepl 5
	if `wherecoef'==5 local leginsidepl 1
    local legtype ring(0) region(lc(gs5) fc(white)) pos(`leginsidepl')
}

local legopts legend(`legtype' `legresize')

// compile labels and legend options
if "`binned'"=="" local leg_obs Observed
if "`binned'"!="" local leg_obs Bin means
if strpos("`fit'","lfit") | "`fit'"=="" local leg_fit Linear
if strpos("`fit'","qfit") local leg_fit Quadratic
if strpos("`fit'","poly") local leg_fit Local polynomial
if strpos("`fit'","lowess") local leg_fit Lowess

foreach ii of numlist 1/4 {
    local n`ii' = `ii'
	if "`mweighted'"!="" local n`ii' = `n`ii''+1
}

if `isthereby'==0 {
	if "`mlabel'"=="" local firstlegel 1 "`leg_obs'" 
	if !strpos("`fit'","ci") {
		local legopts `legopts' legend(order(`firstlegel' `n2' "`leg_fit' fit"))
	}
	else {
		local legopts `legopts' legend(order(`firstlegel' `n3' "`leg_fit' fit" `n2' "`cilvltext' CIs"))
	}
}

if `isthereby'==1 {
	egen distinctby = group(`by')
	sum distinctby
	local maxdistinctby = r(max)
	local coln = 0
	levelsof `by', local(bynum)
	foreach bynum2 in `bynum' {
		if `isbyvarlabeled'==1 local legbyvarl`bynum2' `bylab`bynum2''
		if `isbyvarlabeled'==0 local legbyvarl`bynum2' `by'==`bynum2'
		local coln = `coln'+1
		if `ownpred_plot'==0 {
		    if strpos("`fit'","ci") local ciadd "+`coln'"
			local coln2 = `coln'+`maxdistinctby'`ciadd'
		}
		else {
		    if strpos("`fit'","ci") local ciadd "+`maxdistinctby'"
			local coln2 = `coln'+`maxdistinctby'`ciadd'
		}
		if "`mweighted'"!="" local coln2 = `coln2'+`maxdistinctby'
		if "`mlabel'"=="" local legorder `legorder' `coln' "`legbyvarl`bynum2''" `coln2' ""
		if "`mlabel'"!="" local legorder `legorder' `coln2' "`legbyvarl`bynum2''"
	}
	if "`mlabel'"=="" local legcol 2
	if "`mlabel'"!="" local legcol 1
	local legopts `legopts' legend(order(`legorder') col(`legcol') textfirst)
}


*-------------------------------------------------------------------------------
* overall plot size
*-------------------------------------------------------------------------------

if strpos("`scale'","*") {
	local scale2 = subinstr("`scale'","*","",.)
	local scale = `scale2'
}
if strpos("`xysize'","*") {
	local xysize2 = subinstr("`xysize'","*","",.)
	local xysize = `xysize2'
}

if "`scale'"!="" local plotsize scale(`scale')
if "`xysize'"!="" {
	if `xysize'<=1 {
		local ysize = 100
		local xsize = 100*`xysize'
	}
	if `xysize'>1 {
		local xsize = 100
		local ysize = 100*(1/`xysize')
	}
	local plotsize `plotsize' xsize(`xsize') ysize(`ysize')
}


*-------------------------------------------------------------------------------
* final plot
*-------------------------------------------------------------------------------

// options
local lscatteropts `scheme' `xtitle' `ytitle' `printcoef2' `legopts' `plotsize' `opts' 
if "`bwidth'"!="" local bwidth bw(`bwidth')
if "`jitter'"!="" local jitter jitter(`jitter')

// empty scatter marker plot for correct legend in case of weighted scatter markers
if "`mweighted'"!="" {
	local coln = 0
	foreach bynum2 in `bynum' {
		local coln = `coln'+1
		local coln2 = `coln'
		if `isthereby'==0 local coln2 = `coln'+1
		local sce `sce' (scatter `yplot' `xplot' if `by'==`bynum2' & mi(`yplot'), ``escattermarkers'`coln2'')
	}
}

// compile the main plot
local coln = 0
foreach bynum2 in `bynum' {
	local coln = `coln'+1
	local coln2 = `coln'
	if `isthereby'==0 local coln2 = `coln'+1
	
	* scatter points
	local sc `sc' (scatter `yplot' `xplot' if `by'==`bynum2' `scw', ``mscattermarkers'`coln2'' `markerlab' `jitter')
	
	* fit line & CIs using off-the-shelf fitting
	if `ownpred_plot'==0 {
		local pl `pl' (`fittype' `y' `x' if `by'==`bynum2', `o`coln'' `bwidth' `cilvl')
	}
	
	* own fit lines and CIs from margins
	else {
		if strpos("`fit'","ci") local ci `ci' (rarea cil`bynum2' ciu`bynum2' cix`bynum2', `ciareas`coln'')
		local pl `pl' (line pe`bynum2' cix`bynum2', `o`coln'')
	}
	
}

// draw the final plot
tw `sce' `sc' `ci' `pl', `lscatteropts'




*-------------------------------------------------------------------------------
* finalize
*-------------------------------------------------------------------------------

restore
}
end












