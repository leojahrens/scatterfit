*! version 1.6.3   Leo Ahrens   leo@ahrensmail.de

program define slopefit
version 14.2
	
*-------------------------------------------------------------------------------
* syntax and options
*-------------------------------------------------------------------------------
	
#delimit ;

syntax varlist(min=3 max=3)	[if] [in] [aweight fweight] , [

fit(string)
Method(string)  NQuantiles(numlist max=1) NUnibin(numlist max=1) BINVar(varlist max=1)
INDSLOPESMethod(string) indslopesci
Controls(varlist) Fcontrols(varlist)     
BINARYModel(string)
REGParameters(string) PARPos(string) PARSize(string)
vce(string) Level(numlist max=1)
JITter(numlist max=1)
STANDardize 
LEGINside LEGSize(string)
NOYline
MWeighted MSize(string) MLabel(varlist max=1)
scale(string) XYSize(string)
PLOTScheme(string asis) COLorscheme(string asis) CINTensity(numlist max=1)
opts(string asis)

/* legacy */
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


*-------------------------------------------------------------------------------
* prep
*-------------------------------------------------------------------------------

// suppress output
quietly {

// declare x and y variables
tokenize `varlist'
local y `1'
local x `2'
local z `3'

// weight local
if ("`weight'"!="") local w [`weight'`exp']
if ("`weight'"!="") local weightname = subinstr("`exp'","=","",.)


*-------------------------------------------------------------------------------
* check if options are correct & output errors
*-------------------------------------------------------------------------------

// report that legacy options were specified and changed internally
if "`absorb'"!="" di as error "You specified the legacy option {bf:absorb()}, which is now {bf:fcontrols()}. The setting {bf:fcontrols(`absorb')} is assumed"
if "`covariates'"!="" di as error "You specified the legacy option {bf:covariates()}, which is now {bf:controls()}. The setting {bf:controls(`covariates')} is assumed"
if "`coef'"!="" di as error "You specified the legacy option {bf:coef}, which is now governed by {bf:regparameters()}. The setting {bf:regparameters(coef pval)} is assumed."
if "`coefpos'"!="" | "`coefplace'"!=""  di as error "You specified the legacy option {bf:coefpos()} or {bf:coefplace()}, which is now {bf:parpos()}. The setting {bf:parpos(`coefplace'`coefpos')} is assumed"
if strpos("`regparameters'","beta")  di as error "You specified the legacy option {bf:regparameters(beta)}, which is now {bf:regparameters(coef)}. The setting {bf:regparameters(coef)} is assumed."

// harmonize legacy options
if "`absorb'"!="" & "`fcontrols'"=="" local fcontrols `absorb'
if "`covariates'"!="" & "`controls'"=="" local controls `covariates'
if "`coef'"!="" & "`regparameters'"=="" local regparameters coef pval
if "`coefplace'"!="" & "`parpos'"=="" local parpos `coefplace'
if "`coefpos'"!="" & "`parpos'"=="" local parpos `coefpos'
if strpos("`regparameters'","beta") local regparameters `regparameters' coef

// x and y variables
foreach mm in x z y {
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
if "`binarymodel'"!="" & `yvalcount'!=2 {
	di as error "Logit/probit models require a binary dependent variable."
	exit 498
}

// method
if "`method'"!="" & ("`method'"!="quantiles" & "`method'"!="unibin" & "`method'"!="discrete") {
	di as error "Only {bf:method({it:quantiles})}, {bf:method({it:unibin})}, and {bf:method({it:discrete})} are allowed (if {bf:byvar()} is not specified)."
	exit 498
}
if "`binvar'"!="" & ("`method'"!="" | "`nquantiles'"!="" | "`nunibin'"!="") {
	di as error "The binvar() option cannot be combined with method(), nquantiles(), or nunibin()."
	exit 498
}
if "`nunibin'"!="" & ("`method'"!="unibin" & "`method'"!="") {
	di as error "The nunibin() option requires method(unibin)."
	exit 498
}
if "`nquantiles'"!="" & ("`method'"!="quantiles" & "`method'"!="") {
	di as error "The nquantiles() option requires method(quantiles)."
	exit 498
}
if "`nquantiles'"!="" & "`nunibin'"!="" {
	di as error "nquantiles() and nunibin() cannot be used together."
	exit 498
}
// indslopesmethod
if "`indslopesmethod'"!="" & ("`indslopesmethod'"!="interact" & "`indslopesmethod'"!="stratify") {
	di as error "Only {bf:indslopesmethod({it:stratify})} and {bf:bymethod({it:interact})} are allowed."
	exit 498
}
// fit specification 
if !("`fit'"=="" | "`fit'"=="lfit" | "`fit'"=="lfitci" | "`fit'"=="qfit" | "`fit'"=="qfitci") {
	di as error "fit() must be lfit, lfitci, qfit, or qfitci."
	exit 498
}
// CIs and standard errors
if "`vce'"!="" {
	if !strpos("`fit'","ci") & "`regparameters'"=="" {
		di as error "The vce() option requires that confidence intervals are drawn by fit(lfitci) / fit(qfitci) or that regression parameters are plotted via the regparameters() option."
		exit 498
	}
}
// coefficient print
if "`regparameters'"!="" {
	if strpos("`fit'","qfit") {
		dis as error "Regression parameters can only be plotted for a linear fit."
		exit 498
	}
	if strpos("`regparameters'","sig") & !strpos("`regparameters'","coef") {
		dis as error "{bf:regparameters({it:sig})} requires {bf:regparameters({it:coef})}"
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
drop if mi(`x') | mi(`y') | mi(`z') `covdrop'
if "`if'"!="" keep `if'
if "`in'"!="" keep `in'
keep `x' `y' `z' `controls' `fcontrols' `binvar' `mlabel' `weightname'


*-------------------------------------------------------------------------------
* prep x, y, and z variables
*-------------------------------------------------------------------------------

// check if dependent variable is binary & transform into a dummy if required
qui levelsof `y' `if', local(yval)
local yvalcount: word count `yval'
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

// check if x variable is binary & transform into a dummy if required
qui levelsof `x' `if', local(xval)
local xvalcount: word count `xval'
local binaryx = 0
if "`xval'"=="0 1" local binaryx = 1
capture confirm numeric variable `x'
if `xvalcount'==2 & ("`xval'"!="0 1" | _rc) {
	local binaryx = 1
	local xlab: variable label `x'
	rename `x' old`x'
	egen `x' = group(old`x')
	replace `x' = `x'-1
	lab var `x' "`xlab'"
}

// retrieve names and labels from variables
foreach labb in x z y {
	if !(`binarydv'==1 & "`labb'"=="y") | !(`binaryx'==1 & "`labb'"=="x") {
		local `labb'lab: variable label ``labb''
		if "``labb'lab'"=="" local `labb'lab ``labb''
	}
}

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
		local ylab `y'
	}
	else {
		local ylab Pr(`ylab')
	}
}

if `binaryx'==1 {
	cap confirm variable old`x'
	if _rc {
		local xlabber `x'
	}
	else {
		local xlabber old`x'
	}
	levelsof `xlabber', local(oldxvals)
	foreach kk in `oldxvals' {
		local xlab: label (`xlabber') `kk'
	}
	gen xtitletest = "`xlab'"
	destring xtitletest,replace
	capture confirm numeric variable xtitletest
	if !_rc { 
		local xlab `x'
	}
	else {
		local xlab `xlab'
	}
}

local xtitle xtitle("`zlab'")
local ytitle ytitle("Effect of `xlab'" "on `ylab'")

// standardize x and y
if "`standardize'"!="" {
	if `binarydv'==0 gstats transform (standardize) `y' `x' `w', replace
	if `binarydv'==1 gstats transform (standardize) `x' `w', replace
}


*-------------------------------------------------------------------------------
* generate / specify bin variable
*-------------------------------------------------------------------------------

if "`method'"=="" & "`binvar'"=="" & "`nunibin'"=="" local method quantiles
if "`method'"=="" & "`binvar'"=="" & "`nunibin'"!="" local method unibin

if "`method'"=="quantiles" {
    if "`nquantiles'"=="" local nquantiles 10
	gquantiles `z'_cat = `z' `w', xtile nq(`nquantiles')  
}

else if "`method'"=="unibin" {
    if "`nunibin'"=="" local nunibin 10
	local `nunibin' = `unibin'+1
	gen `z'_cat = .
	local bincount = 0
	sum `z' `w'
	range binrange r(min) r(max) `nunibin'
	foreach bb of numlist 1/`nunibin' {
		local bincount = `bincount'+1
		qui sum binrange if _n==`bb'+1
		local binrange1 = r(mean)
		qui sum binrange if _n==`bb'
		replace `z'_cat = `bb' if `z'>=r(mean) & `z'<`binrange1'
	}
}

else if "`method'"=="discrete" {
    clonevar `z'_cat = `z' 
}

else if "`binvar'"!="" {
    clonevar `z'_cat = `binvar'
}

egen zbin = group(`z'_cat)
levelsof zbin, local(zlvl)
sum zbin 
local zbinlength = r(max)


*-------------------------------------------------------------------------------
* estimate and store regression
*-------------------------------------------------------------------------------

// prep setting 
if "`indslopesmethod'"=="" local indslopesmethod interact 
if "`level'"=="" local level 95
local cilvl level(`level')
local cilvltext `level'%

// factorize 
foreach vv in `fcontrols' {
	local fcontrols_f `fcontrols_f' i.`vv'
}

// prepare specification
local xzfactor c.`x'##i.zbin
local xzcont c.`x'##c.`z'
if strpos("`fit'","qfit") local xzcont `xzcont'##c.`z'

local hdfeabsorb noabsorb
if "`fcontrols'"!="" local hdfeabsorb absorb(`fcontrols')
if "`vce'"!="" local vce vce(`vce')
if "`binarymodel'"=="" local binarymodel logit

// estimate factorized
if "`indslopesmethod'"=="interact" {
	if `binarydv'==0 reghdfe `y' `xzfactor' `controls' `w', `hdfeabsorb' `vce' `cilvl'
	if `binarydv'==1 `binarymodel' `y' `xzfactor' `controls' `fcontrols_f' `w', `vce' `cilvl'
	est sto regm_fact
}
if "`indslopesmethod'"=="stratify" {
	foreach zlvl2 in `zlvl' {
		if `binarydv'==0 reghdfe `y' `xzfactor' `controls' if zbin==`zlvl2' `w', `hdfeabsorb' `vce' `cilvl'
		if `binarydv'==1 `binarymodel' `y' `xzfactor' `controls' if zbin==`zlvl2' `fcontrols_f' `w', `vce' `cilvl'
		est sto regm_fact`zlvl2'
	}
}

// estimate continuous
if `binarydv'==0 {
	reghdfe `y' `xzcont' `controls' `w', `hdfeabsorb' `vce' `cilvl'
	local coef = r(table)[1,3]
	local pval = r(table)[4,3]
}
if `binarydv'==1 `binarymodel' `y' `xzcont' `controls' `fcontrols_f' `w', `vce' `cilvl'
est sto regm_cont


*-------------------------------------------------------------------------------
* gather prediction & CIs from regression
*-------------------------------------------------------------------------------

// gather individual slope coefficients
gen `x'_sl = .
gen `x'_lci = .
gen `x'_uci = .
gen `x'_ci = .

// interact ind slopes method 
if "`indslopesmethod'"=="interact" {
	local count = 0
	est res regm_fact
	margins, dydx(`x') at(zbin=(`zlvl')) post `cilvl'
	foreach zlvl2 in `zlvl' {
		local count = `count'+1
		replace `x'_sl = r(table)[1,`count'] if zbin==`zlvl2'
		replace `x'_lci = r(table)[5,`count'] if zbin==`zlvl2'
		replace `x'_uci = r(table)[6,`count'] if zbin==`zlvl2'
		dis r(table)[6,`count']
	}
}

// stratify ind slopes method
if "`indslopesmethod'"=="stratify" {
	foreach zlvl2 in `zlvl' {
		est res regm_fact`zlvl2'
		margins, dydx(`x') post `cilvl'
			replace `x'_sl = r(table)[1,1] if zbin==`zlvl2'
			replace `x'_lci = r(table)[5,1] if zbin==`zlvl2'
			replace `x'_uci = r(table)[6,1] if zbin==`zlvl2'
	}
}

foreach zlvl2 in `zlvl' {
	sum `x'_lci if zbin==`zlvl2'
	bys zbin: replace `x'_ci = r(min) if zbin==`zlvl2' & _n==1
	sum `x'_uci if zbin==`zlvl2'
	bys zbin: replace `x'_ci = r(max) if zbin==`zlvl2' & _n==2
}

// prepare settings for margins command
local atmean atmeans
if !strpos("`fit'","qfit") 	local margpoints 30
if strpos("`fit'","qfit") 	local margpoints 50
local mcount = 0
sum `z' `w'
range `z'_at r(min) r(max) `margpoints'
levelsof `z'_at, local(margat)

// gather total slope 
gen `z'_sl = .
gen `z'_lci = .
gen `z'_uci = .
gen `z'_p = .

est res regm_cont
margins, dydx(`x') at(`z'=(`margat')) `atmean' `cilvl' post
est sto intmarg 
foreach num of numlist 1/`margpoints' {
    local mcount = `mcount'+1
	replace `z'_lci = r(table)[5,`num'] if _n==`num'
	replace `z'_uci = r(table)[6,`num'] if _n==`num'
	replace `z'_sl = r(table)[1,`num'] if _n==`num'
}

	
*-------------------------------------------------------------------------------
* gather regression parameters
*-------------------------------------------------------------------------------

if "`regparameters'"!="" {
	
	* nobs, r2, adjr2
	est res regm_cont
	local nobs = e(N)
	if `binarydv'==0 {
		local r2 = e(r2)
		local adjr2 = e(r2_a)
	}
	if `binarydv'==1 {
		local r2 = e(r2_p)
		local adjr2 = e(r2_p)
	}
	dis `nobs' `r2'
	*coef, pval, sig
	if `binarydv'==1 {
		sum `z',d
		local zp50 = r(p50)
		local zp502 = `zp50'+1
		est res regm_cont
		margins, dydx(`x') at(`z'=(`zp50' `zp502')) `cilvl' post
		local coef = _b[1._at]-_b[2._at]
		test _b[1._at] = _b[2._at]
		local pval = r(p)
	}

	local siglevel = 99
	if `pval'<.1 & `pval'>.05 local siglevel = .1
	if `pval'<.05 & `pval'>.01 local siglevel = .05
	if `pval'<.01 local siglevel = .01

	// round and store parameters	
	foreach par in nobs r2 adjr2 coef pval {
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
* finalize data points
*-------------------------------------------------------------------------------

gegen `z'_plot = mean(`z') `w', by(zbin)
rename `x'_sl `y'_plot 
if "`mweighted'"!="" gegen scw = count(`y'_plot), by(zbin)
egen tag = tag(zbin) if !mi(`y'_plot)
replace `y'_plot = . if tag!=1

count if !mi(`y'_plot) & !mi(`z'_plot) //& !mi(`by')
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
	
	foreach i of numlist 1/2 {
		local mlines`i' lc(`c`i'') lw(medthick)
		local mlinesci`i' acol(`c`i'o50') alw(none) clc(`c`i'') clw(medthick)
		local ciareas`i' lw(none) fc(`c`i'o30')
		local cicap`i' lw(medium) lc(`c`i'o50') lp(solid)
		local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') mlabc(`c`i'') mlw(thin)
		local efullscatterm`i' `mfullscatterm`i''
		local mscatter75m`i' mfc(`c`i'o50') mlc(`c`i'o75') mlabc(`c`i'o80') mlw(vthin) mlalign(inside)
		local escatter75m`i' `mscatter75m`i''
		foreach h of numlist 25 50  {
			local mscatter`h'm`i' mfc(`c`i'o`h'') mlabc(`c`i'o`h'') mlw(none) mlalign(outside)
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
	local hollowresize = .875

	foreach sizeloc in osize tsize ssize dsize {
	    local e`sizeloc' = ``sizeloc''*`globmsize'
		local m`sizeloc' = ``sizeloc''*`globmsize'`mresize'`mweightedresize'`mbyresize'
	    local he`sizeloc' = `e`sizeloc''*`hollowresize'
		local hm`sizeloc' = `m`sizeloc''*`hollowresize'
	}

	foreach en in e m {
		local m1 m(d) msize(*``en'dsize')
		local m2 m(o) msize(*``en'osize')

		foreach g in `en'fullscatterm `en'scatter25m `en'scatter50m `en'scatter75m { 
			local `g'1 ``g'1' `m1'
			local `g'2 ``g'2' `m2'
			local `g'3 ``g'3' m(t) msize(*``en'tsize')
			local `g'4 ``g'4' m(s) msize(*``en'ssize')
		}
	}

	foreach g in mlines mlinesci {
		local `g'1 ``g'1' lp(solid)
		local `g'2 ``g'2' lp(dash)
	}
}

else {
	local scheme scheme(`plotscheme')
	if "`colorscheme'"!="" {
		foreach i of numlist 1/2 {
			local mlines`i' lc(`c`i'') lw(medthick)
			local mlinesci`i' acol(`c`i'o50') alw(none) clc(`c`i'') clw(medthick)
			local ciareas`i' lw(none) fc(`c`i'o30')
			local cicap`i' lw(medthin) lc(`c`i'o50') lp(solid)
			local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') mlabc(`c`i'') mlw(thin)
			local efullscatterm`i' `mfullscatterm`i''
			local mscatter75m`i' mfc(`c`i'o50') mlc(`c`i'o75') mlabc(`c`i'o80') mlw(vthin) mlalign(inside)
			local escatter75m`i' `mscatter75m`i''
			foreach h of numlist 25 50  {
				local mscatter`h'm`i' mfc(`c`i'o`h'') mlabc(`c`i'o`h'') mlw(none) mlalign(outside)
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
	sum `y'_plot,d
	local ymax = r(max)
	local ymin = r(min)
	local y25 = r(p25)
	local y75 = r(p75)
	count if `y'_plot>`y75'
	local y75larger = r(N)
	count if `y'_plot<`y25'
	local y25smaller = r(N)
	sum `z'_plot,d
	local xmax = r(max)
	local xmin = r(min)
	local xmean = r(mean)
	local x75 = r(p75)
	count if `z'_plot>`x75'
	local x75larger = r(N)
	count if `z'_plot<`x75'
	local x75smaller = r(N)

	count if `y'_plot>`y75' & `z'_plot<`y75'
	if (r(N)/`n')<.03 {
		local textplacey `ymax'
		local textplacex `xmean'
	}
	else {
		count if `y'_plot>`y75' & `z'_plot>`y75'
		if (r(N)/`n')<.03 {
			local textplacey `ymax'
			local textplacex `xmax'
		}
		else {
			count if `y'_plot<`y25' & `z'_plot>`y75'
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
	
	// compile the parameters 
	if strpos("`regparameters'","r2") & !strpos("`regparameters'","adjr2") {
		if `binarydv'==0 local r2par {it:R2}`r2'
		if `binarydv'==1 local r2par {it:Pseudo R2}`r2'
	}
	* adj r2
	if strpos("`regparameters'","adjr2") {
		if `binarydv'==0 local adjr2par `adjr2par' {it:Adj. R2}`adjr2'
		if `binarydv'==1 local adjr2par `adjr2par' {it:Pseudo R2}`r2'
	}
	* observations 
	if strpos("`regparameters'","nobs") local nobspar {it:N}`nobs'
	* coef and sig
	if strpos("`regparameters'","coef") {
		if strpos("`regparameters'","sig") {
			if `siglevel'==.1 local sigstar *
			if `siglevel'==.05 local sigstar **
			if `siglevel'==.01 local sigstar ***
		}
		if `binarydv'==0 local whatcoef {it:{&delta}{&beta}/{&delta}x}
		if `binarydv'==1 local whatcoef {it:{&delta}{&beta}/{&delta}x}
		local coefpar `whatcoef'`coef'`sigstar'
	}
	* p value
	if strpos("`regparameters'","pval") local pvalpar {it:p}`pval'

	// compile 
	foreach ppp in `regparameters' {
		if "`ppp'"!="sig" & strpos("`regparameters'","`ppp'") local rparcol `" `rparcol' "``ppp'par'" "'
	}

	// final text box
	local printcoef2 text(`textplace' `rparcol', ///
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
	*if `wherecoef'==5 local leginsidepl 1
    local legtype ring(0) region(lc(gs5) fc(white)) pos(`leginsidepl')
}

local legopts legend(`legtype' `legresize')

// compile labels and legend options
if strpos("`fit'","lfit") | "`fit'"=="" local leg_fit Linear
if strpos("`fit'","qfit") local leg_fit Quadratic


local addmarkercis = 0
local addemptymark = 0
local add95cis = 0
if "`indslopesci'"!="" local addmarkercis = `zbinlength'
if "`mweighted'"!="" local addemptymark = 1
if strpos("`fit'","ci") local add95cis = 1
foreach i of numlist 1/3 {
	local n`i' = `i'+`addmarkercis'
	local nn`i' = `i'+`addmarkercis'+`addemptymark'
	local nnn`i' = `i'+`addmarkercis'+`addemptymark'+`add95cis'
}
if "`mlabel'"=="" local legindslopes `n1' "Individual slopes"
if "`indslopesci'"!="" local legindsl 1 "95% CIs"
if strpos("`fit'","ci") local legci95 `nn2' "95% CIs"

local legopts `legopts' legend(order(`legindslopes' `legindsl' `nnn2' "`leg_fit' slope fit" `legci95'))
dis `""`legopts'""'

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
if "`bwidth'"!="" local bwidth bw(`bwidth')
if "`jitter'"!="" local jitter jitter(`jitter')
if "`noyline'"=="" local yline yline(0, lc(gs5) lw(thin) lp(solid))
local lscatteropts `scheme' `yline' `xtitle' `ytitle' `printcoef2' `legopts' `plotsize' `opts' 

// empty scatter marker plot for correct legend in case of weighted scatter markers
if "`mweighted'"!="" {
	local sce `sce' (scatter `y'_plot `z'_plot if mi(`y'_plot), ``escattermarkers'2')
}

// cis for individual slopes
*local isl `isl' (line `x'_lci `x'_uci , `cicap2')
foreach zlvl2 in `zlvl' {
	if "`indslopesci'"!="" local isl `isl' (line `x'_ci `z'_plot if zbin==`zlvl2', `cicap2')
}

// individual slopes
local sc `sc' (scatter `y'_plot `z'_plot `scw', ``mscattermarkers'2' `markerlab' `jitter')


// fitted slope
if strpos("`fit'","ci") local ci `ci' (rarea `z'_lci `z'_uci `z'_at, `ciareas1')
local pl `pl' (line `z'_sl `z'_at, `o1')

// draw the final plot
tw `isl' `sce' `sc' `ci' `pl', `lscatteropts'
	

	
	
*-------------------------------------------------------------------------------
* finalize
*-------------------------------------------------------------------------------

restore
}
end













