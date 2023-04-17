*! version 1.3   Leo Ahrens   leo@ahrensmail.de

program define scatterfit
	version 13.1
	
	#delimit ;
	syntax varlist(min=2 max=2 numeric)	[if] [in] [aweight fweight] , 
	[
	fit(string) bw(numlist) polybw(numlist) 
	by(string) 
	BINned DISCrete NQuantiles(numlist) BINVar(varlist)
	Controls(varlist) Fcontrols(varlist) COVariates(varlist) ABSorb(varlist)
	coef COEFPlace(string)
	vce(string)
	JITter(numlist)
	STANDardize 
	opts(string asis)
	plotscheme(string asis)
	colorscheme(string asis)
	] ;
	#delimit cr

	// install dependencies
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

	// harmonize legacy options
	if "`polybw'"!="" & "`bw'"=="" local bw `polybw'
	if "`absorb'"!="" & "`fcontrols'"=="" local fcontrols `absorb'
	if "`covariates'"!="" & "`controls'"=="" local controls `covariates'
	
	// check if everything is specified correctly
	if "`coef'"!="" & "`by'"!="" {
		di as error "The coefficient can only be plotted if the by() option is not specified."
		exit 498
	}
	if ("`bw'"!="" | "`polybw'"!="") & !(strpos("`fit'","poly") | "`fit'"=="lowess") {
		di as error "The bw() or polybw() option requires a local polynomial fit or a lowess smother as the fit line."
		exit 498
	}
	if "`nquantiles'"!="" & ("`discrete'"!="" | "`binned'"=="") {
		di as error "The nquantiles() option requires the option binned, and it cannot be combined with the option discrete."
		exit 498
	}
	if "`coefplace'"!="" & "`coef'"=="" {
		di as error "The coefplace() option requires the coef option."
		exit 498
	}
	if ("`discrete'"!="" | "`binvar'"!="") & "`binned'"=="" {
		di as error "The discrete and binvar() options require the binned option."
		exit 498
	}
	if !("`fit'"=="" | "`fit'"=="lfit" | "`fit'"=="lfitci" | "`fit'"=="qfit" | "`fit'"=="qfitci" | "`fit'"=="poly" | "`fit'"=="polyci" | "`fit'"=="lowess") {
		di as error "fit() must be lfit, lfitci, qfit, qfitci, poly, or polyci."
		exit 498
	}
	if "`vce'"!="" & (!strpos("`fit'","ci") | strpos("`fit'","poly")) {
		di as error "The vce() option requires that confidence intervals are drawn by fit(lfitci) or fit(qfitci). Not available for fit(polyci)."
		exit 498
	}
	if "`controls'`fcontrols'"!="" & strpos("`fit'","ci") & ("`by'"!="" | strpos("`fit'","poly")) {
		di as error "No confidence intervals are drawn by the current version of scatterfit when (a) covariates are specified and (b) any of the following are specified: a polynomial fit, a quadratic fit, or by()."
	}

	quietly {
	preserve

	// declare x and y variables
	tokenize `varlist'
	local y `1'
	local x `2'
		
	// prepare by variable
	local isthereby = 0
	if "`by'"!="" local isthereby = 1
	if "`by'"!="" local byparen by(`by')
	if `isthereby'==0 {
		gen sfitbyvar = 1
		local by sfitbyvar
	}
	capture confirm numeric variable `by'
	if _rc {
		capture which labmask
		if _rc==111 ssc install labutil, replace
		rename `by' oldby
		egen `by' = group(oldby)
		labmask `by', values(oldby)
	}
	levelsof `by', local(bynum)

	// weight local
	if ("`weight'"!="") local w [`weight'`exp']
	if ("`weight'"!="") local weightname = subinstr("`exp'","=","",.)
	
	// clean dataset
	if "`controls'`fcontrols'"!="" | "`binvar'"!="" {
		foreach v of varlist `controls' `fcontrols' `binvar' {
			local covdrop `covdrop' | mi(`v')
		}
	}
	drop if mi(`x') | mi(`y') | mi(`by') `covdrop'
	if "`if'"!="" keep `if'
	if "`in'"!="" keep `in'
	keep `x' `y' `by' `controls' `fcontrols' `binvar' `weightname'
	
	// retrieve labels from variables
	local xlab: variable label `x'
	local xtitle xtitle("`xlab'")
	local ylab: variable label `y'
	local ytitle ytitle("`ylab'")
	if `isthereby'==1 {
		foreach bynum2 in `bynum' {
			local bylab`bynum2': label (`by') `bynum2'
		}
	}

	// generate bin variables
	if "`binned'"!="" {
		local xbin `x'_q `by'
		if "`nquantiles'"=="" local nquantiles 30
		if "`binvar'"=="" {
			if "`discrete'"=="" gquantiles `x'_q = `x' `w', xtile nq(`nquantiles') `byparen'
			if "`discrete'"!="" clonevar `x'_q = `x'
		}
		if "`binvar'"!="" {
			clonevar `x'_q = `binvar'
		}
	}

	// standardize 
	if "`standardize'"!="" gstats transform (standardize) `y' `x' `w', replace

	// estimate model for CIs or gathering parameters
	if !strpos("`fit'","poly") & "`fit'"!="lowess" & `isthereby'==0 & (("`coef'"!="" & !strpos("`fit'","qfit")) | (strpos("`fit'","ci") & ("`controls'`fcontrols'"!="" | "`vce'"!=""))) {
		local hdfeabsorb noabsorb
		if "`fcontrols'"!="" local hdfeabsorb absorb(`fcontrols')
		if "`vce'"!="" & strpos("`fit'","ci") local vce vce(`vce')
		local `x'marg `x'
		if strpos("`fit'","qfit") local `x'marg c.`x'##c.`x'
		reghdfe `y' ``x'marg' `controls' `w', `hdfeabsorb' `vce'
		est sto regmodel
	}
		
	// retrieve linear fit coefficient
	if "`coef'"!="" & !strpos("`fit'","qfit") & !strpos("`fit'","poly") & "`fit'"!="lowess" & `isthereby'==0 {
		est res regmodel
		local beta = _b[`x']
		local pval = 2*normal(-abs(_b[`x']/_se[`x']))
		foreach par in beta pval {
			local smallround`par' = 0
			if ``par''>=10 {
				local `par'round "1"
			}
			else {
				if ``par''>=1 {
					local `par'round ".1"
				}
				else {
					local roundcount = 0
					local `par'string "``par''"
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
			}
			cap if strpos("``par'string'","e") local smallround`par' = 1
			local `par' = round(``par'',``par'round')
			if `smallround`par''==0 {
				local `par' "= ``par''"
			}
			else {
				local `par' "< .00001"
			}
		}
	}

	// gather linear prediction & CIs in presence of controls
	if (strpos("`fit'","lfit") | "`fit'"=="" | strpos("`fit'","qfit")) & strpos("`fit'","ci") & ("`controls'`fcontrols'"!="" | "`vce'"!="") & `isthereby'==0 {
		if "`binned'"!="" gen `x'2 = `x'
		if "`binned'"=="" {
			local hdfeabsorb noabsorb
			if "`fcontrols'"!="" local hdfeabsorb absorb(`fcontrols')
			reghdfe `x' `controls' `w', `hdfeabsorb' res(`x'2)
			sum `x' `w'
			replace `x'2 = `x'2 + r(mean)
		}
		local margpoints 30
		gen cix = .
		gen ciu = .
		gen cil = .
		gen pe = .
		local mcount = 0
		sum `x'2 `w'
		range range r(min) r(max) `margpoints'
		foreach p of numlist 1/`margpoints' {
			local mcount = `mcount'+1
			sum range if _n==`p'
			local ranger`mcount' = r(mean)
			local margat `margat' `ranger`mcount''
		}
		est res regmodel
		margins, at(`x'=(`margat')) atmeans post
		foreach num of numlist 1/`mcount' {
			replace cix = `ranger`num'' if _n==`num'
			replace cil = r(table)[5,`num'] if _n==`num'
			replace ciu = r(table)[6,`num'] if _n==`num'
			replace pe = r(table)[1,`num'] if _n==`num'
		}
	}

	// covariate adjustment of scatter points
	if "`controls'`fcontrols'"!="" {
		local hdfeabsorb noabsorb
		if "`fcontrols'"!="" local hdfeabsorb absorb(`fcontrols')
		
		if "`binned'"=="" {
			foreach v in `y' `x' {
				gen `v'_r = .
				foreach bynum2 in `bynum' {
					reghdfe `v' `controls' if `by'==`bynum2' `w', `hdfeabsorb' res(`v'_r`bynum2')
					sum `v' if `by'==`bynum2' `w'
					replace `v'_r = `v'_r`bynum2' + r(mean) if `by'==`bynum2'
				}
			}
			local y `y'_r 
			local x `x'_r
		}
		
		if "`binned'"!="" {
			gen `y'_r = .
			foreach bynum2 in `bynum' {
				reghdfe `y' `controls' i.`xbin' if `by'==`bynum2' `w', `hdfeabsorb'
				predict `y'_r`bynum2' if e(sample), xb
				if "`controls'"!="" {
					foreach v of varlist `controls' {
						replace `y'_r`bynum2' = `y'_r`bynum2' - _b[`v']*`v' if `by'==`bynum2'
					}
				}
				replace `y'_r = `y'_r`bynum2' if !mi(`y'_r`bynum2')
			}
			sum `y' `w'
			local adjm = r(mean)
			sum `y'_r `w'
			replace `y'_r = `y'_r + (`adjm'-r(mean))
			local y `y'_r
		}
	}
	
	// mean within bins
	if "`binned'"!="" {
		gegen `y'_mean = mean(`y') `w', by(`xbin')
		gegen `x'_mean = mean(`x') `w', by(`xbin')
		egen tag = tag(`xbin') if !mi(`y'_mean)
		replace `y'_mean = . if tag!=1
	}
	
	// specify variable to be plotted depending on binned / non-binned
	if "`binned'"=="" {
		local yplot `y'
		local xplot `x'
	}
	else {
		local yplot `y'_mean
		local xplot `x'_mean
	}

	// overall scheme options
	if "`colorscheme'"=="" {
		local colorscheme Set1
		local colorschemeopts int(1.8)
	}
	if !("`colorscheme'"=="" & "`plotscheme'"!="") {
		colorpalette `colorscheme', `colorschemeopts' nograph local(,prefix(c) nonames)
		foreach i of numlist 25 30 40 50 75 {
			colorpalette `colorscheme', `colorschemeopts' op(`i') nograph local(,prefix(c) suffix(o`i') nonames) 
		}
	}
	if "`plotscheme'"=="" {
		local plotscheme scheme(plotplain) graphregion(lc(white) lw(vthick)) title(,size(medium)) ///
		ysc(lc(black) lw(thin)) ylab(, labs(*1.05) tlc(black) tlw(thin) glc(gs13) glp(solid) glw(thin) gmin gmax) ///
		xsc(lc(black) lw(thin)) xlab(, labs(*1.05) tlc(black) tlw(thin) glc(gs13) glp(solid) glw(thin) gmin gmax)
		
		foreach i of numlist 1/8 {
			local mlines`i' lc(`c`i'') lw(medthick)
			local mlinesci`i' acol(`c`i'o50') alw(none) clc(`c`i'') clw(medthick)
			local ciareas`i' lw(none) fc(`c`i'o30')
			local fullscatterm`i' m(o) mfc(`c`i'o40') mlc(`c`i'') mlw(thin)
			local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') mlw(thin)
			foreach h of numlist 25 50 75 {
				local scatter`h'm`i' m(o) mfc(`c`i'o`h'') mlalign(outside) mlw(none)
				local mscatter`h'm`i' mfc(`c`i'o`h'') mlalign(outside) mlw(none)
			}
		}
		foreach g in mfullscatterm mscatter25m mscatter50m mscatter75m { 
			local `g'1 ``g'1' m(s) 
			local `g'2 ``g'2' m(o) 
			local `g'3 ``g'3' m(d)
			local `g'4 ``g'4' m(t) 
			local `g'5 ``g'5' m(o) 
		}
		foreach g in mlines mlinesci {
			local `g'1 ``g'1' lp(solid)
			local `g'2 ``g'2' lp(dash)
			local `g'3 ``g'3' lp(shortdash)
			local `g'4 ``g'4' lp(shortdash dot)
			local `g'5 ``g'5' lp(dot)
			foreach i of numlist 6/9 {
				local `g'`i' ``g'`i'' lp(solid)
			}
		}
	}
	else {
		local plotscheme scheme(`plotscheme')
		if "`colorscheme'"!="" {
			foreach i of numlist 1/8 {
				local mlines`i' lc(`c`i'')
				local mlinesci`i' acol(`c`i'o50') clc(`c`i'')
				local fullscatterm`i' mfc(`c`i'o50') mlc(`c`i'') 
				local mfullscatterm`i' mfc(`c`i'o50') mlc(`c`i'')
				foreach h of numlist 25 50 75 {
					local scatter`h'm`i' mfc(`c`i'o`h'')
					local mscatter`h'm`i' mfc(`c`i'o`h'')
				}
			}
		}
	}

	// load fit type options
	if "`fit'"=="lfitci" local fittype lfitci
	if "`fit'"=="lfit" | "`fit'"=="" | ("`fit'"=="lfitci" & ("`controls'`fcontrols'"!="" | "`vce'"!="")) local fittype lfit
	if "`fit'"=="qfitci" local fittype qfitci
	if "`fit'"=="qfit" | ("`fit'"=="qfitci" & ("`controls'`fcontrols'"!="" | "`vce'"!="")) local fittype qfit
	if "`fit'"=="polyci" local fittype lpolyci
	if "`fit'"=="poly" | ("`fit'"=="polyci" & "`controls'`fcontrols'"!="") local fittype lpoly
	if "`fit'"=="lowess" local fittype lowess
	
	if "`fit'"=="" | "`fit'"=="lfit" | "`fit'"=="qfit" | "`fit'"=="poly" | "`fit'"=="lowess" {
		foreach ff of numlist 1/9 {
			local o`ff' `mlines`ff''
		}
	}
	if "`fit'"=="lfitci" | "`fit'"=="qfitci" | "`fit'"=="polyci" {
		foreach ff of numlist 1/9 {
			local o`ff' `mlinesci`ff''
		}
	}
	
	// legend options
	local legopts_all legend(region(lc(white)) pos(3) size(*1.05))
	if strpos("`fit'","lfit") | "`fit'"=="" local leg_fit Linear
	if strpos("`fit'","qfit") local leg_fit Quadratic
	if strpos("`fit'","poly") local leg_fit Local polynomial
	if strpos("`fit'","lowess") local leg_fit Lowess
	if `isthereby'==0 {
		if !strpos("`fit'","ci") | (strpos("`fit'","ci") & "`controls'`fcontrols'"!="" & !strpos("`fit'","lfit")) {
			local legopts `legopts_all' legend(order(1 "Observed" 2 "`leg_fit' fit"))
		}
		else {
			local legopts `legopts_all' legend(order(1 "Observed" 3 "`leg_fit' fit" 2 "95% CIs"))
		}
	}
	else {
		egen distinctby = group(`by')
		sum distinctby
		local maxdistinctby = r(max)
		local coln = 0
		levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			local coln = `coln'+1
			if !strpos("`fit'","ci") | (strpos("`fit'","ci") & "`controls'`fcontrols'"!="") {
				local coln2 = `coln'+`maxdistinctby'
			}
			else {
				local coln2 = `coln'+`maxdistinctby'+`coln'
			}
			local legorder `legorder' `coln' "`bylab`bynum2''" `coln2' ""
		}
		local legopts `legopts_all' legend(order(`legorder') col(2) textfirst)
	}
	
	// scatter marker options
	count if !mi(`xplot') & !mi(`yplot') & !mi(`by')
	local n = r(N)
	if `n'<=100 local scattermarkers mfullscatterm
	if `n'>100 & `n'<=300 local scattermarkers mscatter75m
	if `n'>300 & `n'<=2000 local scattermarkers mscatter50m
	if `n'>2000 local scattermarkers mscatter25m
	
	// coefficient text options 
	if (strpos("`fit'","lfit") | "`fit'"=="") & "`coef'"!="" & `isthereby'==0 {
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
				}
				else {
					local textplacey `ymax'
					local textplacex `xmean'
				}
			}
		}
		local textplace `textplacey' `textplacex'
		if "`coefplace'"!="" local textplace `coefplace'
		local printcoef2 text(`textplace' "{it:ß} `beta'" "{it:p} `pval'", placement(center) size(*.85) box fc(white%50) lc(gs5) lw(thin) la(outside) margin(vsmall))
	}

	// generate the plot
	local lscatteropts `plotscheme' `xtitle' `ytitle' `printcoef2' `opts' `legopts'
	if "`bw'"!="" local bw bw(`bw')
	if "`jitter'"!="" local jitter jitter(`jitter')
	local coln = 0
	foreach bynum2 in `bynum' {
		local coln = `coln'+1
		local coln2 = `coln'
		if `isthereby'==0 local coln2 = `coln'+1
		local sc `sc' (scatter `yplot' `xplot' if `by'==`bynum2', ``scattermarkers'`coln2'' `jitter')
		if ("`fit'"=="lfitci" | "`fit'"=="qfitci") & ("`controls'`fcontrols'"!="" | "`vce'"!="") & `isthereby'==0 {
			local ci `ci' (rarea cil ciu cix, `ciareas`coln'')
			local pl (line pe cix, `o`coln'')
		}
		else {
			local pl `pl' (`fittype' `y' `x' if `by'==`bynum2', `o`coln'' `bw')
		}
	}
	
	// draw the final plot
	tw `sc' `ci' `pl', `lscatteropts'
		
	restore
	}
end



