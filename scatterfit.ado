*! version 1.2 Leo Ahrens

program define scatterfit
	syntax varlist(min=2 max=2 numeric), [fit(string) by(string) BINned discrete NQuantiles(numlist) polybw(numlist) COVariates(string) ///
	ABSorb(string) STANDardize coef vce(string) coefplace(string) jitter(numlist) COLorscheme(string) PLOTScheme(string) opts(string asis)]

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
	if "`plotscheme'"!="" capture set scheme plotplain
	if _rc==111 ssc install blindschemes, replace

	// check if everything is specified correctly
	if "`coef'"!="" & "`by'"!="" {
		di as error "The coefficient can only be plotted if the by() option is not specified."
		exit 498
	}
	if "`polybw'"!="" & !strpos("`fit'","poly") {
		di as error "The polybw() options requires a local polynomial fit via fit(poly) or fit(polyci)."
		exit 498
	}
	if "`nquantiles'"!="" & ("`discrete'"!="" | "`binned'"=="") {
		di as error "The nquantiles() option requires the option binned, and it cannot be combined with the option discrete."
		exit 498
	}
	if ("`coefplace'"!="" | "`vce'"!="") & "`coef'"=="" {
		di as error "The coefplace() and vce() options require the coef option."
		exit 498
	}
	if "`discrete'"!="" & "`binned'"=="" {
		di as error "The discrete option requires the binned option."
		exit 498
	}
	if ("`covariates'"!="" | "`absorb'"!="") & strpos("`fit'","ci") {
		di as error "Please note that the displayed confidence intervals are incorrect in the presence of covariates due to not taking the controls into account. This will be fixed in a future version"
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
		rename `by' oldby
		egen `by' = group(oldby)
		labmask `by', values(oldby)
	}
	levelsof `by', local(bynum)

	// clean dataset
	if "`covariates'"!="" {
		foreach v of varlist `covariates' `absorb' {
			local covdrop `covdrop' | mi(`v')
			local covkeep `covkeep' `v'
		}
	}
	drop if mi(`x') | mi(`y') | mi(`by') `covdrop'
	keep `x' `y' `by' `covkeep'

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
		if "`discrete'"=="" {
			if "`nquantiles'"=="" local nquantiles 30
			gquantiles `x'_q = `x', xtile nq(`nquantiles') `byparen'
			local xbin `x'_q `by'
		}
		if "`discrete'"!="" {
		    local xbin `x' `by'
		}
	}

	// standardize 
	if "`standardize'"!="" gstats transform (standardize) `y' `x', replace

	// regression adjustment
	if "`covariates'"!="" | "`absorb'"!="" {
		local hdfeabsorb noabsorb
		if "`absorb'"!="" local hdfeabsorb absorb(`absorb')

		if "`binned'"=="" {
			foreach v in `y' `x' {
				reghdfe `v' `covariates', `hdfeabsorb' res(`v'_r)
				sum `v'
				replace `v'_r = `v'_r + r(mean)
				local `v' `v'_r
			}
		}
		
		if "`binned'"!="" {
			gen `y'_r = .
			foreach bynum2 in `bynum' {
				reghdfe `y' `covariates' i.`xbin' if `by'==`bynum2', `hdfeabsorb'
				predict `y'_r`bynum2' if e(sample), xb
				if "`covariates'"!="" {
					foreach v of varlist `covariates' {
						replace `y'_r`bynum2' = `y'_r`bynum2' - _b[`v']*`v' if `by'==`bynum2'
					}
				}
				replace `y'_r = `y'_r`bynum2' if !mi(`y'_r`bynum2')
			}
			local y `y'_r
		}
	}

	// retrieve linear fit coefficient
	if (strpos("`fit'","lfit") | "`fit'"=="") & "`coef'"!="" & `isthereby'==0 {
		if "`vce'"!="" local vce vce(`vce')
		local hdfeabsorb noabsorb 
		if "`absorb'"!="" local hdfeabsorb absorb(`absorb')
		reghdfe `y' `x' `covariates', `hdfeabsorb' `vce'
		local beta = round(_b[`x'], .01)
		local pval = round(2*normal(-abs(_b[`x']/_se[`x'])),.01)
	}
	
	// mean within bins
	if "`binned'"!="" {
		gegen `y'_mean = mean(`y'), by(`xbin')
		gegen `x'_mean = mean(`x'), by(`xbin')
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
		foreach i of numlist 25 50 75 {
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
			local fullscatterm`i' m(o) mfc(`c`i'o50') mlc(`c`i'') mlw(thin)
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
	if "`fit'"=="lfit" | "`fit'"=="" local fittype lfit
	if "`fit'"=="lfitci" local fittype lfitci
	if "`fit'"=="qfit" local fittype qfit
	if "`fit'"=="qfitci" local fittype qfitci
	if "`fit'"=="poly" local fittype lpoly
	if "`fit'"=="polyci" local fittype lpolyci
	
	if "`fit'"=="" | "`fit'"=="lfit" | "`fit'"=="qfit" | "`fit'"=="poly" {
		foreach ff of numlist 1/9 {
			local o`ff' `mlines`ff''
		}
	}
	if "`fit'"=="lfitci" | "`fit'"=="qfitci" | "`fit'"=="polyci" {
		foreach ff of numlist 1/9 {
			local o`ff' `mlinesci`ff''
		}
	}

	// coefficient text options 
	if (strpos("`fit'","lfit") | "`fit'"=="") & "`coef'"!="" & `isthereby'==0 {
		sum `yplot'
		local textplacey = r(max)+.1*r(max)
		sum `xplot'
		local textplacex = r(mean)
		local textplace2 `textplacey' `textplacex'
		if "`coefplace'"!="" local textplace2 `coefplace'
		local printcoef2 text(`textplace2' "{it:ß} = `beta'" "{it:p} = `pval'", placement(center) size(*.85) box fc(white%50) lc(gs5) lw(thin) la(outside) margin(vsmall))
	}
	
	// legend options
	local legopts_all legend(region(lc(white)) pos(3) size(*1.05))
	if strpos("`fit'","lfit") | "`fit'"=="" local leg_fit Linear
	if strpos("`fit'","qfit") local leg_fit Quadratic
	if strpos("`fit'","poly") local leg_fit Local polynomial
	if `isthereby'==0 {
		if !strpos("`fit'","ci") local legopts `legopts_all' legend(order(1 "Observed" 2 "`leg_fit' fit"))
		if strpos("`fit'","ci") local legopts `legopts_all' legend(order(1 "Observed" 3 "`leg_fit' fit" 2 "95% CIs"))
	}
	else {
		egen distinctby = group(`by')
		sum distinctby
		local maxdistinctby = r(max)
		local coln = 0
		levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			local coln = `coln'+1
			if !strpos("`fit'","ci") local coln2 = `coln'+`maxdistinctby'
			if strpos("`fit'","ci") local coln2 = `coln'+`maxdistinctby'+`coln'
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
	
	// prepare standard settings
	if "`polybw'"!="" local polybw2 bw(`polybw')
	if "`jitter'"!="" local jitter jitter(`jitter')
	
	// generate the plot
	local lscatteropts `plotscheme' `xtitle' `ytitle' `printcoef2' `opts' `legopts'
	local coln = 0
	foreach bynum2 in `bynum' {
		local coln = `coln'+1
		local coln2 = `coln'
		if `isthereby'==0 local coln2 = `coln'+1
		local sc `sc' (scatter `yplot' `xplot' if `by'==`bynum2', ``scattermarkers'`coln2'' `jitter') 
		local pl `pl' (`fittype' `y' `x' if `by'==`bynum2', `o`coln'' `polybw2')
	}
	tw `sc' `pl', `lscatteropts'

				
	restore
	}
end




