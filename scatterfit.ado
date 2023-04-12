*! version 1.0 Leo Ahrens

program scatterfit
	syntax varlist(min=2 max=2 numeric), [fit(string) by(string) BINned discrete NQuantiles(numlist) polybw(numlist) COVariates(string) ///
	ABSorb(string) coef vce(string) coefplace(string) jitter(numlist) COLorscheme(string) PLOTScheme(string) opts(string asis)]

	// install dependencies
	foreach package in reghdfe gtools {
		capture which `package'
		if _rc==111 ssc install `package', replace
		if _rc==111 & "`package'"=="gtools" gtools, upgrade
	}
	capture which colorpalette
	if _rc==111 ssc install palettes, replace
	if _rc==111 ssc install colrspace, replace
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
		di as error "The nquantiles() option, requires the option binned, and it cannot be combined with the option discrete."
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
	
	// prepare dataset
	preserve
	qui drop if mi(`x') | mi(`y')
	if "`by'"!="" qui drop if mi(`by')
	
	// declare x and y variables
	tokenize `varlist'
	local y `1'
	local x `2'
	
	// retrieve labels from variables
	local xlab: variable label `x'
	local xtitle xtitle("`xlab'")
	local ylab: variable label `y'
	local ytitle ytitle("`ylab'")
	if "`by'"!="" {
		qui levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			local bylab`bynum2': label (`by') `bynum2'
		}
	}

	// regression adjustment
	if "`covariates'"!="" | "`absorb'"!="" {
	    local hdfeabsorb noabsorb 
		if "`absorb'"!="" local hdfeabsorb absorb(`absorb')
		qui gstats transform (standardize) `y' `x', replace
		qui reghdfe `y' `covariates', `hdfeabsorb' res(`y'_r)
		qui reghdfe `x' `covariates', `hdfeabsorb' res(`x'_r)
		local y `y'_r
		local x `x'_r
	}
	
	// retrieve linear fit coefficient
	if (strpos("`fit'","lfit") | "`fit'"=="") & "`coef'"!="" & "`by'"=="" {
		if "`vce'"!="" local vce vce(`vce')
		local hdfeabsorb noabsorb 
		if "`absorb'"!="" local hdfeabsorb absorb(`absorb')
		qui reghdfe `y' `x' `covariates', `hdfeabsorb' `vce'
		local beta = round(_b[`x'], .01)
		local pval = round(2*normal(-abs(_b[`x']/_se[`x'])),.01)
	}
	
	// prepare (standard) settings
	if "`nquantiles'"=="" local nquantiles 40
	if "`polybw'"!="" local polybw2 bw(`polybw')
	if "`jitter'"!="" local jitter jitter(`jitter')
	
	// adapt var locals to by / binned settings
	if "`binned'"=="" {
		local yplot `y'
		local xplot `x'
	}
	else {
		local yplot `y'_mean
		local xplot `x'_mean
		local xvarname `x'_q
		if "`discrete'"!="" local xvarname `x'
	}
	
	// generate binned variables
	if "`binned'"!="" {
		if "`by'"=="" {
			qui if "`discrete'"=="" gquantiles `x'_q = `x', xtile nq(`nquantiles')
			qui gegen `y'_mean = mean(`y'), by(`xvarname')
			qui gegen `x'_mean = mean(`x'), by(`xvarname')
			qui egen tag = tag(`xvarname') if !mi(`y'_mean)
			qui replace `y'_mean = . if tag!=1
		}
		if "`by'"!="" {
			qui levelsof `by', local(bynum)
			foreach bynum2 in `bynum' {
				qui if "`discrete'"=="" gquantiles `x'_q`bynum2' = `x' if `by'==`bynum2', xtile nq(`nquantiles')
				qui if "`discrete'"!="" gen `x'`bynum2' = `x' if `by'==`bynum2'
				qui gegen `y'_mean`bynum2' = mean(`y') if `by'==`bynum2', by(`xvarname'`bynum2')
				qui gegen `x'_mean`bynum2' = mean(`x') if `by'==`bynum2', by(`xvarname'`bynum2')
				qui egen tag`bynum2' = tag(`xvarname'`bynum2') if !mi(`y'_mean`bynum2')
				qui replace `y'_mean`bynum2' = . if tag`bynum2'!=1
				
			}
		}
	}
	
	// overall scheme options
	if "`colorscheme'"=="" {
		local colorscheme Set1
		local colorschemeopts int(1.8)
	}
	qui colorpalette `colorscheme', `colorschemeopts' nograph local(,prefix(c) nonames)
	foreach i of numlist 25 50 75 {
		qui colorpalette `colorscheme', `colorschemeopts' op(`i') nograph local(,prefix(c) suffix(o`i') nonames) 
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
	if (strpos("`fit'","lfit") | "`fit'"=="") & "`coef'"!="" & "`by'"=="" {
		qui sum `yplot'
		local textplacey = r(max)+.1*r(max)
		qui sum `xplot'
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
	if "`by'"=="" {
		if !strpos("`fit'","ci") local legopts `legopts_all' legend(order(1 "Observed" 2 "`leg_fit' fit"))
		if strpos("`fit'","ci") local legopts `legopts_all' legend(order(1 "Observed" 3 "`leg_fit' fit" 2 "95% CIs"))
	}
	else {
		qui egen distinctby = group(`by')
		qui sum distinctby
		local maxdistinctby = r(max)
		local coln = 0
		qui levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			local coln = `coln'+1
			if !strpos("`fit'","ci") local coln2 = `coln'+`maxdistinctby'
			if strpos("`fit'","ci") local coln2 = `coln'+`maxdistinctby'+`coln'
			local legorder `legorder' `coln' "`bylab`bynum2''" `coln2' ""
		}
		local legopts `legopts_all' legend(order(`legorder') col(2) textfirst)
	}
	
	// scatter marker options
	if "`by'"=="" | "`binned'"=="" {
		qui count if !mi(`xplot') & !mi(`yplot')
		local n = r(N)
	}
	else {
		local n = 0
		qui levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			qui count if !mi(`xplot'`bynum2') & !mi(`yplot'`bynum2')
			local n = `n'+r(N)
		}
	}
	if `n'<=100 local scattermarkers mfullscatterm
	if `n'>100 & `n'<=300 local scattermarkers mscatter75m
	if `n'>300 & `n'<=2000 local scattermarkers mscatter50m
	if `n'>2000 local scattermarkers mscatter25m

	// generate the plot
	local lscatteropts `plotscheme' `xtitle' `ytitle' `printcoef2' `opts' `legopts'
	
	if "`by'"=="" {
		tw (scatter `yplot' `xplot', ``scattermarkers'2' `jitter') (`fittype' `y' `x', `o1' `polybw2'), `lscatteropts'
	}
	
	if "`by'"!="" {
		local coln = 0
		qui levelsof `by', local(bynum)
		foreach bynum2 in `bynum' {
			if "`binned'"!="" local bynum3 `bynum2'
			local coln = `coln'+1
			local sc `sc' (scatter `yplot'`bynum3' `xplot'`bynum3' if `by'==`bynum2', ``scattermarkers'`coln'' `jitter') 
			local pl `pl' (`fittype' `y' `x' if `by'==`bynum2', `o`coln'' `polybw2')
		}

		tw `sc' `pl', `lscatteropts'
	}
	
	restore
end




