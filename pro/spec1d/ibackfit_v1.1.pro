; THIS VERSION HAS A LOT OF USEFUL CODE, BUT BASICALLY IT SHOULD NOT
; BE USED OVER ISPECLINEFIT.
;
;+
; NAME:
;	IBACKFIT()
;
; PURPOSE:
;	Fit a linear combination of spectral templates to a galaxy
;	spectrum.
;
; CALLING SEQUENCE:
;       backfit = ibackfit(flux,wave,starflux,invvar=,$ 
;          mask=,nback=,ndust=,nmonte=,redindx=,_extra=extra,$
;          /backnegative,/silent)
;
; INPUTS:
;	flux     - input spectrum [NPIX]
;	wave     - corresponding rest wavelength vector [NPIX]
;	starflux - stellar template spectra [NPIX,NSTAR] with the same
;                  wavelength spacing as the data
;
; OPTIONAL INPUTS:
;	invvar        - inverse variance spectrum [NPIX]
;	nback         - number of background polynomial terms to
;                       append to STARFLUX (default: 0)
;       nmonte        - number of Monte Carlo realizations to get the
;                       error in the continuum fit (default: 30) 
;	ndust         - number of reddening terms (= NSTAR or = 0)
;	redindx       - integer array [NSTAR]; stars with the same
;                       REDINDX are constrained to have the same
;                       reddening (default: one reddening for all
;                       templates) 
;       extra         - keywords for K_LAMBDA() to specify the
;                       reddening curve (default O'Donnell 1994)
;	
; KEYWORD PARAMETERS:
;       backnegative - allow the coefficients to go negative
;       silent       - suppress output to STDOUT
;
; OUTPUTS:
;	backfit      - output data structure with the following fields 
;          continuum       - best-fitting continuum
;          continuum_sigma - error in CONTINUUM
;          chi2            - reduced chi2 of the fit, excluding masked
;                            pixels
;          nstar           - number of stellar templates
;          nback           - number of polynomial background terms 
;          maskpix         - masked pixels that were masked - see
;                            EMISSION_MASK()
;          nmonte          - number of Monte Carlo realizations
;          starebv         - best-fitting color excess [NSTAR]
;          redindx         - INPUT
;          starcoeff       - best-fitting coefficients [NSTAR+NBACK] 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURE:
;       We use a non-negative least squares matrix solution to
;       determine an initial guess for the coefficients of the fit.
;       We then improve upon these coefficients using MPFIT() and an
;       optional extinction curve.
;
; COMMENTS:
;	The amount of dust reddening is constrained to be between 0
;	and 1.0 in E(B-V).  Some of the MPFIT code is based on code in
;	D. Schlegel's LINEBACKFIT routine, especially constraining the
;	reddening parameters.
;
; EXAMPLE:
;
; INTERNAL SUPPORT ROUTINES:
;       BVLSBACKFIT(), BACKMODEL(), MPBACKFIT()
;
; PROCEDURES USED:
;	POLY_ARRAY(), BSPLINE_ITERFIT(), MPFITFUN(), EMISSION_MASK(), 
;	SPLOG, K_LAMBDA(), BVLS
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 May 14-15, U of A - major rewrite of
;          existing code
;       jm03jun09uofa - various bug fixes
;       jm03sep18uofa - bug fixes, especially in the treatment of
;                       reddening; added SILENT keyword
;       jm03nov06uofa - added STARCONTINUUM input
;       jm04jan07uofa - bug fix in how STARCONTINUUM is handled
;       jm04jan14uofa - added BACKNEGATIVE keyword
;       jm04mar08uofa - minor bug fix in MPBACKFIT when NDUST=0; added
;                       ALPHA parameter
;       jm04apr29uofa - added POLYCOEFF structure tag
;       jm04may11uofa - pass STEP to MPFIT via the PARINFO structure;
;                       use MPCHI2 rather than my own chi2 value 
;       jm04sep01uofa - various changes to help convergence, thanks to
;                       work by Christy Tremonti; return MPSTATUS and
;                       MPNITER 
;       jm06feb08uofa - do not change the dimensions of STARFLUX if
;                       NBACK>0 
;
; Copyright (C) 2002-2004, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function bvlsbackfit, fitflux, fitsigma, fitstarflux, nmonte=nmonte, silent=silent
; ---------------------------------------------------------------------------    
; BVLS continuum fitting
; ---------------------------------------------------------------------------    

    nstar = (size(fitstarflux,/dimension))[1]
    npix = n_elements(fitflux)
    
    bnd = fltarr(2,nstar) ; boundaries
    bnd[1,*] = 1.0

    if not keyword_set(silent) then begin
       if nmonte gt 0L then $
         splog, 'Monte Carlo fitting the continuum '+strn(nmonte)+' times with BVLS.' else $
         splog, 'Fitting the continuum with BVLS.'
    endif
;   t1 = systime(/seconds)
    for j = 0L, nmonte do begin
        
       Abvls = fitstarflux
       if j eq 0L then Bbvls = fitflux else Bbvls = fitflux + randomn(seed,npix)*fitsigma ; perturbed flux
       
       bvls, Abvls, Bbvls, bnd, coeff, rnorm, nsetp, w, $
         index, ierr, eps=eps, itmax=itmax, iter=iter

       if j eq 0L then eigencoeff = coeff else eigencoeff = [ [eigencoeff], [coeff] ]
       
    endfor 
;   print, systime(/seconds)-t1

return, eigencoeff
end

function mpbackfit, fitflux, fitwave, fitsigma, fitstarflux, redindx=redindx, $
  redcurve=redcurve, coeffguess=coeffguess, nstar=nstar, ndust=ndust, $
  alpha=alpha, nmonte=nmonte, backnegative=backnegative, mpchi2=mpchi2, $
  mpstatus=mpstatus, mpniter=mpniter, silent=silent
; ---------------------------------------------------------------------------    
; MPFIT() continuum fitting    
; ---------------------------------------------------------------------------    

    npix = n_elements(fitflux)
    nalpha = n_elements(alpha)
    
    ntemplate = (size(fitstarflux,/dimension))[1] ; number of fitting templates, including polynomials
    nback = ntemplate - nstar                     ; number of polynomials
    nparams = ntemplate + ndust + nalpha          ; total number of parameters 
    
; compute the template normalization factors and normalize the data 

    normstar = dblarr(ntemplate)
    
    for istar = 0L, ntemplate-1L do begin

       get_element, fitwave, djs_mean(fitwave)+[-25.0,+25.0], xx
       normstar[istar] = median(fitstarflux[xx[0]:xx[1],istar])
       fitstarflux[*,istar] = fitstarflux[*,istar] / normstar[istar]
       
    endfor

    normflux = median(fitflux)
    fitflux = fitflux / normflux
    fitsigma = fitsigma / normflux

; initialize the fitting parameters.  constrain the stellar
; contribution to be positive and also constrain the range of the
; reddening

    parinfo = {$
      parname:     'Coeff', $
      value:          0.0D, $
      fixed:            0L, $
      limited:     [0L,0L], $
      limits:  [0.0D,0.0D], $
      step:            0.0, $ ; fix(0) = automatic
      relstep:         0.0, $ ; fix(0) = automatic
      mpside:            0, $ ; automatic
      mpmaxstep:         0, $ ; automatic
      tied:             '', $
      mpprint:           0}
    parinfo = replicate(parinfo,nparams)

    parinfo[0:ntemplate-1L].step = 0.0      ; template step size 
    parinfo[0:ntemplate-1L].relstep = 0.05  ; template step size 

    if (ndust gt 0L) then begin ; note that NALPHA = NDUST

       parinfo[ntemplate].mpprint = 1
       parinfo[ntemplate:ntemplate+ndust-1L].parname = 'E(B-V)'
       parinfo[ntemplate:ntemplate+ndust-1L].step = 0.005         ; E(B-V) step size 
       parinfo[ntemplate:ntemplate+ndust-1L].relstep = 0.0        ; E(B-V) relative step size 
       parinfo[ntemplate:ntemplate+ndust-1L].limited = [1L,1L]    ; lower/upper constraint on E(B-V)
       parinfo[ntemplate:ntemplate+ndust-1L].limits = [0.0D,6.0D] ; lower/upper limits on E(B-V)

       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].parname = 'alpha'
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].step = 0.001          ; step size 
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].relstep = 0.0         ; relative step size 
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].fixed = alpha.fixed
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].value = alpha.value
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].limited = alpha.limited
       parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].limits = alpha.limits
       
    endif

    if (not keyword_set(backnegative)) then $
      parinfo[0:nstar-1L].limited[0] = 1L else $ ; the stellar contribution must be positive
      parinfo[0:ntemplate-1L].limited[0] = 0L

    if (n_elements(coeffguess) eq 0L) then begin

       parinfo[0L:ntemplate-1L].value = 1.0 ; template coefficients (**do not use 0.0**)
       if (ndust gt 0L) then parinfo[ntemplate:ntemplate+ndust-1L].value = 0.15 ; E(B-V)

    endif else begin

;      parinfo.step = 0.0    ; restore the MPFIT defaults
;      parinfo.relstep = 0.0

       parinfo[0L:ntemplate-1L].value = coeffguess[0L:ntemplate-1L] * normstar / normflux ; template coefficients
;      parinfo[0L:ntemplate-1L].value = coeffguess[0L:ntemplate-1L] ; template coefficients

       if (ndust gt 0L) then begin

          parinfo[ntemplate:ntemplate+ndust-1L].value = coeffguess[ntemplate:ntemplate+ndust-1L] ; E(B-V)
          parinfo[ntemplate+ndust:ntemplate+ndust+nalpha-1L].value = $
            coeffguess[ntemplate+ndust:ntemplate+ndust+nalpha-1L]

       endif

;      struct_print, parinfo
       
    endelse

; initialize the reddening parameters
    
    if (ndust gt 0L) then begin 

; tie reddening values to one another if requested
    
       allindx = redindx[uniq(redindx,sort(redindx))]
       for i = 0L, n_elements(allindx)-1L do begin
          ii = where(redindx eq allindx[i],count)
          if (count gt 1L) then for j = 1L, count-1L do $
            parinfo[ntemplate+ii[j]].tied = string(ntemplate+ii[0],format='("P[",I0,"]")')
       endfor 

; tie the values of ALPHA to one another if requested

       allindx = alpha[uniq(alpha.indx,sort(alpha.indx))].indx
       for i = 0L, n_elements(allindx)-1L do begin
          ii = where(alpha.indx eq allindx[i],count)
          if (count gt 1L) then for j = 1L, count-1L do $
            parinfo[ntemplate+ndust+ii[j]].tied = string(ntemplate+ndust+ii[0],format='("P[",I0,"]")')
       endfor 
       
    endif 

; additional arguments

    functargs = {wave: fitwave, redcurve: redcurve, nstar: nstar, $
      nback: nback, ndust: ndust}

;   struct_print, parinfo
    
; Monte Carlo the fit to find the error in the continuum

    if not keyword_set(silent) then begin
       if nmonte gt 0L then $
         splog, 'Monte Carlo fitting the continuum '+strn(nmonte)+' times with MPFIT.' else $
         splog, 'Fitting the continuum with MPFIT.'
    endif

;   t1 = systime(/seconds)
    oflux = fitflux 
    for j = 0L, nmonte do begin

       if j gt 0L then fitflux = oflux + randomn(seed,npix)*fitsigma ; perturbed flux

       if (n_elements(coeffguess) eq 0L) then begin
       
;         struct_print, parinfo

          coeff = mpfitfun('backmodel',fitstarflux,fitflux,fitsigma,$
            /autoderivative,bestnorm=bestnorm,covar=covar,dof=dof,$
            ftol=ftol,functargs=functargs,gtol=gtol,maxiter=30,yfit=yfit,$
            nfev=nfev,nfree=nfree,niter=mpniter,npegged=npegged,$
            parinfo=parinfo,perror=perror,quiet=1,resdamp=0,$
            status=mpstatus,xtol=xtol)
 
          splog, 'INIT  MPFIT: nfev=', nfev, ', niter=', mpniter, ', status=', mpstatus

; only set the reddening and alpha values here, not the template values          

          if (mpstatus gt 0L) then begin
             if (ndust gt 0L) or (nalpha gt 0L) then parinfo[ntemplate:ntemplate+ndust+nalpha-1L].value = $
               coeff[ntemplate:ntemplate+ndust+nalpha-1L]
             template_coeffs = coeff[0L:ntemplate-1L]
             zero = where(template_coeffs eq 0.0,nzero)
             if (nzero ne 0L) then template_coeffs[zero] = 0.1
             parinfo[0L:ntemplate-1L].value = template_coeffs
          endif

;         if (mpstatus gt 0L) then parinfo.value = coeff

       endif 
          
; use better starting values; also hack the STEP structure tag to make
; MPFIT figure out an appropriate value          
          
       parinfo.step = 0.0
       parinfo.relstep = 0.0

;      struct_print, parinfo

       coeff = mpfitfun('backmodel',fitstarflux,fitflux,fitsigma,$
         /autoderivative,bestnorm=bestnorm,covar=covar,dof=dof,$
         ftol=ftol,functargs=functargs,gtol=gtol,maxiter=maxiter,$
         nfev=nfev,nfree=nfree,niter=mpniter,npegged=npegged,yfit=yfit,$
         parinfo=parinfo,perror=perror,quiet=1,resdamp=0,$
         status=mpstatus,xtol=xtol)
 
;      niceprint, coeff

       mpchi2 = sqrt(bestnorm/dof)
       splog, 'FINAL MPFIT: nfev=', nfev, ', niter=', mpniter, ', status=', mpstatus
       if (mpstatus eq 5L) then splog, 'WARNING: Maximum number of iterations reached: ', mpniter

       if j eq 0L then eigencoeff = coeff else eigencoeff = [ [eigencoeff], [coeff] ]

    endfor 
;   print, systime(/seconds)-t1

    eigencoeff[0L:ntemplate-1L] = eigencoeff[0L:ntemplate-1L] * normflux / normstar

return, eigencoeff
end

function ibackfit, flux, wave, starflux1, invvar=invvar, mask=mask, $
  nback=nback, coeffguess=coeffguess, alpha=alpha, ndust=ndust, nmonte=nmonte, $
  redindx=redindx, _extra=extra, backnegative=backnegative, silent=silent

    light = 2.99792458D5

    if n_elements(flux) eq 0L then begin
       print, 'Syntax = backfit = ibackfit(flux,wave,starflux,invvar=,mask=$'
       print, '   nback=,coeffguess=,alpha=,ndust=,nmonte=,redindx=,$'
       print, '   _extra=extra,/backnegative,/silent)'
       return, -1
    endif 
       
    ndim = size(flux,/n_dimension) 
    if ndim ne 1L then message, 'FLUX vector is not one dimensional.' 
    fsize = size(flux,/dimension)
    npix = fsize[0]

    if (n_elements(wave) eq 0L) or (n_elements(wave) ne npix) then $
      message, 'Either WAVE is not defined or FLUX and WAVE have incompatible dimensions.'

    if n_elements(invvar) eq 0L then invvar = replicate(1.0,npix) else $
      if n_elements(invvar) ne npix then message, 'FLUX and INVVAR have incompatible dimensions.'

    if n_elements(mask) eq 0L then mask = replicate(1,npix) else $
      if n_elements(mask) ne npix then message, 'FLUX and MASK have incompatible dimensions.'

    if n_elements(nback) eq 0L then nback = 0L
    if n_elements(nmonte) eq 0L then nmonte = 0L

    stardim = size(starflux1,/n_dimension)
    starsize = size(starflux1,/dimension)
    if starsize[0] ne npix then message, 'FLUX and STARFLUX have incompatible dimensions.'
    if stardim eq 1L then nstar = 1L else nstar = starsize[1]

    if n_elements(redindx) eq 0L then redindx = lonarr(nstar) else $
      if n_elements(redindx) ne nstar then message, 'STARFLUX and REDINDX have incompatible dimensions.'

    if n_elements(ndust) eq 0L then ndust = nstar
    nalpha = n_elements(alpha)

; reddening curve    
    
    if (ndust eq 0L) then redcurve = fltarr(npix) else $
      redcurve = k_lambda(wave,_extra=extra)

; add polynomial templates to the stellar templates and normalize;
; exclude the "constant" polynomial (order=0)

    if (nback gt 0L) then begin

       polyflux = poly_array(npix,nback);*djs_median(flux)

;      polyflux = dblarr(npix,nback)
;      xpix = findgen(npix)/npix
;      for i = 0L, nback-1L do polyflux[*,i] = xpix^(i+1)
       
       starflux = [ [starflux1],[polyflux] ]       

    endif else starflux = starflux1 ; jm06feb08uofa

    ntemplate = nstar + nback ; number of stellar plus polynomial templates
    
    fitinvvar = invvar*mask
    good = where(fitinvvar gt float(0),ngood,comp=maskpix,ncomp=nmask)

    fitflux = flux[good]
    fitwave = wave[good]
    fitsigma = 1.0/sqrt(fitinvvar[good])

    fitstarflux = starflux[good,*]

; fit the continuum.  find the best-fitting solution with BVLS and, in
; general, Monte Carlo this solution to obtain the error in the
; coefficients.  then input this guess into MPFIT to refine the
; solution and to solve for the reddening
    
; BVLS

;   bvlscoeff = bvlsbackfit(fitflux,fitsigma,fitstarflux,nmonte=nmonte,silent=silent)
;   if nmonte gt 0L then coeffguess = bvlscoeff[*,0] else coeffguess = bvlscoeff
;   bvlsfit = starflux # bvlscoeff
    
; MPFIT    

    mpfitcoeff = mpbackfit(fitflux,fitwave,fitsigma,fitstarflux,redindx=redindx,$
      redcurve=redcurve[good],coeffguess=coeffguess,nstar=nstar,ndust=ndust,$
      alpha=alpha,nmonte=0,backnegative=backnegative,mpchi2=mpchi2,mpstatus=mpstatus,$
      mpniter=mpniter,silent=silent)

    continuum = backmodel(starflux,mpfitcoeff,wave=wave,$
      redcurve=redcurve,nstar=nstar,nback=nback,ndust=ndust)

; if there are stellar templates *and* polynomial templates then only
; store the stellar coefficients

    if (ndust gt 0L) then begin
       staralpha = mpfitcoeff[ntemplate+ndust:ntemplate+ndust+nalpha-1L] ; alpha
       reddening = staralpha*mpfitcoeff[ntemplate:ntemplate+ndust-1L]    ; E(B-V)
    endif else begin
       staralpha = dblarr(nstar)+1.0D
       reddening = dblarr(nstar)
    endelse

    if (nback gt 0L) then $
      polycoeff = mpfitcoeff[nstar:ntemplate-1L] else $
      polycoeff = dblarr(nstar)

    if (nstar gt 0L) and (nback gt 0L) then begin
       starcoeff = mpfitcoeff[0L:nstar-1L] ; final coefficients
    endif else begin
       starcoeff = mpfitcoeff[0L:ntemplate-1L]
    endelse

; evaluate the error in the continuum if NMONTE > 0

    if nmonte gt 0L then begin
    
       continuum_array = continuum*0.0 # (fltarr(nmonte)+1)
       for k = 0L, nmonte-1L do continuum_array[*,k] = backmodel(starflux,$
         [bvlscoeff[*,k+1],reddening],wave=wave,redcurve=redcurve,nstar=nstar,$
         nback=nback,ndust=ndust)
       continuum_sigma = continuum*0.0

       for i = 0L, npix-1L do continuum_sigma[i] = stddev(continuum_array[i,*])
       
    endif else continuum_sigma = continuum*0.0
    
;   djs_plot, wave, flux, xsty=3, ysty=3, ps=10
;   for j = 0L, nmask-1L do plots, wave[maskpix[j]], flux[maskpix[j]], ps=4, color=djs_icolor('red')
;   djs_oplot, wave, continuum, color='green', thick=2.0
;   cc = get_kbrd(1)
;   djs_plot, wave, flux-continuum, xsty=3, ysty=3, ps=10

    snr = djs_median(sqrt(fitinvvar[good])*flux[good])
    chi2 = total(fitinvvar[good]*(flux[good]-continuum[good])^2.0) / $
      (ngood-ntemplate-1L)
;   splog, chi2, mpchi2
    chi2 = mpchi2
    if (not keyword_set(silent)) then splog, 'Reduced chi2 of the continuum fit = ', chi2

    if (ndust gt 0L) and (not keyword_set(silent)) then $
      splog, 'Best-fitting E(B-V) = ', reddening[uniq(reddening)]
    
    backfit = {$
      log_continuum:  float(continuum), $
      continuum_chi2: float(chi2), $
      continuum_snr:  float(snr), $
      nstar:          long(nstar), $
      nback:          long(nback), $
      nmonte:         long(nmonte), $
      redindx:        long(redindx), $
      starebv:        float(reddening), $
      starcoeff:      double(starcoeff), $ ; [M_sun/cm2]
      staralpha:      float(staralpha), $
      polycoeff:      double(polycoeff), $ ; [M_sun/cm2]
      mpfitcoeff:     double(mpfitcoeff),$
      mpfit_status:   mpstatus, $
      mpfit_niter:    mpniter}

;  djs_plot, wave, flux, ps=10, xsty=3, ysty=3
;  ploterror, wave, flux, 1.0/sqrt(invvar), ps=10, xsty=3, ysty=3
;  djs_oplot, wave, continuum, color='red', ps=10

return, backfit
end    
