;+
; NAME:
;   IBACKFIT()
;
; PURPOSE:
;   Fit a non-negative linear combination of spectral templates to
;   a galaxy spectrum.
;
; INPUTS:
;   wave         - *observed* wavelength vector [NPIX]
;   flux         - *observed* spectrum (erg/s/cm2/A) [NPIX]  
;   templateflux - template models [NPIX,NTEMPLATE]
;
; OPTIONAL INPUTS:
;   zobj        - 
;   templateres - rest-frame FWHM resolution of the templates (default
;                 0.0 A)
;   specres     - *observed* FWHM instrumental resolution
;
;   invvar    - inverse variance spectrum [NPIX]
;   zsnrcross - the default behavior is to tweak the redshift by
;               cross-correlating the data with the best-fitting
;               model, but only if the median S/N>ZSNRCROSS; to
;               suppress the cross-correlation set ZSNRCROSS to a
;               negative number (default 2.0)

;   extra     - keywords for K_LAMBDA() to specify the
;               reddening curve (default O'Donnell 1994) 
;   
; KEYWORD PARAMETERS:
;       backnegative - allow the coefficients to go negative
;       silent       - suppress output to STDOUT
;
; OUTPUTS:
;   backfit      - output data structure with the following fields 
;          continuum       - best-fitting continuum
;          continuum_sigma - error in CONTINUUM
;          chi2            - reduced chi2 of the fit, excluding masked
;                            pixels
;          ntemplate           - number of stellar templates
;          nback           - number of polynomial background terms 
;          maskpix         - masked pixels that were masked - see
;                            EMISSION_MASK()
;          nmonte          - number of Monte Carlo realizations
;          starebv         - best-fitting color excess [NTEMPLATE]
;          redindx         - INPUT
;          starcoeff       - best-fitting coefficients [NTEMPLATE+NBACK] 
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
;   The amount of dust reddening is constrained to be between 0
;   and 1.0 in E(B-V).  Some of the MPFIT code is based on code in
;   D. Schlegel's LINEBACKFIT routine, especially constraining the
;   reddening parameters.
;
; EXAMPLE:
;
; INTERNAL SUPPORT ROUTINES:
;       MPBACKFIT()
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 May 14-15, U of A - major rewrite of
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
;       jm06feb08uofa - do not change the dimensions of TEMPLATEFLUX if
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

function mpbackfit, flux, templateflux, invvar, redcurve=redcurve, $
  mpchi2=mpchi2, mpstatus=mpstatus, mpniter=mpniter, silent=silent
; ---------------------------------------------------------------------------    
; MPFIT() continuum fitting    
; ---------------------------------------------------------------------------    

    npix = n_elements(flux)
    ntemplate = (size(templateflux,/dimension))[1]
    
; initialize the fitting parameters

    parinfo = {$
      parname:     'Coeff', $
      value:          0.0D, $
      fixed:            0L, $
      limited:     [0L,0L], $
      limits:  [0.0D,0.0D], $
      tied:             '', $
      mpprint:           0}
    parinfo = replicate(parinfo,ntemplate+1L)

; template initial parameters

    parinfo[0:ntemplate-1L].limited[0] = 1L ; enforce positivity
    parinfo[0L:ntemplate-1L].value = 0.5    ; coefficients (**do not use 0.0**)

    parinfo[ntemplate].mpprint = 1
    parinfo[ntemplate].parname = 'E(B-V)'
    parinfo[ntemplate].limited = [1L,1L]    ; lower/upper constraint on E(B-V)
    parinfo[ntemplate].limits = [0.0D,6.0D] ; lower/upper limits on E(B-V)
    parinfo[ntemplate].value = 0.1          ; initial guess
;   struct_print, parinfo

    functargs = {redcurve: redcurve}
    coeff = mpfitfun('ibackmodel',templateflux,flux,weights=invvar,$
      /autoderivative,bestnorm=bestnorm,covar=covar,dof=dof,ftol=ftol,$
      functargs=functargs,gtol=gtol,maxiter=maxiter,nfev=nfev,nfree=nfree,$
      niter=mpniter,npegged=npegged,resdamp=0,yfit=model,parinfo=parinfo,$
      perror=perror,quiet=1,status=mpstatus,xtol=xtol)
    if (not keyword_set(silent)) then splog, 'MPFIT: nfev=', $
      nfev, ', niter=', mpniter, ', status=', mpstatus

    if (dof le 0.0) then mpchi2 = -1.0 else mpchi2 = sqrt(bestnorm/dof)
;   mpchi2 = total(invvar*(flux-model)^2.0)/float((ngood-(ntemplate+1L)-1L))
    
return, coeff
end

function ibackfit, wave, flux1, invvar1, templatewave, templateflux, $
  zobj=zobj, templateres=templateres, specres=specres, vdisp=vdisp, $
  normwave=normwave, normdwave=normdwave, zsnrcross=zsnrcross, $
  vmaxtweak=vmaxtweak, vminwave=vminwave, vmaxwave=vmaxwave, niterfit=niterfit, $
  continuum_minwave=continuum_minwave, continuum_maxwave=continuum_maxwave, $
  pspacefactor=pspacefactor, silent=silent

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    npix = n_elements(wave)
    if (npix eq 0L) or (n_params() lt 5) then begin
       doc_library, 'ibackfit'
       return, -1
    endif

; template parameters

    ntemplatedim = size(templateflux,/n_dim)
    templatedim = size(templateflux,/dim)
    ntemplate = templatedim[1]
    if (ntemplatedim eq 3) then Znindx = templatedim[2] else Znindx = 1

    if (n_elements(vdisp) eq 0L) then vdisp = 150.0 ; [km/s]
    if (n_elements(templateres) eq 0L) then templateres = 0.0 ; bad!
    if (n_elements(specres) eq 0L) then specres = 0.0 ; bad!

    if (n_elements(niterfit) eq 0L) then niterfit = 2L
    if (n_elements(zsnrcross) eq 0L) then zsnrcross = 2.0
    if (n_elements(pspacefactor) eq 0L) then pspacefactor = 0.03
    if (n_elements(vmaxtweak) eq 0L) then vmaxtweak = 500.0 ; [km/s]

; define some wavelength parameters    

    goodwave = where((wave gt 0.0) and (invvar1 gt 0.0) and $
      (finite(invvar1) eq 1),ngoodwave)
    if (ngoodwave eq 0L) then begin
       splog, 'This is very, very bad'
       return, -1
    endif
    
    if (n_elements(continuum_minwave) eq 0L) then continuum_minwave = min(wave[goodwave])
    if (n_elements(continuum_maxwave) eq 0L) then continuum_maxwave = max(wave[goodwave])
    if (n_elements(normwave) eq 0L) then normwave = mean(wave[goodwave])
    if (n_elements(normdwave) eq 0L) then normdwave = 50.0 ; normalization width [Angstrom]

    if (n_elements(vminwave) eq 0L) then vminwave = min(wave[goodwave])
    if (n_elements(vmaxwave) eq 0L) then vmaxwave = max(wave[goodwave])

; begin main iteration loop

    for iZ = 0L, Znindx-1L do begin ; loop on metallicity

       if (not keyword_set(silent)) then splog, $
         format='("Fitting metallicity template '+$
         '",I1,"/",I1)', iZ+1, Znindx

       templateflux1 = reform(templateflux[*,*,iZ])
       zabs = zobj & zabs_err = 0.0

       logwave = alog10(wave)
       pixsize= abs(alog(10.0)*(logwave[npix-1]-$ ; mean pixel size [log-A]
         logwave[0])/double(npix)) 

; normalize the data

       inrange = where((wave ge normwave-normdwave/2.0) and $
         (wave le normwave+normdwave/2.0) and (invvar1 gt 0.0),$
         ninrange)
       if (ninrange eq 0L) then message, 'This is very, very bad'

       fluxnorm = abs(median(flux1[inrange])) ; note absolute value!
       flux = flux1/fluxnorm
       invvar = invvar1*fluxnorm^2.0

; initialize the continuum mask
       
       mask = flux*0.0+1.0
       mask = mask * ((wave ge continuum_minwave) and (wave le continuum_maxwave))

       for iter = 0L, niterfit-1L do begin 

; compute the *observed frame* instrumental velocity dispersion, and
; the total velocity dispersion; note that in general the instrumental
; dispersion varies with wavelength, but we ignore that here and
; choose the value at NORMWAVE

          instr_vdisp = interpol(light*sqrt(specres^2.0-$ ; [sigma, km/s]
            ((1.0+zabs)*templateres)^2.0)/fwhm2sig/wave,wave,$
            normwave)>0.0
          vdisp_tot = sqrt(vdisp^2.0+instr_vdisp^2.0) ; [km/s]

; redshift the models, resample onto the logarithmic wavelength vector
; of the data, and then convolve by to the total velocity dispersion;
; finally, normalize each template at a single wavelength NORMWAVE
; (see, e.g., Cid-Fernandes et al. 2004,2005,2007)

          ztemplatewave = templatewave*(1.0+zabs)
          ztemplateflux = dblarr(n_elements(wave),ntemplate)
          ztemplatenorm = dblarr(ntemplate)
          for ii = 0L, ntemplate-1L do begin
             linterp, alog10(ztemplatewave), templateflux1[*,ii], logwave, ztemplateflux1
             ztemplateflux[*,ii] = k_smooth(logwave,ztemplateflux1,vdisp_tot)
             ztemplatenorm[ii] = interpol(ztemplateflux[*,ii],logwave,alog10(normwave))
             ztemplateflux[*,ii] = ztemplateflux[*,ii]/ztemplatenorm[ii]
          endfor

; define the reddening curve (in the rest-frame!); rebin the reddening
; vector to maximize vector operations in IBACKMODEL()

          redcurve = k_lambda(wave/(1.0+zabs),/charlot)
          redcurve = rebin(reform(redcurve,npix,1),npix,ntemplate)

; add pixels affected by telluric absorption, sky-subtraction
; residuals, and emission lines to the bad pixel mask and then fit! 
          mask = mask and iemission_mask(wave,z=zobj,$
            /telluric,/sky,/nebular,/qso)
          good = where((invvar*mask gt 0.0) and (finite(invvar) eq 1),ngood)
          if (ngood eq 0L) then message, 'Definitely bad'

          t0 = systime(1)
          mpfitcoeff = mpbackfit(flux[good],ztemplateflux[good,*],$
            invvar[good],redcurve=redcurve[good,*],mpchi2=mpchi2,$
            mpstatus=mpstatus,mpniter=mpniter,silent=silent)
          if (not keyword_set(silent)) then splog, $
            format='("CPU time for continuum fitting on '+$
            'iteration ",I1,"/",I1,": ",G0," seconds")', iter+1, $
           niterfit, systime(1)-t0
          if (not keyword_set(silent)) then splog, 'Reduced chi2 '+$
            'of the continuum fit = ', mpchi2

; special case of all zero coefficients; set the reddening to zero
; (sometimes it returns a non-zero value); then rebuild the
; best-fitting continuum model
          
          if (total(mpfitcoeff[0:ntemplate-1]) eq 0.0) then begin
             splog, 'All template coefficients zero!'
             mpfitcoeff[0:ntemplate] = 0.0D
             continuum = wave*0.0 ; note float()
          endif else begin
             continuum = fluxnorm*ibackmodel(ztemplateflux,$
               mpfitcoeff,redcurve=redcurve)
          endelse

          ebv = mpfitcoeff[ntemplate]
          if (not keyword_set(silent)) then splog, 'Best-fitting E(B-V) = ', ebv

          continuum_coeff = mpfitcoeff[0L:ntemplate-1L]*$ ; [M_sun/cm^2]
            fluxnorm/ztemplatenorm
          mpfitcoeff[0L:ntemplate-1L] = continuum_coeff

; reject deviant pixels and add them to the mask

          rejmask = mask*0.0+1.0
          djs_iterstat, flux1[good]-continuum[good], sigrej=3.0, mask=rejmask1
          rejmask[good] = rejmask1
          mask = mask and rejmask
          
; improve the redshift guess by cross-correlating the best-fitting
; continuum spectrum with the emission-line subtracted data spectrum
; but only for high signal-to-noise spectra

          good = where((invvar*mask gt 0.0) and (finite(invvar) eq 1) and $
            (wave ge vminwave) and (wave le vmaxwave),ngood)
          medsnr = median(flux1[good]*sqrt(invvar1[good])) > 0.0 ; note!
          if (zsnrcross gt 0.0) and (medsnr gt zsnrcross) and $
            (iter lt niterfit-1L) then begin          

             pmax = floor(alog10(1.0+vmaxtweak/light)/pixsize)
             zans = zcompute(flux1[good],invvar1[good],continuum[good],$
               pspace=pspacefactor*pmax,pmin=-pmax,pmax=pmax,doplot=0,/silent)
             if (zans.dof gt 0) and (zans.z_err gt 0.0) then begin
                ztweak = zabs + 10.0^(pixsize*zans.z)-1.0
                ztweak_err = abs(alog(10.0)*pixsize*zans.z_err*(1.0+ztweak))
                if (not keyword_set(silent)) then splog, 'ZCROSSCOR converged: z = '+$
                  string(zabs,format='(F0.0)')+' --> '+string(ztweak,format='(F0.0)')
                zabs = ztweak & zabs_err = ztweak_err
             endif else $
               if (not keyword_set(silent)) then splog, 'ZCROSSCOR did not converge'
          endif

       endfor ; close iteration loop

; pack the continuum-fitting results into a structure

       backfit1 = {$
         template_bestindx:           0, $
         z_abs:             float(zabs), $
         z_abs_err:         float(zabs_err), $
;        mpfitcoeff:        double(mpfitcoeff),$
         instr_vdisp:       float(instr_vdisp), $
;        add_vdisp:         float(add_vdisp), $
         continuum_vdisp:   float(vdisp), $
         continuum_minwave: float(continuum_minwave), $
         continuum_maxwave: float(continuum_maxwave), $
         continuum_chi2:    float(mpchi2), $
         continuum_snr:     float(medsnr), $
         continuum_ebv:     float(ebv), $
         continuum_coeff:   double(continuum_coeff), $ ; [M_sun/cm2]
         mpfit_status:      fix(mpstatus), $
         mpfit_niter:       fix(mpniter), $
         continuum:         float(continuum)}

       if (iZ eq 0L) then backfit = backfit1 else $
         backfit = [backfit,backfit1]

    endfor  ; close metallicity loop

; pick the best one    
    
    if (Znindx gt 1L) then begin
       minchi2 = min(backfit.continuum_chi2,Zbest)
       if (not keyword_set(silent)) then splog, 'Selecting best-fitting '+$
         'metallicity template '+string(Zbest+1,format='(I0)')
       backfit = backfit[Zbest]
       backfit.template_bestindx = Zbest
    endif

return, backfit
end    
