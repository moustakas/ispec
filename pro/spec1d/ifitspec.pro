;+
; NAME:
;       IFITSPEC
;
; PURPOSE:
;	Fully fit a galaxy spectrum. 
;
; INPUTS:
;       flux      - observed spectrum [erg/s/cm2/A] [NPIX]
;       wave      - observed wavelength [Angstrom] [NPIX] 
;       eigenflux - eigen-template flux vector [NEIGENPIX,NSTAR] or if
;                   ZMULTI=1 then an [NEIGENPIX,NSTAR,ZNINDX] vector 
;       eigenwave - wavelength vector for EIGENFLUX [NEIGENPIX]
;       eigenres  - scalar spectral resolution for EIGENFLUX [FWHM,
;                   Angstrom]  
;       specres   - observed scalar instrumental resolution for FLUX
;                   [FWHM, Angstrom] 
;       linepars  - emission line structure array from READ_LINEPARS()
;                   with the following structure tags [NLINE]
;          line   - emission line name (string)
;          wave   - central wavelength (air wavelength) [Angstrom]
;          zindex - see ILINEFIT() documentation
;          windex - see ILINEFIT() documentation
;          findex - see ILINEFIT() documentation
;          fvalue - see ILINEFIT() documentation
;
; OPTIONAL INPUTS:
;       invvar    - inverse variance spectrum corresponding to FLUX
;                   [NPIX]; if not provided then the output errors
;                   will not be correct
;       starvdisp - rest-frame stellar velocity dispersion (default
;                   100.0 km/s)
;       zobj      - approximate galaxy redshift (default 0.0) 
;       zline     - emission-line redshift guess (default ZOBJ)
;       snrcut    - compute upper limits on lines with signal to noise
;                   less than SNRCUT (default 1.0)
;       vmaxshift - if ZCROSSCOR=1 then do not allow the continuum
;                   redshift to change by more than +/-VMAXSHIFT
;                   [km/s] (default 300.0, see PROCEDURE)
;       vminwave  - if ZCROSSCOR=1 then minimum wavelength to consider
;                   for cross-correlation (default WAVE[0])
;       vmaxwave  - if ZCROSSCOR=1 then maximum wavelength to consider
;                   for cross-correlation (default WAVE[NPIX-1]) 
;       niter     - number of continuum fitting iterations (default 1) 
;       zsnrcross - only use ZCROSSCOR is the median S/N is >
;                   ZSNRCROSS (default 10)
;       maskwidth - mask all pixels within MASKWIDTH Angstroms of the
;                   each emission-line center (default 20 Angstrom)
;       zindex    - optionally overwrite LINEPARS.ZINDEX
;       windex    - optionally overwrite LINEPARS.WINDEX
;       findex    - optionally overwrite LINEPARS.FINDEX
;       fvalue    - optionally overwrite LINEPARS.FVALUE
;	extra     - extra parameters for IBACKFIT() and ILINEFIT()  
;	
; KEYWORD PARAMETERS:
;       debug          - generate diagnostic plots at various stages
;                        of the fitting (useful for debugging) 
;       zcrosscor      - cross-correlate the best-fitting continuum
;                        spectrum with FLUX to improve the
;                        absorption-line redshift (see PROCEDURE and
;                        ZSNRCROSS)
;       combine_blends - combine unresolved emission-line doublets in
;                        LINEFIT into a single emission line
;       Zmulti         - find the best-fitting metallicity template
;                        set (see EIGENFLUX)
;       no_medcontinuum - do not median-subtract the residuals
;
; OUTPUTS:
;       All outputs are optional.  For example, the fitting can be
;       carried out and examined with the keyword DEBUG.
;
; OPTIONAL OUTPUTS:
;       backfit     - see IBACKFIT()
;       linefit     - see ILINEFIT()
;       balmerabs   - see IBALMERABS()
;       indices     - see SPECTRAL_INDICES()
;       speclinefit - fitted emission-line spectrum [NPIX]
;       starflux    - identical to EIGENFLUX but broadened and
;                     wavelength resampled to match FLUX and WAVE
;                     [NPIX,NSTAR]
;       Zbestindx   - if ZMULTI=1 then return the best-fitting index
;                     corresponding to the third dimension of
;                     EIGENFLUX 
;
; PROCEDURE:
;       The essential pieces of this routine are a spectrum, an
;       arbitrary number of template spectra, and a data structure
;       (LINEPARS) specifying the emission lines to fit.  IFITSPEC
;       matches the wavelength range and spectral resolution of the
;       templates to the spectrum [IMATCH_TEMPLATES()] and finds the
;       best-fitting linear combination of templates for the spectrum
;       [IBACKFIT()] based on an initial guess of the redshift of the
;       galaxy.
;
;       The main iteration loop proceeds in the following way: [1]
;       Match the templates (in wavelength spacing and resolution) to
;       the data [IMATCH_TEMPLATES()]; [2] Mask pixels near emission
;       lines [EMISSION_MASK()]; [3] Fit the continuum [IBACKFIT()];
;       [4] Subtract the continuum spectrum and fit the pure
;       emission-line spectrum [ILINEFIT()]; [5] If ZCROSSCOR=1 then
;       improve the redshift of the continuum by cross-correlating the
;       best-fitting continuum spectrum with the emission-line
;       subtracted data spectrum [IM_ZTWEAK()]; [6] Steps [1-5] are
;       iterated a second time to improve the fitting; [7] Measure
;       color and Lick indices in the input spectrum
;       [SPECTRAL_INDICES()]; [8] Measure the Balmer absorption lines
;       from the best-fitting continuum spectrum [IBALMERABS()]; [9]
;       Measure Balmer- and forbidden-line equivalent widths using the
;       results of the simultaneous Gaussian fitting to constrain the
;       line profiles.
;
;       Upper limits on emission lines are also computed in 
;       IUPPER_LIMITS().  
;
;       In the emission-line fitting we identify three cases:
;       [1] well-measured line, [2] upper limits, and [3] un-measured
;       lines.  The error codes for cases [1] and [2] are detailed
;       below:
;
;       Case [2]: Upper limits are computed on lines that were either
;       dropped from the fit by MPFIT [see ILINEFIT()], or they have
;       S/N < SNRCUT.  They are identified by the following error
;       codes and quantities:
;
;          LINESIGMA         =  either tied to a well-detected line,
;                               or equal to LINERES at line center 
;          LINESIGMA_ERR     = -3.0
;
;          LINEAREA          =  computed in IUPPER_LIMITS()
;          LINEAREA_ERR      = -3.0
;
;          LINEBOX           =  not meaningful 
;          LINEBOX_ERR       =  not meaningful 
;
;          LINECONTLEVEL     =  computed in IUPPER_LIMITS()
;          LINECONTLEVEL_ERR = -3.0
;
;          LINEEW_AREA       =  based on LINEAREA and LINECONTLEVEL 
;          LINEEW_AREA_ERR   = -3.0
;
;          LINEEW_BOX        =  not meaningful
;          LINEEW_BOX_ERR    =  not meaningful
;
;          LINENPIX          =  not meaningful
;          LINEDOF           =  not meaningful
;          LINECHI2          = -3.0
;
;       Case [3]: Un-measured lines are either completely or partially
;       out of the wavelength range of the spectrum.  They have the
;       following (equivalent) error codes:
;
;          LINEZ             =  ZOBJ
;          LINEZ_ERR         = -2.0
;          LINESIGMA         =  0.0
;          LINESIGMA_ERR     = -2.0
;          LINEAREA          =  0.0
;          LINEAREA_ERR      = -2.0
;          LINEBOX           =  0.0
;          LINEBOX_ERR       = -2.0
;          LINECONTLEVEL     =  0.0
;          LINECONTLEVEL_ERR = -2.0
;          LINEEW_AREA       =  0.0
;          LINEEW_AREA_ERR   = -2.0
;          LINEEW_BOX        =  0.0
;          LINEEW_BOX_ERR    = -2.0
;          LINENPIX          =    0
;          LINEDOF           =  0.0
;          LINECHI2          = -2.0
;
; COMMENTS:
;       For best results the pixel size in Angstrom per pixel of the
;       galaxy and the eigen-templates should be constant (assumed
;       true).  
;
;       Currently both EIGENRES and SPECRES must be scalars.
;
; TODO:
;       [1] Allow EIGENRES and SPECRES to be a vectors.
;       [2] Give the option of not fitting the emission-line spectrum, 
;           or, vice-versa, not fitting the continuum.
;
; EXAMPLES:
;
; PROCEDURES USED:
;       SPLOG, EMISSION_MASK(), IBACKFIT(), IM_ZTWEAK(),
;       STRUCT_ADDTAGS(), ILINEFIT(), SPECTRAL_INDICES(),
;       IUPPER_LIMITS(), IBALMERABS(), ILOCAL_CONTINUUM(), FILL_EW(),
;       REMOVE, DJS_PLOT, DJS_OPLOT, LEGEND, DJS_MEAN(),
;       ICOMBINE_BLENDS(), BALMERINDX(), IM_STRUCT_TRIMTAGS
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002-2003, U of A - developed
;       jm03jul28uofa - added IFITSPEC_COMBINE_BLENDS() sub-routine
;                       and COMBINE_BLENDS keyword to properly treat
;                       the [O II] doublet 
;       jm03nov26uofa - many modifications and simplifications,
;                       especially the way in which equivalent widths
;                       are computed
;       jm04jan04uofa - added NITER optional input
;       jm04mar09uofa - added ZMULTI and ZBESTINDX keywords; cleaned
;                       up and streamlined the code and documentation;
;                       compute upper limits on lines with LINESIGMA=0 
;       jm04mar29uofa - use the *rest-frame* spectral resolution!
;       jm04apr07uofa - added ZSNRCROSS optional input
;       jm04may09uofa - added STARVDISP optional input; SPECRES must
;                       be a scalar; use new routine 
;                       NEW_IMATCH_TEMPLATES(); fit in log-lambda
;       jm04may13uofa - on the final iteration mask emission lines
;                       more liberally (e.g., all the high-order
;                       Balmer lines)
;       jm04jun09uofa - bug fix when two Balmer lines cannot be
;                       measured 
;       jm04sep13uofa - bug fix when the galaxy redshift is
;                       sufficiently high that no Balmer absorption
;                       lines are measured
;       jm04oct26uofa - bug fix: do not divide the fluxes by (1+z)
;                       since ILINEFIT() is called with rest-frame
;                       quantities 
;       jm04nov30uofa - error check if the redshift is so high that no
;                       emission lines fall in the spectral range
;       jm05jan20uofa - include spectral index measurements made off
;                       the best-fitting continuum spectrum
;       jm05jan28uofa - bug fix when matching Balmer emission lines to
;                       the Balmer absorption lines
;       jm06dec06nyu  - added ZLINE optional input
;       jm07apr10nyu  - always median smooth the continuum-subtracted
;                       spectrum, before fitting for the emission
;                       lines (needs to be tested!)
;       jm07aug30nyu  - added MEDSMOOTH_WINDOW and SMOOTH_WINDOW
;                       optional inputs
;       jm07sep06nyu  - added MASK_BACKFIT additional input
;       jm07sep27nyu  - added NO_MEDCONTINUUM keyword
;
; Copyright (C) 2002-2007, John Moustakas
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

pro ifitspec, flux, wave, eigenflux, eigenwave, eigenres=eigenres, specres=specres, $
  linepars=linepars, invvar=invvar, starvdisp=starvdisp, zobj=zobj, zline=zline, snrcut=snrcut, $
  vmaxshift=vmaxshift, vminwave=vminwave, vmaxwave=vmaxwave, niter=niter, $
  zsnrcross=zsnrcross, maskwidth=maskwidth, zindex=zindex, windex=windex, findex=findex, $
  fvalue=fvalue, _extra=extra, backfit=backfit, linefit=linefit, balmerabs=balmerabs, $
  indices=indices, speclinefit=speclinefit, starflux=starflux, Zbestindx=Zbestindx, $
  debug=debug, zcrosscor=zcrosscor, combine_blends=combine_blends, Zmulti=Zmulti, $
  force_chisq_one=force_chisq_one, medsmooth_window=medsmooth_window, smooth_window=smooth_window, $
  mask_backfit=mask_backfit, no_medcontinuum=no_medcontinuum

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

; defaults and error checking

    nwave = n_elements(wave)
    npix = n_elements(flux)
    neigenwave = n_elements(eigenwave)

    if (nwave eq 0L) or (npix eq 0L) or (neigenwave eq 0L) or (n_elements(eigenflux) eq 0L) then begin
       doc_library, 'ifitspec'
       return
    endif

    if (nwave ne npix) then begin
       splog, 'WAVE and FLUX do not have the same number of elements.'
       return
    endif

    eigen_ndim = size(eigenflux,/n_dimension)
    if keyword_set(Zmulti) then begin
       if (eigen_ndim ne 3L) then begin
          splog, 'EIGENFLUX must be a three-dimensional array.'
          return
       endif
    endif else begin
       if (eigen_ndim ne 2L) then begin
          splog, 'EIGENFLUX must be a two-dimensional array.'
          return
       endif 
    endelse

    eigendim = size(eigenflux,/dimension)
    nstar = eigendim[1]
    if keyword_set(Zmulti) then Znindx = eigendim[2] else Znindx = 1L

    if (neigenwave ne eigendim[0]) then begin
       splog, 'EIGENWAVE and EIGENFLUX do not have the same number of elements.'
       return
    endif

    nline = n_elements(linepars)
    if (nline eq 0L) then begin
       splog, 'Structure LINEPARS must be passed.'
       return
    endif
    
    if (n_elements(zobj) eq 0L) then zobj = float(0.0) else zobj = float(zobj)
    if (n_elements(zline) eq 0L) then zline = zobj
    if (n_elements(starvdisp) eq 0L) then starvdisp = 100.0 ; [km/s]
    
    nivar = n_elements(invvar)
    if nivar eq 0L then begin
       splog, 'WARNING:  INVVAR not passed - errors will be incorrect.'
       invvar = flux*0D0+1.0/(djs_median(flux)/50.0)^2.0 ; "intelligent" default INVVAR
    endif else begin
       if nivar ne npix then begin
          splog, 'FLUX and INVVAR do not have the same number of elements.'
          return
       endif
    endelse

    ferr = 1.0/sqrt(invvar)

    if (n_elements(mask_backfit) eq 0L) then mask_backfit = bytarr(npix)+1B else begin
       if (n_elements(mask_backfit) ne npix) then begin
          splog, 'FLUX and MASK_BACKFIT do not have the same number of elements.'
          return
       endif
    endelse

; error checking on the spectral resolution vectors

    nspecres = n_elements(specres)
    if (nspecres eq 0L) then begin
       splog, 'SPECRES must be specified.'
       return
    endif

    if (nspecres ne 1L) then begin
       splog, 'SPECRES must be a scalar.'
       return
    endif

;   if (nspecres ne npix) then begin
;      splog, 'FLUX and SPECRES do not have the same number of elements.'
;      return
;   endif

    neigenres = n_elements(eigenres)
    if (neigenres eq 0L) then begin
       splog, 'EIGENRES must be specified.'
       return
    endif

    if (neigenres ne 1L) then begin
       splog, 'EIGENRES must be a scalar.'
       return
    endif

    if n_elements(niter) eq 0L then niter = 2L
    niter = niter>1L

    if (n_elements(zsnrcross) eq 0L) then zsnrcross = 10.0
    if (n_elements(snrcut) eq 0L) then snrcut = 1.0
    if (n_elements(maskwidth) eq 0L) then maskwidth = 20.0         ; [Angstrom]
    if (n_elements(medsmooth_window) eq 0L) then medsmooth_window = 150L
    if (n_elements(smooth_window) eq 0L) then smooth_window = 50L

    if n_elements(vmaxshift) eq 0L then vmaxshift = 300.0   ; [km/s]
    if n_elements(vminwave) eq 0L then vminwave = min(wave) ; [Angstrom]
    if n_elements(vmaxwave) eq 0L then vmaxwave = max(wave) ; [Angstrom]

; ---------------------------------------------------------------------------
; initialize variables
; ---------------------------------------------------------------------------

; initialize variables for the iterative fitting

    newmaskwidth = maskwidth

    ziter = float(zobj)              ; absorption-line redshift guess [updated by IM_ZTWEAK()]
    zguess1 = replicate(zline,nline) ; emission-line redshift [overwritten by ILINEFIT()]
    zabs1 = float(zobj)              ; default absorption-line redshift and error
    zabs1_err = -999.0

    dwave = wave[1]-wave[0]       ; observed pixel size [Angstrom]
    coeff1 = 2D-4                 ; desired pixel size, ~60 km/s [log-Angstrom]
    logwave = alog10(double(wave))
    
    sigma_vdisp = starvdisp/coeff1/light ; stellar velocity dispersion [pixel]

; ---------------------------------------------------------------------------    
; begin main iteration loop: measure everything in the rest frame 
; ---------------------------------------------------------------------------    

;   if keyword_set(force_chisq_one) then niter = niter + 1L ; NOTE!

    for iZ = 0L, Znindx-1L do begin ; loop on metallicity

       if keyword_set(Zmulti) then $
         splog, format='("Fitting metallicity template '+$
           '",I1,"/",I1,".")', iZ+1, Znindx

;      chisq_one_factor = 1.0   ; see below
    
       eigenflux1 = reform(eigenflux[*,*,iZ])

       for iter = 0L, niter-1L do begin 

; shift everything to the rest frame; note that FLUX is in units of
; erg/s/cm2/A, and therefore it must be multiplied by (1+z); also
; compute the instrumental Gaussian broadening width in pixels 

          splog, 'Shifting to the rest frame-of-reference.'

          logrestwave = logwave - alog10(1D0+ziter)    ; rest wavelength [log-Angstrom]
          restwave = 10^logrestwave                    ; rest wavelength [Angstrom]
          restflux = flux*(1D0+ziter)                  ; rest flux [erg/s/cm2/A]
          restferr = ferr*(1D0+ziter)                  ; rest flux error [erg/s/cm2/A]
          restinvvar = 1D0/restferr^2.0                ; rest inverse variance
          restspecres = specres/(1D0+ziter)            ; rest spectral resolution [Angstrom]
          restdwave = restwave[1]-restwave[0]          ; rest pixel size [Angstrom/pixel]
;         restferr = chisq_one_factor*ferr*(1D0+ziter) ; rest flux error [erg/s/cm2/A]

          sigma_instr = sqrt((restspecres^2.0 - eigenres^2.0)>0.0)/fwhm2sig/restdwave ; [pixel]

; resample the spectrum such that each pixel has a constant width in
; km/s equal to ~60 km/s; COMBINE1FIBER is relegated here (use,
; instead, INTERPOL)
          
          coeff0 = logrestwave[0]
          lognpix = ceil(1D0 + (max(logrestwave)-coeff0)/coeff1)
          logvconstwave = coeff0 + coeff1*findgen(lognpix)
          vconstwave = 10^logvconstwave

          restfluxvconst = interpol(restflux,restwave,vconstwave)
          restinvvarvconst = sqrt(interpol(restinvvar,restwave,vconstwave)^2.0) ; prevents negative INVVAR!
          restferrvconst = 1.0/sqrt(restinvvarvconst)
          restmaskbackfitvconst = (interpol(mask_backfit,restwave,vconstwave) ge 1)
          
;;        inrestinvvar = restinvvar 
;;        combine1fiber, logrestwave, restflux, inrestinvvar, newloglam=logvconstwave, $
;;          newflux=restfluxvconst, newivar=restinvvarvconst
;;        restferrvconst = 1.0/sqrt(restinvvarvconst) ; RESTINVVARVCONST may have zeros

; mask emission-line wavelengths and pixels affected by telluric
; absorption; to get the telluric features right we need to pass WAVE
; and ZITER to EMISSION_MASK() rather than RESTWAVE and 0.0  

          splog, 'Masking +/-'+strtrim(string(newmaskwidth,format='(F12.1)'),2)+' Angstroms.'
          restmaskvconst = emission_mask(vconstwave*(1+ziter),z=ziter,$
            width=newmaskwidth,bluewidth=0.5*newmaskwidth,/telluric,/qsomask,$
            tellmask=tellmask)
;           spectrum=restfluxvconst,sigrej_cosmic=3.0,/cosmic,/qsomask)

;         restmaskvconst = emission_mask(vconstwave*(1+ziter),z=ziter,$
;           width=newmaskwidth,/telluric,bluemask=iter eq (niter-1L))
;         w = where(restmaskvconst eq 0)
;         djs_plot, vconstwave, restfluxvconst, ps=10, xsty=3, ysty=3
;         djs_oplot, vconstwave[w], restfluxvconst[w], ps=10, color='purple'
          
; prepare the eigen-templates for fitting; first, broaden to the
; instrumental resolution of the data; next, resample onto the rest
; wavelength vector of the galaxy (having constant velocity pixel
; size); finally, broaden by the stellar velocity dispersion 
          
          splog, 'Re-sampling and broadening the fitting templates.'
          starfluxvconst1 = restfluxvconst # replicate(1.0,nstar)

          for istar = 0L, nstar-1L do begin

             eigenflux1[*,istar] = gconvolve(eigenflux1[*,istar],sigma_instr,/edge_truncate)

             vbroadflux = interpol(reform(eigenflux1[*,istar]),eigenwave,vconstwave)
;            combine1fiber, alog10(eigenwave), reform(eigenflux1[*,istar]), $
;              newloglam=logvconstwave, newflux=vbroadflux

             vbroadflux = gconvolve(vbroadflux,sigma_vdisp,/edge_truncate)
             starfluxvconst1[*,istar] = vbroadflux

          endfor

; fit the continuum

          t0 = systime(1)
          backfit1 = ibackfit(restfluxvconst,vconstwave,starfluxvconst1,$
            invvar=restinvvarvconst,mask=restmaskvconst*restmaskbackfitvconst,$
            coeffguess=coeffguess,_extra=extra)
          splog, format='("CPU time for continuum fitting on iteration ",I1,"/",I1,": ",G0," '+$
            'seconds.")', iter+1, niter, systime(1)-t0

;         coeffguess = backfit1.mpfitcoeff ; <-- no longer needed! jm04aug31uofa

;         if keyword_set(force_chisq_one) then begin
;            chisq_one_factor = sqrt(backfit1.continuum_chi2)
;            splog, 'Increasing the errors by factor = '+string(chisq_one_factor,format='(G0.0)')
;         endif
          
; improve the redshift guess by cross-correlating the best-fitting
; continuum spectrum with the emission-line subtracted data spectrum
; but only for high signal-to-noise spectra

          medsnr = backfit1.continuum_snr
;         medsnr = djs_median(sqrt(restinvvarvconst/chisq_one_factor^2.0)*restflux)

          if keyword_set(zcrosscor) and (iter lt niter-1L) then begin

             if (medsnr gt zsnrcross) then begin

                mkhdr, hdr, restfluxvconst
                sxaddpar, hdr, 'COEFF0', coeff0
                sxaddpar, hdr, 'COEFF1', coeff1

                if keyword_set(debug) and (iZ eq 0L) then window, 1, xsize=500, ysize=400
                zupdate = zfind(restfluxvconst,restinvvarvconst*restmaskvconst,hdr=hdr,starloglam0=coeff0,$
                  starflux=reform(backfit1.log_continuum,lognpix,1),zmin=-vmaxshift/light,$
                  zmax=+vmaxshift/light,doplot=0)

                if (zupdate.z_err gt 0.0) then begin
                   zabs1 = float(ziter + zupdate.z) ; absorption-line redshift and error
                   zabs1_err = float(zupdate.z_err)
                   splog, 'ZCROSSCOR converged: z = '+string(ziter,format='(F0.0)')+' --> '+$
                     string(zabs1,format='(F0.0)')
                   ziter = float(zabs1)
                endif else begin
                   splog, 'ZCROSSCOR did not converge.'
                endelse

;               zupdate = im_ztweak(restflux,restwave,backfit1.continuum,restwave,$
;                 specivar=mask*restinvvar,vmin=-vmaxshift,vmax=+vmaxshift,$
;                 minwave=vminwave,maxwave=vmaxwave,/silent,doplot=keyword_set(debug))
;               if (zupdate.errflag eq 0L) then begin
;                  zabs1 = float(ziter + zupdate.zshift) ; absorption-line redshift and error
;                  zabs1_err = float(zupdate.zshift_err)
;                  splog, 'ZCROSSCOR converged: z = '+string(ziter,format='(F0.0)')+' --> '+$
;                    string(zabs1,format='(F0.0)')
;                  ziter = float(zabs1)
;               endif else splog, 'ZCROSSCOR did not converge.'

             endif else splog, 'ZCROSSCOR - Continuum S/N < '+string(zsnrcross,format='(I0)')

          endif

       endfor ; close iteration loop

       if (iZ eq 0L) then begin
          
          starfluxvconst = starfluxvconst1
          zabs = zabs1
          zabs_err = zabs1_err
          backfit = backfit1
          
       endif else begin

          starfluxvconst = [ [ [starfluxvconst] ], [ [starfluxvconst1] ] ]
          zabs = [zabs,zabs1]
          zabs_err = [zabs_err,zabs1_err]
          backfit = [ [backfit], [backfit1] ]

       endelse

;      niter = 1L ; only iterate once for subsequent metallicity templates

;      delvarx, coeffguess
       
    endfor  ; close metallicity loop

    if keyword_set(Zmulti) then begin

       backfit = reform(backfit)       
       
       minchi2 = min(backfit.continuum_chi2,Zbestindx)
       splog, 'Selecting best-fitting metallicity template '+string(Zbestindx+1,format='(I0)')+'.'
              
       backfit = backfit[Zbestindx]
       starfluxvconst = reform(starfluxvconst[*,*,Zbestindx])
       zabs = zabs[Zbestindx]
       zabs_err = zabs_err[Zbestindx]

    endif

; finally resample the best-fitting continuum and continuum templates
; onto the wavelength vector with constant pixel size in Angstrom
; rather than km/s

    restwave = wave / (1.0+zabs)
    
    newcontinuum = interpol(backfit.log_continuum,vconstwave,restwave)
;   combine1fiber, logvconstwave, backfit.log_continuum, $
;     newloglam=logrestwave, newflux=newcontinuum
    
    backfit = struct_addtags({continuum: newcontinuum, medcontinuum: newcontinuum*0.0},$
      struct_trimtags(backfit,except='LOG_CONTINUUM'))

    starflux = restflux # replicate(1.0,nstar)
    for istar = 0L, nstar-1L do begin

       newstarflux = interpol(reform(starfluxvconst[*,istar]),vconstwave,restwave)
;      combine1fiber, logvconstwave, reform(starfluxvconst[*,istar]), $
;        newloglam=logrestwave, newflux=newstarflux
       starflux[*,istar] = newstarflux

    endfor
    
; ---------------------------------------------------------------------------    
; median-smooth the continuum subtracted spectrum in preparation for
; fitting the emission lines
; ---------------------------------------------------------------------------    
       
    espectrum = restflux - backfit.continuum 
    espectrum_invvar = restinvvar

    if keyword_set(no_medcontinuum) then begin

       medcontinuum = backfit.continuum*0.0

    endif else begin
       
       medcontinuum = smooth(medsmooth(espectrum,medsmooth_window),$ ; median correction
         smooth_window,/nan,/edge_truncate) 
;      medcontinuum = smooth(im_medsmooth(espectrum,250,$$
;        weight=espectrum_invvar),100,/nan,/edge_truncate)

;      moremask = telluric_mask(wave)
;      medcontinuum = medcontinuum*moremask

    endelse
    
    espectrum = espectrum-medcontinuum
    backfit.medcontinuum = medcontinuum

    normconst = max(espectrum)
    espectrum = espectrum/normconst
    espectrum_invvar = espectrum_invvar*normconst^2.0

; ---------------------------------------------------------------------------    
; fit the emission-line spectrum
; ---------------------------------------------------------------------------    

    if (n_elements(zindex) eq nline) then zindex1 = zindex else zindex1 = linepars.zindex
    if (n_elements(windex) eq nline) then windex1 = windex else windex1 = linepars.windex
    if (n_elements(findex) eq nline) then findex1 = findex else findex1 = linepars.findex
    if (n_elements(fvalue) eq nline) then fvalue1 = fvalue else fvalue1 = linepars.fvalue
;   struct_print, linepars

; use 70 km/s as an initial guess for the intrinsic Gaussian line
; width for all the emission lines (the result is not very sensitive
; to this guess, but helps speed up convergence)

    sigguess = replicate(70.0,nline) / light / alog(10.0) ; log-Angstrom units

; convert the Gaussian sigma line width of the emission lines due to
; instrumental broadening to km/s

    lineres = replicate(restspecres,nline)                ; rest-frame FWHM resolution [Angstrom]
;   lineres = interpol(specres,wave,linepars.wave)
    elineres = light * lineres / linepars.wave / fwhm2sig ; rest-frame sigma width [km/s]

; to speed up the line fitting, only pass pixels that are near
; emission lines (within 2.35*10000 km/s for QSO emission lines, or
; within 2.35*1000 km/s for other emission lines); thanks to C.
; Tremonti for this idea; jm07sep16nyu: drop emission lines whose
; central wavelengths isn't at least 2*LINERES within either edge of
; the spectrum, which if you don't do sometimes wreaks havoc with
; MPFIT

    igood = where(espectrum_invvar gt 0.0)
    for iline = 0L, nline-1L do begin
       if (strmatch(windex[iline],'*qso*',/fold)) then $
         wavewindow = fwhm2sig*linepars[iline].wave*1E4/light else $ ; QSO line
         wavewindow = fwhm2sig*linepars[iline].wave*1E3/light        ; other emission line
       if (linepars[iline].wave gt (min(restwave[igood])+2*lineres[iline])) and $
         (linepars[iline].wave lt (max(restwave[igood])-2*lineres[iline])) then begin
          eindx1 = where((restwave gt linepars[iline].wave-wavewindow) and $
            (restwave lt linepars[iline].wave+wavewindow),neindx1)
       endif else neindx1 = 0L
       if (neindx1 ne 0L) then if (n_elements(eindx) eq 0L) then $
         eindx = eindx1 else eindx = cmset_op(eindx,'OR',eindx1)
    endfor

    if (n_elements(eindx) eq 0L) then eindx = lindgen(npix) ; jm04nov30uofa

    t0 = systime(1)
    splog, 'Fitting the emission-line spectrum.'
;   eindx = eindx[30:n_elements(eindx)-1L]
    linefit = ilinefit(espectrum[eindx],restwave[eindx],linepars.wave,elineres,$
      invvar=espectrum_invvar[eindx],linename=linepars.line,zindex=zindex1,windex=windex1,$
      findex=findex1,fvalue=fvalue1,zguess=1E-6,vmaxshift=vmaxshift,sigguess=sigguess,$
      specfit=specfit,linefit_chi2=linefit_chi2,linefit_dof=linefit_dof,$
      linefit_niter=linefit_niter,linefit_status=linefit_status,_extra=extra)
    speclinefit = espectrum*0.0
    speclinefit[eindx] = specfit

;   niceprint, linefit.linesigma, linefit.linesigma_instr
;   print, linefit[4].linearea/linefit[3].linearea, get_ebv(linefit[4].linearea/linefit[3].linearea,/hbeta)
;   djs_plot, restwave, espectrum, ps=10, xsty=3, ysty=3, xr=[4800,5100]
;   djs_oplot, restwave, speclinefit, ps=10, color='yellow'    
    splog, 'CPU time for line fitting '+string(systime(1)-t0,format='(G0.0)')+' seconds.'

    goodlines = where((linefit.linearea_err gt 0.0),ngoodlines,comp=badlines,ncomp=nbadlines)
    if (nbadlines ne 0L) then linefit[badlines].linez = zabs

    if (ngoodlines ne 0L) then begin
       linefit[goodlines].linez = linefit[goodlines].linez + zabs
       restlinewave = wave / (1D0+djs_median(linefit[goodlines].linez) )
    endif else begin
       linefit.linez = zabs
       restlinewave = restwave
    endelse
    
; put back the normalization constant into the emission-line fits and
; the emission-line fluxes; multiply the continuum coefficients by the
; normalization constant

    backfit = struct_addtags(struct_trimtags(backfit,except='mpfitcoeff'),$
      {z_abs: zabs, z_abs_err: zabs_err, linefit_chi2: linefit_chi2, $
      linefit_dof: linefit_dof, linefit_niter: linefit_niter, $
      linefit_status: linefit_status})

    espectrum = espectrum*normconst
    espectrum_invvar = espectrum_invvar/normconst^2.0

    speclinefit = speclinefit*normconst ; rest [erg/s/cm2/A]
    
    goodfit = where((linefit.linearea_err gt 0.0),ngoodfit)
    if ngoodfit ne 0L then begin
       linefit[goodfit].linearea = linefit[goodfit].linearea*normconst         ; [erg/s/cm2]
       linefit[goodfit].linearea_err = linefit[goodfit].linearea_err*normconst ; [erg/s/cm2]
    endif

; box fluxes

    goodfit = where((linefit.linebox_err gt 0.0),ngoodfit)
    if ngoodfit ne 0L then begin
       linefit[goodfit].linebox = linefit[goodfit].linebox*normconst         ; [erg/s/cm2]
       linefit[goodfit].linebox_err = linefit[goodfit].linebox_err*normconst ; [erg/s/cm2]
    endif

; combine the [O II] doublet (not generalized!); go back to the
; spectrum and compute the LINEBOX measurement properly

    ilinefit = icombine_blends(linefit)

    oii_1 = where(strmatch(linefit.linename,'*OII_3727_1*',/fold),noii_1)
    oii_2 = where(strmatch(linefit.linename,'*OII_3727_2*',/fold),noii_2)
    if (noii_1 eq 1L) and (noii_2 eq 1L) then linefit = [linefit[oii_1],linefit[oii_2],ilinefit]

    oii = where(strmatch(linefit.linename,'OII_3727',/fold),noii)
    if (noii ne 0L) then if (linefit[oii].linearea_err ne -2.0) then begin

       dlam = restwave[1]-restwave[0]
       sigwidth = linefit[oii].linesigma_total/light*linefit[oii].linewave
       indx = where((restwave ge linefit[oii].linewave - 3.0*sigwidth) and $ ; +/- 3-sigma
         (restwave le linefit[oii].linewave + 3.0*sigwidth),nindx)
       
       lineflux = espectrum[indx]
       lineivar = espectrum_invvar[indx]*tellmask[indx]
       denom = total(lineivar+(lineivar eq 0.0))
       
       linefit[oii].linebox = dlam*nindx*total(lineivar*lineflux)/denom
       linefit[oii].linebox_err = dlam*float(nindx) / sqrt(denom)
       
;      print, linefit[oii].linearea, linefit[oii].linebox, total(lineflux)*dlam, $
;        linefit[oii].linebox_err, linefit[oii].linearea_err

    endif
    
; ---------------------------------------------------------------------------    
; measure the local continuum
; ---------------------------------------------------------------------------    

; compute upper limits and measure the local continuum for all the
; lines in range; this way we can replace the measured line-fluxed and
; EWs with the upper limits in post-analysis as necessary; also
; construct the EWs

    ulimit = where(linefit.linearea_err gt -2.0,nulimit)
    if (nulimit ne 0L) then begin
       ulinefit = linefit[ulimit]
       ulinefit = iupper_limits(restwave,restflux-speclinefit,ulinefit,$
         snrcut=1.0,glosigma=5.0,ghisigma=20.0,telluric=1,$
         debug=debug)
       linefit[ulimit].linecontlevel     = ulinefit.linecontlevel
       linefit[ulimit].linecontlevel_err = ulinefit.linecontlevel_err
       linefit[ulimit].linelimit         = ulinefit.linearea
    endif

    linefit = fill_ew(linefit)

; ---------------------------------------------------------------------------    
; measure rest-frame spectral indices off the raw data, the
; emission-line subtracted spectrum, and the best-fitting continuum
; spectrum; do not include the error in the continuum
; ---------------------------------------------------------------------------    

    splog, 'Measuring spectral indices.'
    indices = spectral_indices(restwave,restflux-speclinefit,ferr=restferr,debug=debug,/silent)
    raw_indices = spectral_indices(restwave,restflux,ferr=restferr,debug=debug,/silent)
    model_indices = spectral_indices(restwave,backfit.continuum,ferr=restferr,/silent)

    newindices = [indices.indices,strtrim(model_indices.indices,2)+'_model',$
      strtrim(raw_indices.indices,2)+'_raw']
    
    indices = struct_addtags({indices: newindices},struct_addtags(struct_trimtags(indices,$
      except='INDICES'),struct_trimtags(im_struct_trimtags(raw_indices,select=$
      tag_names(raw_indices),newtags=tag_names(raw_indices)+'_raw'),except='INDICES_RAW')))
    indices = struct_addtags({indices: newindices},struct_addtags(struct_trimtags(indices,$
      except='INDICES'),struct_trimtags(im_struct_trimtags(model_indices,select=$
      tag_names(model_indices),newtags=tag_names(model_indices)+'_model'),except='INDICES_MODEL')))

; ---------------------------------------------------------------------------    
; measure the Balmer absorption lines from the best-fitting continuum
; spectrum.  identify all the measured Balmer emission lines and
; compute their Gaussian sigma widths [Angstrom] based on their
; measured emission-line widths.  we do not include any extra error
; from the continuum fit; you must be measuring the emission lines to
; get the absorption line measurement (because we need the line-sigma)
; ---------------------------------------------------------------------------    

    allbalmerindx = balmerindx(linefit.linewave,/allbalmer)
    if (allbalmerindx[0] eq -1L) then begin

       balmerabs = ibalmerabs(restwave,backfit.continuum,/return_balmerabs) ; initialize this data structure          

    endif else begin

       balmerwaves = linefit[allbalmerindx].linewave ; [Angstrom]
;      balmersigma = linefit[allbalmerindx].linesigma*balmerwaves/light
       balmersigma = linefit[allbalmerindx].linesigma_total*balmerwaves/light ; [Angstrom]
       
       splog, 'Measuring Balmer absorption line equivalent widths for ['+$
         strjoin(linefit[allbalmerindx].linename,', ')+'].'

       balmerabs = ibalmerabs(restwave,backfit.continuum,ferr=restferr,$
         balmerwaves=balmerwaves,balmersigma=balmersigma,gsigma=2.0,$
         debug=debug,/silent)

    endelse
    
; ---------------------------------------------------------------------------
; generate a debugging plot of the continuum and emission-line fits
; ---------------------------------------------------------------------------

    if keyword_set(debug) then begin

       if !d.window ne 4L then window, 4, xs=600, ys=400
       scale = 1E15
       djs_plot, restwave, scale*restflux, xsty=3, ysty=3, ps=10, charsize=1.5, $
         charthick=2.0, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
         ytitle='Normalized Flux', position=[0.13,0.45,0.96,0.96], $
         yrange=scale*minmax(restflux-speclinefit)*[1.0,1.5]
       djs_oplot, restwave, backfit.continuum+speclinefit, ps=10, color='red'
       djs_oplot, restwave, scale*backfit.continuum, ps=10, color='purple', thick=2.0
       legend, ['Data','Best Fit','Continuum','Residuals'], charsize=1.5, $
         charthick=2.0, color=djs_icolor(['default','red','purple','cyan']), $
         /right, /top, box=0, thick=2.0, line=[0,0,0,0]

       residuals = 100*(restflux-(backfit.continuum+speclinefit))/restflux
       yrange = (djs_median(residuals)+djsig(residuals,sigrej=3.0)*[-5,+5])
       djs_plot, wave, residuals, ps=10, color='cyan', $
         position=[0.13,0.12,0.96,0.45], ytitle='Residuals [%]', $
         xtitle='Rest Wavelength [\AA]', $
         /noerase, xsty=3, ysty=3, charsize=1.5, charthick=2.0, $
         xthick=2.0, ythick=2.0, yrange=yrange
       djs_oplot, !x.crange, [0,0], line=0, thick=2.0
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
       
    endif

;   struct_print, struct_trimtags(linefit,select=['*name','*area*','*area_ew*','*sigma*','*linez*'])

return
end
