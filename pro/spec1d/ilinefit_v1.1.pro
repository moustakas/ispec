;+
; NAME:
;       ILINEFIT()
;
; PURPOSE:
;       Fit an emission-line spectrum.
;
; INPUTS:
;	flux       - emission-line spectrum (erg/s/cm2/A) [NPIX] 
;       wave       - wavelength array (Angstroms), assumed constant
;                    linear dispersion in Angstrom per pixel [NPIX]   
;       linelambda - rest-frame wavelengths of the emission lines of
;                    interest (Angstroms) [NLINE] 
;       lineres    - Gaussian sigma spectral resolution line width
;                    (km/s) [NLINE] 
;
; OPTIONAL INPUTS:
;       invvar     - inverse variance spectrum [NPIX].  if not
;                    provided then the emission line flux errors will
;                    be incorrect and the fitting may even fail
;       linename   - String name(s) of line(s) to be copied to output
;                    structure.
;       zindex     - lines with the same ZINDEX are constrained to
;                    have the same redshift; default to a unique list
;                    of indices [NLINE]
;       windex     - Lines with the same WINDEX are constrained to
;                    have the same width [NLINE]; default to a unique 
;                    list of indices
;       findex     - Lines with the same FINDEX are constrained to
;                    have the same flux ratio as input in FVALUE;
;                    defult to a unique list of indices [NLINE]
;       fvalue     - If FINDEX is specified, then constrain the lines 
;                    to have the same flux ratios as in FVALUE
;                    [NLINE]; default to all values of 1.  These
;                    values are also used as the initial guesses for 
;                    the line strengths
;       sigmin     - minimum intrinsic emission-line width [km/s]
;                    (default 0); see COMMENTS
;       sigmax     - maximum intrinsic emission-line width [km/s]
;                    (default 500) 
;       qso_sigmax - maximum intrinsic emission-line width for QSO
;                    emission lines (e.g., Mg II 2800), as identified
;                    by string matching on WINDEX (default 10000 km/s) 
;                    5000)  
;       vmaxshift  - maximum shift in the emission-line center [km/s]
;                    (default 1000) - THIS PARAMETER HAS BEEN
;                                     RELEGATED (SEE BELOW)
;       background - background vector(s) to fit simultaneously with
;                    the lines, where the scaling of each vector is
;                    fit [NPIX,NBACK]. The redshift of these
;                    background vectors is not fit, but rather we
;                    maintain the one-to-one correspondence of each
;                    BACKGROUND pixel to each FLUX pixel.  The initial
;                    guess for the scaling of each background vector
;                    is unity.
;       sigguess   - initial guess for *intrinsic* sigmas of all lines
;                    can be a scalar or a vector with one entry per
;                    line [log-Angstrom] (default 100 km/s)
;       zguess     - redshift guess (needs to be pretty close, or set
;                    to zero and pass the rest wavelength vector); can
;                    be a scalar or a vector with one entry per
;                    emission line 
;
; KEYWORD PARAMETERS:
;       positive_emission - restrict the emission lines to be positive
;                           (default is to allow lines to be positive
;                           or negative)
;
; OUTPUTS:
;       linefit - output structure with result of line fits [NLINE]  
;          LINENAME          - string name of line copied from input 
;                              parameter by the same name, or else 
;                              constructed from the rounded-down 
;                              wavelength of each line in air  
;          LINEWAVE          - rest-frame wavelength [Angstrom] (copy
;                              of LINELAMBDA) 
;          LINEZ             - emission-line redshift
;          LINEZ_ERR         - error in LINEZ 
;          LINESIGMA         - intrinsic Gaussian sigma width [km/s] 
;          LINESIGMA_ERR     - error in LINESIGMA 
;          LINESIGMA_INSTR   - instrumental Gaussian sigma width
;                              [km/s] 
;          LINESIGMA_TOTAL   - instrumental plus intrinsic Gaussian 
;                              sigma width [km/s] 
;          LINEAREA          - Gaussian area [erg/s/cm2] 
;          LINEAREA_ERR      - error in LINEAREA 
;          LINEBOX           - total flux within +/-3 sigma of the
;                              line center; box-car flux estimate; the
;                              dispersion must be constant [erg/s/cm2]   
;          LINEBOX_ERR       - error in LINEBOX
;          LINECONTLEVEL     - continuum level at line center
;                              [erg/s/cm2/A]; if the line center is  
;                              outside the wavelength range, then
;                              return the nearest value (either the
;                              first or last value) 
;          LINECONTLEVEL_ERR - error in LINECONTLEVEL, or -1L if the
;                              line center is outside the wavelength
;                              range 
;          LINEEW_AREA       - emission-line equivalent width based on
;                              LINEAREA [Angstroms]
;          LINEEW_AREA_ERR   - error in LINEEW_AREA
;          LINEEW_BOX        - emission-line equivalent width based on
;                              LINEBOX [Angstroms]
;          LINEEW_BOX_ERR    - error in LINEEW_BOX
;          LINENPIX          - number of pixels within +/- 3 sigma of
;                              the line center that have INVVAR > 0.
;          LINEDOF           - LINENPIX minus the number of terms fit
;                              for that line, which could be fractional
;                              (if one parameter is fixed between N lines,
;                              then we say only 1/N-th of that parameter
;                              is fit in each of those lines).  This
;                              quantity can be zero or negative. 
;          LINECHI2          - chi^2 for all points within +/- 3 sigma
;                              of the line center; -1L if no such
;                              points 
;
; OPTIONAL OUTPUTS:
;       specfit - fitted spectrum including lines and background terms
;                 [NPIX] 
;       bfit    - fitted spectrum including only background terms
;                 [NPIX] 
;       bterms  - coefficients for background terms [NBACK] 
;
; COMMENTS:
;       This routine fits multiple emission lines assuming Gaussian 
;       line-profiles with arbitrary constrains among lines regarding 
;       their redshift, intrinsic velocity width, or amplitude or 
;       total flux.  The fitting is done in log-lambda units according 
;       to some really cool math that David Schlegel worked out. 
;       Optionally, multiple background terms (e.g., some linear 
;       combination of polynomials) can also be fitted.  All 
;       wavelengths should be in Angstroms.  If a line was dropped 
;       from the fit (for example, if there were no points to fit), 
;       then LINEAREA=0 and LINEAREA_ERR=-1L.  Also, if LINENPIX=0 for 
;       a line, then that line is also removed from the fit.  In 
;       practice, SIGMIN is not used because  occasionally MPFIT will
;       unexplicably drop a line from the fit  if one of the parameter
;       boundaries is set to zero.  
;
;       Possible bug:  Do not use lines with no points to fit in the 
;       computation of degrees of freedom for other lines. 
;
;       LINECONTLEVEL has been deprecated in this version of 
;       ILINEFIT()  in favor of the algorithm implemented in 
;       IFITSPEC().  Consequently, equivalent widths are not
;       computed. 
;
; EXAMPLES:
;
; PROCEDURES USED:
;       ICREATE_LINEFIT(), MPFITFUN(), MANYGAUSS(), ONEGAUSS(), SPLOG 
;
; MODIFICATION HISTORY:
;       05-Feb-2002  Written by D. Schlegel, Princeton
;       J. Moustakas, 2002 March 15, U of A, incorporated into
;          iSPEC1d; major changes implemented 
;       jm04jan06uofa - added LINESIGMA_TOTAL
;       jm04apr28uofa - added VMAXSHIFT keyword
;       jm05jan12uofa - code streamlined, documentation cleaned up 
;       jm06feb16uofa - added QSO_SIGMAX optional input
;       jm07mar04nyu  - compute the weighted LINEBOX
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

function ilinefit, flux, wave, linelambda, lineres, invvar=invvar, $
  linename=linename, zindex=zindex, windex=windex, findex=findex, $
  fvalue=fvalue, sigmin=sigmin, sigmax=sigmax, qso_sigmax=qso_sigmax, $
  vmaxshift=vmaxshift, background=background, zguess=zguess1, $
  sigguess=sigguess1, specfit=specfit, bfit=bfit, bterms=bterms, $
  linefit_chi2=linefit_chi2, linefit_dof=linefit_dof, linefit_niter=linefit_niter, $
  linefit_status=linefit_status, positive_emission=positive_emission, silent=silent

    light = 2.99792458D5 ; speed of light [km/s]

    npix = n_elements(flux)
    nwave = n_elements(wave)
    nline = n_elements(linelambda)

    if (npix eq 0L) or (nwave eq 0L) or (nline eq 0L) then begin
       doc_library, 'ilinefit'
       return, -1
    endif
    
    if (nwave ne npix) then begin
       splog, 'WAVE and FLUX do not have the same number of elements.'
       return, -1L
    endif
    
    nivar = n_elements(invvar)
    if (nivar eq 0L) then invvar = make_array(npix,value=1.0,/float) else begin
       if nivar ne npix then begin
          splog, 'FLUX and INVVAR do not have the same number of elements.'
          return, -1L
       endif
    endelse

    ndim = size(background,/n_dimension)
    if (ndim eq 0L) then begin
       background = flux*0.0
       nback = 0L
    endif else begin
       dims = size(background,/dimension)
       if (ndim eq 1) then begin
          nback = 1L
          if n_elements(background) ne npix then begin
             splog, 'FLUX and BACKGROUND do not have the same number of elements.'
             return, -1L
          endif
       endif
       if (ndim eq 2) then begin
          nback = dims[1]
          if n_elements(background[*,0]) ne npix then begin
             splog, 'FLUX and BACKGROUND[*,0] do not have the same number of elements.'
             return, -1L
          endif
       endif
    endelse 

; convert the spectral resolution from km/s to log-Angstrom units

    nlineres = n_elements(lineres)
    if (nlineres eq 0L) then begin
       splog, 'LINERES must be provided.'
       return, -1L
    endif

    if (nlineres ne nline) then begin
       splog, 'LINERES must have the same number of elements as LINELAMBDA.'
       return, -1L
    endif

    sigmares = lineres / light / alog(10.0) ; [log-Angstrom]

; set defaults
    
    if (NOT keyword_set(zindex)) then zindex = lindgen(nline)
    if (NOT keyword_set(windex)) then windex = lindgen(nline)
    if (NOT keyword_set(findex)) then findex = lindgen(nline)
    if (NOT keyword_set(fvalue)) then fvalue = fltarr(nline) + 1.0

    if (keyword_set(zguess1)) then begin
       if (n_elements(zguess1) EQ 1) then begin
          zguess = replicate(zguess1[0], nline)
       endif else if (n_elements(zguess1) EQ nline) then begin
          zguess = zguess1
       endif else begin
          splog, 'Wrong number of elements for ZGUESS.'
          return, -1L
       endelse
    endif else begin
       zguess = replicate(0D,nline)
    endelse

;   splog, '###########################################################################'
;   print, zguess
;   splog, '###########################################################################'

    if (keyword_set(sigguess1)) then begin
       if (n_elements(sigguess1) EQ 1) then begin
          sigguess = replicate(sigguess1[0], nline)
       endif else if (n_elements(sigguess1) EQ nline) then begin
          sigguess = sigguess1
       endif else begin
          splog, 'Wrong number of elements for SIGGUESS.'
          return, -1L
       endelse
    endif else begin
       sigguess = replicate(100.0/light/alog(10.0),nline) ; default sigma width [log-Angstrom]
    endelse 
    
    loglam = alog10(wave) ; log-Angstrom

; initialize the LINEFIT output structure

    linefit = icreate_linefit(nline)

    linefit.linewave = linelambda
    if (n_elements(linename) ne nline) then linename = strtrim(string(long(linelambda),format='(I0)'),2)
    linefit.linename = linename

; initialize the structure to be passed to MPFIT 

    parinfo = replicate({parname: '', value:0.0D, fixed: 0L, limited:[0L,0L], $
      tied:'', limits:[0.0D,0.0D]}, nline*3+nback)
    functargs = {nline: nline, nback: nback, loglam: loglam, $
      sigmares: sigmares, background: background}

    parinfo[lindgen(nline)*3].parname = linename
    
; set the initial guesses of the fitting parameters; if we set SIGMIN
; equal to zero then occasionally MPFIT drops the line from the fit;
; therefore set it to a negative number 

    if (n_elements(sigmin) eq 0L) then sigmin = -50.0D ; 0.0D ; constrain the sigma widths [km/s]
    if (n_elements(sigmax) eq 0L) then sigmax = 500.0D
    if (n_elements(qso_sigmax) eq 0L) then qso_sigmax = 10000.0D
    if (n_elements(vmaxshift) eq 0L) then vmaxshift = 1000.0 ; [km/s] ; maximum shift in line center
    
    for iline=0, nline-1 do begin

       parinfo[0+iline*3].value = fvalue[iline]*1D-2 ; the factor of 0.001 make FVALUE a better initial guess
       parinfo[1+iline*3].value = alog10(linelambda[iline]*(1+zguess[iline])) ; restrict the line-width
       parinfo[2+iline*3].value = sigguess[iline]

       if keyword_set(positive_emission) then $ ; restrict the sign of the amplitude
         parinfo[0+iline*3].limited[0] = 1L else $
         parinfo[0+iline*3].limited[0] = 0L

; test code jm06feb26uofa       

;      parinfo[0+iline*3].limited = [0L,1L]
;      parinfo[0+iline*3].limits[1] = 1.0

; constrain the sigma widths between [SIGMIN,SIGMAX] (log-Angstrom);
; jm04mar10uofa: do not constrain SIGMIN (occasionally drops lines
; from the fit if SIGMIN=0 (unknown why)

       parinfo[2+iline*3].limited = [0,1] ; <-- do not constrain SIGMIN
;      parinfo[2+iline*3].limits[0] = sigmin/alog(10.0)/light
       parinfo[2+iline*3].limits[1] = sigmax/alog(10.0)/light
       if (strmatch(windex[iline],'*qso*',/fold)) then parinfo[2+iline*3].limits[1] = qso_sigmax/alog(10.0)/light

; do not allow the line center to move more than +/-VMAXSHIFT km/s!
; THIS PARAMETER IS NOW RELEGATED (jm07sep16nyu); after a horrible
; amount of crap, it appears that strange/non-linear MPFIT behavior
; can appear when this parameter is restricted, especially for lines
; that fall near the edges of the spectrum (e.g., with AGES); I just
; need to be sure to check for Z_LINE values that are significantly
; different from Z_OBJ

       parinfo[1+iline*3].limited = 1B
       parinfo[1+iline*3].limits[0] = alog10(linelambda[iline]*(1.0+zguess[iline])-linelambda[iline]*vmaxshift/light)
       parinfo[1+iline*3].limits[1] = alog10(linelambda[iline]*(1.0+zguess[iline])+linelambda[iline]*vmaxshift/light)
       
    endfor

; set the initial guess for the initial background level to be the
; median flux

    if (nback ne 0L) then parinfo[nline*3+0].value = median(flux)
    
;; set the initial guess for the scaling of each background vector to unity.
;;   for iback=0, nback-1 do begin
;;      parinfo[nline*3+iback].value = 1.0
;;   endfor

; make a list of the number of fitting terms per line.  If a parameter
; is constrained between N lines, then we say each of those lines is
; only fitting 1/N-th of that parameter

    nfitterms = fltarr(nline)

; apply constraints to peak flux values (not integrated flux in a
; line) 

    allindex = findex[uniq(findex,sort(findex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(findex EQ allindex[iall],ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          fratio = fvalue[ii[jj]] / fvalue[ii[0]]
          parinfo[0+ii[jj]*3].tied = string(fratio,0+ii[0]*3,format='(E12.5," * P[",I0,"]")')
;         parinfo[0+ii[jj]*3].tied = string(fratio,0+ii[0]*3,format='(E12.5," * P(",I0,")")')
       endfor 
    endfor 

; apply constraints to couple redshifts

    allindex = zindex[uniq(zindex,sort(zindex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(zindex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          lamshift = alog10(linelambda[ii[jj]] / linelambda[ii[0]])
          parinfo[1+ii[jj]*3].tied = string(lamshift,1+ii[0]*3,format='(E12.5," + P[",I0,"]")')
;         parinfo[1+ii[jj]*3].tied = string(lamshift,1+ii[0]*3,format='(E12.5," + P(",I0,")")')
       endfor 
    endfor 

; apply constraints to couple widths

    allindex = windex[ uniq(windex,sort(windex)) ]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(windex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          parinfo[2+ii[jj]*3].tied = string(2+ii[0]*3,format='("P[",I0,"]")')
;         parinfo[2+ii[jj]*3].tied = string(2+ii[0]*3,format='("P(",I0,")")')
       endfor
    endfor 

;   struct_print, parinfo

; do the fit!

    specfit = fltarr(npix)

    igood = where((invvar gt 0),ngood)
    linefit_status = 0
    if (ngood GT 0) then begin
       lfit = mpfitfun('manygauss',igood,double(flux[igood]),1.0/sqrt(double(invvar[igood])),$
         parinfo=parinfo,covar=covar,perror=perror,yfit=specfit1,functargs=functargs,$
         nfev=nfev,niter=linefit_niter,status=linefit_status,quiet=1,bestnorm=linefit_chi2,$
         dof=linefit_dof,/autoderivative)
       if (not keyword_set(silent)) then splog, 'MPFIT nfev=', nfev, ' niter=', $
         linefit_niter, ' status=', linefit_status, ' chi2=', linefit_chi2, ' dof=', linefit_dof
       if (linefit_status EQ 5) then splog, 'Warning: Maximum number of iterations reached: ', linefit_niter
       specfit[igood] = specfit1
;      dfpsclose
;      djs_plot, wave[igood], flux[igood], ps=10;, xr=[6400,6900]
;      djs_oplot, wave[igood], specfit1, ps=10, color='red'
;;     plot, wave, manygauss(lindgen(npix),parinfo.value,_extra=functargs)    
;      stop
    endif

;   niceprint, string(parinfo.parname,format='(A10)'), string(parinfo.tied,format='(A25)'), lfit, perror
    
    if ((ngood eq 0L) or (linefit_status eq 0L)) then begin
       splog, 'Too few points to fit ', ngood
       nparam = 3 * nline + nback
       lfit = fltarr(nparam)
       perror = fltarr(nparam)
       covar = fltarr(nparam,nparam)
    endif

; convert -0.0 to 0.0

;   negzero = where_negzero(perror,negcount)
;   if negcount ne 0L then perror[negzero] = 0.0
    
; for parameters that are fixed, assign them the same errors as those
; parameters to which they are tied; for the flux area, scale those
; errors according to the ratio of the flux levels

    allindex = findex[uniq(findex,sort(findex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(findex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+0] = perror[ii[0]*3+0]*fvalue[ii]/fvalue[ii[0]]
    endfor

    allindex = zindex[uniq(zindex,sort(zindex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(zindex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+1] = perror[ii[0]*3+1]
    endfor

    allindex = windex[uniq(windex,sort(windex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(windex EQ allindex[iall], ct)
       if (ct GT 1) then perror[ii*3+2] = perror[ii[0]*3+2]
    endfor

; construct the line-measure outputs (and their errors); set the
; redshift to zero if the parameter is zero 

    linefit.linearea      = lfit[lindgen(nline)*3+0]*alog(10)*10.0^lfit[lindgen(nline)*3+1]
    linefit.linez         = (10.0^lfit[lindgen(nline)*3+1]/linelambda-1)*(lfit[lindgen(nline)*3+1] gt 0)
    linefit.linesigma     = abs(lfit[lindgen(nline)*3+2])*alog(10.)*light

    linefit.linearea_err  = perror[lindgen(nline)*3+0]*alog(10)*10.0^lfit[lindgen(nline)*3+1]
    linefit.linez_err     = perror[lindgen(nline)*3+1]*alog(10)*(linefit.linez+1)
    linefit.linesigma_err = perror[lindgen(nline)*3+2]*alog(10)*light

; if a line was dropped from the fit by MPFIT (I'm not sure under what
; circumstances this occurs) then set the LINEAREA to 0 and the
; LINEAREA_ERR to -1L.  we will compute an upper limit on these lines
; in IFITSPEC(); jm07mar03nyu: also set LINESIGMA to zero because
; sometimes a non-zero value is returned

    ibad = where(perror[lindgen(nline)*3+0] LE 0)
    if (ibad[0] NE -1) then begin
       linefit[ibad].linearea = 0.0
       linefit[ibad].linearea_err = -1.0
       linefit[ibad].linesigma = 0.0
       linefit[ibad].linesigma_err = -1.0
    endif

; assign a linewidth of zero to lines whose width is negative
; (unphysically narrow lines); compute upper limits on these lines in
; IFITSPEC(); jm04mar22uofa - actually, do not mess with these lines;
; sometimes MPFIT does not need the intrinsic line width and fits just
; fine with the instrumental width, for some reason, which is
; sufficient 
    
    ibad = where(lfit[lindgen(nline)*3+2] le 0.0)
;   ibad = where(lfit[lindgen(nline)*3+2] le SIGMIN/alog(10.0)/light)
    if (ibad[0] NE -1) then begin
;      lfit[ibad*3+2] = sigmares[ibad]
       linefit[ibad].linesigma = abs(linefit[ibad].linesigma)
       linefit[ibad].linesigma_err = abs(linefit[ibad].linesigma_err)
;      linefit[ibad].linesigma = 0.0
;      linefit[ibad].linesigma_err = -1.0
    endif

; fill the output structure with the instrumental and the total
; (instrumental plus intrinsic) line width

    linefit.linesigma_instr = lineres
    linefit.linesigma_total = sqrt(linefit.linesigma^2.0 + lineres^2.0)

; find the background levels only

    if (nback EQ 0) then begin
       bterms = 0
       bfit = fltarr(npix)
       berr = bfit*0.0
    endif else begin
       bterms = lfit[nline*3:nline*3+nback-1]
       bcovar = covar[nline*3:nline*3+nback-1,nline*3:nline*3+nback-1]

; the following two methods for evaluating bfit are equivalent:

;      bfit = manygauss(lindgen(npix), bterms, nline=0, nback=nback, $
;      loglam=loglam, background=background)
       bfit = bterms ## background

       berr = fltarr(npix)
       for ipix=0, npix-1 do $
         berr[ipix] = sqrt( (transpose(background[ipix,*]) $
         # bcovar # transpose(background[ipix,*])) )
    endelse

; For each line, determine the background level at the line center and
; the number of pixels and chi^2 of each line fit

; DEPRECATED

;;    logmin = min(loglam)
;;    logmax = max(loglam)
;;    for iline=0, nline-1 do begin
;;       if (lfit[iline*3+1] LT logmin) then begin
;;; Case where the line center is blueward of the entire spectrum.
;;          linefit[iline].linecontlevel = 0.0      ; bfit[0]
;;          linefit[iline].linecontlevel_err = -1.0
;;       endif else if (lfit[iline*3+1] GT logmax) then begin
;;; Case where the line center is redward of the entire spectrum.
;;          linefit[iline].linecontlevel = 0.0      ; bfit[npix-1]
;;          linefit[iline].linecontlevel_err = -1.0
;;       endif else begin
;;; Select the nearest pixel for evaluating the background
;;; level at this line center.
;;          junk = min(abs(loglam - lfit[iline*3+1]), ipix)
;;          linefit[iline].linecontlevel = bfit[ipix]
;;          linefit[iline].linecontlevel_err = berr[ipix]
;;       endelse
;;    endfor
    
; find the pixels that are within +/- 3 sigma of the line center; note 
; that if the line is very (unphysically) narrow, it is possible to
; have no lines within this domain; reject those fits;  lines that
; floor against SIGMIN or SIGMAX will have very large errors, so don't 
; flag those cases as special rejection; also compute the box-car area
; within +/-3 sigma of the line-center; assume the dispersion is
; constant 

    dlam = wave[1]-wave[0] ; [Angstrom/pixel]
;   var = 1.0/invvar       ; variance
    
    for iline=0, nline-1 do begin

       if (strmatch(windex[iline],'*qso*',/fold)) then sigfactor = 1.0 else sigfactor = 3.0
       sigwidth = sqrt(sigmares[0]^2.0 + lfit[iline*3+2]^2.0) ; total line width
       
       indx = where(loglam ge lfit[iline*3+1] - sigfactor*sigwidth $ ; +/- 3-sigma
         and loglam le lfit[iline*3+1] + sigfactor*sigwidth, nindx)
;      indx = where(loglam ge lfit[iline*3+1] - 3*lfit[iline*3+2] $ ; +/- 3-sigma
;        and loglam le lfit[iline*3+1] + 3*lfit[iline*3+2] )

       if (indx[0] NE -1) then begin

          linefit[iline].linenpix = total(invvar[indx] GT 0)
          linefit[iline].linedof = linefit[iline].linenpix-nfitterms[iline]

          if (linefit[iline].linedof GT 0) then begin
             linefit[iline].linechi2 = total( (flux[indx]-specfit[indx])^2*invvar[indx])
          endif 

; box-car flux measurement and error (jm02apr18uofa); measure the
; total flux contained within +/- 3-sigma of line center.  this
; measurement is meaningless for blended lines

          denom = total(invvar[indx]+(invvar[indx] eq 0.0))

          linefit[iline].linebox = dlam*nindx*total(invvar[indx]*(flux[indx]-bfit[indx]))/denom ; [erg/s/cm2]
          linefit[iline].linebox_err = dlam*float(nindx) / sqrt(denom)                          ; [erg/s/cm2]

;         linefit[iline].linebox = total(flux[indx]-bfit[indx])*dlam                ; [erg/s/cm2]
;         linefit[iline].linebox_err = sqrt(total(var[indx] + berr[indx]^2.0))*dlam ; [erg/s/cm2]

; select the nearest pixel for evaluating the background
; level at this line center.  DEPRECATED

;         junk = min(abs(loglam-lfit[iline*3+1]),ipix)
;         linefit[iline].linecontlevel = bfit[ipix]
;         linefit[iline].linecontlevel_err = berr[ipix]

       endif 

; special-case rejection.  set LINEAREA=0L, LINEAREA_ERR=-2L if there
; are no data points within the line-fitting region.  also reject if
; the edges of the line encroach on the endpoints of the wavelength
; range.  also set the fitted redshift to the redshift guess
; (jm03mar10uofa) 

       if (linefit[iline].linenpix EQ 0) or (indx[0] eq 0L) or $
         (indx[n_elements(indx)-1L] eq (npix-1L)) then begin

;         struct_print, struct_trimtags(linefit,select=['*name','*area*','*area_ew*','*sigma*','*linez*','*npix*'])
          
          linefit[iline].linez             = zguess[iline]
          linefit[iline].linez_err         = -2.0
          linefit[iline].linesigma         =  0.0
          linefit[iline].linesigma_err     = -2.0
;         linefit[iline].linesigma_instr   =  0.0
;         linefit[iline].linesigma_total   =  0.0
          linefit[iline].linearea          =  0.0
          linefit[iline].linearea_err      = -2.0
          linefit[iline].linebox           =  0.0
          linefit[iline].linebox_err       = -2.0
          linefit[iline].lineew_area       =  0.0
          linefit[iline].lineew_area_err   = -2.0
          linefit[iline].lineew_box        =  0.0
          linefit[iline].lineew_box_err    = -2.0
          linefit[iline].linecontlevel     =  0.0
          linefit[iline].linecontlevel_err = -2.0
          linefit[iline].linenpix          =  0
          linefit[iline].linedof           =  0.0
          linefit[iline].linechi2          = -2.0

; set these line-fit coefficients equal to zero for when we
; re-evaluate SPECFIT

          lfit[iline*3+0] = 0
          
       endif

    endfor

; re-evaluate such that we get the functional fit at the rejected
; wavelengths too.  This also re-evaluates the fit for lines that have
; been rejected and removed due to no valid data points within the
; line-fitting region. 

    if (arg_present(specfit)) then specfit = manygauss(lindgen(npix),lfit,_extra=functargs)

return, linefit
end
