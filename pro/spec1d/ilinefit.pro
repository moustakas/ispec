;+
; NAME:
;   ILINEFIT()
;
; PURPOSE:
;   Fit an emission-line spectrum.
;
; INPUTS:
;   wave     - rest wavelength array (A) [NPIX]   
;   flux     - rest emission-line spectrum (erg/s/cm2/A) [NPIX]
;   invvar   - rest inverse variance spectrum 1/(erg/s/cm2/A)^2 [NPIX] 
;   linepars - emission-line parameter structure (see ISPECLINEFIT and
;     READ_LINEPARS) [NLINE]
;
; OPTIONAL INPUTS:
;   instr_lineres - Gaussian sigma instrumental resolution line width
;                   (default 0.0 km/s) [NLINE] 
;   zguess     - redshift guess (default 0); can be a scalar or a vector
;                with one entry per emission line
;   sigguess   - initial guess for *intrinsic* sigmas of all lines
;                can be a scalar or a vector with one entry per
;                line (default 70 km/s)
;   sigmin     - minimum intrinsic line-width (default 1.0 km/s)
;   sigmax     - maximum intrinsic line-width (default 500 km/s)
;   qso_sigmax - maximum intrinsic line-width for emission lines
;                identified as QSO lines (default 10000 km/s)
;   vlinemaxtweak - maximum shift in the emission-line center
;                   (default 1000 km/s) 
;   zindex   - override the values of ZINDEX in LINEPARS [NLINE] 
;   windex   - override the values of WINDEX in LINEPARS [NLINE] 
;   findex   - override the values of FINDEX in LINEPARS [NLINE] 
;   fvalue   - override the values of FVALUE in LINEPARS [NLINE] 
;
; KEYWORD PARAMETERS:
;   positive_emission - restrict the emission-line amplitudes to be
;     positive 
;
; OUTPUTS:
;   linefit - output structure with result of line fits [NLINE]  
;      LINENAME      - string name of line copied from input 
;              parameter by the same name, or else 
;              constructed from the rounded-down 
;              wavelength of each line in air  
;      LINEWAVE      - rest-frame wavelength [Angstrom] (copy
;              of LINEWAVE) 
;      LINEZ         - emission-line redshift
;      LINEZ_ERR     - error in LINEZ 
;      LINESIGMA     - intrinsic Gaussian sigma width [km/s] 
;      LINESIGMA_ERR     - error in LINESIGMA 
;      LINESIGMA_INSTR   - instrumental Gaussian sigma width
;              [km/s] 
;      LINESIGMA_TOTAL   - instrumental plus intrinsic Gaussian 
;              sigma width [km/s] 
;      LINEAREA      - Gaussian area [erg/s/cm2] 
;      LINEAREA_ERR      - error in LINEAREA 
;      LINEBOX       - total flux within +/-3 sigma of the
;              line center; box-car flux estimate; the
;              dispersion must be constant [erg/s/cm2]   
;      LINEBOX_ERR   - error in LINEBOX
;      LINECONTLEVEL     - continuum level at line center
;              [erg/s/cm2/A]; if the line center is  
;              outside the wavelength range, then
;              return the nearest value (either the
;              first or last value) 
;      LINECONTLEVEL_ERR - error in LINECONTLEVEL, or -1L if the
;              line center is outside the wavelength
;              range 
;      LINEEW_AREA   - emission-line equivalent width based on
;              LINEAREA [Angstroms]
;      LINEEW_AREA_ERR   - error in LINEEW_AREA
;      LINEEW_BOX    - emission-line equivalent width based on
;              LINEBOX [Angstroms]
;      LINEEW_BOX_ERR    - error in LINEEW_BOX
;      LINENPIX      - number of pixels within +/- 3 sigma of
;              the line center that have INVVAR > 0.
;      LINEDOF       - LINENPIX minus the number of terms fit
;              for that line, which could be fractional
;              (if one parameter is fixed between N lines,
;              then we say only 1/N-th of that parameter
;              is fit in each of those lines).  This
;              quantity can be zero or negative. 
;      LINECHI2      - chi^2 for all points within +/- 3 sigma
;              of the line center; -1L if no such
;              points 
;
; OPTIONAL OUTPUTS:
;   speclinefit - fitted spectrum including lines and background terms
;         [NPIX] 
;   bfit    - fitted spectrum including only background terms
;         [NPIX] 
;   bterms  - coefficients for background terms [NBACK] 
;
; COMMENTS:
;   This routine fits multiple emission lines assuming Gaussian 
;   line-profiles with arbitrary constrains among lines regarding 
;   their redshift, intrinsic velocity width, or amplitude or 
;   total flux.  The fitting is done in log-lambda units according 
;   to some really cool math that David Schlegel worked out. 
;   Optionally, multiple background terms (e.g., some linear 
;   combination of polynomials) can also be fitted.  All 
;   wavelengths should be in Angstroms.  If a line was dropped 
;   from the fit (for example, if there were no points to fit), 
;   then LINEAREA=0 and LINEAREA_ERR=-1L.  Also, if LINENPIX=0 for 
;   a line, then that line is also removed from the fit.  In 
;   practice, SIGMIN is not used because  occasionally MPFIT will
;   unexplicably drop a line from the fit  if one of the parameter
;   boundaries is set to zero.  
;
;   Possible bug:  Do not use lines with no points to fit in the 
;   computation of degrees of freedom for other lines. 
;
;   LINECONTLEVEL has been deprecated in this version of 
;   ILINEFIT()  in favor of the algorithm implemented in 
;   IFITSPEC().  Consequently, equivalent widths are not
;   computed. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   05-Feb-2002  Written by D. Schlegel, Princeton
;   J. Moustakas, 2002 March 15, U of A, incorporated into
;      iSPEC1d; major changes implemented 
;   jm04jan06uofa - added LINESIGMA_TOTAL
;   jm04apr28uofa - added VMAXSHIFT keyword
;   jm05jan12uofa - code streamlined, documentation cleaned up 
;   jm06feb16uofa - added QSO_SIGMAX optional input
;   jm07mar04nyu  - compute the weighted LINEBOX
;   jm08oct27nyu  - iterate the line-fitting twice
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

function ilinefit, wave, flux1, invvar1, linepars, instr_lineres=instr_lineres, $
  zguess=zguess1, sigguess=sigguess1, ampguess=ampguess1, sigmin=sigmin, $
  sigmax=sigmax, vlinemaxtweak=vlinemaxtweak, zindex=zindex, windex=windex, $
  findex=findex, fvalue=fvalue, speclinefit=speclinefit, linefit_chi2=linefit_chi2, $
  linefit_niter=linefit_niter, linefit_status=linefit_status, $
  positive_emission=positive_emission, silent=silent

    light = 2.99792458D5 ; speed of light [km/s]

    npix = n_elements(wave)
    nflux = n_elements(flux1)
    ninvvar = n_elements(invvar1)
    nline = n_elements(linepars)

    if (npix eq 0L) or (nflux eq 0L) or (ninvvar eq 0L) or $
      (nline eq 0L) then begin
       doc_library, 'ilinefit'
       return, -1
    endif
    
    if (nflux ne npix) or (ninvvar ne npix) then begin
       splog, 'Dimensions of WAVE, FLUX, and INVVAR must match'
       return, -1L
    endif
    
    if (n_elements(instr_lineres) eq 0L) then instr_lineres = fltarr(nline)
    if (n_elements(instr_lineres) ne nline) then begin
       splog, 'Dimensions of INSTR_LINERES and LINEPARS must match'
       return, -1L
    endif

    if (n_elements(sigmin) eq 0L) then sigmin =   1.0 ; [km/s]
    if (n_elements(sigmax) eq 0L) then sigmax = 500.0 ; [km/s]
    if (n_elements(vlinemaxtweak) eq 0L) then vlinemaxtweak = 500.0 ; [km/s]
    
    if (n_elements(zindex) eq 0L) then zindex = linepars.zindex
    if (n_elements(windex) eq 0L) then windex = linepars.windex
    if (n_elements(findex) eq 0L) then findex = linepars.findex
    if (n_elements(fvalue) eq 0L) then fvalue = linepars.fvalue

    if (n_elements(zguess1) ne 0L) then begin
       if (n_elements(zguess1) eq 1L) then begin
          zguess = replicate(zguess1[0], nline)
       endif else if (n_elements(zguess1) eq nline) then begin
          zguess = zguess1
       endif else begin
          splog, 'Dimensions of ZGUESS and LINEPARS must match'
          return, -1L
       endelse
    endif else zguess = replicate(0.0,nline)

    if (n_elements(sigguess1) ne 0L) then begin
       if (n_elements(sigguess1) eq 1L) then begin
          sigguess = replicate(sigguess1[0],nline)
       endif else if (n_elements(sigguess1) eq nline) then begin
          sigguess = sigguess1
       endif else begin
          splog, 'Dimensions of SIGGUESS and LINEPARS must match'
          return, -1L
       endelse
    endif else sigguess = replicate(70.0,nline) ; [km/s]
    sigguess_loglam = sigguess/light/alog(10.0) ; log-lambda [Angstrom]
    
    if (n_elements(ampguess1) ne 0L) then begin
       if (n_elements(ampguess1) eq 1L) then begin
          ampguess = replicate(ampguess1[0],nline)
       endif else if (n_elements(ampguess1) eq nline) then begin
          ampguess = ampguess1
       endif else begin
          splog, 'Dimensions of AMPGUESS and LINEPARS must match'
          return, -1L
       endelse
    endif

; normalize the input spectrum

    fluxnorm = abs(djs_median(flux1[where(invvar1 gt 0.0)])) ; note absolute value!
    flux = flux1/fluxnorm
    invvar = invvar1*fluxnorm^2.0
    
; initialize the output LINEFIT structure

    linefit = icreate_linefit(nline)
    linefit.linename = linepars.line
    linefit.linewave = linepars.wave ; rest frame [Angstrom]

; initialize the MPFIT parameter structure

    parinfo = {$
      parname:         '', $
      value:         0.0D, $
      fixed:            0, $
      limited:      [0,0], $
      tied:            '', $
      limits:   [0.0D,0.0D]}
    parinfo = replicate(parinfo,3*nline)
    parinfo[lindgen(nline)*3].parname = linepars.line

; set the initial guesses of the fitting parameters

    linewave = linepars.wave
    ampindex = 0+lindgen(nline)*3
    lamindex = 1+lindgen(nline)*3
    sigindex = 2+lindgen(nline)*3

; amplitude; don't take out the sqrt(2*!pi)*sigma factor (it has to
; match what's in IMULTIGAUSS), otherwise DOOM!

    if (n_elements(ampguess) eq 0L) then begin
       linterp, wave, smooth(flux,10), linewave*(1.0+zguess), $
         ampguess, missing=1E-2
    endif
;   parinfo[ampindex].value = fvalue
    parinfo[ampindex].value = 1.0
;   parinfo[ampindex].value = (sqrt(2.0*!pi)*sigguess_loglam*ampguess) > 1E-2 ; note!

    if keyword_set(positive_emission) then $
      parinfo[ampindex].limited[0] = 1 else $
      parinfo[ampindex].limited[0] = 0

; line-center (redshift); do not allow the line to shift more than
; +/-VLINEMAXTWEAK km/s from ZGUESS

    zminshift = zguess-vlinemaxtweak/light
    zmaxshift = zguess+vlinemaxtweak/light
    
    parinfo[lamindex].value = alog10(linewave*(1+zguess)) ; log-lambda [Angstrom]
    parinfo[lamindex].limited = [1,1]
    parinfo[lamindex].limits[0] = alog10(linewave*(1+zminshift)) ; log-lambda [Angstrom]
    parinfo[lamindex].limits[1] = alog10(linewave*(1+zmaxshift)) ; log-lambda [Angstrom]

; line-width; restrict the width in the interval [SIGMIN,SIGMAX]
       
    parinfo[sigindex].value = sigguess_loglam 
    parinfo[sigindex].limited = [1,1]
    parinfo[sigindex].limits = [sigmin,sigmax]/light/alog(10.0) ; log-lambda [Angstrom]
    
; apply constraints to peak flux values

    nfitterms = fltarr(nline)

    allindex = findex[uniq(findex,sort(findex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(findex EQ allindex[iall],ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          fratio = fvalue[ii[jj]] / fvalue[ii[0]]
          parinfo[0+ii[jj]*3].tied = string(fratio,0+ii[0]*3,format='(E12.5," * P[",I0,"]")')
       endfor 
    endfor 

; apply redshifts constraints

    allindex = zindex[uniq(zindex,sort(zindex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(zindex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          lamshift = alog10(linewave[ii[jj]] / linewave[ii[0]])
          parinfo[1+ii[jj]*3].tied = string(lamshift,1+ii[0]*3,format='(E12.5," + P[",I0,"]")')
       endfor 
    endfor 

; apply line-width constraints

    allindex = windex[uniq(windex,sort(windex))]
    for iall=0, n_elements(allindex)-1 do begin
       ii = where(windex EQ allindex[iall], ct)
       nfitterms[ii] = nfitterms[ii] + 1.0/ct
       if (ct gt 1) then for jj=1, ct-1 do begin
          parinfo[2+ii[jj]*3].tied = string(2+ii[0]*3,format='("P[",I0,"]")')
       endfor
    endfor 

;   struct_print, parinfo

; to speed up the line fitting (and to reduce the formal chi2 values
; of the line-fits), only give MPFIT pixels that are near the emission
; lines of interest; also drop emission lines whose central
; wavelengths aren't at least 2*INSTR_LINERES away from the edges of
; the spectrum

    restwave = wave/(1.0+zguess[0])
    emask = wave*0.0

    for iline = 0L, nline-1L do begin
       thissig = sqrt(instr_lineres[iline]^2.0+sigmax^2.0)
       if (linewave[iline] gt (min(restwave)+2*linewave[iline]*thissig/light)) and $
         (linewave[iline] lt (max(restwave)-2*linewave[iline]*thissig/light)) then begin
          emask = emask or (abs(light*(restwave-linewave[iline])/$
            linewave[iline]) lt ((10.0*thissig)>2000.0)<20000.0)
       endif 
    endfor 

;   emask = wave*0.0+1.0
    
; do the fit!

    sigmares = instr_lineres/light/alog(10.0) ; [log-Angstrom]
    functargs = {nline: nline, sigmares: sigmares}

    linefit_status = 0
    speclinefit = fltarr(npix)

    good = where((invvar*emask gt 0.0) and (finite(invvar) eq 1),ngood)

    niterfit = 2L
    for iter = 0L, niterfit-1L do begin 

       if (ngood GT 0) then begin
          t0 = systime(1)
          lfit = mpfitfun('imultigauss',alog10(wave[good])*1.0D,flux[good]*1.0D,$
            weights=invvar[good]*1.0D,parinfo=parinfo,covar=covar,perror=perror,$
            yfit=speclinefit1,functargs=functargs,nfev=nfev,niter=linefit_niter,$
            status=linefit_status,quiet=1,bestnorm=linefit_chi2,dof=linefit_dof,$
            /autoderivative,npegged=npegged)
          if (not keyword_set(silent)) then splog, 'MPFIT nfev=', nfev, $
            ' niter=', linefit_niter, ' status=', $
            linefit_status, ' chi2=', linefit_chi2, ' dof=', linefit_dof
          if (linefit_status eq 5) then splog, 'Warning: Maximum number '+$
            'of iterations reached: ', linefit_niter
          if (linefit_dof le 0.0) then linefit_chi2 = -1.0 else $
            linefit_chi2 = sqrt(linefit_chi2/linefit_dof)
          speclinefit[good] = speclinefit1
       endif

       if (ngood eq 0L) then begin
          splog, 'Too few points to fit ', ngood
          nparam = 3 * nline
          lfit = fltarr(nparam)
          perror = fltarr(nparam)
          covar = fltarr(nparam,nparam)
          linefit_chi2 = -1.0 & linefit_niter = -1 & linefit_status = -1 
          continue
       endif else begin
          if (linefit_status eq 0L) then begin
             if (not keyword_set(silent)) then splog, $
               format='("All lines dropped from fit on '+$
               'iteration ",I1,"/",I1,": ",G0," seconds")', iter+1, $
               niterfit, systime(1)-t0
          endif else begin
             if (not keyword_set(silent)) then splog, $
               format='("CPU time for emission-line fitting on '+$
               'iteration ",I1,"/",I1,": ",G0," seconds")', iter+1, $
               niterfit, systime(1)-t0
          endelse
          parinfo.value = lfit  ; improve the initial guess!
          parinfo[sigindex].tied = ''
          parinfo[lamindex].tied = ''
       endelse 

    endfor ; close the iteration loop
       
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

    linefit.linearea      = lfit[ampindex]*alog(10.0)*10.0^lfit[lamindex]
    linefit.linez         = (10.0^lfit[lamindex]/linewave-1.0)*(lfit[lamindex] gt 0)
    linefit.linesigma     = abs(lfit[sigindex])*alog(10.0)*light

    linefit.linearea_err  = perror[ampindex]*alog(10.0)*10.0^lfit[lamindex]
    linefit.linez_err     = perror[lamindex]*alog(10.0)*(linefit.linez+1.0)
    linefit.linesigma_err = perror[sigindex]*alog(10.0)*light

    linefit.linesigma_instr = instr_lineres
    linefit.linesigma_total = sqrt(linefit.linesigma^2.0 + instr_lineres^2.0)

; flag lines that were dropped from the fit by MPFIT (not sure under
; what circumstances this occurs)

    bad = where((perror[ampindex] le 0),nbad)
    if (nbad ne 0L) then begin
       linefit[bad].linearea = 0.0
       linefit[bad].linearea_err = -1.0
       linefit[bad].linesigma = 0.0
       linefit[bad].linesigma_err = -1.0
    endif

; compute chi2    
    
    for iline=0, nline-1 do begin
       sigwidth = sqrt(sigmares[iline]^2.0 + lfit[iline*3+2]^2.0) ; total line-width
       indx = where((alog10(wave) ge lfit[iline*3+1] - 3.0*sigwidth) and $ ; +/- 3-sigma
         (alog10(wave) le lfit[iline*3+1] + 3.0*sigwidth),nindx)

       if (nindx ne 0L) then begin
          linefit[iline].linenpix = total(invvar[indx] gt 0)
          linefit[iline].linedof = linefit[iline].linenpix-nfitterms[iline]
          if (linefit[iline].linedof GT 0) then begin
             linefit[iline].linechi2 = total((flux[indx]-speclinefit[indx])^2.0*invvar[indx])
          endif 
       endif

; if the line was not fitted (outside the wavelength range) or is
; right on the edge of the spectrum then reject it using a different
; error code        

;      fix this!!!!!!
       
       if (linefit[iline].linenpix eq 0) or (indx[0] eq 0L) or $
         (indx[n_elements(indx)-1L] eq (npix-1L)) then begin

          linefit[iline].linez             = zguess[iline]
          linefit[iline].linez_err         = -2.0
          linefit[iline].linesigma         =  0.0
          linefit[iline].linesigma_err     = -2.0
          linefit[iline].linearea          =  0.0
          linefit[iline].linearea_err      = -2.0
          linefit[iline].lineew_area       =  0.0
          linefit[iline].lineew_area_err   = -2.0
          linefit[iline].linecontlevel     =  0.0
          linefit[iline].linecontlevel_err = -2.0
          linefit[iline].linenpix          =  0
          linefit[iline].linedof           =  0.0
          linefit[iline].linechi2          = -2.0

       endif

    endfor

; set the amplitudes of reject lines equal to zero for when we
; re-evaluate SPECLINEFIT; for all other lines scale by the
; normalization constant

    keep = where((linefit.linearea_err gt 0.0),nkeep,$
      comp=rej,ncomp=nrej)
    if (nrej ne 0L) then lfit[ampindex[rej]] = 0.0
    if (nkeep ne 0L) then begin
       linefit[keep].linearea = fluxnorm*linefit[keep].linearea
       linefit[keep].linearea_err = fluxnorm*linefit[keep].linearea_err
    endif
    
; re-evaluate SPECLINEFIT and return

    if (arg_present(speclinefit)) then speclinefit = fluxnorm*$
      imultigauss(alog10(wave),lfit,_extra=functargs)

return, linefit
end
