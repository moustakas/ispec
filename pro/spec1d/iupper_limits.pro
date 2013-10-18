;+
; NAME:
;       IUPPER_LIMITS()
;
; PURPOSE:
;       Compute emission-line upper limits.
;
; CALLING SEQUENCE:
;       linefit = iupper_limits(wave,flux,linefit,snrcut=,$
;          glosigma=,ghisigma,/telluric,/debug)
;
; INPUTS:
;       wave     - rest wavelength vector [A]
;       flux     - emission-line subtracted data spectrum [erg/s/cm2/A]
;       linefit  - linefit structure from IFITSPEC() for lines
;                  requiring upper limits
;
; OPTIONAL INPUTS:
;       snrcut     - S/N cut applied to the emission lines (default
;                    1.0)  
;       glosigma   - compute the noise properties of the continuum in
;                    the interval +/- [GLOSIGMA,GHISIGMA] * LINESIGMA
;                    Angstrom on either side of line center, where
;                    LINESIGMA is the Gaussian width of the line, but
;                    see COMMENTS (default 5.0)
;       ghisigma   - see GLOSIGMA (default 15.0)
;
; KEYWORD PARAMETERS:
;       telluric - mask telluric absorption
;       debug    - generate a plot for debugging
;
; OUTPUTS:
;       linefit  - the following fields are modified and updated 
;          LINESIGMA_ERR
;          LINEAREA
;          LINEAREA_ERR
;          LINEBOX
;          LINEBOX_ERR
;          LINECONTLEVEL
;          LINECONTLEVEL_ERR
;          LINECHI2
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Upper limits are computed by measuring the statistical
;       variance in the continuum surrounding the emission line.
;       Forbidden and Balmer emission lines are treated differently.
;       The forbidden line continuum is defined as the interval
;       [LINEWAVE-GHISIGMA*LINESIGMA:LINEWAVE-GLOSIGMA*LINESIGMA] and
;       [LINEWAVE+GHISIGMA*LINESIGMA:LINEWAVE+GLOSIGMA*LINESIGMA].
;       The variance in the continuum around a Balmer line is measured
;       in the continuum windows defined in BALMER_DATA().  
; 
;       The upper limit is defined as
;       SNRCUT*LINECONTLEVEL_ERR*SQRT(2*!PI)*LINESIGMA.  
;
;       Also see ILOCAL_CONTINUUM().  
;
; INTERNAL SUPPORT ROUTINES:
;       IUPPER_COMPUTE_BOUNDARIES
;
; PROCEDURES USED:
;       REPSTR(), GET_ELEMENT, DJS_ITERSTAT, DJS_PLOT,
;       TELLURIC_MASK(), BALMER_MASK()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 June 11, U of A, excised from IFITSPEC() 
;       jm03nov26uofa - various modifications to reflect changes in
;                       IFITSPEC(); added DEBUG keyword and improved
;                       the documentation
;       jm04jan04uofa - treat Balmer and forbidden emission lines
;                       separately; added TELLURIC keyword
;       jm04mar01uofa - bug fix: use LINESIGMA_TOTAL when computing the
;                       total area of the line
;       jm04sep03uofa - use the actual error in the continuum rather
;                       than setting it equal to -3.0
;       jm07mar04nyu  - added SCALE optional input; compute the MEAN
;                       continuum level, not the median; also compute
;                       the ERROR in the MEAN continuum, not the
;                       standard deviation
;
; Copyright (C) 2003-2004, 2007, John Moustakas
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

pro iupper_compute_boundaries, wave, linefit, linesigma=linesigma, $
  linewave=linewave, glosigma=glosigma, ghisigma=ghisigma, xxlo=xxlo, $
  xxhi=xxhi

    light = 2.99792458D5            ; speed of light [km/s]

    linewave = linefit.linewave
;   linesigma = linewave*linefit.linesigma/light       ; [Angstrom]
    linesigma = linewave*linefit.linesigma_total/light ; [Angstrom]

    if (n_elements(glosigma) ne 0L) and (n_elements(ghisigma) ne 0L) then begin
    
       lo = glosigma*linesigma
       hi = ghisigma*linesigma
;      lo = (glosigma*linesigma)>5.0
;      hi = (ghisigma*linesigma)>20.0

       lolowave = linewave - hi
       lohiwave = linewave - lo
       get_element, wave, [lolowave,lohiwave], xxlo

       hilowave = linewave + lo
       hihiwave = linewave + hi
       get_element, wave, [hilowave,hihiwave], xxhi

    endif
       
return
end    
    
function iupper_limits, wave, flux, linefit, snrcut=snrcut, glosigma=glosigma, $
  ghisigma=ghisigma, telluric=telluric, scale=scale, debug=debug

    npix = n_elements(wave)
    nflux = n_elements(flux)
    nline = n_elements(linefit)

    if (npix eq 0L) or (nflux eq 0L) or (nline eq 0L) then begin
       print, 'Syntax - linefit = iupper_limits(wave,flux,linefit,$'
       print, '   snrcut=,glosigma=,ghisigma=,/telluric,/debug)'
       return, -1L
    endif

    if (npix ne nflux) then begin
       print, 'Dimensions of WAVE and FLUX do not agree.'
       return, -1L
    endif

    if n_elements(snrcut) eq 0L then snrcut = 1.0
    if n_elements(glosigma) eq 0L then glosigma = 3.0
    if n_elements(ghisigma) eq 0L then ghisigma = 10.0

    tmask = telluric_mask(wave)           ; telluric features
    bmask = balmer_mask(wave,bdata=bdata) ; Balmer absorption lines
    if keyword_set(telluric) then mask = ((tmask + bmask) eq 2B) else mask = (bmask eq 1B)
    badmask = where(mask eq 0B,nbadmask)
    
    for j = 0L, nline-1L do begin
       
       factor = 1.0

       linename = strupcase(repstr(linefit[j].linename,'_',' '))
       
; if this emission line is a Balmer line then compute the noise in the
; continuum using the continuum windows defined in BALMER_DATA()

       match = where(strmatch(bdata.line,'*'+strtrim(linefit[j].linename,2)+'*',/fold) eq 1B,nmatch)
       if (nmatch eq 1L) then begin

          get_element, wave, bdata[match].llimit + bdata[match].lwidth*[-1,+1]/2.0, xxlo
          get_element, wave, bdata[match].ulimit + bdata[match].uwidth*[-1,+1]/2.0, xxhi
          localmask = [mask[xxlo[0]:xxlo[1]],mask[xxhi[0]:xxhi[1]]]
          good = where(localmask,ngood)

          iupper_compute_boundaries, wave, linefit[j], $
            linesigma=linesigma, linewave=linewave

       endif else begin

; compute the Gaussian width in [Angstrom]; define the boundaries over
; which we will determine the continuum; if there are fewer than 5
; pixels with which to compute the continuum then multiply GLOSIGMA
; and GHISIGMA by FACTOR until the criterion is satisfied.  this will
; occur if the line is extremely narrow (tens of km/s) or for
; forbidden lines that are near Balmer lines, whose pixels are
; typically masked

          ngood = npix
          while ((ngood lt 5L) or (ngood eq npix)) and (factor le 10.0) do begin

             iupper_compute_boundaries, wave, linefit[j], linesigma=linesigma, $
               linewave=linewave, glosigma=factor*glosigma, ghisigma=factor*ghisigma, $
               xxlo=xxlo, xxhi=xxhi

             localmask = [mask[xxlo[0]:xxlo[1]],mask[xxhi[0]:xxhi[1]]]
             good = where(localmask,ngood) ; check for masked pixels

             factor = factor + 1.0

             if (ngood lt 5L) then $
               splog, linename+': Insufficient continuum pixels: GLOSIGMA, GHISIGMA, FACTOR = '+$
               string(factor,format='(I0)')
;            if (factor gt 10.0) then message, 'Problem here!'
          
          endwhile

       endelse 

       if (factor le 10.0) then begin ; this should never/rarely happen
          
          localflux = [flux[xxlo[0]:xxlo[1]],flux[xxhi[0]:xxhi[1]]]
          localwave = [wave[xxlo[0]:xxlo[1]],wave[xxhi[0]:xxhi[1]]]
          
          djs_iterstat, localflux[good], sigrej=3.0, mean=linecontlevel, $ ; median=linecontlevel, $
            sigma=linecontlevel_err, mask=lmask
          lpts = where(lmask,npts)
          linearea = snrcut*linecontlevel_err*sqrt(2.0*!pi)*linesigma

; fill the output structure
          
          linefit[j].linesigma_err = -3.0
          linefit[j].linearea = linearea
          linefit[j].linearea_err = -3.0
          linefit[j].linebox = 0.0
          linefit[j].linebox_err = -3.0
          linefit[j].linecontlevel = linecontlevel
          linefit[j].linecontlevel_err = linecontlevel_err/sqrt(npts>1L)
;         linefit[j].linecontlevel_err = -3.0
          linefit[j].linechi2 = -3.0

          if keyword_set(debug) then begin

             if !d.window ne 8 then window, 8, xs=550, ys=450

             if (n_elements(scale) eq 0L) then scale = 1E17
             ytitle = 'Flux Density [10^{-17} '+flam_units()+']'

             djs_plot, wave[(xxlo[0]-20)>0L:(xxhi[1]+20)<(npix-1)], $
               scale*flux[(xxlo[0]-20)>0L:(xxhi[1]+20)<(npix-1)], xsty=3, ysty=3, $
               color='yellow', charsize=1.5, charthick=2.0, ps=10, $
               xtitle='Wavelength [\AA]', ytitle=ytitle, $
               title=linename+' Upper Limit', thick=3.0, xthick=2.0, ythick=2.0
             polyfill, [wave[xxlo[0]],wave[xxlo[1]],wave[xxlo[1]],wave[xxlo[0]]], $ ; blue
               [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
               /line_fill, orientation=45, color=djs_icolor('blue'), thick=2.0
             polyfill, [wave[xxhi[0]],wave[xxhi[1]],wave[xxhi[1]],wave[xxhi[0]]], $ ; red
               [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
               /line_fill, orientation=45, color=djs_icolor('red'), thick=2.0
             djs_oplot, !x.crange, scale*linecontlevel*[1,1], line=0, thick=3.0, color='cyan'
             djs_oplot, !x.crange, scale*(linecontlevel*[1,1]-linecontlevel_err), line=2, thick=3.0, color='cyan'
             djs_oplot, !x.crange, scale*(linecontlevel*[1,1]+linecontlevel_err), line=2, thick=3.0, color='cyan'
             djs_oplot, linewave*[1,1], !y.crange, line=0, thick=3.0, color='green'
             if badmask[0] ne -1L then djs_oplot, wave[badmask], scale*flux[badmask], ps=4, $
               color='white', syms=1.5
             
             if j eq 0L then splog, 'Press any key to continue.'
             cc = get_kbrd(1)
             
          endif

       endif
          
    endfor

return, linefit
end
