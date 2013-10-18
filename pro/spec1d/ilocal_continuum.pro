;+
; NAME:
;       ILOCAL_CONTINUUM()
;
; PURPOSE:
;       Measure the local continuum around an emission line. 
;
; CALLING SEQUENCE:
;       linefit = ilocal_continuum(wave,flux,speclinefit,$
;          linename=,linewave=,linesigmakms=,glosigma=,$
;          ghisigma,/telluric,/debug,/postscript)
;
; INPUTS:
;       wave         - rest wavelength vector [A] 
;       flux         - data spectrum [erg/s/cm2/A] 
;       speclinefit  - fit to the emission-line spectrum from
;                      ILINEFIT() 
;       linename     - line name
;       linewave     - line wavelength [A]
;       linesigmakms - Gaussian line width [km/s]
;
; OPTIONAL INPUTS:
;       glosigma   - compute the noise properties of the continuum in
;                    the interval +/- [GLOSIGMA,GHISIGMA] * LINESIGMA
;                    Angstrom on either side of line center, where
;                    LINESIGMA is the Gaussian width of the line
;                    (default 5.0)
;       ghisigma   - see GLOSIGMA (default 15.0)
;
; KEYWORD PARAMETERS:
;       telluric   - mask telluric absorption
;       debug      - generate a plot useful for debugging
;       postscript - if set then do not open a window and suppress
;                    waiting for a keystroke between plots
;
; OUTPUTS:
;       ilocal  - properties of the local continuum structure with the
;                 following fiels:
;          LINECONTLEVEL     - local continuum
;          LINECONTLEVEL_ERR - local continuum error
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; INTERNAL SUPPORT ROUTINES:
;       ILOCAL_COMPUTE_BOUNDARIES
;
; PROCEDURES USED:
;       REPSTR(), GET_ELEMENT, DJS_ITERSTAT, DJS_PLOT,
;       TELLURIC_MASK(), BALMER_MASK(), FLAM_UNITS(), TEXTOIDL(),
;       INICE_LINENAME()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 November 26, U of A, written
;       jm03dec1uofa - minor updates and modifications; added
;                      POSTSCRIPT keyword
;       jm04jan6uofa - better error checking when encountering
;                      extremely narrow emission lines near masked
;                      Balmer lines
;       jm04feb04uofa - bug fix in the way the continuum error was
;                       being calculated
;       jm04apr22uofa - plotting improvements
;
; Copyright (C) 2003-2004, John Moustakas
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

pro ilocal_compute_boundaries, wave, linesigma=linesigma, linewave=linewave, $
  glosigma=glosigma, ghisigma=ghisigma, xxlo=xxlo, xxhi=xxhi

    lolowave = linewave - ghisigma*linesigma ; lower bound of lower window
    lohiwave = linewave - glosigma*linesigma ; upper bound of lower window
    get_element, wave, [lolowave,lohiwave], xxlo

    hilowave = linewave + glosigma*linesigma ; lower bound of upper window
    hihiwave = linewave + ghisigma*linesigma ; upper bound of upper window
    get_element, wave, [hilowave,hihiwave], xxhi

return
end    
    
function ilocal_continuum, wave, flux, speclinefit, linename=linename, $
  linewave=linewave, linesigmakms=linesigmakms, glosigma=glosigma, $
  ghisigma=ghisigma, telluric=telluric, debug=debug, postscript=postscript

    npix = n_elements(wave)
    nflux = n_elements(flux)
    nspeclinefit = n_elements(speclinefit)
    nline = n_elements(linename)

    if (npix eq 0L) or (nflux eq 0L) or (nspeclinefit eq 0L) or (nline eq 0L) then begin
       print, 'Syntax - linefit = ilocal_continuum(wave,flux,speclinefit,$'
       print, '   linename=,linewave=,linesigmakms=,glosigma=,ghisigma=,$'
       print, '   /telluric,/debug,/postscript)'
       return, -1L
    endif

    if (npix ne nflux) or (npix ne nspeclinefit) then begin
       print, 'Dimensions of WAVE, FLUX, and SPECLINEFIT do not agree.'
       return, -1L
    endif

    if (nline ne n_elements(linewave)) and (nline ne n_elements(linesigmakms)) then begin
       print, 'Dimensions of LINENAME, LINEWAVE, and LINESIGMAKMS do not agree.'
       return, -1L
    endif
    
    if n_elements(glosigma) eq 0L then glosigma = 3.0
    if n_elements(ghisigma) eq 0L then ghisigma = 10.0

    light = 2.99792458D5        ; speed of light [km/s]

    tmask = telluric_mask(wave)           ; telluric features
    bmask = balmer_mask(wave,bdata=bdata) ; Balmer absorption lines
    if keyword_set(telluric) then mask = ((tmask + bmask) eq 2B) else mask = (bmask eq 1B)
    badmask = where(mask eq 0B,nbadmask,comp=goodmask,ncomp=ngoodmask)

    efreeflux = flux - speclinefit ; remove the emission lines

; initialize the output structure

    ilocal = {linecontlevel: 0.0, linecontlevel_err: 0.0}
    ilocal = replicate(ilocal,nline)
    
    for j = 0L, nline-1L do begin

       factor = 1.0
       
       linewavej = linewave[j]
       linenamej = inice_linename(linename[j])
;      linenamej = strupcase(repstr(linename[j],'_',' '))
       linesigmaj = linewavej*linesigmakms[j]/light ; [Angstrom]
       
; if this emission line is a Balmer line then compute the noise in the
; continuum using the continuum windows defined in BALMER_DATA()

       match = where(strmatch(bdata.line,'*'+strtrim(linename[j],2)+'*',/fold) eq 1B,nmatch)
       if (nmatch eq 1L) then begin

          get_element, wave, bdata[match].llimit + bdata[match].lwidth*[-1,+1]/2.0, xxlo
          get_element, wave, bdata[match].ulimit + bdata[match].uwidth*[-1,+1]/2.0, xxhi
          localmask = [mask[xxlo[0]:xxlo[1]],mask[xxhi[0]:xxhi[1]]]
          good = where(localmask,ngood)

       endif else begin

; compute the Gaussian width in [Angstrom]; define the boundaries over
; which we will determine the continuum; if there are fewer than 5
; pixels with which to compute the continuum then multiply GLOSIGMA
; and GHISIGMA by FACTOR until the criterion is satisfied.  this will
; occur if the line is extremely narrow (tens of km/s) or for
; forbidden lines that are near Balmer lines, whose pixels are
; typically masked

          ngood = npix
          while (ngood lt 5L) or (ngood eq npix) do begin

             ilocal_compute_boundaries, wave, linesigma=linesigmaj, linewave=linewavej, $
               glosigma=factor*glosigma, ghisigma=factor*ghisigma, xxlo=xxlo, xxhi=xxhi

             localmask = [mask[xxlo[0]:xxlo[1]],mask[xxhi[0]:xxhi[1]]]
             good = where(localmask,ngood) ; check for masked pixels
          
             factor = factor + 1.0

             if (ngood lt 5L) then $
               splog, linenamej+': Insufficient continuum pixels: GLOSIGMA, GHISIGMA, FACTOR = '+$
               string(factor,format='(I0)')

             if (factor gt 10.0) then message, 'Problem here!'

          endwhile

       endelse
          
       localflux = [efreeflux[xxlo[0]:xxlo[1]],efreeflux[xxhi[0]:xxhi[1]]]
       localwave = [wave[xxlo[0]:xxlo[1]],wave[xxhi[0]:xxhi[1]]]

; compute the error in the continuum as the sigma-clipped rms of the
; continuum pixels divided by the square root of the number of pixels
; contributing to the rms
       
       djs_iterstat, localflux[good], sigrej=3.0, mask=cmask, median=linecontlevel, sigma=csigma
       linecontlevel_err = csigma / sqrt(n_elements(where(cmask)))

; fill the output structure
       
       ilocal[j].linecontlevel = linecontlevel
       ilocal[j].linecontlevel_err = linecontlevel_err

       if keyword_set(postscript) then $
         colorlist = ['black','dark blue','dark red','purple','red','dark green'] else $
         colorlist = ['yellow','blue','red','white','cyan','green']

       if keyword_set(debug) or keyword_set(postscript) then begin

          if keyword_set(postscript) then begin
             tthick = 8.0
          endif else begin
             if !d.window ne 8L then window, 8, xs=550, ys=450
             tthick = 2.0
          endelse

          scale = 1E15
          ytitle = 'Flux Density [10^{-15} '+flam_units()+']'

          pixoff = 15L
          xrange = fltarr(2)
          xrange[0] = interpol(wave,findgen(npix),xxlo[0]-pixoff)
          xrange[1] = interpol(wave,findgen(npix),xxhi[1]+pixoff)

          yrange = scale*[min(flux[(xxlo[0]-pixoff)>0L:(xxhi[1]+pixoff)<(npix-1)]),$
            max(flux[xxlo[1]:xxhi[0]])]*[1.0,1.1]
          
          djs_plot, wave[(xxlo[0]-pixoff)>0L:(xxhi[1]+pixoff)<(npix-1)], $
            scale*flux[(xxlo[0]-pixoff)>0L:(xxhi[1]+pixoff)<(npix-1)], xsty=3, ysty=3, $
            color=colorlist[0], charsize=1.5, charthick=tthick, ps=10, thick=tthick, $
            ytitle='', xtickname=replicate(' ',10), xthick=tthick, ythick=tthick, $
            position=[0.16,0.45,0.95,0.93], xrange=xrange, yrange=yrange
          legend, textoidl(linenamej), /left, /top, box=0, charsize=2.0, charthick=tthick

          polyfill, [wave[xxlo[0]],wave[xxlo[1]],wave[xxlo[1]],wave[xxlo[0]]], $ ; blue
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=0, color=djs_icolor(colorlist[1]), thick=1.0, spacing=0.04, linestyle=1
          polyfill, [wave[xxlo[0]],wave[xxlo[1]],wave[xxlo[1]],wave[xxlo[0]]], $ ; blue
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=90, color=djs_icolor(colorlist[1]), thick=1.0, spacing=0.04, linestyle=1

          polyfill, [wave[xxhi[0]],wave[xxhi[1]],wave[xxhi[1]],wave[xxhi[0]]], $ ; red
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=0, color=djs_icolor(colorlist[2]), thick=1.0, spacing=0.04, linestyle=1
          polyfill, [wave[xxhi[0]],wave[xxhi[1]],wave[xxhi[1]],wave[xxhi[0]]], $ ; red
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=90, color=djs_icolor(colorlist[2]), thick=1.0, spacing=0.04, linestyle=1

;         djs_oplot, wave, scale*(linecontlevel+speclinefit), ps=10, line=0, $
;           color=colorlist[2], thick=tthick

          djs_oplot, linewavej*[1,1], !y.crange, line=0, color=colorlist[5], thick=tthick+3
;         djs_oplot, lowave*[1,1], !y.crange, line=2, color='green', thick=tthick
;         djs_oplot, hiwave*[1,1], !y.crange, line=2, color='green', thick=tthick
          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]), line=0, color=colorlist[4], thick=tthick
          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]+linecontlevel_err), $
            line=2, color=colorlist[4], thick=tthick
          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]-linecontlevel_err), $
            line=2, color=colorlist[4], thick=tthick

          if (badmask[0] ne -1L) then begin

             maskwave = wave[badmask]
             maskflux = scale*flux[badmask]
             show = where((maskwave gt !x.crange[0]) and (maskwave lt !x.crange[1]),nshow)
             if (nshow ne 0L) then begin
                djs_oplot, maskwave[show], maskflux[show], line=0, thick=tthick+3, $
                  color=colorlist[3], ps=10
             endif

          endif

;         if badmask[0] ne -1L then djs_oplot, wave[badmask], scale*flux[badmask], ps=4, $
;           color=colorlist[3], syms=1.5

          xyouts, 0.08, 0.52, textoidl(ytitle), /normal, align=0.5, orientation=90, $
            charsize=1.5, charthick=tthick
          
; emission-line subtracted plot          

          djs_plot, wave[(xxlo[0]-pixoff)>0L:(xxhi[1]+pixoff)<(npix-1)], $
            scale*efreeflux[(xxlo[0]-pixoff)>0L:(xxhi[1]+pixoff)<(npix-1)], xsty=3, ysty=3, $
            color=colorlist[0], charsize=1.5, charthick=tthick, ps=10, /noerase, $
            ytitle='', xtitle='Rest Wavelength [\AA]', position=[0.16,0.13,0.95,0.45], $
            yrange=scale*(max(abs(localflux[good]-linecontlevel))*[-1.5,1.5]+linecontlevel), $
            xthick=tthick, ythick=tthick, thick=tthick, xrange=xrange
          djs_oplot, linewavej*[1,1], !y.crange, line=0, color=colorlist[5], thick=tthick+3

          polyfill, [wave[xxlo[0]],wave[xxlo[1]],wave[xxlo[1]],wave[xxlo[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=0, color=djs_icolor(colorlist[1]), thick=1.0, spacing=0.04, linestyle=1
          polyfill, [wave[xxlo[0]],wave[xxlo[1]],wave[xxlo[1]],wave[xxlo[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=90, color=djs_icolor(colorlist[1]), thick=1.0, spacing=0.04, linestyle=1

          polyfill, [wave[xxhi[0]],wave[xxhi[1]],wave[xxhi[1]],wave[xxhi[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=0, color=djs_icolor(colorlist[2]), thick=1.0, spacing=0.04, linestyle=1
          polyfill, [wave[xxhi[0]],wave[xxhi[1]],wave[xxhi[1]],wave[xxhi[0]]], $
            [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
            /line_fill, orientation=90, color=djs_icolor(colorlist[2]), thick=1.0, spacing=0.04, linestyle=1

          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]), line=0, color=colorlist[4], thick=tthick

          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]+linecontlevel_err), $
            line=2, color=colorlist[4], thick=tthick
          djs_oplot, !x.crange, scale*(linecontlevel*[1,1]-linecontlevel_err), $
            line=2, color=colorlist[4], thick=tthick

          if (badmask[0] ne -1L) then begin

             maskwave = wave[badmask]
             maskflux = scale*efreeflux[badmask]
             show = where((maskwave gt !x.crange[0]) and (maskwave lt !x.crange[1]),nshow)
             if (nshow ne 0L) then begin
                djs_oplot, maskwave[show], maskflux[show], line=0, thick=tthick+3, $
                  color=colorlist[3], ps=10
             endif

          endif

;         if badmask[0] ne -1L then djs_oplot, wave[badmask], scale*efreeflux[badmask], ps=4, $
;           color=colorlist[3], syms=1.5

          if not keyword_set(postscript) then begin
             if (nline gt 1L) then begin
                if (j eq 0L) then splog, 'Press any key to continue.'
                cc = get_kbrd(1)
             endif
          endif
          
       endif

    endfor

return, ilocal
end
