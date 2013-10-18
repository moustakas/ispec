;+
; NAME:
;       IFIGURE_SPECLINEFIT
;
; PURPOSE:
;       Generate a QA plot of the spectral fitting results.  
;
; CALLING SEQUENCE:
;       ifigure_speclinefit, specdata, specfit, starflux, $
;          /postscript, _extra=extra
;
; INPUTS:
;       specdata - see ISPECLINEFIT()()
;       specfit  - see ISPECLINEFIT()()
;       starflux - see ISPECLINEFIT()()
;
; OPTIONAL INPUTS:
;       extra - optional inputs for K_LAMBDA()
;
; KEYWORD PARAMETERS:
;       postscript - if generating postscript then set this keyword to
;                    prevent initializing IDL windows
;
; OUTPUTS:
;       Three plots are generated.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;      SPLOG, EMISSION_MASK(), MRDFITS(), PAGEMAKER, DJS_ITERSTAT,
;      DJS_PLOT, DJS_OPLOT, GET_ELEMENT, LEGEND, STRUCT_TRIMTAGS(),
;      PLOTSYM, IM_STATS(), INICE_LINENAME(), TEXTOIDL(), IM_WINDOW 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 12, U of A - based on an earlier code 
;       jm04apr28uofa - minor improvements; overplot polynomials;
;          added EXTRA keyword
;       jm05jul26uofa - better ranges on the emission-line plot
;       jm07apr12nyu  - special case of no templates used 
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

pro ifigure_speclinefit, specdata1, specfit, starflux, postscript=postscript, $
  _extra=extra

    if (n_elements(specdata1) eq 0L) or (n_elements(specfit) eq 0L) then begin
       print, 'Syntax - ifigure_speclinefit, specdata, specfit, $'
       print, '   starflux, /postscript, _extra=extra'
       return
    endif

    specdata = specdata1
    
    light = 2.99792458D5        ; speed of light [km/s]

; define some convenient internal variables

    galaxy = strcompress(specdata.galaxy,/remove)
    specfile = strcompress(specdata.specfile,/remove)
    if strmatch(specfile,'*ms*') then $
      specfile = strmid(specfile,0,strpos(specfile,'.ms.fits')) else $
      specfile = strmid(specfile,0,strpos(specfile,'.fits'))

    z_obj = specdata.z_abs
    z_line = specdata.z_line
    
;   elinewave = specfit[*,0]/(1+z_line) 
    elinewave = reform(specfit[*,0])*(1.0+z_obj)/(1.0+z_line) ; emission line rest wavelength vector
;   elinewave = specfit[*,0]
    wave = specfit[*,0]        ; continuum rest wavelength vector
    flux = specfit[*,1]        ; rest flux
    continuum = specfit[*,2]   ; rest continuum
    speclinefit = specfit[*,3] 
    medcontinuum = specfit[*,4]
    espectrum = flux-continuum-medcontinuum
    npix = n_elements(wave)
    ntemplate = specdata.ntemplate
;   lineres = specdata.specres

; do not plot the broad and narrow components, if they exist
    
;   linekeep = where((strmatch(specdata.linename,'*narrow*',/fold) eq 0B) and $
;     (strmatch(specdata.linename,'*broad*',/fold) eq 0B),nlinekeep)
;   if (nlinekeep ne 0L) then begin
;      specdata = struct_addtags({linename: specdata.linename[linekeep]},struct_trimtags(specdata,except='LINENAME'))
;   endif
    nline = n_elements(specdata.linename)
    
    mask = emission_mask(elinewave*(1+z_line),z=z_line,width=20.0,$
      bad=maskpix,good=keeppix,/telluric,/bluemask,bluewidth=10.0)

    charsize = 1.1
    plotsym, 0, 0.3, /fill

    if keyword_set(postscript) then begin
       postthick1 = 1.0
       postthick2 = 1.0
       postthick3 = 1.0
       clear = 1L
    endif else begin
       postthick1 = 2.0
       postthick2 = 1.0
       postthick3 = 1.0
       clear = 0L
    endelse
    
; ---------------------------------------------------------------------------    
; page 1 postscript    
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=2, position=position, /normal, $
      xmargin=[0.8,0.2], ymargin=[0.4,1.0], yspace=0.0

    if (not keyword_set(postscript)) then im_window, 0, xratio=0.7, yratio=0.5

    stats = im_stats(flux,sigrej=5.0)
    yrange = [-0.15,1.2]*stats.maxrej
;   yrange = stats.median_rej+[-7,+8]*stats.sigma_rej
;   yrange[0] = yrange[0]>(0.2*min(continuum))
;   yrange[1] = yrange[1]<(1.4*max(continuum))
    xrange = minmax(wave)
    norm = stats.median

    djs_plot, wave, flux/norm, xsty=3, ysty=3, ps=10, xthick=postthick1, xrange=xrange, $
      ythick=postthick1, charsize=charsize, charthick=postthick1, yrange=yrange/norm, thick=postthick3, $
      position=position[*,0], xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      ytitle='Relative Flux', color='grey'
    djs_oplot, wave, (continuum+medcontinuum)/norm, line=0, ps=10, thick=postthick2+1.0, color='red'
    djs_oplot, wave, continuum/norm, line=0, ps=10, thick=postthick2, color='yellow'
;   djs_oplot, wave, smooth(espectrum,6)/norm, ps=10, thick=postthick2, color='grey'
;   djs_oplot, wave, espectrum/norm, ps=10, thick=postthick2, color='grey'
;   djs_oplot, wave, speclinefit/norm, ps=10, thick=postthick2
    djs_oplot, !x.crange, [0,0], thick=pthick, color='blue'

; figure out the indices for contiguously masked regions

    nmaskpix = n_elements(maskpix)
    shiftpix = maskpix-shift(maskpix,1)
    shiftpix[0L] = 1L
    lopoints = [0,where(shiftpix ne 1L)]
    hipoints = [where(shiftpix ne 1L)-1L,nmaskpix-1L]
    nregions = n_elements(lopoints)

    for i = 0L, nregions-1L do begin

       lo = maskpix[lopoints[i]]
       hi = maskpix[hipoints[i]]
       
;      polyfill, [wave[lo],wave[hi],wave[hi],wave[lo]], $
;        [!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]], /data, $
;        color=djs_icolor('purple'), linestyle=1, /line_fill, orientation=90,$
;        spacing=0.01;, thick=3.0

    endfor
    
; information legend

    label = [$
      textoidl('\chi^2 = '+strtrim(string(specdata.continuum_chi2,format='(F7.2)'),2)),$
      'S/N = '+strtrim(string(specdata.continuum_snr,format='(I4)'),2),$
      textoidl('z_{abs} = '+strtrim(string(specdata.z_abs,format='(F12.4)'),2))]
    legend, label, /left, /top, charsize=charsize, charthick=postthick1, box=0, clear=clear

    label = [strupcase(specdata.template_imf)+' IMF',$
      'Z = '+string(specdata.continuum_Z,format='(F5.3)')]
    legend, label, /right, /top, charsize=charsize, charthick=postthick1, box=0, clear=clear

; now overplot the templates

    mass_fraction = 100*specdata.continuum_coeff/total(specdata.continuum_coeff)

; compute the mass-to-light ratio at the midpoint of the wavelength
; range of the spectrum

    midwave = fix(djs_mean(wave))
    strmidwave = string(midwave,format='(I0)')+' \AA'
    
    ebv = specdata.continuum_ebv

;;  template_ml = dblarr(ntemplate)
;;  mlwave = [3500.0,4400.0,5500.0,6500.0] ; we have M/L ratios at these wavelengths
;;  big_template_ml = [ [specdata.template_ml_u], [specdata.template_ml_b], $
;;    [specdata.template_ml_v], [specdata.template_ml_r] ]
;;  for i = 0L, ntemplate-1L do template_ml[i] = interpol(big_template_ml[i,*],mlwave,midwave)

    kl = k_lambda(wave,_extra=extra)
    klmidwave = k_lambda(midwave,_extra=extra)
    
    fraction = specdata.continuum_coeff
;;  fraction = specdata.continuum_coeff/template_ml            ; [L_sun/cm2]
    if (total(fraction) gt 0.0) then begin
       light_fraction = 100*fraction/total(fraction)
       good = where((light_fraction gt 1.0),ngood)
       if (ngood ne 0L) then begin
          plotflux = dblarr(npix,ngood)
          for it = 0L, ngood-1L do plotflux[*,it] = light_fraction[good[it]] * $ ; <-- NOT GENERAL
            starflux[*,good[it]] * 10^(-0.4*ebv[it]*(kl-klmidwave)) / interpol(starflux[*,good[it]],wave,midwave)
          age = string(round(specdata.template_age[good]/1E6),format='(I5)')
;;        mfraction = string(mass_fraction[good],format='(F5.1)')
          lfraction = string(light_fraction[good],format='(F5.1)')
          ebv = string(specdata.continuum_ebv[good],format='(F5.2)')
          nage = string(ngood,format='(I0)')
       endif
    endif else ngood = 0L

    if (total(fraction) le 0.0) or (ngood eq 0L) then begin
       age = '...' & mfraction = '...' & lfraction = '...' & ebv = '...' & nage = '1'
    endif
       
    polygood = where(specdata.continuum_polycoeff ne 0.0,npoly)
    if (npoly ne 0L) then begin
       polyflux = poly_array(npix,npoly)
       for ip = 0L, npoly-1L do begin
          polyflux[*,ip] = polyflux[*,ip]*specdata.continuum_polycoeff[ip]
          polyfrac = 100/interpol(continuum,wave,midwave)
          polyflux[*,ip] = polyfrac*polyflux[*,ip]
       endfor          
    endif

    if (n_elements(plotflux) eq 0L) then yrange = [0,1] else yrange = [min(plotflux),max(plotflux)*1.3]
   
    ndim = size(plotflux,/n_dimension)
    dims = size(plotflux,/dimension)
    case ndim of
       0L: nplot = 0L
       1L: nplot = 1L
       2L: nplot = dims[1]
       else: message, 'Problem here.'
    endcase
;   if (ndim eq 1L) then nplot = 1L else nplot = dims[1]

    djs_plot, [0], [0], /nodata, xrange=xrange, /noerase, charsize=charsize, $
      charthick=postthick1, thick=postthick3, xthick=postthick1, ythick=postthick1, position=position[*,1], $
      color='blue', xtitle='Rest Wavelength (\AA)', xsty=3, ysty=3, $
      yrange=yrange, ps=10, yminor=2, ytitle='Light Contribution at '+strmidwave+' (%)'; ytickname=replicate(' ',10), 
    for it = 0L, nplot-1L do djs_oplot, wave, plotflux[*,it], color='blue', thick=postthick3, ps=10

; overplot the polynomials

    if (npoly ne 0L) then djs_oplot, wave, total(reform(polyflux,npix,npoly),2), color='dark green', thick=postthick3
;   if (npoly ne 0L) then for ip = 0L, npoly-1L do djs_oplot, wave, $
;     polyflux[*,ip], color='dark green', thick=postthick2
    
; fractional contribution legend    
    
;   label = ['Age   '+strjoin(age,' '),'Mass  '+strjoin(mfraction,' '),'E(B-V)'+strjoin(ebv,' ')]
    
    label = [$
      'Age    '+string(age,format='('+nage+'A6)')+' Myr',$
      'Light  '+string(lfraction,format='('+nage+'A6)')+'   %',$
;;    'Mass   '+string(mfraction,format='('+nage+'A6)')+'   %',$
      'E(B-V) '+string(ebv,format='('+nage+'A6)')+' mag' $
      ]

;   xyouts, 0.1, 0.48, label[0], /normal, charsize=charsize, charthick=postthick1
;   xyouts, 0.1, 0.46, label[1], /normal, charsize=charsize, charthick=postthick1
;   xyouts, 0.1, 0.44, label[2], /normal, charsize=charsize, charthick=postthick1
;   xyouts, 0.1, 0.42, label[3], /normal, charsize=charsize, charthick=postthick1
    
    legend, label, /left, /top, charsize=charsize, charthick=postthick1, box=0, $
      corners=corn, clear=clear

; title    
    
    refpos = reform(position[*,0])
    xpos = (refpos[2]-refpos[0])/2.0+refpos[0]
    ypos = refpos[3]*1.01
    
    xyouts, xpos, ypos, galaxy+' ['+specfile+'] Continuum', /normal, charsize=1.5, charthick=postthick1, align=0.5

; ---------------------------------------------------------------------------
; page 2 postscript
; ---------------------------------------------------------------------------

    pagemaker, nx=4, ny=ceil(nline/4.0), position=position, /normal, $
      xmargin=[0.6,0.2], ymargin=[0.4,0.8], yspace=0.0, xspace=0.0

    if (not keyword_set(postscript)) then im_window, 3, xratio=0.35, /square

    box = 100.0
    tags = tag_names(specdata[0])
    for k = 0L, nline-1L do begin

       if k eq 0L then noerase = 0L else noerase = 1L

       line = inice_linename(specdata.linename[k])

; emission-line wavelength
       
       match = where(specdata.linename[k]+'_WAVE' eq tags)
       match2 = where(specdata.linename[k]+'_LINEZ' eq tags)
       match3 = where(specdata.linename[k]+'_CHI2' eq tags)

;      linewave = specdata.(match)*(1+z_obj)/(1+z_line)
       linewave = specdata.(match)
       linez = (specdata.(match2))[0]
;      linewave = specdata.(match)*(1.0+(specdata.(match2))[0])
;      linewave = specdata.(match)*(1.0+specdata.z_line)

; total sigma line-width       

       match = where(specdata.linename[k]+'_SIGMA_TOTAL' eq tags)
       lineres = linewave*specdata.(match)/light ; [Angstrom]

       leftbox  = linewave-((15.0*lineres)<130.0)
       rightbox = linewave+((15.0*lineres)<130.0)
       get_element, elinewave, [leftbox,rightbox], xx

       leftbox_zoom  = linewave-1.0*lineres
       rightbox_zoom = linewave+1.0*lineres
       get_element, elinewave, [leftbox_zoom,rightbox_zoom], xx_zoom

       yrange = fltarr(2)
       yrange[1] = max(espectrum[xx_zoom[0]:xx_zoom[1]])*1.8
;      yrange[0] = 0.0 ; min(espectrum[xx[0]:xx[1]])
       yrange[0] = -0.1*yrange[1]

       match = where(specdata.linename[k] eq tags)
       errflag = specdata.(match)[1]
       
       if errflag eq -2.0 then begin

          upper = 'Not Measured'
          label = [line,upper]

          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
            xthick=postthick1, ythick=postthick1, thick=postthick2, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange

       endif else begin
          
          djs_plot, elinewave[xx[0]:xx[1]], espectrum[xx[0]:xx[1]], xsty=3, ysty=3, $
            xthick=postthick1, ythick=postthick1, thick=postthick2, noerase=noerase, ps=10, $
            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange, color='grey'
          djs_oplot, elinewave[xx[0]:xx[1]], speclinefit[xx[0]:xx[1]], thick=postthick2, $
            color='red', ps=10
          djs_oplot, linewave*[1,1], !y.crange, line=2, thick=postthick2, color='navy'
;         djs_oplot, (linewave*(1.0+z_line)/(1.0+linez))*[1,1], !y.crange, line=2, thick=postthick2, color='navy'
          djs_oplot, !x.crange, [0,0], line=0, thick=postthick2, color='navy'

; if this is a Balmer line then plot the extent of the Balmer
; absorption window

          balmermatch = where(strmatch(specdata.babs_linename,'*'+specdata.linename[k]+'*',/fold) eq 1B,nbmatch)
          if (nbmatch ne 0L) then begin

             match1 = where('BABS_'+strupcase(specdata.babs_linename[balmermatch[0]])+'_LLINE' eq tags)
             match2 = where('BABS_'+strupcase(specdata.babs_linename[balmermatch[0]])+'_ULINE' eq tags)

             djs_oplot, specdata.(match1)*[1,1], !y.crange, line=0, thick=postthick2, color='dark green'
             djs_oplot, specdata.(match2)*[1,1], !y.crange, line=0, thick=postthick2, color='dark green'
          
          endif

          case errflag of
             -3.0: label = [line,'Upper Limit']
             else: begin
                snr = 'S/N = '+string((specdata.(match))[0]/(specdata.(match))[1],format='(I0)')
                chi2 = '\chi^{2} = '+strtrim(string(specdata.(match3),format='(F12.1)'),2)
                label = [line,snr+', '+chi2]
             endelse
          endcase
       
       endelse

       legend, textoidl(label), /left, /top, box=0, charsize=1.0, $
         charthick=postthick1, clear=keyword_set(postscript)

    endfor

; title    
    
    xyouts, xpos, ypos, galaxy+' Line Fluxes', /normal, charsize=1.5, charthick=postthick1, align=0.5
;   xyouts, xpos, ypos, galaxy+' ['+specfile+'] Line Fluxes', /normal, charsize=1.5, charthick=postthick1, align=0.5
    
;; ---------------------------------------------------------------------------
;; page 3 postscript
;; ---------------------------------------------------------------------------
;
;    pagemaker, nx=4, ny=ceil(nline/4.0), position=position, /normal, $
;      xmargin=[0.6,0.2], ymargin=[0.4,0.8], yspace=0.0, xspace=0.0
;
;    if (not keyword_set(postscript)) then im_window, 1, xratio=0.35, /square
;
;    minibox = 50.0
;    tags = tag_names(specdata[0])
;
;; plot the emission lines
;
;    for k = 0L, nline-1L do begin
;
;       if k eq 0L then noerase = 0L else noerase = 1L
;
;       line = inice_linename(specdata.linename[k])
;       
;; emission-line wavelength
;       
;       match = where(specdata.linename[k]+'_WAVE' eq tags)
;;      match2 = where(specdata.linename[k]+'_LINEZ' eq tags)
;;      linewave = specdata.(match)*(1+z_obj)/(1+z_line)
;       linewave = specdata.(match)
;;      linewave = specdata.(match)*(1.0+(specdata.(match2))[0])
;;      linewave = specdata.(match)*(1.0+specdata.z_line)
;
;       leftbox  = linewave-20.0*lineres
;       rightbox = linewave+20.0*lineres
;       get_element, elinewave, [leftbox,rightbox], xx
;
;; error flag       
;
;       match = where(specdata.linename[k]+'_EW' eq tags)
;       errflag = specdata.(match)[1]
;
;       if (errflag eq -2.0) then begin
;
;          upper = 'Not Measured'
;          label = [line,upper]
;
;          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
;            xthick=postthick1, ythick=postthick1, thick=postthick2, noerase=noerase, ps=10, $
;            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;            position=position[*,k], yrange=yrange
;
;       endif else begin
;
;; continuum measurement
;       
;          cmatch = where(specdata.linename[k]+'_CONTINUUM' eq tags)
;          line_continuum = (specdata.(cmatch))[0]
;          line_continuum_err = (specdata.(cmatch))[1]
;
;          djs_iterstat, flux[xx[0]:xx[1]], sigma=csig, sigrej=2.5
;;         djs_iterstat, continuum[ww[0]:ww[1]], sigma=csig, sigrej=5.0
;          yrange = line_continuum+[-5,+4]*csig
;
;          djs_plot, elinewave[xx[0]:xx[1]], flux[xx[0]:xx[1]], xsty=3, ysty=3, $
;            xthick=postthick1, ythick=postthick1, thick=postthick2, noerase=noerase, ps=10, $
;            xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;            position=position[*,k], yrange=yrange, color='grey'
;
;          djs_oplot, elinewave[xx[0]:xx[1]], continuum[xx[0]:xx[1]]+speclinefit[xx[0]:xx[1]], $
;            color='red', ps=10, thick=postthick2
;;         djs_oplot, wave[xx[0]:xx[1]], continuum[xx[0]:xx[1]], color='green', ps=10, thick=postthick2
;          djs_oplot, !x.crange, line_continuum*[1,1], line=0, thick=postthick2, color='navy'
;          djs_oplot, !x.crange, line_continuum*[1,1]-line_continuum_err, line=2, thick=postthick2, color='navy'
;          djs_oplot, !x.crange, line_continuum*[1,1]+line_continuum_err, $
;            line=2, thick=postthick2, color='navy'
;          djs_oplot, linewave*[1,1], !y.crange, line=2, thick=postthick2, color='navy'
;
;; if this is a Balmer line then plot the extent of the Balmer
;; absorption window
;
;          balmermatch = where(strmatch(specdata.babs_linename,'*'+specdata.linename[k]+'*',/fold) eq 1B,nbmatch)
;          if (nbmatch ne 0L) then begin
;
;             match1 = where('BABS_'+strupcase(specdata.babs_linename[balmermatch[0]])+'_LLINE' eq tags)
;             match2 = where('BABS_'+strupcase(specdata.babs_linename[balmermatch[0]])+'_ULINE' eq tags)
;
;             djs_oplot, specdata.(match1)*[1,1], !y.crange, line=0, thick=postthick2, color='dark green'
;             djs_oplot, specdata.(match2)*[1,1], !y.crange, line=0, thick=postthick2, color='dark green'
;          
;          endif
;
;          case errflag of
;             -3.0: label = [line,'Upper Limit']
;             else: begin
;                snr = 'S/N = '+string((specdata.(match))[0]/(specdata.(match))[1],format='(I0)')
;                label = [line,snr]
;             endelse
;          endcase
;       
;       endelse
;          
;       legend, label, /left, /bottom, box=0, charsize=1.0, $
;         charthick=postthick1, clear=keyword_set(postscript)
;
;    endfor
;
;; title    
;    
;    xyouts, xpos, ypos, galaxy+' EWs', /normal, charsize=1.5, charthick=postthick1, align=0.5

return
end    
