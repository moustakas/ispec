;+
; NAME:
;   QAPLOT_ISPECLINEFIT
;
; PURPOSE:
;   Generate a QA plot for ISPECLINEFIT.
;
; INPUTS:
;   specdata - see ISPECLINEFIT()()
;   specfit  - see ISPECLINEFIT()()
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   postscript - if generating postscript then set this keyword to
;                prevent initializing IDL windows
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 May 12, U of A - based on an earlier code 
;   jm04apr28uofa - minor improvements; overplot polynomials;
;      added EXTRA keyword
;   jm05jul26uofa - better ranges on the emission-line plot
;   jm07apr12nyu  - special case of no templates used 
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

pro qaplot_ispeclinefit, specdata, specfit, templateinfo, $
  unfluxed=unfluxed, postscript=postscript

    if (n_elements(specdata) eq 0L) or (n_elements(specfit) eq 0L) then begin
       doc_library, 'qaplot_ispeclinefit'
       return
    endif

    light = 2.99792458D5 ; speed of light [km/s]

;    if keyword_set(postscript) then begin
;       postthick1 = 2.5 ; 3.0
;       postthick2 = 1.0
;       postthick3 = 1.5 ; 2.0
;       postthick4 = 3.5 ; 4.0
;;      clear = 1L
;    endif else begin
;       postthick1 = 1.0
;       postthick2 = 1.0
;       postthick3 = 1.0
;       postthick4 = 1.0
;       clear = 0L
;    endelse
    
    charsize1 = 1.4
    colors = fsc_color(['medium gray','firebrick',$
      'goldenrod','firebrick'],[10,11,12,13])

; define some convenient internal variables

    invvar = reform(specfit[*,5])
    good = where(invvar gt 0.0)
    invvar = invvar[good]

    wave = reform(specfit[good,0])
    flux = reform(specfit[good,1])
    continuum = reform(specfit[good,2])
    speclinefit = reform(specfit[good,3])
    smooth_continuum = reform(specfit[good,4])
    lineflux = flux - continuum - smooth_continuum

    xtitle1 = 'Rest Wavelength ('+cgsymbol('angstrom')+')'
    power = ceil(abs(alog10(median(flux))))
    if keyword_set(unfluxed) then begin
       scale = 10.0^(-power)
       ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} ADU)'
    endif else begin
       scale = 10.0^power
       ytitle1 = 'Flux (10^{-'+string(power,format='(I0)')+'} '+flam_units()+')'
    endelse
    
; ---------------------------------------------------------------------------    
; page1 - continuum plot
; ---------------------------------------------------------------------------    
    
    pagemaker, nx=1, ny=2, position=position, /normal, $
      xmargin=[0.8,0.2], ymargin=[0.4,0.8], yspace=0.7

    xrange = [min(wave),mean(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /xsty, /ysty, xthick=postthick1, ythick=postthick1, $
      xrange=xrange, yrange=scale*yrange, charsize=charsize1, charthick=postthick1, $
      position=position[*,0], title=strtrim(repstr(specdata.galaxy,'_',' '),2)+$
      ', z_{abs} = '+strtrim(string(specdata.z_abs,format='(F12.5)'),2)+$
      ' (z_{line} = '+strtrim(string(specdata.z_line,format='(F12.5)'),2)+', '+$
      ' z_{obj} = '+strtrim(string(specdata.z_obj,format='(F12.5)'),2)+')', $
      ytitle=ytitle1
    djs_oplot, wave, scale*flux, ps=10, thick=postthick2, color=colors[0]
    djs_oplot, wave, scale*continuum, ps=10, thick=postthick3, color=colors[1]
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      thick=postthick3, color=colors[2]

    if keyword_set(unfluxed) then begin
       leftlabel = ['S/N = '+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)]
    endif else begin
       leftlabel = [$
         textoidl('\chi^2 = '+strtrim(string(specdata.continuum_chi2,format='(F12.2)'),2)),$
         'S/N = '+strtrim(string(specdata.continuum_snr,format='(F12.1)'),2)]
       rightlabel = [$
;        strupcase(specdata.template_imf),$
         'Z = '+string(templateinfo.template_z,format='(F5.3)'),$
         'E(B-V) = '+strtrim(string(specdata.continuum_ebv,format='(F12.3)'),2)]
       im_legend, rightlabel, /right, /top, box=0, clear=clear, $
         charsize=charsize1, charthick=postthick1, margin=0
    endelse
    im_legend, leftlabel, /left, /top, box=0, clear=clear, charsize=charsize1, $
      charthick=postthick1, margin=0

    xrange = [mean(wave),max(wave)]
    get_element, wave, xrange, xx
    stats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
    yrange = [-0.12,1.2]*stats.maxrej

    djs_plot, [0], [0], /nodata, /noerase, /xsty, /ysty, xthick=postthick1, ythick=postthick1, $
      xrange=xrange, yrange=scale*yrange, charsize=charsize1, charthick=postthick1, $
      position=position[*,1], ytitle=ytitle1, xtitle=xtitle1
    djs_oplot, wave, scale*flux, ps=10, thick=postthick2, color=colors[0]
    djs_oplot, wave, scale*continuum, ps=10, thick=postthick3, color=colors[1]
    djs_oplot, wave, scale*(continuum+smooth_continuum), ps=10, $
      thick=postthick3, color=colors[2]
    
 ; ---------------------------------------------------------------------------
 ; page 2 - emission lines
 ; ---------------------------------------------------------------------------

    nline = n_elements(specdata.linename)
    
    pagemaker, nx=4, ny=ceil(nline/4.0), position=position, /normal, $
      xmargin=[0.6,0.2], ymargin=[0.4,0.8], yspace=0.0, xspace=0.0
    
    box = 100.0
    tags = tag_names(specdata[0])
    for k = 0L, nline-1L do begin
       
       if k eq 0L then noerase = 0L else noerase = 1L
       
       line = inice_linename(specdata.linename[k])
 
; emission-line wavelength
       
       match = where(specdata.linename[k]+'_WAVE' eq tags)
       match2 = where(specdata.linename[k]+'_LINEZ' eq tags)
       match3 = where(specdata.linename[k]+'_CHI2' eq tags)
       
       linez = (specdata.(match2))[0]
       linewave = specdata.(match)*(1+(linez-specdata.z_obj)) ; adjust for the z differences
;      linewave = specdata.(match)*(1+(linez-specdata.z_abs)) ; adjust for the z differences
       
; total sigma line-width       
       
       match4 = where(specdata.linename[k]+'_SIGMA' eq tags)
       match5 = where(specdata.linename[k]+'_SIGMA_TOTAL' eq tags)
       lineres = linewave*specdata.(match5)/light ; [Angstrom]
       
       leftbox  = linewave-((15.0*lineres)<130.0)
       rightbox = linewave+((15.0*lineres)<130.0)
       get_element, wave, [leftbox,rightbox], xx

       leftbox_zoom  = linewave-1.0*lineres
       rightbox_zoom = linewave+1.0*lineres
       get_element, wave, [leftbox_zoom,rightbox_zoom], xx_zoom
       
       match = where(specdata.linename[k] eq tags)
       errflag = specdata.(match)[1]
       
       if errflag eq -2.0 then begin
          
          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick1, $
            ythick=postthick1, noerase=noerase, ps=10, xtickname=replicate(' ',10), $
            ytickname=replicate(' ',10), position=position[*,k], yrange=yrange
          im_legend, textoidl([line,'Not Measured']), /left, /top, box=0, charsize=1.2, $
            charthick=postthick1, clear=0, margin=0
          
       endif else begin
          
;         stats = im_stats(lineflux[xx[0]:xx[1]],sigrej=3.0)
          yrange = fltarr(2)
          yrange[1] = (max(abs(lineflux[xx_zoom[0]:xx_zoom[1]]))*1.8);>(3.0*stats.sigma_rej)
          yrange[0] = -0.1*yrange[1]
;         yrange[0] = (-0.1*yrange[1])<(-1.5*stats.sigma_rej)

          djs_plot, wave[xx[0]:xx[1]], lineflux[xx[0]:xx[1]], xsty=3, ysty=3, $
            xthick=postthick1, ythick=postthick1, thick=postthick3, noerase=noerase, $
            ps=10, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
            position=position[*,k], yrange=yrange, color=colors[0]
          djs_oplot, wave[xx[0]:xx[1]], speclinefit[xx[0]:xx[1]], thick=postthick1, $
            color=colors[3], ps=10
          djs_oplot, linewave*[1,1], !y.crange, line=5, thick=postthick2

          case errflag of
             -3.0: begin
                im_legend, textoidl([line,'Upper Limit']), /left, /top, box=0, charsize=1.2, $
                  charthick=postthick1, clear=0, margin=0
             end
             else: begin
                snr = 'S/N = '+string((specdata.(match))[0]/(specdata.(match))[1],format='(I0)')
                chi2 = '\chi^{2} = '+strtrim(string(specdata.(match3),format='(F12.1)'),2)
                sigma = '\sigma = '+strtrim(string((specdata.(match5))[0],format='(F12.1)'),2);+' km/s'
;               leftlabel = [line,snr+', '+chi2]
                im_legend, textoidl([line,snr]), /left, /top, box=0, charsize=1.2, $
                  charthick=postthick1, clear=0, margin=0
                im_legend, textoidl([sigma,chi2]), /right, /top, box=0, charsize=1.2, $
                  charthick=postthick1, clear=0, margin=0
             endelse
          endcase
          
       endelse

    endfor 

return
end    
