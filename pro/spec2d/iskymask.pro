;+
; NAME:
;	ISKYMASK()
;
; PURPOSE:
;	Mask and interpolate over bad sky residuals.
;
; CALLING SEQUENCE:
;       newspec = iskymask(spec,sigspec,wave,mask=,nsig=,$
;          skywaves=,/doplot)
;
; INPUTS:
;	spec    - one dimensional spectrum to repair
;       sigspec - corresponding error spectrum
;	wave    - corresponding wavelength vector
;
; OPTIONAL INPUTS:
;	mask     - corresponding bad pixel mask (updated on output) 
;       skywaves - wavelengths of night sky lines to repair
;       nsig     - flag NSIG outliers (default 2.0)
;	
; KEYWORD PARAMETERS:
;       doplot - generate a plot of the masked lines
;
; OUTPUTS:
;       newspec - cleaned one-dimensional spectrum
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The adopted algorithm interpolates over masked pixels by
;       fitting a line around the masked regions and adding artificial
;       noise.  The masked pixels are added to the bad pixel mask and
;       assigned a bit number according to IMASK_BITS(). 
;
; EXAMPLE:
;
; PROCEDURES USED:
;	GET_ELEMENT, DJS_MASKINTERP(), DJS_PLOT, DJS_OPLOT, LEGEND,
;	IM_POLY_ITER, IMASK_BITS() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 April 23, U of A
;       jm02oct17uofa - added noise to the bad sky pixels for visual
;                       purposes 
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm05jun22uofa - slight cosmetic improvements; use IM_WINDOW to
;                       launch a resolution independent window; only
;                       mask 5577 and 6300 
;       jm05aug02uofa - added NSIG and SKYWAVES optional inputs 
;
; Copyright (C) 2002-2003, 2005, John Moustakas
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

function iskymask, spec, sigspec, wave, mask=mask, skywaves=skywaves, $
  nsig=nsig, doplot=doplot

    if n_elements(spec) eq 0L then begin
       print, 'Syntax - newspec = iskymask(spec,sigspec,wave,mask=,$'
       print, '   skywaves=,nsig=,/doplot)'
       return, -1
    endif

    npix = n_elements(spec)

    if (n_elements(nsig) eq 0L) then nsig = 2.0
    
    if (n_elements(skywaves) eq 0L) then begin
       skywaves = [5577.345,6300.32] ; sky wavelengths
;      skywaves = [5577.345,5769.5982,5889.950,6300.32,6363.81]
    endif
       
    inrange = where((skywaves gt min(wave)) and (skywaves lt max(wave)),nrange)
    if nrange ne 0L then skywaves = skywaves[inrange] else return, spec
    nskywaves = n_elements(skywaves)

    dwave = median(wave-shift(wave,1)) ; median dispersion
    
; technique 1 - uses statistics to replace bad pixels

;    newspec = spec
;    for i = 0L, nskywaves-1L do begin
;
;       skymask = bytarr(npix)   ; bad pixel mask
;
;       get_element, wave, skywaves[i]+[-50,+50], x12
;       djs_iterstat, spec[x12[0]:x12[1]], sigma=specsig, median=specmed, sigrej=2.0
;
;       get_element, wave, skywaves, skypix
;       skymask[skypix] = 1B
;
;       skymask = smooth(float(skymask),7) gt 0B 
;       badpix = where(skymask eq 1B,nbadpix,comp=goodpix,ncomp=ngoodpix)
;
;       for j = 0L, nbadpix-1L do newspec[badpix[j]] = specmed+randomn(seed,1)*3*specsig ; perturb by 3-sigma normally
;       
;    endfor

; technique 2 - doesn't work well with just linear interpolation
;               around 5577
    
    skymask = intarr(npix) ; bad pixel mask
;   get_element, wave, skywaves, skypix
;   skymask[skypix] = 1B

    newspec = spec
    
; grow the mask by 2 pixels on either side
    
;   skymask = smooth(float(skymask),7) gt 0B 
;   badpix = where(skymask eq 1B,nbadpix,comp=goodpix,ncomp=ngoodpix)

; technique 2a:
    
;   newspec = djs_maskinterp(spec,skymask,wave)

; technique 2b: iteratively fit a line around the contaminated pixels
; and replace with a simple model+noise (2.5-sigma)

    for i = 0L, nskywaves-1L do begin

       get_element, wave, skywaves[i]+dwave*10.0*[-1.0,+1.0], x12
       x = wave[x12[0]:x12[1]]
       y = spec[x12[0]:x12[1]]
       sig = sigspec[x12[0]:x12[1]]
;      m = skymask[x12[0]:x12[1]]

       if (n_elements(x) gt 3L) then begin

          djs_iterstat, y, sigrej=2.0, mask=good
          w = where(good,nw)
          if (nw gt 5L) then $
            coeff = robust_poly_fit(x[w],y[w],1,yfit,nsig) else $
            coeff = robust_poly_fit(x,y,1,yfit,nsig)
            
          yfit = poly(x,coeff)
          djs_iterstat, y-yfit, sigrej=nsig, mask=good
          
;         im_poly_iter, x, y, 1, nsig, yfit, coeff=coeff, good=good

          replace = where(good eq 0B,nreplace)
          
          if (nreplace ne 0L) then begin
             skymask[x12[0]+replace] = imask_bits(/skymask)
             newspec[x12[0]+replace] = yfit[replace] + randomn(seed,nreplace)*1.0*djs_median(sig[good])
          endif

       endif
          
;      plot, x, y, ps=10, xsty=3, ysty=3, charsize=1.5, charthick=2.0, $
;        xthick=2.0, ythick=2.0, thick=2.0, xtitle='Wavelength', $
;        ytitle='Flux/Counts', xmargin=[12,3]
;      djs_oplot, x, yfit, line=0, thick=2.0, color='dark green'
;      djs_oplot, wave, newspec, ps=10, color='red'
;      djs_oplot, skywaves[i]*[1,1], !y.crange, line=2, thick=2
;      cc = get_kbrd(1)
       
;      coeff = func_fit(x,y,2,invvar=(x*0.0+1.0)*(m eq 0B),yfit=yfit)
;      w = where(m eq 1B,nw)

    endfor

    goodpix = where(skymask eq fix(0),ngoodpix,comp=badpix,ncomp=nbadpix)
    
    if keyword_set(doplot) then begin

       inrange = where((wave gt min(skywaves)-30*dwave) and $
         (wave lt max(skywaves)+30*dwave))
       xrange = minmax(wave[inrange])
;      yrange = median(spec[inrange])+djsig(spec[inrange],sigrej=3.0)*[-3.0,5.0]
       yrange = minmax(spec[inrange])

       norm = max(yrange)
       
       im_window, 3, xratio=0.25, /square
       djs_plot, wave, spec/norm, ps=10, xsty=3, ysty=3, xrange=xrange, $
         yrange=yrange/norm, charsize=1.2, charthick=2.0, xthick=2.0, $
         ythick=2.0, xtitle='Wavelength (\AA)', ytitle='Relative Flux'
       if (nbadpix ne 0L) then djs_oplot, wave[badpix], spec[badpix]/norm, $
         ps=2, syms=2, color='light blue'
       djs_oplot, wave, newspec/norm, ps=10, color='green'
       legend, ['Masked pixels'], /right, /top, /box, charsize=1.3, $
         charthick=2.0, color=djs_icolor('light blue'), psym=2, thick=2.0
       
    endif
       
; identify these replaced pixels as sky-subtraction residuals in the
; bad pixel mask

    if (n_elements(mask) ne 0L) then if (nbadpix ne 0L) then $
      mask[badpix] = mask[badpix] + fix(2^5) ; 0 is good
    
return, newspec
end


