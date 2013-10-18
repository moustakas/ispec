;+
; NAME:
;	ITRACEORDERS
;
; PURPOSE:
;	Trace the spatial centroid of a spectrum as a function of
;	wavelength or column number. 
;
; CALLING SEQUENCE:
;       itracespec, image, header, tracetype=, refrow=, refcol=, $
;          sbox=, traceorder=, trace_sigrej=, trace_mincol=, $
;          trace_maxcol=, traceinfo=, tracewindow=, /gauss, /doplot
;
; INPUTS:
;	image       - input two-dimensional spectrum
;	header      - corresponding FITS image header 
;
; OPTIONAL INPUTS:
;	refcol       - identify orders starting from this reference
;                      column (default to the middle column)
;       norder       - number of echellete orders to identify (default
;                      20) 
;       order_sigma  - significance threshold for each order (default
;                      5.0) 
;
;	sbox         - if REFROW is given, then only trace within
;                      plus-or-minus SBOX pixels (default 20)
;	traceorder   - order of the polynomial fit to the trace for
;                      'default' (default linear)
;	trace_sigrej - rejection threshold for points when fitting the
;                      trace 
;	trace_mincol - fit the trace starting at this minimum (first)
;                      column (default 0)
;	trace_maxcol - fit the trace ending at this minimum (last)
;                      column (default NCOLS-1)
;
; KEYWORD PARAMETERS:
;	gauss       - instead of flux-weighted centroiding i 'default'
;                     use gaussian-weighted centroiding
;	doplot      - plot the trace and the fit
;
; OPTIONAL OUTPUTS:
;	traceinfo   - structure containing the functional form of the
;                     trace
;	tracewindow - window number where the trace is drawn
;
; COMMENTS:
;	The default tracing method uses max-weighted centering as an
;	initial guess and flux-weighted centroiding as the final
;	trace.  If max-weighted centroiding is requested then no
;	function is fitted to the trace and the coefficient of the fit
;	(TCOEFF) is -1.
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2007 January 11, NYU - written
;
; Copyright (C) 2007, John Moustakas
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

pro itraceorders, image, refcol=refcol, norder=norder, order_sigma=order_sigma, $
  doplot=doplot

;  sbox=sbox, traceorder=traceorder, trace_sigrej=trace_sigrej, $
;  trace_mincol=trace_mincol, trace_maxcol=trace_maxcol, traceinfo=traceinfo, $
;  tracewindow=tracewindow, gauss=gauss, doplot=doplot

    if (n_elements(image) eq 0L) then begin
       doc_library, 'itraceorders'
       return
    endif

    imsize = size(image,/dimension)
    ncols = imsize[0] & nrows = imsize[1]
    midrow = nrows/2L & midcol = ncols/2L
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

    if (n_elements(trace_sigrej) eq 0L) then trace_sigrej = 3.0
    if (n_elements(trace_mincol) eq 0L) then trace_mincol = 0L
    if (n_elements(trace_maxcol) eq 0L) then trace_maxcol = ncols-1L
    if (n_elements(traceorder) eq 0L) then traceorder = 3L
    if (n_elements(refcol) eq 0L) then refcol = midcol
    if (n_elements(norder) eq 0L) then norder = 20L

; identify orders using peak-finding

    splog, 'Reference column '+string(refcol,format='(I0)')+'.'
    
    xpeak = find_npeaks(image[refcol,*],colaxis,nfind=norder,ypeak=ypeak,npeak=npeak)
    good = where(ypeak gt median(image[refcol,*])+5.0*djsig(image[refcol,*]),ngood)
    splog, 'Identified '+string(ngood,format='(I0)')+'/'+string(norder,format='(I0)')+$
      ' significant orders.'
    if (ngood eq 0L) then return else begin
       xpeak = xpeak[good]
       ypeak = ypeak[good]
       npeak = ngood
    endelse

    if keyword_set(doplot) then begin
       im_window, 0, yratio=0.6, xratio=0.7
       djs_plot, image[refcol,*], xsty=3, ysty=3, xtitle='Row', ytitle='Intensity', $
         charsize=1.8, charthick=2.0, xthick=2.0, ythick=2.0
       for ii = 0, npeak-1L do djs_oplot, xpeak[ii]*[1,1], !y.crange, color='red'
    endif

    
    
stop    
    
    y1 = 0L
    y2 = nrows-1L

stop    

; initial guess: max-weighted centroiding    

    center = fltarr(ncols)
    for j = 0L, ncols-1L do begin ; maximum
       ymax = max(image[j,y1:y2],cenmax)
       center[j] = cenmax
    endfor
    
    cen = center

    case tracemeth of

       0L: begin

          if keyword_set(gauss) then $
            for i = 0L, 2L do cen = trace_gweight(transpose(image),cen,colaxis,sigma=sigma) else $
            for i = 0L, 2L do cen = trace_fweight(transpose(image[*,y1:y2]),cen,colaxis)
          cen = cen + y1 ; add the starting ROW back in

          djs_iterstat, cen, sigrej=trace_sigrej, sigma=censig, $
            median=cenmed, mask=cenmask
          cengood = where(cenmask eq 1B,ngood)

          if (ngood eq 0L) then message, 'No good points to fit to trace.'
          
          xfit = colaxis[trace_mincol>0L:trace_maxcol<(ncols-1L)]
          yfit = cen[trace_mincol>0L:trace_maxcol<(ncols-1L)]
          
;         tcoeff = poly_fit(xfit[cengood],yfit[cengood],traceorder)
          tcoeff = robust_poly_fit(xfit[cengood],yfit[cengood],traceorder)
          center = poly(colaxis,tcoeff)

       end

       1L: begin

          cen = cen+y1
          cenmed = median(cen)
          censig = 1.0
          
       end
          
    endcase

; plot the trace
    
    airmass = sxpar(header,'AIRMASS',count=nairmass)
    if (nairmass eq 0L) then airmass = 'Unknown' else $
      airmass = strtrim(string(airmass,format='(F12.3)'),2)
    if keyword_set(doplot) then begin
       
       title = sxpar(header,'OBJECT')

       tracewindow = 1

       im_window, tracewindow, xratio=0.25, /square
;      window, tracewindow, xs=250, ys=250
;      window, tracewindow, xs=450, ys=450
       djs_plot, wave, cen, ps=3, xsty=3, ysty=11, xtit='Wavelength or Column Number', ytit='Row Number', $
         charthick=2.0, charsize=1.2, title=title, yrange=cenmed+censig*[-5.0,5.0], $
         xthick=2.0, ythick=2.0, xmargin=[7,7], xcharsize=1.0;, xticks=4
       axis, /yaxis, ythick=2.0, charthick=2.0, xsty=3, $
         charsize=1.0, ytitle='Deviation (arcsec)', yrange=(!y.crange-cen[midcol])*pscale
       djs_oplot, wave, center, line=0, thick=2.0
       legend, 'Airmass = '+strtrim(string(airmass,format='(F12.3)'),2), $
         /left, /top, /box, charsize=1.3, charthick=2.0, thick=2.0

;      splog, 'Press any key to continue.' & cc = get_kbrd(1)
;      stay = '' & read, stay, prompt='ITRACESPEC: Press ENTER to continue: '; & wdelete, 1
       
    endif
    
    traceinfo = create_struct('npix', ncols, 'airmass', airmass, 'trace', $
      center, 'ncoeff', traceorder+1, 'coeff', reform(tcoeff))

return
end
    
