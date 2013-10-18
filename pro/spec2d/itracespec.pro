;+
; NAME:
;	ITRACESPEC
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
;       tracetype    - trace the spectrum using either 'default' or
;                      'max' (but see GAUSS keyword, below, and the
;                      COMMENTS)
;	refrow       - center the trace on REFROW (default to the
;                      middle row)
;	refcol       - if REFROW is not set, then look for the
;                      reference row starting at this column (default
;                      to the middle column)
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
; PROCEDURES USED:
;	SXPAR(), SPLOG, TRACE_FWEIGHT(), TRACE_GWEIGHT(),
;	DJS_ITERSTAT, DJS_PLOT, DJS_OPLOT, LEGEND, IM_WINDOW,
;	ROBUST_POLY_FIT() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 24, U of A
;       jm02oct16uofa - cleaned up and improved documentation 
;       jm02oct23uofa - cleaned up the tracing options and added the
;                       GAUSS keyword
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03dec22uofa - added TRACE_SIGREJ keyword
;       jm05jun20uofa - call IM_WINDOW; improved QA plot
;       jm05jul21uofa - use ROBUST_POLY_FIT() in favor of POLY_FIT();
;                       added TRACE_MINCOL and TRACE_MAXCOL optional
;                       inputs; removed internal routine FIT_TRACE() 
;       jm07jun20nyu  - improved the debugging plot; switched to using
;                       IM_POLY_ITER rather than ROBUST_POLY_FIT,
;                       since the latter crashes when TRACEORDER gets
;                       large 
;
; Copyright (C) 2002-2003, 2005, 2007, John Moustakas
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

pro itracespec, image, header, tracetype=tracetype, refrow=refrow, $
  refcol=refcol, sbox=sbox, traceorder=traceorder, trace_sigrej=trace_sigrej, $
  trace_mincol=trace_mincol, trace_maxcol=trace_maxcol, traceinfo=traceinfo, $
  tracewindow=tracewindow, gauss=gauss, doplot=doplot

    if (n_elements(image) eq 0L) then begin
       print, 'Syntax - itracespec, image, header, tracetype=, refrow=, $'
       print, '   refcol=, sbox=, traceorder=, trace_sigrej=, trace_mincol=, $'
       print, '   traceinfo=, tracewindow=, /gauss, /doplot'
       return
    endif

    if (n_elements(trace_sigrej) eq 0L) then trace_sigrej = 3.0
    
    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    midrow = nrows/2L
    midcol = ncols/2L
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

    if (n_elements(trace_mincol) eq 0L) then trace_mincol = 0L
    if (n_elements(trace_maxcol) eq 0L) then trace_maxcol = ncols-1L

    if long(sxpar(header,'CRVAL1')) ne 0L then $ 
      wave = make_wave(header) else $ ; wavelength vector
      wave = colaxis

    pscale = sxpar(header,'CD2_2')
    
    tracemeth = 0L ; default trace method (max+flux-weighted centroiding)
    if (n_elements(traceorder) eq 0L) then traceorder = 1L
    if (traceorder le 0L) then traceorder = 1L
    
    if n_elements(tracetype) ne 0L then begin
       case tracetype of
          'default':            ; default (tracemeth = 0L)
          'max': begin
             tracemeth = 1L
             tcoeff = [-1.0]
             traceorder = -1
          end
          else: splog, 'Unknown trace type...selecting DEFAULT.'
       endcase
    endif
    
    if (n_elements(refrow) ne 0L) then begin
       
       if n_elements(sbox) eq 0L then sbox = 20L ; default search window size when tracing (radius)
       
       refrowcen = (0L > refrow) < (nrows-1L)
       y1 = (refrowcen-sbox)>0L
       y2 = (refrowcen+sbox)<(nrows-1L)
       
    endif else begin

       if n_elements(refcol) eq 0L then refcol = midcol
       
       y1 = 0L
       y2 = nrows-1L
;      refrow = nrows/2L 
       junk = max(image[refcol,*],refrow) ; make the reference row at REFCOL
       
    endelse

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
          im_poly_iter, xfit[cengood], yfit[cengood], traceorder, coeff=tcoeff, nsig=7.0
;         tcoeff = robust_poly_fit(xfit[cengood],yfit[cengood],traceorder)
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

       im_window, tracewindow, xratio=0.3, /square
;      window, tracewindow, xs=250, ys=250
;      window, tracewindow, xs=450, ys=450

       pagemaker, nx=1, ny=2, position=pos, xmargin=[1.2,1.2], $
         ymargin=[0.8,1.3], yspace=0.0, /normal

       stats = im_stats(cen,sigrej=3.0)
       djs_plot, wave, cen, ps=3, xsty=3, ysty=11, xtitle='', ytitle='Row Number', $
         xtickname=replicate(' ',10), $
         charthick=2.0, charsize=1.2, title=title, yrange=[stats.min,stats.max], $;cenmed+censig*[-5.0,5.0], $
         xthick=2.0, ythick=2.0, xmargin=[7,7], xcharsize=1.0, position=pos[*,0];, xticks=4
       axis, /yaxis, ythick=2.0, charthick=2.0, xsty=3, $
         charsize=1.0, ytitle='Deviation (arcsec)', yrange=(!y.crange-cen[midcol])*pscale
       djs_oplot, wave, center, line=0, thick=2.0, color='green'
       legend, 'Airmass = '+strtrim(string(airmass,format='(F12.3)'),2), $
         /right, /top, /box, charsize=1.3, charthick=2.0, thick=2.0

       resid_stats = im_stats(abs(cen-center),sigrej=3.0)
       djs_plot, wave, cen-center, ps=3, xsty=3, ysty=11, xtitle='Wavelength or Column Number', $
         ytitle='Residuals (pixel)', charthick=2.0, charsize=1.2, yrange=resid_stats.maxrej*[-1,1], $
         xthick=2.0, ythick=2.0, xmargin=[7,7], xcharsize=1.0, position=pos[*,1], /noerase;, xticks=4
       axis, /yaxis, ythick=2.0, charthick=2.0, xsty=3, $
         charsize=1.0, ytitle='Deviation (arcsec)', yrange=!y.crange*pscale
       djs_oplot, !x.crange, [0,0], color='green'

;      splog, 'Press any key to continue.' & cc = get_kbrd(1)
;      stay = '' & read, stay, prompt='ITRACESPEC: Press ENTER to continue: '; & wdelete, 1
       
    endif
    
    traceinfo = create_struct('npix', ncols, 'airmass', airmass, 'trace_data', cen, $
      'trace', center, 'ncoeff', traceorder+1, 'coeff', reform(tcoeff))

return
end
