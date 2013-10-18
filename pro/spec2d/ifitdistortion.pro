;+
; NAME:
;	IFITDISTORTION
;
; PURPOSE:
;	Model the spatial distortion in two-dimensional spectra by
;	tracing the bright spectral features (such as the continua of
;	standard stars). 
;
; CALLING SEQUENCE:
;       ifitdistortion, speclist, datapath=, npad=, tracename=, $
;          tinfo=, finaltrace=, /doplot, _extra=extra
;
; INPUTS:
;	speclist - list of images from which to model the distortion 
;
; OPTIONAL INPUTS:
;	datapath  - path to the data and the output trace structure 
;	npad      - pad the distortion map by 2*NPAD pixels in the
;                   spatial dimension
;	tracename - file name of the output data structure containing
;                   the trace model, e.g. trace_98mar21.idlsave (if
;                   set then the trace structure is written, otherwise
;                   the model is simply plotted) 
;	extra     - keywords for ITRACESPEC
;
; OUTPUTS:
;	finaltrace - final trace structure (written if TRACENAME is
;                    defined) 
;
; OPTIONAL OUTPUTS:
;	tinfo      - trace data structure (fitting parameters)
;
; KEYWORD PARAMETERS: 
;       doplot     - generate a plot of the traces
;
; COMMENTS:
;	This routine takes a list of two-dimensional spectra and
;	traces the spatial profile as a function of column number or 
;	wavelength.  The order of the trace is arbitrary.  The traces
;	are plotted and the best fit coefficients are determined
;	statistically (iterative rejection + the objects at the
;	highest airmass are rejected).  A padded two-dimensional
;	distortion map is also created.
;
; TODO:  
;       Trace several apertures on the CCD (e.g., if a standard
;	star was observed at multiple on the CCD during a single
;	exposure). 
;
;       Generalize the distortion map-making to incorporate a varying
;       amount of distortion with spatial position.
;    
; PROCEDURES USED:
;	CWD(), DJS_ITERSTAT(), RD2DSPEC(), ITRACESPEC, CMSAVE,
;	STRN(), DJS_PLOT, DJS_OPLOT, SPLOG, ICLEANUP
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 December 27, U of A
;	jm02jan29uofa - modified and documented
;	jm02feb19uofa - adopted from IFITTRACE; generalized to an
;                       arbitrary order distortion 
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm07jun20nyu  - improved the debugging plot; do not reject any
;                       stars based on airmass; let the user
;                       customized the TRACELIST file
;
; Copyright (C) 2001-2003, 2007, John Moustakas
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

pro ifitdistortion, speclist, datapath=datapath, npad=npad, tracename=tracename, $
  tinfo=tinfo, finaltrace=finaltrace, doplot=doplot, _extra=extra

    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - ifitdistortion, speclist, datapath=, npad=, tracename=, $'
       print, '   tinfo=, finaltrace=, /doplot, _extra=extra'
       return
    endif
    
    if not keyword_set(datapath) then datapath = cwd()
    if not keyword_set(npad) then npad = 10L ; pixels
    
; read the trace list and trace each image
    
    tcube = rd2dspec(speclist,datapath=datapath)

    ncols = tcube[0].naxis1
    nrows = tcube[0].naxis2
    refcol = ncols/2.0
    colaxis = findgen(ncols)
    
    for i = 0L, nspec-1L do begin
       image = tcube[i].image
       header = *tcube[i].header
       itracespec, image, header, _extra=extra, traceinfo=info, doplot=0
       if (i eq 0L) then tinfo = info else tinfo = [ [tinfo], [info] ]
;      cc = get_kbrd(1)       
    endfor
    tinfo = reform(tinfo)

; sort by airmass and reject the trace at the highest airmass.  but
; always use at least two standards; jm07jun20nyu - DON'T REJECT ANY! 
    
;   sortair = sort(tinfo.airmass) 
;   tinfo = tinfo[sortair]
;   if nspec gt 2L then begin
;      tinfo = tinfo[0:nspec-2L]
;      nspec = nspec-1L
;   endif

; now go through each trac and normalize along the middle column    
    
    normtrace = fltarr(nspec)
    for it = 0L, nspec-1L do normtrace[it] = interpol(tinfo[it].trace,colaxis,refcol)
    for it = 0L, nspec-1L do begin
       tinfo[it].trace = tinfo[it].trace-normtrace[it]+djs_median(normtrace)
       tinfo[it].trace_data = tinfo[it].trace_data-normtrace[it]+djs_median(normtrace)
    endfor

; OLD CODE: iteratively derive the best coefficients

;   ncoeff = tinfo[0].ncoeff
;   coeff = fltarr(ncoeff)
;   for j = 0L, ncoeff-1L do begin
;      djs_iterstat, tinfo.coeff[j], sigrej=3.0, median=med
;      coeff[j] = med
;   endfor
;
;   trace = poly(colaxis,coeff) ; best trace

; NEW CODE: derive the best coefficients by fitting to all the traces
; simultaneously!

    im_poly_iter, rebin(reform(colaxis,ncols,1),ncols,nspec)-refcol, $
      tinfo.trace_data, tinfo[0].ncoeff-1.0, sigma=tracefit_sig, $
      coeff=coeff, nsig=5.0
;   coeff = robust_poly_fit(rebin(reform(colaxis,ncols,1),ncols,nspec)-refcol,$
;     tinfo.trace_data,tinfo[0].ncoeff-1.0,junk,tracefit_sig) ; TRACEORDER
;   coeff[0] = coeff[0] - coeff[1]*refcol

    trace = poly(colaxis-refcol,coeff) ; best trace

; generate the distortion map.  we initialize the map by computing
; indices from 0 to NROWS-1 for every column, then modify the indices
; by the distortion modeled by the best trace.  this distortion map
; can then be rectified by TRIANGULATE and TRIGRID (see ICALIBRATE).
; pad the map by NPAD rows on either end of the CCD; pivot about the
; middle column
    
; TODO:  generalize to incorporate a varying amount of distortion with
; spatial position
    
    pivot = interpol(trace,colaxis,refcol)
    tcurve = trace - pivot
    
    dmap = (fltarr(ncols)+1.0) # (findgen(nrows+2*npad)-npad)
    for k = 0L, 2*npad+nrows-1L do dmap[*,k] = dmap[*,k] - tcurve
    
; --------------------------------------------------
; test that the distortion can be removed

;   xyccd = findgen(ncols) # (fltarr(nrows+2*npad)+1) ; column number
;   triangulate, xyccd, dmap, tr
;   goodmap = trigrid(xyccd,dmap,dmap,tr,[1,1],[0,-npad,ncols-1,npad+nrows-1])
;   goodmap = goodmap[*,npad:npad+nrows-1]
;   plot, dmap[*,npad+60], xsty=3, ysty=3 & oplot, goodmap[*,60], line=2, thick=2
    
; --------------------------------------------------
    
; convert the trace slope to an clockwise angle.  the angle should be
; in the range [-90,90]

;   deltay = double(meantrace[npix-1L]-meantrace[0])
;   deltax = double(npix)
;   theta = atan(deltay/deltax)*!radeg ; [degrees]
    
; plot the individual traces and the best trace

    if keyword_set(doplot) then begin
    
       im_window, 0, xratio=0.45, /square
       pagemaker, nx=1, ny=2, position=pos, xmargin=[1.1,0.1], $
         ymargin=[0.8,1.0], yspace=0.0, /normal

       djs_plot, [0], [0], /nodata, xrange=minmax(colaxis), yrange=minmax(tinfo.trace), xsty=3, ysty=3, $
         xthick=2.0, ythick=2.0, charsize=1.6, charthick=2.0, xtitle='', xtickname=replicate(' ',10), $
         ytitle='Row Number', title='Spectral Tracing', position=pos[*,0]
       for j = 0L, nspec-1L do oplot, colaxis, tinfo[j].trace_data, ps=3
;      for j = 0L, nspec-1L do oplot, colaxis, tinfo[j].trace, line=2, thick=2.0
       djs_oplot, colaxis, trace, line=0, thick=2.0, color='green'

       djs_plot, [0], [0], /nodata, xrange=minmax(colaxis), yrange=tracefit_sig*[-5,5], $
         xsty=3, ysty=3, line=2, thick=2.0, xthick=2.0, ythick=2.0, charsize=1.6, $
         charthick=2.0, xtitle='Column Number', ytitle='Residuals (pixel)', /noerase, $
         position=pos[*,1]
       for j = 0L, nspec-1L do djs_oplot, colaxis, tinfo[j].trace_data-trace, ps=3;, color='red'
;      for j = 0L, nspec-1L do djs_oplot, colaxis, tinfo[j].trace-trace, line=2, thick=2.0
       djs_oplot, !x.crange, [0,0], color='green', thick=1.0

       for j = 0L, nspec-1L do print, speclist[j]+'  ', tinfo[j].airmass, $
         min(tinfo[j].trace-trace), max(tinfo[j].trace-trace) ;tinfo[j].coeff, $

    endif
       
; write the best trace and the distortion map if there were no errors 

    finaltrace = create_struct('speclist', speclist, 'datapath', datapath, $
      'trace', trace, 'ncols', ncols, 'nrows', nrows, $
      'coeff', coeff, 'npad', npad, 'dmap', dmap)
    
    if (n_elements(tracename) ne 0L) and (size(tracename,/type) eq 7L) and $
      (strcompress(tracename,/remove) ne '') then begin
       print, 'Writing '+datapath+tracename+'.'
       cmsave, filename=datapath+tracename, finaltrace, names=['tracefit']
    endif else begin
       splog, 'Invalid TRACENAME file '+tracename+'.'
    endelse
    
    icleanup, tcube

return
end
