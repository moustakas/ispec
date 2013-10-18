;+
; NAME:
;       I2DSTARFIT()
;
; PURPOSE:
;       Model and subtract a point source (e.g., a star) from a 2D
;       spectrum.  For example, this routine can be used to subtract a
;       contaminating foreground star from a extended background
;       source like a galaxy.
;
; CALLING SEQUENCE:
;       starimage = i2dstarfit(image,sigmamap,star_refrow=,$
;          frac_rows=,frac_cols=,fit_fraction=,trace_searchbox=,$
;          ngauss_terms=,center_order=,sigma_order=,/debug,/doplot)
;
; INPUTS:
;	image    - input 2D spectrum
;       sigmamap - corresponding 2D error map
;       header   - corresponding FITS header
;
; OPTIONAL INPUTS:
;       star_refrow     - reference row to help center the trace; this 
;                         input may be an array, indicating several
;                         stars should be subtracted [NSTAR]
;       frac_rows       - fraction of (spatial) rows to include in the 
;                         Gaussian fit, centered on the trace center
;                         at each position (default 0.1)
;       frac_cols       - fraction of (dispersion) columns to average 
;                         before Gaussian fitting (default 0.02) 
;       fit_fraction    - generate average spatial profiles every  
;                         FIT_FRACTION columns (default 0.2)
;       trace_searchbox - see SBOX in ITRACESPEC (default 5)
;       ngauss_terms    - number of terms in the Gaussian function
;                         (default 5 = Gaussian+line)
;       center_order    - order of the polynomial fit to the Gaussian
;                         center (default 2)
;       sigma_order     - order of the polynomial fit to the Gaussian
;                         sigma (default 2)
;
; KEYWORD PARAMETERS:
;       moffat - use a MOFFAT profile rather than a Gaussian for more
;                flexibility; NGAUSS_TERMS is set to be at least 6 in
;                this case
;       debug  - make a plot for debugging
;       doplot - generate a QA plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SPLOG, ITRACESPEC, MPFITPEAK(), DJS_PLOT, DJS_OPLOT,
;       IM_WINDOW, PAGEMAKER
;
; COMMENTS:
;       Two iterations are required to properly model the point
;       source.  First, we generate average spatial profiles with
;       spacing equal to FIT_FRACTION.  Each spatial profile is the
;       average of 2*FRAC_COLS columns and a cutout of FRAC_ROWS
;       spatial rows.  Gaussian profiles are fitted and low-order
;       functions are fitted to the resulting center and sigma-width,
;       which are assumed to be slowly varying functions of
;       wavelength.  Next, a Gaussian plus line are fitted to every
;       column, independently of every other column, to determine the
;       peak amplitude.  On this second iteration the sigma width and
;       center are fixed.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 20, U of A
;       jm05jul25uofa - added CENTER_ORDER and SIGMA_ORDER optional
;                       inputs 
;       jm07jun20nyu  - added MOFFAT keyword
;
; Copyright (C) 2005, 2007, John Moustakas
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

function i2dstarfit, image, sigmamap, header, star_refrow=star_refrow, $
  frac_rows=frac_rows, frac_cols=frac_cols, fit_fraction=fit_fraction, $
  trace_searchbox=trace_searchbox, ngauss_terms=ngauss_terms, $
  center_order=center_order, sigma_order=sigma_order, moffat=moffat, $
  debug=debug, doplot=doplot, _extra=extra

    if (n_elements(image) eq 0L) or (n_elements(sigmamap) eq 0L) then begin
       print, 'Syntax - starimage = i2dstarfit(image,sigmamap,header,$'
       print, '   star_refrow=,frac_rows=,frac_cols=,fit_fraction=,$'
       print, '   trace_searchbox=,ngauss_terms=,center_order=,$'
       print, '   sigma_order=,/moffat,/debug,/doplot)'
       return, image
    endif
    
; call this routine recursively     
    
    nstar = n_elements(star_refrow)
    if (nstar gt 1L) then begin
       image1 = image
       starimage = image*0.0
       for istar = 0L, nstar-1L do begin
          starimage1 = i2dstarfit(image1,sigmamap,header,star_refrow=star_refrow[istar],$
            frac_rows=frac_rows,frac_cols=frac_cols,fit_fraction=fit_fraction,$
            trace_searchbox=trace_searchbox,ngauss_terms=ngauss_terms,$
            center_order=center_order, sigma_order=sigma_order,moffat=moffat,$
            debug=debug,doplot=doplot)
          image1 = image - starimage1
          starimage = starimage + starimage1
       endfor
       return, starimage
    endif

    if (n_elements(frac_rows) eq 0L) then frac_rows = 0.1
    if (n_elements(frac_cols) eq 0L) then frac_cols = 0.02
    if (n_elements(fit_fraction) eq 0L) then fit_fraction = 0.2
    if (n_elements(trace_searchbox) eq 0L) then trace_searchbox = 5L
    if (n_elements(ngauss_terms) eq 0L) then ngauss_terms = 5L ; Gaussian + line
    if (n_elements(center_order) eq 0L) then center_order = 2L
    if (n_elements(sigma_order) eq 0L) then sigma_order = 2L

    gaussian = 1L
    if keyword_set(moffat) then begin
       ngauss_terms = ngauss_terms > 6L ; NOTE!
       gaussian = 0L
    endif

    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

    nrows_med = long(frac_rows*nrows)
    ncols_med = long(frac_cols*ncols)

    plotscale = 1E15

; construct the initial trace    

    itracespec, image, header, refrow=star_refrow, refcol=refcol, $
      traceinfo=traceinfo, sbox=trace_searchbox, doplot=debug, _extra=extra

    tracecenter = long(traceinfo.trace)

    splog, 'Fitting star near row '+$
      string(tracecenter[ncols/2L],format='(I0.0)')+'.'

; construct an average spatial profile every FIT_FRACTION columns 
    
    nsample = long(fit_fraction*ncols)>1L
    colspace = ncols/nsample
    colpos = fix(findgen(nsample)*colspace+colspace/2.0)

; initialize the first-iteration and final fitting results
    
    starinfo1 = {$
      colpos:       0L, $
      center:      0.0, $
      peak:        0.0, $
      sigma:       0.0}
    starinfo1 = replicate(starinfo1,nsample)
    starinfo1.colpos = colpos

    starinfo = {$
      colaxis:     0.0, $
      peak:        0.0, $
      peak_err:    0.0, $
      peak_fixed:  0.0, $
      fitsigma:    0.0, $
      fitcenter:   0.0}
    starinfo = replicate(starinfo,ncols)
    starinfo.colaxis = colaxis

; loop on each average position

    parinfo = replicate({value: 1.0, limited: [0,0], limits: [0.0,0.0]},ngauss_terms)

    parinfo[0].value = max(image[ncols/2.0,*]) ; initial guess at the peak amplitude
;   parinfo[0].value = 1.0/plotscale
    parinfo[1].value = nrows_med     ; center

    for isample = 0L, nsample-1L do begin

; zoom in on the "star"       
       
       minrow = (tracecenter[colpos[isample]]-nrows_med)>0L
       maxrow = (tracecenter[colpos[isample]]+nrows_med)<(nrows-1L)

; generate the average profile
       
       mincol = (colpos[isample]-ncols_med)>0L
       maxcol = (colpos[isample]+ncols_med)<(ncols-1L)

       subimage = image[mincol:maxcol,minrow:maxrow]
       suberror = sigmamap[mincol:maxcol,minrow:maxrow]
       naverage = (size(subimage,/dimension))[0]

       starprofile = total(subimage,1)/float(naverage)
       stareprofile = sqrt(total(suberror^2,1)/float(naverage))
       gaussaxis = findgen(n_elements(starprofile))

       gfit = mpfitpeak(gaussaxis,starprofile,a,error=stareprofile,$
         nterms=ngauss_terms,/positive,moffat=moffat,gaussian=gaussian)

       starinfo1[isample].center = a[1]-nrows_med+tracecenter[colpos[isample]]
       starinfo1[isample].peak = a[0] > 0.0 ; NOTE!
       starinfo1[isample].sigma = a[2]

;      if keyword_set(debug) then begin
;         if (isample eq 0L) then im_window, 5, xratio=0.4, /square
;         djs_plot, gaussaxis, starprofile, xsty=3, ysty=3, ps=10, $
;           xthick=2.0, ythick=2.0, thick=2.0, charsize=1.3, charthick=2.0, $
;           xtitle='Row', ytitle='Flux', title='Column '+string(colpos[isample],format='(I0)')
;         djs_oplot, gaussaxis, gfit, color='red', ps=10
;         wait, 0.1
;;        cc = get_kbrd(1)
;      endif
       
    endfor

; fit low-order functions to the Gaussian centroid and sigma width 

    nsig = 2.5
    
;   center_coeff = robust_poly_fit(starinfo1.colpos,starinfo1.center,$
;     center_order,fitcenter1,nsig)
    im_poly_iter, starinfo1.colpos,starinfo1.center, center_order, fitcenter1, $
      coeff=center_coeff, nsig=nsig
    starinfo.fitcenter = poly(starinfo.colaxis,center_coeff)

;   sigma_coeff = robust_poly_fit(starinfo1.colpos,starinfo1.sigma,$
;     sigma_order,fitsigma1,nsig)
    im_poly_iter, starinfo1.colpos,starinfo1.sigma, sigma_order, fitsigma1, $
      coeff=sigma_coeff, nsig=nsig
    starinfo.fitsigma = poly(starinfo.colaxis,sigma_coeff)

; now that the sigma-width and center are well-constrained, loop on
; every column and fit just to the peak amplitude (plus linear term,
; to account for the underlying galaxy)

    parinfo = replicate({value: 1.0, fixed: 0B},ngauss_terms)
;   parinfo[0].value = 1.0/plotscale
    parinfo[0].value = max(image[ncols/2.0,*])    
    for icol = 0L, ncols-1L do begin

       minrow = (tracecenter[icol]-nrows_med)>0L
       maxrow = (tracecenter[icol]+nrows_med)<(nrows-1L)

       starprofile = djs_median(image[(icol-2L)>0L:(icol+2L)<(ncols-1L),minrow:maxrow],1)
;      starprofile = image[icol,minrow:maxrow]
       stareprofile = sigmamap[icol,minrow:maxrow]
       gaussaxis = findgen(n_elements(starprofile))

; constrained Gaussian fit       

       parinfo[1L:2L].fixed = 1B
       
       parinfo[1].value = starinfo[icol].fitcenter-tracecenter[icol]+nrows_med
       parinfo[2].value = starinfo[icol].fitsigma
       
       gfit = mpfitpeak(gaussaxis,starprofile,a,error=stareprofile,$
         nterms=ngauss_terms,/positive,parinfo=parinfo,$
         perror=perror,estimates=parinfo.value,moffat=moffat,gaussian=gaussian)

       starinfo[icol].peak = a[0]
       starinfo[icol].peak_err = perror[0]
;      print, a[1]-nrows_med+tracecenter[icol], starinfo[icol].fitcenter, $
;        a[2], starinfo[icol].fitsigma

       if keyword_set(debug) then begin

          if (icol eq 0L) then im_window, 7, xratio=0.4, /square
          djs_plot, gaussaxis, starprofile, xsty=3, ysty=3, ps=10, $
            xthick=2.0, ythick=2.0, thick=2.0, charsize=1.3, charthick=2.0, $
            xtitle='Row', ytitle='Flux', title='Column '+string(icol,format='(I0)')
          djs_oplot, gaussaxis, gfit, color='red', ps=10
          wait, 0.05
;         cc = get_kbrd(1)

       endif
       
    endfor

; to ensure continuity interpolate over 3-sigma "peak flux" outliers

    medpeak = median(starinfo.peak,10)
    djs_iterstat, starinfo.peak-medpeak, sigrej=4.0, mask=peak_mask

    peak_outliers = where(peak_mask eq 0L,noutliers)
    
    peak_fixed = djs_maskinterp(starinfo.peak,(peak_mask eq 0L),colaxis)
    starinfo.peak_fixed = peak_fixed
    
; now generate the 2D image by evaluating the appropriate Gaussian at
; each column

    starimage = image*0.0
    for icol = 0L, ncols-1L do starimage[icol,*] = starinfo[icol].peak_fixed * $
      exp(-0.5*(rowaxis-starinfo[icol].fitcenter)^2/starinfo[icol].fitsigma^2)
    
; generate a QA plot    
    
    if keyword_set(doplot) then begin

       im_window, 4, xratio=0.7, yratio=0.6
       plotsym, 0, 0.5
       
       pagemaker, nx=1, ny=3, xmargin=[1.1,0.3], ymargin=[0.3,1.3], $
         xspace=0.0, yspace=0.0, position=pos, /normal

       yrange = plotscale*(median(starinfo.peak_fixed)+djsig(starinfo.peak_fixed,sigrej=2.5)*[-5,5])
;      yrange = plotscale*minmax(starinfo.peak_fixed)       
       djs_plot, starinfo.colaxis, plotscale*starinfo.peak_fixed, xsty=3, ysty=3, ps=8, $
         xthick=2.0, ythick=2.0, thick=2.0, charsize=1.3, charthick=2.0, $
         xtickname=replicate(' ',10), xtitle='', ytitle='Relative Peak Flux', $
         position=pos[*,0], yminor=3, yrange=yrange
       if (noutliers ne 0L) then $
         djs_oplot, starinfo[peak_outliers].colaxis, plotscale*starinfo[peak_outliers].peak, $
           psym=7, syms=2.0, color=djs_icolor('light blue')

       yrange = minmax(starinfo.fitcenter)
       djs_plot, starinfo1.colpos, starinfo1.center, xsty=3, ysty=3, ps=8, $
         xthick=2.0, ythick=2.0, thick=2.0, charsize=1.3, charthick=2.0, $
         xtickname=replicate(' ',10), xtitle='', ytitle='Center [pixel]', $
         position=pos[*,1], yminor=3, /noerase, yrange=yrange
       djs_oplot, starinfo.colaxis, starinfo.fitcenter, line=0, thick=2.0, color='green'

       yrange = minmax(starinfo.fitsigma)
       djs_plot, starinfo1.colpos, starinfo1.sigma, xsty=3, ysty=3, ps=8, $
         xthick=2.0, ythick=2.0, thick=2.0, charsize=1.3, charthick=2.0, $
         xtitle='Column Number', ytitle='Sigma [pixel]', $
         position=pos[*,2], /noerase, yminor=3, yrange=yrange
       djs_oplot, starinfo.colaxis, starinfo.fitsigma, line=0, thick=2.0, color='red'
              
    endif

return, starimage
end
    
