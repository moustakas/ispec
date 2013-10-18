;+
; NAME:
;	IILLUMFLAT()
;
; PURPOSE:
;	Generate an illumination correction flat field.
;
; CALLING SEQUENCE:
;       illumflat = iillumflat(simage,flat,illumfile,datapath=,$
;          illumname=,illumbins=,norder_illum=,norder_space=, $
;          illumheader=,/doplot,/write,/checkillum,_extra=extra)
;
; INPUTS:
;	simage    - twilight sky flat
;       flat      - dome flat field output from IMASTERFLAT()
;       illumfile - file name corresponding to SIMAGE
;
; OPTIONAL INPUTS:
;	datapath     - write the illumination flat and the QA plot to
;                      this path 
;       flatname     - name of the output illumination flat (default 
;                     'illumflat.fits')
;	illumbins    - number of columns to median when creating the
;                      illumination function (default 41L)
;	norder_illum - order of the fit to the illumination function
;                      at each "median" column (along rows) (default
;                      3)
;	norder_space - order of the fit in the spatial dimension
;                      (along columns) (default 2)
;       extra        - fitting/rejection parameters for
;                      BSPLINE_ITERFIT()
;
; KEYWORD PARAMETERS:
;       doplot     - generate a QA plot
;       write      - write ILLUMFLAT and the QA plot to DATAPATH 
;       checkillum - generate a plot of the illumination function fit
;                    at each "median" column
;
; OUTPUTS:
;	illumflat - two-dimensional illumination flat
;
; OPTIONAL OUTPUTS:
;       illumheader - corresponding header for ILLUMFLAT
;
; COMMENTS:
;
; PROCEDURES USED:
;	DJS_MEDIAN(), BSPLINE_ITERFIT(), FUNC_FIT(), MKHDR, LEGEND, 
;	DJS_ITERSTAT, ODD(), DJS_PLOT, DJS_OPLOT, SXADDPAR, SXADDHIST, 
;	SPLOG, DFPSPLOT, DFPSCLOSE, ISPEC_VERSION(), WRITEFITS,
;	PLOTHIST, POLY_ITER, DJS_ICOLOR(), TEXTOIDL(), IM_WINDOW
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August/October, U of A - written
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm03apr18uofa - added NORDER_SPACE parameter ISPEC v1.0.0
;       jm05jun17uofa - prettier histogram plot
;       jm05jun21uofa - use IM_WINDOW to spawn monitor resolution 
;                       independent windows
;       jm07jun19nyu  - use IM_POLY_ITER instead of POLY_ITER (better
;                       error catching)
;
; Copyright (C) 2001, 2003, 2005, 2007, John Moustakas
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

function iillumflat, simage, flat, illumfile, datapath=datapath, illumname=illumname, $
  illumbins=illumbins, norder_illum=norder_illum, norder_space=norder_space, $
  illumheader=illumheader, doplot=doplot, write=write, checkillum=checkillum, $
  _extra=extra

    if (n_elements(simage) eq 0L) or (n_elements(flat) eq 0L) then begin
       print, 'Syntax - illumflat = iillumflat(simage,flat,illumfile,datapath=,$'
       print, '   illumname=,illumbins=,norder_illum=,norder_space=, $'
       print, '   illumheader=,/doplot,/write,/checkillum,_extra=extra)'
       return, -1L
    endif

    if n_elements(illumname) eq 0L then illumname = 'illumflat.fits'
    if n_elements(norder_illum) eq 0L then norder_illum = 3L
    if n_elements(norder_space) eq 0L then norder_space = 2L
    if n_elements(illumbins) eq 0L then illumbins = 41L ; 11L

    imsize = size(simage)
    ncols = imsize[1]
    nrows = imsize[2]
    midrow = nrows/2L           ; middle row

    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

; flat-field the sky flat
    
    sflat = simage / flat
    
; use the sky image to construct the illumination correction.  take
; the median of ILLUMBINS columns and fit in the row dimension
    
    ntot = ncols/illumbins
    splog, 'Fitting '+strn(ntot)+' median columns.'

    if odd(illumbins) then midcol = illumbins/2L+1L else midcol = illumbins/2L
    
    flataxis = fltarr(ntot)
    medfit = fltarr(ntot,nrows)
 
    bkptres = (max(rowaxis)-min(rowaxis))/norder_illum ; break point resolution
    bkpt = findgen(norder_illum)*bkptres+min(rowaxis)

    if keyword_set(checkillum) then im_window, 3, xratio=0.3, /square
    for i = 0L, ntot-1L do begin

       skycols = lindgen(illumbins)+i*illumbins ; column indices
       medsky = djs_median(sflat[skycols,*],1)
       
       sset = bspline_iterfit(rowaxis,medsky,nord=3,yfit=fit,bkpt=bkpt,/silent)
;      result = func_fit(rowaxis,medsky,norder_illum+1,function_name='poly',yfit=fit)

       medfit[i,*] = fit / fit[midrow] ; normalize to unity at the center of the slit
       flataxis[i] = skycols[midcol]

       if keyword_set(checkillum) then begin

          norm = max(medsky)
          djs_plot, rowaxis, medsky/norm, xsty=3, ysty=3, ps=10, charthick=2.0, $
            charsize=1.3, xthick=2.0, ythick=2.0, xtitle='Row Number', $
            ytitle='Normalized Counts', thick=2.0, yrange=minmax(medsky/norm)*[1.0,1.01]
          djs_oplot, rowaxis, fit/norm, color='green', thick=3.0
          legend, 'Median column '+strn(i,format='(I0)')+'/'+strn(ntot,format='(I0)'), $
            /left, /top, box=0, charsize=1.3, charthick=2.0
          cc = get_kbrd(1)

       endif
       
    endfor

; generate the two dimensional illumination function by fitting a
; function with iterative rejection along the column direction at
; every row 
    
    illumflat = sflat*0.0
    for i = 0L, nrows-1L do begin

;      sset = bspline_iterfit(flataxis,reform(medfit[*,i]),nord=2L,yfit=yfit,bkpt=flataxis,/silent)
;      illumfit = bspline_valu(colaxis,sset)

       im_poly_iter, flataxis, reform(medfit[*,i]), norder_space, yfit, coeff=result, nsig=3.0
;      poly_iter, flataxis, reform(medfit[*,i]), norder_space, 3.0, yfit, coeff=result
;      result = func_fit(flataxis,reform(medfit[*,i]),3,function_name='poly',yfit=yfit)
       illumfit = poly(colaxis,result)
       illumflat[*,i] = illumfit

       if keyword_set(checkillum) then begin

          djs_plot, flataxis, medfit[*,i], xsty=3, ysty=3, charthick=2.0, $
            charsize=1.3, xthick=2.0, ythick=2.0, xtitle='Column Number', $
            ytitle='Normalized Counts', thick=2.0, ps=4
          djs_oplot, flataxis, yfit, color='green', thick=3.0
          legend, 'Row '+strn(i,format='(I0)')+'/'+strn(nrows,format='(I0)'), $
            /left, /top, box=0, charsize=1.3, charthick=2.0
          cc = get_kbrd(1)

       endif
       
    endfor

    if keyword_set(checkillum) then wdelete, !d.window
    
; normalize the illumination function and create the master flat and
; corresponding header

    illumflat = illumflat / mean(illumflat)

    mkhdr, illumheader, illumflat ; basic FITS header

    sxaddpar, illumheader, 'DATE', im_today(), 'Illumination flat field created'
    sxaddhist, 'ISPEC '+ispec_version(), illumheader
    sxaddhist, 'Sky flat image '+illumfile, illumheader

    if keyword_set(write) then begin

       splog, 'Writing '+datapath+illumname
       writefits, datapath+illumname, illumflat, illumheader
       
    endif
    
    masterflat = flat * illumflat
    masterflat = masterflat / mean(masterflat)
    
; compute the percent to which the image is flat
    
    sflattened = simage / masterflat
    for k = 0L, ncols-1L do sflattened[k,*] = sflattened[k,*]/mean(sflattened[k,*])
    djs_iterstat, sflattened, sigma=sfsigma, mean=sfmean, median=sfmedian, sigrej=4.5
    splog, 'Sky image is flat to '+strn(sfsigma*100.0,format='(F10.3)')+' percent.'

    xrange = sfmedian+sfsigma*[-5,+10]

    stats = strtrim(string(sfmean,format='(F12.3)'),2)+' \pm '+$
      strtrim(string(sfsigma,format='(F12.3)'),2)+$
      ' ('+strtrim(string(sfmedian,format='(F12.3)'),2)+')'

    if keyword_set(doplot) or keyword_set(write) then begin
    
       qaplotname = 'qaplot_illum_'+strmid(illumname,0,strpos(illumname,'.fits'))+'.ps'

       if keyword_set(doplot) then iter = 1L else iter = 0L
       
       for i = 0L, iter do begin

          if (i eq 0L) then begin
             dfpsplot, datapath+qaplotname, /square, /color, /isolatin1
             postthick = 8.0
          endif else begin
             im_window, 1, xratio=0.3, /square
             postthick = 2.0
          endelse

          plothist, sflattened, bin=sfsigma, xbin, ybin, /noplot
          plothist, sflattened, bin=sfsigma, charsize=1.3, charthick=postthick, $
            ytitle='Number', xtitle='Counts', title='Flattened Sky Image', $
            thick=postthick, xrange=xrange, yrange=minmax(ybin)*[1.0,1.15], $
            subtitle='['+illumfile+']', ymargin=[5,2], /fill, /fline, $
            fcolor=djs_icolor('grey'), forientation=45, fspacing=0.05, $
            xthick=postthick, ythick=postthick, color=djs_icolor('')
          legend, textoidl(stats), /right, /top, box=0, charsize=1.3, charthick=postthick

          if i eq 0L then begin
             dfpsclose
             !p.font = -1
          endif
          
       endfor

    endif 

;   flathist = fltarr(ncols)
;   for k = 0L, ncols-1L do flathist[k] = (max(sflattened[k,*])-min(sflattened[k,*]))/stddev(sflattened[k,*])

; ----------------------------------------------------------------------
; this code fits to each column of the flattened sky image to create
; the sky flat, but it doesn't work because it fits to the spectral
; lines too much
; ----------------------------------------------------------------------

; fit to each column of the flattened sky image.  ccdxy is an
; [ncols,nrows] index array 

;    ccdxy = replicate(1.0,ncols) # findgen(nrows) 
;    xy2traceset, transpose(ccdxy), transpose(sflat), tset, invvar=transpose(sflatinvvar), $
;      yfit=skyfit, maxiter=maxiter, func=func, outmask=skymask, /silent
;
;    skyflat = transpose(skyfit) / mean(skyfit) ; normalize to the mean
;    skymask = transpose(skymask)
;
;    masterflat = flat * skyflat ; multiply the dome flat and the sky flat together

;   skytest = simage / masterflat
;   print, mean(skytest), median(skytest), stddev(skytest)

; ----------------------------------------------------------------------
    
return, illumflat
end
