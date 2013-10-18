;+
; NAME:
;	IOVERSCAN
;
; PURPOSE:
;	Subtract the overscan region and trim the data.
;
; CALLING SEQUENCE:
;       ioverscan, image, header, overscan=, norder_overscan=, $
;          bkpt_overscan=, /checkoverscan, /spine_overscan
;
; INPUTS:
;	image    - two dimensional raw image
;	header   - FITS header
;
; OPTIONAL INPUTS:
;	overscan        - [x1,x2,y1,y2] defining the overscan region 
;                         (zero-indexed)
;	norder_overscan - order of the fit to the overscan region
;                         (default 2)
;       bkpt_overscan   - explicit breakpoints to be used for spline
;                         (only relevant if SPLINE_OVERSCAN is set)
;	
; KEYWORD PARAMETERS:
;       checkoverscan     - generate a plot of the overscan region and fit
;       spline_overscan   - fit a spline instead of a legendre function
;       constant_overscan - use a median value for the overscan (no fit)
;
; OUTPUTS:
;	image  - (modified)
;	header - (modified)
;
; COMMENTS:
;	Default trim and overscan regions are taken from the data
;	headers.  The overscan section of any flat fields is set to
;	norder_overscan+1.  The overscan region is determined to be
;	along the columns or rows by assuming that the overscan axis
;	is the longer of the overscan region dimensions.
;
; PROCEDURES USED:
;	FUNC_FIT(), DJS_MEDIAN(), SXDELPAR, SXADDPAR,
;	BSPLINE_ITERFIT()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August, U of A
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm04jun21uofa - additional error checking; automatically
;                       determine whether the overscan region is along
;                       rows or along columns
;       jm04sep21uofa - if NORDER_OVERSCAN is negative then do not
;                       subtract the overscan vector
;       jm05jun17uofa - removed QFLAG input; simplify the OVERSCAN
;                       header keyword (no longer follows the IRAF
;                       convention) 
;       arm05jul12uofa - SPLINE_OVERSCAN and BKPT_OVERSCAN added
;       jm06jun28uofa - CONSTANT_OVERSCAN keyword added
;
; Copyright (C) 2001, 2003-2005, John Moustakas
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
    
pro ioverscan, image, header, overscan=overscan, norder_overscan=norder_overscan, $
  bkpt_overscan=bkpt_overscan, checkoverscan=checkoverscan, constant_overscan=constant_overscan, $
  spline_overscan=spline_overscan

    if (n_elements(norder_overscan) eq 0L) then norder_overscan = 1L

    if (norder_overscan lt 0L) then begin
       sxaddhist, 'No overscan subtraction performed.', header
       return
    endif
    
    origimage = image           ; original data
    image = origimage*0.0

    imsize = size(origimage,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]

    colaxis = findgen(ncols)
    rowaxis = findgen(nrows)
    
; pad the default overscan region by two columns    

    if keyword_set(overscan) then overscan = long(overscan) else begin 

       biassec = sxpar(header,'BIASSEC',count=bcount)
       if (bcount eq 0L) then begin
          splog, 'No overscan region specified.'
          return
       endif
       
       overscan = long(strsplit(biassec,'[:,] ',/extract))-1L

       overscan[0] = overscan[0] + 2L
       overscan[1] = overscan[1] - 2L

    endelse

    oscanstr = '['+strjoin(strcompress(overscan[0:1],/remove),',')+':'+$
      strjoin(strcompress(overscan[2:3],/remove),',')+']'
    
    x1 = overscan[0] > 0L
    x2 = overscan[1] < (ncols-1L)

    y1 = overscan[2] > 0L
    y2 = overscan[3] < (nrows-1L)

;amd060622
;stop

; figure out if the overscan is along rows or along columns; do this
; by assuming that the range of the overscan region will always be
; along the smallest dimension

    alongrows = 0L
    alongcols = 0L
    
    if (x2-x1) lt (y2-y1) then alongcols = 1L else alongrows = 1L
    
    imtype = sxpar(header,'imagetyp',count=icount)
    if (icount eq 0L) then imtype = '' else imtype = strlowcase(strcompress((imtype),/remove))

    object = sxpar(header,'object',count=ocount)
    if (ocount eq 0L) then object = '' else object = strlowcase(strcompress((object),/remove))

    if strmatch(object,'*flat*') then nord = norder_overscan + 1L else nord = norder_overscan
    
; median of the overscan region

    ny = n_elements(rowaxis[y1:y2])
    nx = n_elements(colaxis[x1:x2])

    if alongcols then begin

       scanvector = djs_median(origimage[x1:x2,y1:y2],1L)
       scanerror = scanvector*0.0
       for i = 0L, ny-1L do scanerror[i] = djsig(origimage[x1:x2,y1+i],sigrej=3.0)
    
; fit and subtract the overscan

       if keyword_set(constant_overscan) then begin

           fit = scanvector*0.0+median(scanvector)

       endif else begin
           
           if KEYWORD_SET(spline_overscan) then begin
               if (N_ELEMENTS(bkpt_overscan) eq 0L) then begin
                   bkptres = (y2-y1)/nord ; break point resolution
                   bkpt_overscan = findgen(nord)*bkptres+min(rowaxis)
               endif
               sset = bspline_iterfit(rowaxis[y1:y2], scanvector, nord=3, $
                                      yfit=fit,bkpt=bkpt_overscan, /silent)
           endif else coef = func_fit(rowaxis[y1:y2],scanvector,nord, $
                                      function_name='flegendre',yfit=fit)

       endelse 
           
       for k = 0L, ncols-1L do image[k,y1:y2] = origimage[k,y1:y2] - fit

    endif else begin

       scanvector = djs_median(origimage[x1:x2,y1:y2],2L)
       scanerror = scanvector*0.0
       for i = 0L, nx-1L do scanerror[i] = djsig(origimage[x1+i,y1:y2],sigrej=3.0)

       if keyword_set(constant_overscan) then begin

           fit = scanvector*0.0+median(scanvector)

       endif else begin
       
           if KEYWORD_SET(spline_overscan) then begin
               if N_ELEMENTS(bkpt_overscan) eq 0 then begin
                   bkptres = (x2-x1)/nord ; break point resolution
                   bkpt_overscan = findgen(nord)*bkptres+min(colaxis)
               endif
               sset = bspline_iterfit(colaxis[x1:x2], scanvector, nord=3, $
                                      yfit=fit, bkpt=bkpt_overscan,/silent)
           endif else coef = func_fit(colaxis[x1:x2],scanvector,nord, $
                                      function_name='flegendre',yfit=fit)

       endelse
           
       for k = 0L, nrows-1L do image[x1:x2,k] = origimage[x1:x2,k] - fit

    endelse
       
    residuals = 100.0*(scanvector-fit)/fit
    yrange = median(scanvector)+djsig(scanvector,sigrej=3.0)*[-3.0,3.0]
;   yrange = [min(scanvector-scanerror),max(scanvector+scanerror)]*[0.999,1.001]
    
    if keyword_set(checkoverscan) then begin

       if (!d.window ne 3L) then window, 3, xs=400, ys=400

       if alongcols then begin

          djs_plot, rowaxis[y1:y2], scanvector, xsty=3, ysty=3, ps=3, charsize=1.5, $
            charthick=2.0, xtitle='Row Number', ytitle='Counts', yrange=yrange, $
            title=strupcase(object)+' ('+strupcase(imtype)+')', color='green'
          djs_oplot, rowaxis[y1:y2], fit, thick=3.0, color='red'
          djs_oplot, rowaxis[y1:y2], scanvector+scanerror, ps=3, color='yellow'
          djs_oplot, rowaxis[y1:y2], scanvector-scanerror, ps=3, color='yellow'

       endif else begin

          djs_plot, colaxis[x1:x2], scanvector, xsty=3, ysty=3, ps=3, charsize=1.5, $
            charthick=2.0, xtitle='Column Number', ytitle='Counts', yrange=yrange, $
            title=strupcase(object)+' ('+strupcase(imtype)+')', color='green'
          djs_oplot, colaxis[x1:x2], fit, thick=3.0, color='red'
          djs_oplot, colaxis[x1:x2], scanvector+scanerror, ps=4, color='yellow'
          djs_oplot, colaxis[x1:x2], scanvector-scanerror, ps=4, color='yellow'

       endelse
          
       legend, oscanstr, /left, /top, box=0, charthick=2.0, charsize=1.5
       cc = get_kbrd(1)

    endif
       
; update the header, remembering that IRAF indices start at one

    sxdelpar, header, 'BIASSEC'
    sxaddpar, header, 'OVERSCAN', '['+string(x1+1,format='(I0)')+':'+string(x2+1,format='(I0)')+','+$
      string(y1+1,format='(I0)')+':'+string(y2+1,format='(I0)')+']', ' overscan region (mean='+$
      strtrim(string(djs_mean(scanvector),format='(F12.2)'),2)+')', before='HISTORY'

;   sxaddpar, header, 'OVERSCAN', im_today()+' Overscan section '+$
;     '['+string(x1+1,format='(I0)')+':'+string(x2+1,format='(I0)')+','+$
;     string(y1+1,format='(I0)')+':'+string(y2+1,format='(I0)')+'] with mean='+$
;     strtrim(string(djs_mean(scanvector),format='(F12.2)'),2), before='HISTORY'
    
return
end
