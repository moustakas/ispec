;+
; NAME:
;	IMASTERFLAT()
;
; PURPOSE:
;	Generate a flat field that corrects for pixel-to-pixel
;	variations.
;
; CALLING SEQUENCE:
;	masterflat = imasterflat(dimage,header,datapath=,flatname=,$
;          nmed_flat=,bsorder_flat=,nord_flat=,flatheader=,/doplot,$
;          /write,_extra=extra)
;
; INPUTS:
;	dimage   - two-dimensional dome flat
;       header   - corresponding FITS header for DIMAGE
;
; OPTIONAL INPUTS:
;	datapath  - write the flat field and the QA plot to this path
;       flatname  - name of the output flat field (default
;                   'masterflat.fits')
;	nmed_flat - number of dome flat rows to median about the middle
;                   row before fitting the spectral response (default
;                   10)
;
; PARAMETERS FOR BSPLINE_ITERFIT():
;       bsorder_flat - b-spline breakpoint interval, analogous to the
;                      IRAF cubic spline (if nord_flat=3)
;	nord_flat    - order of the B-spline fit along the dispersion 
;                      dimension to the dome flat (default 3) 
;	extra        - additional fitting/rejection parameters
;
; KEYWORD PARAMETERS:
;       doplot - generate a QA plot
;       write  - write MASTERFLAT and the QA plot to DATAPATH
;
; OUTPUTS:
;	masterflat - two-dimensional flat field with a mean of 1
;
; OPTIONAL OUTPUTS:
;       flatheader - master flat field header (based on HEADER)
;
; COMMENTS:
;	Since modeling the dome flat is very insensitive to cosmic
;	rays, the bad pixel masks are not used; perhaps a TODO item.
; 
; PROCEDURES USED:
;	DJS_MEDIAN(), BSPLINE_ITERFIT(), SXADDPAR, DJS_ITERSTAT,
;	DFPSPLOT, DFPSCLOSE, DJS_PLOT, WRITEFITS, SPLOG, LEGEND,
;	SXADDHIST, ISPEC_VERSION(), WRITEFITS, IM_WINDOW
;
; TODO:
;       Make the fitting interactive.
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August/October, U of A - written
;	jm01oct23uofa, added BSORDER optional input
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm05jun21uofa - use IM_WINDOW to spawn monitor resolution 
;                       independent windows
;
; Copyright (C) 2001, 2003, 2005, John Moustakas
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

function imasterflat, dimage, header, datapath=datapath, flatname=flatname, $
  nmed_flat=nmed_flat, bsorder_flat=bsorder_flat, nord_flat=nord_flat, $
  flatheader=flatheader, doplot=doplot, write=write, _extra=extra
    
    if not keyword_set(flatname) then flatname = 'masterflat.fits'

    if n_elements(nmed_flat) eq 0L then nmed_flat = 10L
    if n_elements(nord_flat) eq 0L then nord_flat = 3L
    if n_elements(bsorder_flat) eq 0L then bsorder_flat = 50.0 else $
      bsorder_flat = float(bsorder_flat)

    ndims = size(dimage,/n_dimension)

    imsize = size(dimage,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]

; call this routine recursively if there are multiple dome flats
    
    if ndims eq 3L then begin

       nflat = imsize[2]
       masterflat = dimage*0.0
       
       for k = 0L, nflat-1L do begin

          masterflat1 = imasterflat(dimage[*,*,k],header[*,k],datapath=datapath,$
            flatname=flatname[k],nmed_flat=nmed_flat,bsorder_flat=bsorder_flat,$
            nord_flat=nord_flat,doplot=doplot,write=write,flatheader=flatheader1)

          masterflat[*,*,k] = masterflat1
          if k eq 0L then flatheader = flatheader1 else $
            flatheader = [ [flatheader], [flatheader1] ]
          
       endfor
       return, masterflat

    endif

    midrow = nrows/2L ; middle row
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

; median filter about MIDROW +/- NMED_FLAT to improve fitting (should we
; box-car smooth as well?!?)
    
    dflat = djs_median(dimage[*,midrow-nmed_flat/2L:midrow+nmed_flat/2L],2)           ; median response

; iteratively fit a b-spline to the median spectral response function
; with break points every BKSPACE pixels

    bkpt = 0L
    sset = bspline_iterfit(colaxis,dflat,everyn=ncols/bsorder_flat,nord=nord_flat,$
      /eachgroup,yfit=dflatfit,bkpt=bkpt,outmask=outmask,_extra=extra,/silent)

; percent residuals of the fit

    good = where(outmask eq 1B,ngood) ; there better be at least one good point! 
    if ngood ne 0L then begin
       resid = dflat[good]-dflatfit[good]
       dresid = 100.0*resid/dflatfit[good]
    endif else message, 'No good points in the response fit!'
    
; divide the dome flat by the response and multiply by the slit
; function

    flat = dimage*0.0
    for i = 0L, nrows-1L do flat[*,i] = dimage[*,i] / dflatfit

    masterflat = flat / mean(flat)    ; normalize to the mean

;; derive the slit function
;; ----------------------------------------------------------------------
;    
;    slitflat = flat*0.0
;    for k = 0L, ncols-1L do slitflat[k,*] = flat[k,*] / mean(flat[k,*])
;
;    sfunction = total(slitflat,1)/ncols
;    for j = 0L, ncols-1L do flat[j,*] = flat[j,*] / sfunction
;
;; ----------------------------------------------------------------------
    
    flatheader = header
    sxaddpar, flatheader, 'DATE', im_today(), 'Flat field created'
    sxaddhist, 'ISPEC '+ispec_version(), flatheader
    
;   mkhdr, flatheader, masterflat ; basic FITS header
;   sxaddpar, flatheader, 'DATE', im_today(), 'Flat field created'
;   sxaddpar, flatheader, 'DOMEFILE', domefile

; write out     

    if keyword_set(write) then begin
       
       splog, 'Writing '+datapath+flatname
       writefits, datapath+flatname, masterflat, flatheader

    endif
       
; quality assurance plot

    if keyword_set(doplot) or keyword_set(write) then begin
    
       qaplotname = 'qaplot_response_'+strmid(flatname,0,strpos(flatname,'.fits'))+'.ps'

       if keyword_set(doplot) then iter = 1L else iter = 0L
       
       for i = 0L, iter do begin

          if (i eq 0L) then begin
             dfpsplot, datapath+qaplotname, /square, /color, /isolatin1
             postthick = 5.0
          endif else begin
             im_window, 0, xratio=0.5, /square
             postthick = 2.0
          endelse

;         title = repstr(strjoin(strsplit(qaplotname,'_',/extract),' '),'.PS','')
          
          djs_plot, colaxis, dflat, thick=postthick, line=0, color='red', xsty=3, ysty=11, $
            xtickname=replicate(' ',10), position=[0.16,0.35,0.85,0.93], $
            charthick=postthick, charsize=1.3, ytitle='Response [counts]', xthick=postthick, $
            ythick=postthick, title='Response Function Fit'
          oplot, colaxis, dflatfit, thick=postthick, line=0
          axis, /yaxis, yrange=minmax(dflat), ysty=3, ytitle='Response [counts]', $
            charsize=1.3, charthick=postthick, ythick=postthick
          
          legend, ['b-spline '+strn(long(bsorder_flat)),'order '+strn(nord_flat)], $
            /left, /top, box=0, charsize=1.3, charthick=postthick, spacing=1.7

          djs_iterstat, dresid, sigrej=5.0, mean=rmean, sigma=rsig, median=rmedian
          yrange = 10*rsig*[-1,1]

          stats = strtrim(string(rmean,format='(F12.3)'),2)+' \pm '+strtrim(string(rsig,format='(F12.3)'),2)+$
            ' ('+strtrim(string(rmedian,format='(F12.3)'),2)+')'
          legend, textoidl(stats), /right, /bottom, box=0, charsize=1.3, charthick=postthick
          
          djs_plot, colaxis, dresid, yrange=yrange, $
            xsty=3, ysty=3, position=[0.16,0.15,0.85,0.35], charthick=postthick, charsize=1.3, $
            thick=postthick, xtit='Column Number', ytitle='Residuals [%]', /noerase, yminor=2, $
            xthick=postthick, ythick=postthick, color='red', xrange=minmax(colaxis), $
            subtitle='['+repstr(flatname,'.fits','')+']'
          oplot, [!x.crange[0],!x.crange[1]], [0,0], line=0, thick=postthick
          
          if (i eq 0L) then dfpsclose
          
       endfor

    endif 
    
return, masterflat
end
