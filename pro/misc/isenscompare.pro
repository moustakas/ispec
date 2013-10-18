;+
; NAME:
;       ISENSCOMPARE
;
; PURPOSE:
;       Compare sensitivity functions from different nights or runs.
;
; CALLING SEQUENCE:
;       isenscompare, senslist, datapath=, meansens=, title=, $
;          label=, psname=, /postscript
;
; INPUTS:
;       senslist - FITS file list of sensitivity functions to compare
;
; OPTIONAL INPUTS:
;       datapath - data path to SENSLIST (default to PWD)
;       meansens - compare each sensitivity function in SENSLIST to
;                  the mean function, MEANSENS; if MEANSENS is not
;                  passed then the mean sensitivity function is the
;                  average of the functions in SENSLIST
;       title    - plot title
;       label    - optional label for each plot (e.g., the date)
;       psname   - name for the postscript output
;
; KEYWORD PARAMETERS:
;       postscript - generate postscript output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The sensitivity functions need to have the same number of
;       elements. 
;
; PROCEDURES USED:
;       IM_FITS_CUBE(), MAKE_WAVE, DJS_PLOT, LEGEND, DFPSPLOT,
;       DFPSCLOSE, REPSTR(), IRDSENSFUNC(), IPARSE_QALOGFILE() 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 Nov 27, U of A
;       jm03feb02uofa - include the seeing measurements in the
;                       sensitivity function headers 
;       jm05jun23uofa - obsolete internal routine MEAN_SEEING()
;                       removed 
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

pro isenscompare, senslist, datapath=datapath, meansens=meansens, title=title, $
  label=label, psname=psname, postscript=postscript

; will need to make datapath a vector?
    
    nsens = n_elements(senslist)
    if nsens eq 0L then begin
       print, 'Syntax - isenscompare, senslist, datapath=, meansens=, title=, $'
       print, '   label=, psname=, /postscript'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()

    if n_elements(label) eq 0L then label = reform(replicate('',nsens),1,nsens) else begin

       ndim = size(label,/n_dimension)
       case ndim of

          1L: begin
             if n_elements(label) ne nsens then begin
                splog, 'LABEL and SENSLIST do not have the same number of elements.' 
                return
             endif else label = reform(label,1,nsens)
          end
          2L: begin
             if n_elements(label[0,*]) ne nsens then begin
                splog, 'LABEL and SENSLIST do not have the same number of elements.'
                return
             endif
          end

          else: begin
             splog, 'LABEL has too many dimensions!'
             return
          endelse

       endcase
             
    endelse
    
    if n_elements(psname) eq 0L then psname = 'qaplot_sens_compare.ps'

    pushd, datapath

; grey-shifted sensitivity functions and log lists

    greysenslist = repstr(senslist,'.fits','_grey.fits')
    loglist = repstr('qalog_'+senslist,'.fits','.log')
    
; read in the individual sensitivity functions

;   senscube = irdsensfunc(senslist)
;   senswave = senscube[0].wave
    
; if MEANSENS was not passed then generate a mean sensitivity curve,
; otherwise read in MEANSENS
    
    if (n_elements(meansens) eq 0L) then begin

       splog, 'Not supported yet!'
       return
;      meanfunc = 2.5*alog10(total(10^(0.4*senscube),2)/nsens) 

    endif else begin ; remove the zero point shift from the mean curve
       
       meansenscube = irdsensfunc(meansens)
       zptshift = sxpar(meansenscube.header,'ZPTSHIFT')
       meansenscube.sens = meansenscube.sens - zptshift
       meansenscube.greysens = meansenscube.greysens - zptshift
       
    endelse

; generate plots comparing the sensitivity functions

    if keyword_set(postscript) then begin
       splog, 'Writing '+psname+'.'
       dfpsplot, psname, /color, /square
       postthick = 5.0       
    endif else begin
       window, 0, xs=500, ys=500
       postthick = 2.0
    endelse
    
    if nsens gt 1L then nx = 2L else nx = 1L
    pagemaker, nx=nx, ny=round(nsens/2.0), position=position, /normal, $
      xspace=0.7, xmargin=[1.5,0.3], ymargin=[1.0,1.5], yspace=0.0

    if nsens eq 1L then begin
       position = reform(position,4,1)
       ttitle = ''
       if n_elements(title) ne 0L then ttitle = title
       charsize = 1.8
       tcharsize = 1.8
       lcharsize = 1.8
    endif else begin
       charsize = 1.1
       tcharsize = 1.1
       lcharsize = 1.1
    endelse
    
    for i = 0L, nsens-1L do begin

; read in the each sensitivity function

       senscube = irdsensfunc(senslist[i])
       senswave = senscube.wave
       
       diff = 100*(10^(0.4*(senscube.sens-meansenscube.sens))-1) ; percent
       diffmag = senscube.sens-meansenscube.sens                 ; magnitude
;      diffmagerr = senscube.senserr
       diffmagerr = sqrt(senscube.senserr^2.0 + meansenscube.senserr^2.0)
       
       yrange = max(abs(minmax(diffmag)))*[-1.1,+1.1]

       if odd(i) eq 0L then ytitle = 'Observed - Mean (mag)'
       
       if ((i eq (nsens-2L)>0L) or (i eq (nsens-1L)>0L)) then $
         xtitle = 'Wavelength (\AA)' else $
         xtickname = replicate(' ',10)
       
       plot, senswave, diffmag, xsty=3, ysty=3, line=0, $
;      ploterror, senswave, diffmag, diffmagerr, xsty=3, ysty=3, line=0, ycharsize=ycharsize, $
         charsize=charsize, charthick=postthick, xthick=postthick, ythick=postthick, $
         thick=postthick, yrange=yrange, /noerase, position=position[*,i], $
         xtickname=xtickname, ytitle=ytitle, xtitle=xtitle, title=ttitle
       oplot, !x.crange, [0,0], line=2, thick=postthick
       legend, label[*,i], /right, /top, box=0, charsize=lcharsize, charthick=postthick
;      legend, 'Observed', /left, /top, box=0, charsize=1.2, charthick=postthick

       photo = iparse_qalogfile(loglist[i],stats=stats)
       nstds = n_elements(photo)

       djs_iterstat, photo.seeing, sigma=sigseeing, mean=meanseeing
       djs_iterstat, 100.0*photo.energy, sigma=sigenergy, mean=meanenergy

       legstats = ['Seeing = '+string(meanseeing,format='(F4.2)')+'"',$
         'Energy = '+string(meanenergy,format='(F5.1)')+'%',$
         'N = '+strn(nstds)]

       legend, legstats, /right, /bottom, box=0, charsize=lcharsize, charthick=postthick
       
       delvarx, xtitle, ytitle, xtickname
       
    endfor

    if (n_elements(title) ne 0L) and (nsens gt 1L) then $
      xyouts, position[2,0]+(position[0,1]-position[2,0])/2.0, position[3,0]*1.01, $
      title, /normal, charsize=tcharsize, align=0.5, charthick=postthick
    
    if keyword_set(postscript) then dfpsclose
    
    icleanup, senscube
    icleanup, meansenscube

    popd
    
return
end    
