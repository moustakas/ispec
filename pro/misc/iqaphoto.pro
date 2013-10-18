;+
; NAME:
;       IQAPHOTO
;
; PURPOSE:
;       Analyze the ISENSFUNC() output to determine the photometric
;       quality of the data.
;
; CALLING SEQUENCE:
;       iqaphoto, loglist, datapath=, title=, label=, psname=, /postscript
;
; INPUTS:
;       loglist - file list of ISENSFUNC() ".log" files to compare
;                 [LOGCOUNT]
;
; OPTIONAL INPUTS:
;       datapath - path for I/O (default current directory)
;       title    - plot title (e.g., '2000 April') (default none)
;       label    - string array of labels corresponding to LOGLIST
;                  [NCOMMENTS,LOGCOUNT] (default none)
;       psname   - name of the output postscript file name (default
;                  'qaplot_photo.ps') 
;
; KEYWORD PARAMETERS:
;       postscript - generate postscript output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; PROCEDURES USED:
;       READCOL, SPLOG, DFPSPLOT, DFPSCLOSE, PAGEMAKER,
;       DJS_READILINES(), TEXTOIDL(), LEGEND, DJS_ITERSTAT, ODD(),
;       CWD(), IPARSE_QALOGFILE()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 April 15, U of A, ISPEC v.1.0.0
;       jm05jun23uofa - documentation and postscript output improved 
;
; Copyright (C) 2003, 2005, John Moustakas
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

pro iqaphoto, loglist, datapath=datapath, title=title, label=label, $
  psname=psname, postscript=postscript

    logcount = n_elements(loglist)
    if logcount eq 0L then begin
       print, 'Syntax - iqaphoto, loglist, datapath=, title=, label=, $'
       print, '   psname=, /postscript'
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()

    if n_elements(label) eq 0L then label = reform(replicate('',logcount),1,logcount) else begin

       ndim = size(label,/n_dimension)
       case ndim of

          1L: begin
             if n_elements(label) ne logcount then begin
                splog, 'LABEL and LOGLIST do not have the same number of elements.' 
                return
             endif else label = reform(label,1,logcount)
          end
          2L: begin
             if n_elements(label[0,*]) ne logcount then begin
                splog, 'LABEL and LOGLIST do not have the same number of elements.'
                return
             endif
          end

          else: begin
             splog, 'LABEL has too many dimensions!'
             return
          endelse

       endcase
             
    endelse

    if n_elements(psname) eq 0L then psname = 'qaplot_photo.ps'

    pushd, datapath

    residinfo = {logfile: '', mean: 0.0, median: 0.0, min: 0.0, max: 0.0, sigma: 0.0}
    residinfo = replicate(residinfo,logcount)
    residinfo.logfile = loglist
    
    if keyword_set(postscript) then begin
       splog, 'Writing '+psname+'.'
       dfpsplot, psname, /color, /isolatin1;, /square
       postthick = 5.0
    endif else begin
       postthick = 2.0
    endelse

    xcharsize = 1.0
    ycharsize = 1.3

    pagemaker, nx=2<logcount, ny=round(logcount/2.0), position=position, /normal, $
      xspace=0.5, xmargin=[1.0,0.1], ymargin=[0.5,0.8], yspace=0.5

; --------------------------------------------------
; mean sensitivity grey shift histogram
; --------------------------------------------------
       
    if not keyword_set(postscript) then window, 0, xs=500, ys=500

    for i = 0L, logcount-1L do begin
  
       if odd(i) eq 0L then ytitle = 'Number'
       
       if ((i eq (logcount-2L)>0L) or (i eq (logcount-1L)>0L)) then $
         xtitle = 'Mean Sensitivity Grey Shift (mag)'

       if (i gt 0L) then noerase = 1L

       photo = iparse_qalogfile(loglist[i],stats=stats)
       nstds = n_elements(photo)
       
       djs_iterstat, photo.greyshift, sigma=gsig, mean=gmean, $
         median=gmed, sigrej=5.0
       djs_iterstat, photo.seeing, sigma=sigseeing, mean=meanseeing
       djs_iterstat, 100.0*photo.energy, sigma=sigenergy, mean=meanenergy

       legstats = textoidl(['\Delta = '+strn(gmean,format='(F6.3)'),$
         '\Delta_{Med} = '+strn(gmed,format='(F6.3)'),$
         '\sigma_{\Delta} = '+strn(gsig,format='(F6.3)'),$
         'Seeing = '+string(meanseeing,format='(F4.2)')+'"',$
         'Energy = '+string(meanenergy,format='(F5.1)')+'%',$
;        'Seeing = '+string(meanseeing,format='(F4.2)')+' +/- '+string(sigseeing,format='(F4.2)'),$
;        'Energy = '+string(meanenergy,format='(F5.1)')+' +/- '+string(sigenergy,format='(F3.1)')+'%',$
         'N = '+strn(nstds)])

       if nstds gt 1L then begin
       
          plothist, photo.greyshift, bin=gsig, xbin, ybin, /noplot
          if n_elements(xbin) eq 1L then begin
             xbin = xbin[0]
             ybin = ybin[0]
          endif

          xrange = minmax(xbin)+[-2*gsig,+2*gsig]
          yrange = minmax(ybin)*[0.0,2]

       endif else begin

          xrange = gmean+[-0.05,0.05]
          yrange = [0,1.2]
          
       endelse
          
       plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xsty=3, ysty=1, $
         xthick=postthick, ythick=postthick, thick=postthick, xcharsize=xcharsize, ycharsize=ycharsize, $
         charthick=postthick, xtickname=xtickname, xtitle=xtitle, ytitle=ytitle, $
         position=position[*,i], noerase=noerase
       if nstds gt 1L then $
         plothist, photo.greyshift, bin=gsig, /overplot, $
         thick=postthick, /fline, /fill, forientation=45 else $
         oplot, [gmean,gmean], [0,1.0], line=0, thick=postthick
         
       legend, [label[*,i]], /left, /top, box=0, charsize=1.2, charthick=postthick
       legend, legstats, /right, /top, charsize=1.2, charthick=postthick, box=0

       delvarx, xtitle, ytitle, xtickname, noerase, xrange, yrange
       
    endfor

    if n_elements(title) ne 0L then $
      xyouts, position[2,0]+(position[0,1]-position[2,0])/2.0, charthick=postthick, $
      position[3,0]*1.01, title, /normal, charsize=1.5, align=0.5

; --------------------------------------------------
; zero point shift histogram
; --------------------------------------------------
       
    if not keyword_set(postscript) then window, 2, xs=500, ys=500

    for i = 0L, logcount-1L do begin
  
       if odd(i) eq 0L then ytitle = 'Number'
       
       if ((i eq (logcount-2L)>0L) or (i eq (logcount-1L)>0L)) then $
         xtitle = 'Sensitivity Zero Point (mag)'

       if (i gt 0L) then noerase = 1L

       photo = iparse_qalogfile(loglist[i],stats=stats)
       nstds = n_elements(photo)
       
       djs_iterstat, photo.senszero, sigma=ssig, mean=smean, $
         median=smed, sigrej=5.0
       djs_iterstat, photo.seeing, sigma=sigseeing, mean=meanseeing
       djs_iterstat, 100.0*photo.energy, sigma=sigenergy, mean=meanenergy

       legstats = textoidl(['\Delta = '+strn(smean,format='(F6.3)'),$
         '\Delta_{Med} = '+strn(smed,format='(F6.3)'),$
         '\Delta_{Max} = '+strn(max(photo.senszero),format='(F6.3)'),$
         '\sigma_{\Delta} = '+strn(ssig,format='(F6.3)'),$
         'Seeing = '+string(meanseeing,format='(F4.2)')+'"',$
         'Energy = '+string(meanenergy,format='(F5.1)')+'%',$
;        'Seeing = '+string(meanseeing,format='(F4.2)')+' +/- '+string(sigseeing,format='(F4.2)'),$
;        'Energy = '+string(meanenergy,format='(F5.1)')+' +/- '+string(sigenergy,format='(F3.1)')+'%',$
         'N = '+strn(nstds)])

       if nstds gt 1L then begin

          plothist, photo.senszero, xbin, ybin, bin=ssig, /noplot              
          if n_elements(xbin) eq 1L then begin
             xbin = xbin[0]
             ybin = ybin[0]
          endif
       
          xrange = minmax(xbin)+[-2*ssig,+3*ssig]
          yrange = minmax(ybin)*[0.0,2.5]
       
       endif else begin

          xrange = smean+[-0.05,0.05]
          yrange = [0,1.2]
          
       endelse
          
       plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xsty=3, ysty=1, $
         xthick=postthick, ythick=postthick, thick=postthick, xcharsize=xcharsize, ycharsize=ycharsize, $
         charthick=postthick, xtickname=xtickname, xtitle=xtitle, ytitle=ytitle, $
         position=position[*,i], noerase=noerase
       if nstds gt 1L then $
         plothist, photo.senszero, bin=ssig, /overplot, thick=postthick, /fline, $
         /fill, forientation=45 else $
         oplot, [smean,smean], [0,1.0], line=0, thick=postthick
         
       legend, [label[*,i]], /left, /top, box=0, charsize=1.2, charthick=postthick
       legend, legstats, /right, /top, charsize=1.2, charthick=postthick, box=0

       delvarx, xtitle, ytitle, xtickname, noerase, xrange, yrange
       
    endfor

    if n_elements(title) ne 0L then $
      xyouts, position[2,0]+(position[0,1]-position[2,0])/2.0, charthick=postthick, $
      position[3,0]*1.01, title, /normal, charsize=1.5, align=0.5
       
    if keyword_set(postscript) then dfpsclose

;   struct_print, stats
    
    popd
    
return
end
