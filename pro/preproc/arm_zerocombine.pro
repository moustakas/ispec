;+
; NAME: ARM_ZEROCOMBINE
;       
; CATEGORY: astronomy
;
; PURPOSE: create master bias frame (aka, zero frame)
;
; CALLING SEQUENCE: ARM_ZEROCOMBINE, infile, outfile, [ext0=, next=, $
;                      nhigh=, nlow=, psfile=, memsaver=, /silent, /gzip]
;
; INPUTS:
;   infile  - string array of bias frame file names
;   outfile - string file name for master bias frame output
;       
; OPTIONAL INPUTS:
;   ext0   - first FITS extension to be read (default=0)
;   next   - number of FITS extensions to be read (default=1)
;   nhigh  - number of high values to discard (default=1)
;   nlow   - number of low values to discard (default=1)
;   psfile - string name for diagnostic plots postscript file
;   memsaver - if non-zero, individual frames are read, stored in
;              memory and evaluated several columns at a time rather
;              than as entire images - this has no impact on the
;              master bias generated, but may result in a net gain in
;              speed if execution of the routine is slowed by memory
;              limitations (despite the corresponding increase in
;              file-reading operations); MEMSAVER should be a
;              decimal number in the range (0,1] indicating the
;              fraction of columns to be grouped (eg, if
;              MEMSAVER=0.25, then 25% of the columns are stored in
;              memory and evaluated at a given time, requiring each
;              frame to be read 4 times); default=1.0
;
; KEYWORDS:
;   silent - suppress status messages
;   gzip   - compress output FITS and postscript files
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED: ARM_HISTOMAKER, ARM_PLOTCONFIG, ARM_GETINDEX(),
;                  SXADDPAR, SXADDHIST, MWRFITS, MRDFITS(), DJS_ITERSTAT
;
; COMMENTS: The master bias frame is constructed by discarding the
;           largest NHIGH values and smallest NLOW values for each
;           pixel and averaging the remaining values.  For discrete
;           values this is preferable to median combining.  In
;           principle, NHIGH should equal NLOW in order to avoid
;           biasing the mean of the distribution; however, in practice
;           this is not likely to be of much consequence.  Generally
;           speaking, the default of 1 for NHIGH and NLOW is adequate
;           to elliminate all potential cosmic rays.
;
;           Warning messages are printed to screen if one particular
;           bias frame is responsible for more high or low discarded
;           pixels than expected from simple statistical
;           considerations.  This is meant to draw the user's
;           attention to files which are not bias frames and might
;           have been inadvertantly included.  Note, however, that
;           true bias frames may trigger this warning if the mean bias
;           level varies amongst bias frames. 
; 
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2004 Oct 20
;    multiple FITS extensions supported, ARM, 2005 Jun 06
;    DISPLAY keyword elliminated, ARM, 2005 Jul 11
;    cleaner error-induced exits, ARM, 2005 Jul 11
;    EXTO parameter added, _EXTRA removed, bug fixed, ARM, 2005 Jul 19
;    minor plotting changes, ARM, 2005 Jul 19
;    MRDFITS called with UNSIGNED keyword, ARM, 2005 Jul 20
;    MEMSAVER parameter added; ARM, 2005 Jul 22
;    faster due to taking mean w/o rejection; ARM, 2005 Jul 22
;    SILENT and GZIP keywords added; ARM, 2005 Jul 22
;    MWRFITS called with CREATE keyword; ARM, 2005 Jul 25
;    EXT indexing bug fixed; ARM, 2005 Jul 26
;
; BUG REPORT: Please report any bugs to Andrew R. Marble.
;
; Copyright (C) 2004, 2005, Andrew R. Marble
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

pro arm_zerocombine_plot, mean, stdev, avgs, histodata, nrows, ncols, $
                          outfile, thick, exp, sig, nfile, nhigh, nlow

; specify plotting preferences
       
    p_orig = !p
    x_orig = !x
    y_orig = !y

    !p.charsize  = 2
    !p.charthick = thick
    !p.thick     = thick
    !x.thick     = thick
    !y.thick     = thick

; determine desired plotting coordinates

    xmargin1 = [1.25, 0.6]
    xmargin  = [1.25, 0.75]
    ymargin  = [1.5, 1.0]
    yspace   = [1.7, 0.75]
    yspace1  = [2.2, 1.5]

    ARM_PLOTCONFIG, nx=2, ny=6, coords=coords2, $
      xmargin=xmargin, ymargin=ymargin, yspace=[0.1,yspace[0],0.1,yspace[1],0]
    ARM_PLOTCONFIG, nx=2, ny=3, coords=coords, $
      xmargin=xmargin1, ymargin=ymargin, yspace=yspace1

    y = [0.915, 0.60, 0.57, 0.255]
    x = [0.05, 0.235, 0.32, 0.68, 0.765, 0.95]
    
    SHADE_SURF, mean, pos=coords[*,0], xst=7, yst=7, charsize=4
    PLOT,  TOTAL(mean, 1)/ncols, linestyle=0, thick=2, pos=coords2[*,1], yst=3, xst=3, $
      xtickname=REPLICATE(' ', 30), psym=10
    AXIS, xaxis=1, xsty=1, xtitle='row number (pixels)'
    PLOT,  TOTAL(mean, 2)/nrows, linestyle=0, thick=2, pos=coords2[*,3], yst=3, xst=3, $
      xtitle='column number (pixels)', psym=10, /noerase

    SHADE_SURF, stdev, pos=coords[*,2], xst=7, yst=7, charsize=4
    PLOT,  TOTAL(stdev, 1)/ncols, linestyle=0, thick=2, pos=coords2[*,5], yst=3, xst=3, $
      xtickname=REPLICATE(' ', 30), psym=10
    AXIS, xaxis=1, xsty=1, xtitle='row number (pixels)'
    PLOT,  TOTAL(stdev, 2)/nrows, linestyle=0, thick=2, pos=coords2[*,7], yst=3, xst=3, $
      xtitle='column number (pixels)', psym=10, /noerase
     
    PLOT, histodata.x, avgs, psym=10, pos=coords[*,4], xst=3, yst=3, $
      xtitle='frame sequence number', ytitle='mean bias level'
    
    PLOT, histodata.x, /nodata, pos=coords[*,5], $
      xtitle='frame sequence number!C(solid/dotted lines denote min/max values)!C'+ $
      textoidl('[shaded: \pm3\sigma expectation (based on random)]'), $
      ytitle='min/max pixel dist. (%) ', xst=3, yst=3, $
      yr=1d2/(ncols*nrows)*[MIN([histodata.y1,histodata.y2])<(exp-4*sig), $
                            MAX([histodata.y1,histodata.y2])>(exp+4*sig)]
    LOADCT, 0
    POLYFILL, [0,[1,1]*MAX(histodata.x),0,0], (exp+[-1,-1,+1,+1,-1]*3*sig)/nrows/ncols*1d2, color=200
    OPLOT, histodata.x, 1d2*histodata.y1/(ncols*nrows), linestyle=0, psym=10
    OPLOT, histodata.x, 1d2*histodata.y2/(ncols*nrows), linestyle=1, psym=10

; plot boxes

    PLOT, [0,0], /nodata, /noerase, position=[0,0,1,1], xst=5, yst=5

    OPLOT, [x[0], x[5]], [y[1], y[1]]
    OPLOT, [x[0], x[0]], [y[1], y[0]]
    OPLOT, [x[5], x[5]], [y[1], y[0]]
    OPLOT, [x[0], x[2]], [y[0], y[0]]
    OPLOT, [x[3], x[5]], [y[0], y[0]]

    OPLOT, [x[0], x[5]], [y[3], y[3]]
    OPLOT, [x[0], x[0]], [y[3], y[2]]
    OPLOT, [x[5], x[5]], [y[3], y[2]]
    OPLOT, [x[0], x[1]], [y[2], y[2]]
    OPLOT, [x[4], x[5]], [y[2], y[2]]

; annotate plots

    XYOUTS, 0.5, 0.95, align=0.5, /normal, charsize=1.25, charthick=1.5*thick, $
      outfile+' ('+STRTRIM(STRING(nfile,f='(i)'),2)+ $
      ' bias frames, NLOW='+STRTRIM(STRING(nlow,f='(i)'),2)+ $ 
      ', NHIGH='+STRTRIM(STRING(nhigh,f='(i)'),2)+')'
    XYOUTS, 0.5, y[0], align=0.5, /normal, charsize=1.5, $
      'mean of unrejected pixels'
    XYOUTS, 0.5, y[2],  align=0.5, /normal, charsize=1.5, $
      'standard deviation of unrejected pixels'
    XYOUTS, 0.05, 0.03, align=0, /normal, charsize=1.0, $
      systime()

; restore original plotting preferences

    !p = p_orig
    !x = x_orig
    !y = y_orig

 end

pro arm_zerocombine, infile, outfile, trim=trim, next=next, nhigh=nhigh, nlow=nlow, $
                     ext0=ext0, psfile=psfile, memsaver=memsaver, $
                     silent=silent, gzip=gzip 

; defaults

    if N_ELEMENTS(nhigh)    eq 0L then nhigh = 1
    if N_ELEMENTS(nlow)     eq 0L then nlow  = 1
    if N_ELEMENTS(next)     eq 0L then next  = 1    
    if N_ELEMENTS(ext0)     eq 0L then ext0  = 0
    if N_ELEMENTS(memsaver) eq 0L then memsaver = 1.0

; error checking

    if memsaver eq 0.0 or memsaver gt 1.0 then begin
       PRINT, 'ARM_ZEROCOMBINE: invalid MEMSAVER value... using default.'
       memsaver = 1.0
    endif

    nfile = N_ELEMENTS(infile)
    if nhigh + nlow ge nfile then begin
       PRINT, 'ARM_ZEROCOMBINE: All values rejected!'
       return
    endif

    bad = WHERE(FILE_TEST(infile) eq 0L, count)
    if count gt 0L then begin
       PRINT, 'ARM_ZEROCOMBINE: File(s) '+STRJOIN(infile[bad],', ')+' could not be found!'
       return
    endif

    ntrim = n_elements(trim) ; jm06nov06nyu
    if (ntrim ne 0L) then begin
       if (ntrim ne 4L) then begin
          print, 'TRIM must be a 4-element array.'
          return
       endif
    endif

; initialization

    avgs = FLTARR(nfile)

    if N_ELEMENTS(psfile) gt 0 then ARM_PLOTCONFIG, psfile=psfile, /writeover

; read and store indivial frames

    len = STRLEN(STRING(nfile, f='(i0)'))
    fmt = '(i'+STRING(len,f='(i0)')+'.'+STRING(len,f='(i0)')+')'

    for ext = ext0, ext0+next-1L do begin

       for i = 0L, nfile-1 do begin
          
          if not KEYWORD_SET(silent) then $
            PRINT, format='("ARM_ZEROCOMBINE: Reading file ",'+fmt+'," of ",'+fmt+','+ $
            '" (FITS extension '+STRTRIM(STRING(ext, f='(i)'),2)+')",a1,$)',i+1, nfile, $
            STRING(13b)
          
          im = MRDFITS(infile[i], ext, hdr, /silent, /unsigned) * 1d0
          
          if (ntrim ne 0L) then im = im[trim[0]:trim[1],trim[2]:trim[3]] ; jm06nov06nyu

          DJS_ITERSTAT, im, median=mn;, mean=mn ; jm06nov06nyu
          avgs[i] = mn

          im2 = im[0:ROUND(memsaver*N_ELEMENTS(im[*,0]))-1,*]
          if i eq 0L then cube = im2 else cube = [[[cube]],[[im2]]]

       endfor
       
       PRINT
       
; combine the frames
       
       ncols = N_ELEMENTS(im[*,0])
       nrows = N_ELEMENTS(im[0,*])
       
       mean  = DBLARR(ncols, nrows)
       stdev = DBLARR(ncols, nrows)
       
       min = LONARR(ncols, nrows)
       max = LONARR(ncols, nrows)
       
       len1 = STRLEN(STRING(ncols, f='(i0)'))
       fmt1 = 'i'+STRING(len1,f='(i0)')+'.'+STRING(len1,f='(i0)')
       len2 = STRLEN(STRING(nrows, f='(i0)'))
       fmt2 = 'i'+STRING(len2,f='(i0)')+'.'+STRING(len2,f='(i0)')
       
       for l = 0, CEIL(1.0/memsaver)-1 do begin
          
          offset = ROUND(l*memsaver*ncols) > 0

          if l ne 0 then begin
             for k = 0L, nfile-1 do begin                
                if not KEYWORD_SET(silent) then $
                  PRINT, format='("ARM_ZEROCOMBINE: Reading file ",'+fmt+'," of ",'+fmt+','+ $
                  '" (FITS extension '+STRTRIM(STRING(ext, f='(i)'),2)+')",a1,$)',k+1, nfile, $
                  STRING(13b)
                im = MRDFITS(infile[k], ext, /silent, /unsigned) * 1d0 
                im = im[ROUND(l*memsaver*ncols):ROUND((l+1)*memsaver*ncols-1)<(ncols-1),*]
                if k eq 0L then cube = im else cube = [[[cube]],[[im]]]
             endfor
             PRINT
          endif
             
          for i = ROUND(l*memsaver*ncols),ROUND((l+1)*memsaver*ncols-1)<(ncols-1) do begin
             
             if not KEYWORD_SET(silent) then $
               PRINT, format='("ARM_ZEROCOMBINE: Evaluating pixel [",'+fmt1+',",",'+fmt2+',"] of [",'+ $
               fmt1+',",",'+fmt2+',"]",a1,$)', i, nrows-1, ncols-1, nrows-1, STRING(13b)
             
             for j = 0, nrows-1 do begin
                
                vals = cube[i-offset,j,*]
                                   
; get indices of min/max values; in the case of degenerate min/max
; values, pick randomly to avoid a biased distribution...
    
                whmin = ARM_GETINDEX(vals, MIN(vals), nmatch=nwhmin)
                whmax = ARM_GETINDEX(vals, MAX(vals), nmatch=nwhmax)
                if nwhmin eq 1 then min[i,j] = whmin[0] else $
                  min[i,j] = whmin[FLOOR(RANDOMU(seed,1)*nwhmin)]
                if nwhmax eq 1 then max[i,j] = whmax[0] else $
                  max[i,j] = whmax[FLOOR(RANDOMU(seed,1)*nwhmax)]
                
                ord  = SORT(vals)
                vals = (vals[ord])[nlow:nfile-1-nhigh]
                
                mean[i,j]  = TOTAL(vals) / (nfile-nlow-nhigh)
                if N_ELEMENTS(vals) gt 1 then stdev[i,j] = STDDEV(vals) else stdev[i,j] = 0

             endfor
             
          endfor
          
          PRINT
       
       endfor
      
; check for anomalies

       ARM_HISTOMAKER, min, max, histodata=histodata, min=-0.5, /quiet, bin=1, /noplot
       exp = ncols * nrows / nfile
       sig = SQRT(exp)
       dummy = WHERE(histodata.y2 gt (exp + 3*sig), nmaxdev)
       dummy = WHERE(histodata.y1 gt (exp + 3*sig), nmindev)
       
; write to postscript if desired
       
       if N_ELEMENTS(psfile) gt 0 then ARM_ZEROCOMBINE_PLOT, $
         mean, stdev, avgs, histodata, nrows, ncols, outfile, 2, exp, sig, nfile, nhigh, nlow
       
; output warning message if necessary
       
       if nmaxdev gt 0 then begin
          maxdev = ARM_GETINDEX(histodata.y2, MAX(histodata.y2))
          PRINT, 'ARM_ZEROCOMBINE: Note - '+STRTRIM(STRING(1d2*histodata.y2[maxdev[0]]/(ncols*nrows), $
             f='(i)'),2)+'% of highest values come from '+infile[maxdev]+'.'
       endif
       if nmindev gt 0 then begin
          mindev = ARM_GETINDEX(histodata.y1, MAX(histodata.y1))
          PRINT, 'ARM_ZEROCOMBINE: Note - '+STRTRIM(STRING(1d2*histodata.y1[mindev[0]]/(ncols*nrows), $
            f='(i)'),2)+'% of lowest values come from '+infile[mindev]+'.'
       endif
       
       if ext eq 0 then begin
          
; update the header
          
          h = hdr
          for i = 1L, nfile do SXADDPAR, h, 'ZERO_'+STRING(i,format=fmt), $
            infile[i-1], 'individual bias frame FITS file', before='HISTORY'
          
          now = STRMID(SYSTIME(),20)+' '+STRMID(SYSTIME(),4,12)
          SXADDHIST, "'Master bias frame created from ZERO_* files ("+now+")'", h 
          
; output the combined image
          
          MWRFITS, FLOAT(mean), outfile, h, /create

       endif else MWRFITS, FLOAT(mean), outfile
       
    endfor
    
    if N_ELEMENTS(psfile) gt 0 then begin
       DEVICE, /close
       SET_PLOT, 'x'
    endif

    if KEYWORD_SET(gzip) then begin
       SPAWN, 'gzip -f '+outfile
       SPAWN, 'gzip -f '+psfile
    endif
       
 end
 
 









