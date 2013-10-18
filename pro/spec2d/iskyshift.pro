;+
; NAME:
;       ISKYSHIFT()
;
; PURPOSE:
;       Improve the wavelength solution using the night sky spectrum.
;
; CALLING SEQUENCE:
;       skyshift = iskyshift(sky,header,skyspecfile=)
;
; INPUTS:
;       sky    - one dimensional sky spectrum [NCOLS]
;       header - corresponding FITS header
;
; OPTIONAL INPUTS:
;	skylinefile - name of the sky spectrum file to read in and
;                     cross-correlate with the observed sky spectrum 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       skyshift - pixel shift
;       header   - (modified) with the new wavelength coordinates and
;                  a history note 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       If SKYSPECFILE is passed then use that spectrum and
;       cross-correlation to improve the wavelength solution,
;       otherwise use the sky lines in ${ISPEC_DIR}/etc/skylines.dat.
;       if the sky lines file does not exist, or if all the lines in
;       that file are outside the wavelength range then do not apply a 
;       shift. 
;
; EXAMPLE:
;
; PROCEDURES USED:
;       SPLOG, MRDFITS(), COMBINE1FIBER, FINDCHI2MIN, READCOL,
;       MPFITPEAK(), SXADDHIST, SXADDPAR, IM_TODAY(), SXPAR(),
;       TRACE_FWEIGHT(), TRACE_GWEIGHT()
;
; DATA FILES:
;       ${ISPEC_DIR}/etc/skylines.dat
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 16, U of A  excised from IEXTRACT
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm05jun21uofa - improved centroiding that relies on more than
;                       just the 5577 sky line; use iterative flux-
;                       and gaussian-weighted centroiding rather than
;                       a gaussian fit to the sky lines
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

function iskyshift, sky, header, skyspecfile=skyspecfile

    if (n_elements(sky) eq 0L) then begin
       print, 'Syntax - skyshift = iskyshift(sky,header,skyspecfile=)'
       return, 0
    endif
    
    if (n_elements(skyspecfile) eq 0L) then skyspecfile = ''

    wave = make_wave(header,cd1_1=dwave)    ; wavelength vector

    ncols = n_elements(sky)
    colaxis = findgen(ncols)
    
    etcpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    
    if (file_test(etcpath+skyspecfile,/regular) eq 1L) then begin

; use cross-correlation to find the shift

       splog, 'Reading '+etcpath+skyspecfile+'.'
       skyfits = mrdfits(etcpath+skyspecfile,0,skyhead,/silent)
       invvar = skyfits*0.0+1.0

       combine1fiber, alog10(make_wave(skyhead)), skyfits, invvar, $
         newloglam=alog10(wave), newflux=newflux
       skyfits = newflux

       skyfits = skyfits/max(skyfits)
       
       maxshift = 10.0          ; maximum allowable shift [pixels]
       nsamp = 1.0              ; oversampling factor
       npix = maxshift*nsamp
       
       lags = lindgen(2*npix)-npix
       nlags = n_elements(lags)
       chi2array = fltarr(nlags)

;      skynorm = interpol(sky,wave,5500.0) ; normalize the sky spectrum
       skynorm = max(sky)

       for j = 0L, nlags-1L do chi2array[j] = total((skyfits-shift(sky/skynorm,lags[j]))^2.0)
;      minchi2 = min(chi2array,lagbest) ; minimum chi2
;      plot, lags, chi2array, ysty=3, ps=10, xsty=3

; fit a parabolla to the minimum          
       
       findchi2min, lags, chi2array, minchi2, skyshift, doplot=doplot
       splog, 'Applying a '+strn(skyshift*dwave,format='(F6.2)')+$
         ' Angstrom shift based on cross-correlation.'

;      djs_iterstat, sky, median=skymed, sigma=skysig
;      djs_plot, wave, skynorm*skyfits, ps=10, /xsty, /ysty, $
;        yrange=skymed+[-3*skysig,10*skysig], xr=[5000,6500]
;      djs_oplot, wave, shift(sky,skyshift), ps=10, color='green'
       
    endif else begin            

; instead of cross-correlation use flux-weighted centroiding of
; individual night sky lines   
       
       if (file_test(etcpath+'skylines.dat',/regular) eq 1L) then begin ; sky line file exists

          splog, 'Reading '+etcpath+'skylines.dat.'
          readcol, etcpath+'skylines.dat', skywaves, /silent, comment='#', format='D'
          good = where((skywaves gt min(wave)) and (skywaves lt max(wave)),nskylines)

          if (nskylines eq 0L) then begin ; no sky lines in range

             splog, 'No sky lines in the wavelength range . . . no shift applied.'
             skyshift = 0.0

          endif else begin      ; sky lines in range

             skywaves = skywaves[good] ; [Angstrom]
             skypredict = (skywaves-min(wave))/dwave ; predicted positions [pixels]
             
; iteratively find the centers of the night sky lines 
             
             skypixels = trace_fweight(sky,skypredict,replicate(0,nskylines),$
               radius=3.0,invvar=1.0/sky)

             for iiter = 0L, 1L do begin
                medianshift = median(skypredict-skypixels)
                skypixels = trace_fweight(sky,skypredict-medianshift,$
                  replicate(0,nskylines),radius=3.0,invvar=1.0/sky)
             endfor

             skypixels = trace_gweight(sky,skypixels,replicate(0,nskylines),$
               sigma=1.0,invvar=1.0/sky)

             djs_iterstat, skypredict-skypixels, median=skyshift, mask=skymask, sigrej=2.0
             junk = where(skymask eq 1L,nskylines)
             if (nskylines eq 0L) then skyshift = 0.0
             
;;           for j = 0L, nskylines-1L do begin
;;
;;              local = where((colaxis gt skypredict[j]-10) and (colaxis lt skypredict[j]+10),nlocal)
;;              plot, wave[local], sky[local], ps=10, xsty=3, ysty=3
;;              cc = get_kbrd(1)
;;               
;;           endfor

;            nterm = 4L
;            
;            skypixels = fltarr(nskylines)
;            for j = 0L, nskylines-1L do begin
;
;               local = where((colaxis gt skypredict[j]-50) and (colaxis lt skypredict[j]+50),nlocal)
;               if (nlocal ne 0L) then begin
;                  gsfit = mpfitpeak(colaxis[local],sky[local],a,nterm=nterm,$
;                    /gaussian,/positive,yfit=yfit)
;                  skypixels[j] = a[1]
;               endif
;               
;            endfor

;            skycenters = min(wave)+skypixels*dwave ; [Angstrom]

             skyshift = fix(skyshift*100.0)/100.0 ; two significant digits

          endelse 

          splog, 'Applying a '+strtrim(string(skyshift*dwave,format='(F12.2)'),2)+$
            ' Angstrom shift based on '+string(nskylines,format='(I0)')+' line(s).'
          
       endif else begin         ; sky line file does not exist

          splog, 'Sky line file '+etcpath+'skylines.dat not found . . . no sky shift applied.'
          skyshift = 0.0
          
       endelse

    endelse 

; update the header and add a history note

    wave = wave+skyshift*dwave  ; update the wavelength vector

    minwave = min(wave)
    sxaddpar, header, 'CRVAL1', float(minwave)
    sxaddpar, header, 'SKYSHIFT', float(skyshift*dwave), $
      ' sky-line wavelength shift [Angstrom]', before='HISTORY'

    if (abs(skyshift) gt 0.0) then sxaddhist, "'A sky line shift of "+$
      strtrim(string(skyshift,format='(F12.2)'),2)+" pixels was applied "+im_today()+"'", header

return, skyshift
end    
