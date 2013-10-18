;+
; NAME:
;	IRDSENSFUNC()
;
; PURPOSE:
;	Read one or more sensitivity functions.
;
; CALLING SEQUENCE:
;       senscube = irdsensfunc(senslist,datapath=,/silent)
;
; INPUTS:
;	senslist    - sensitivity function file list (non
;                     grey-shifted)
;
; OPTIONAL INPUTS:
;	datapath  - path to the data to read in
;	
; KEYWORD PARAMETERS:
;	silent    - suppress output to STDOUT
;
; OUTPUTS:
;       senscube  - output data structure
;
; COMMENTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;	MRDFITS(), MAKE_WAVE(), CWD()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2003 April 24, U of A, written, ISPEC v1.0.0
;       jm03dec24uofa - added support for reading in the telluric
;                       spectrum in the second FITS extension
;       jm05jun23uofa - the sensitivity functions no longer handle the
;                       telluric spectra
;       jm05jun29uofa - also read in the header for the grey-shifted
;                       sensitivity function
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

function irdsensfunc, senslist, datapath=datapath, silent=silent

    nsens = n_elements(senslist)
    if (nsens eq 0L) then begin
       splog, 'Syntax - senscube = irdsensfunc(senslist,datapath=,/silent)'
       return, -1L
    endif
    
    if (size(senslist[0],/type) ne 7L) or (strmatch(senslist[0],'*.fits*') eq 0B) then begin
       splog, 'SENSLIST list must be type string FITS files.'
       return, -1L
    endif
    
    if (n_elements(datapath) eq 0L) then datapath = cwd()
    
    senslist = strcompress(senslist,/remove)
    
; check to see if each element in SENSLIST is the observed or the
; grey-shifted sensitivity function

    slist = strarr(nsens)
    gslist = strarr(nsens)
    
    for j = 0L, nsens-1L do begin
       
       if strmatch(senslist[j],'*grey*') eq 1B then begin
          gslist[j] = senslist[j]
          slist[j] = repstr(senslist[j],'_grey.fits','.fits')
       endif else begin
          slist[j] = senslist[j]
          gslist[j] = repstr(senslist[j],'.fits','_grey.fits')
       endelse

    endfor

;   greysenslist = repstr(senslist,'.fits','_grey.fits')
    
; read the first extension of the first sensitivity function to
; initialize the arrays 

    if file_test(datapath+slist[0],/regular) eq 0L then begin
       splog, 'Sensitivity function '+datapath+slist[0]+' not found.'
       return, -1L
    endif else begin
       if not keyword_set(silent) then splog, 'Reading '+slist[0]+'.'
       sens = mrdfits(datapath+slist[0],0,header,/silent)
    endelse

    if (size(sens,/n_dimension) ne 1L) then begin
       splog, 'Sensitivity function '+strn(slist[0])+' is not one-dimensional!'
       return, -1L
    endif
    
    if file_test(datapath+gslist[0],/regular) eq 0L then begin
       splog, 'Grey-shifted sensitivity function '+datapath+gslist[0]+' not found.'
       return, -1L
    endif else begin
       if not keyword_set(silent) then splog, 'Reading '+gslist[0]+'.'
       greysens = mrdfits(datapath+gslist[0],0,greyheader,/silent)
    endelse
    
    senssize = size(sens,/dimension)
    npix = senssize[0]

    if (nsens eq 1L) then begin
       h = strarr(n_elements(header))
       greyh = strarr(n_elements(greyheader))
    endif else begin
       h = ptr_new()
       greyh = ptr_new()
    endelse

    senscube = {$
      sensname:      '',           $ ; sensitivity function name
      greysensname:  '',           $ ; grey-shifted sensitivity function name
      datapath:      '',           $ ; data path
      header:        h,            $ ; header
      greyheader:    greyh,        $ ; grey-shifted sensitivity function header
      npix:          0L,           $ ; number of pixels
      sens:          fltarr(npix), $ ; sensitivity function
      senserr:       fltarr(npix), $ ; sensitivity function error
      greysens:      fltarr(npix), $ ; grey-shifted sensitivity function
      greysenserr:   fltarr(npix), $ ; grey-shifted sensitivity function error
      wave:          fltarr(npix)}   ; wavelength vector
    if nsens gt 1L then senscube = replicate(senscube,nsens)

    senscube[0].sensname = slist[0]
    senscube[0].greysensname = gslist[0]
    senscube[0].datapath = datapath
    senscube[0].npix = npix
    senscube[0].sens = sens
    senscube[0].senserr = mrdfits(datapath+slist[0],1,/silent)
    senscube[0].greysens = greysens
    senscube[0].greysenserr = mrdfits(datapath+gslist[0],1,/silent)
    senscube[0].wave = make_wave(header)
    
    if (nsens eq 1L) then begin
       senscube[0].header = header 
       senscube[0].greyheader = greyheader 
    endif else begin
       senscube[0].header = ptr_new(header)
       senscube[0].greyheader = ptr_new(greyheader)
       ptr_free, h, greyh
    endelse

; read each spectrum

    for i = 1L, nsens-1L do begin
       
       if file_test(datapath+slist[i],/regular) eq 0L then begin
          splog, 'Sensitivity function '+datapath+slist[i]+' not found.'
          icleanup, senscube
          return, -1L
       endif else begin
          if not keyword_set(silent) then splog, 'Reading ', slist[i]+'.'
          sens = mrdfits(datapath+slist[i],0,header,/silent)
       endelse

       szt = size(sens,/dimension)
       if (size(sens,/n_dimension) ne 1L) then begin
          splog, 'Sensitivity function '+strn(slist[i])+' is not one-dimensional!'
          icleanup, senscube
          return, -1L
       endif else if (senssize[0] ne szt[0]) then begin
          splog, 'Sensitivity functions are not the same dimension!' ; <-- could generalize this
          icleanup, senscube
          return, -1L
       endif
       
       if file_test(datapath+gslist[i],/regular) eq 0L then begin
          splog, 'Grey-shifted sensitivity function '+datapath+gslist[i]+' not found.'
          return, -1L
       endif else begin
          if not keyword_set(silent) then splog, 'Reading ', gslist[i]+'.'
          greysens = mrdfits(datapath+gslist[i],0,greyheader,/silent)
       endelse

       senscube[i].sensname = slist[i]
       senscube[i].greysensname = gslist[i]
       senscube[i].datapath = datapath
       senscube[i].header = ptr_new(header)
       senscube[i].greyheader = ptr_new(greyheader)
       senscube[i].npix = npix
       senscube[i].sens = sens
       senscube[i].senserr = mrdfits(datapath+slist[i],1,/silent)
       senscube[i].greysens = greysens
       senscube[i].greysenserr = mrdfits(datapath+gslist[i],1,/silent)
       senscube[i].wave = make_wave(header)

    endfor

return, senscube
end
