;+
; NAME:
;	RD1DSPEC()
;
; PURPOSE:
;	Read ISPEC one dimensional spectra into a structure. 
;
; INPUTS:
;	speclist  - FITS file list
;
; OPTIONAL INPUTS:
;	datapath  - I/O path
;       width     - see IM_NORMALIZE() (default 50 Angstrom)
;	
; KEYWORD PARAMETERS:
;	silent    - suppress message output to STDOUT
;       normalize - normalize by the mean flux around the mean
;                   wavelength +/- WIDTH Angstroms
;
; OUTPUTS:
;       speccube  - output data structure
;          specname - spectrum name
;          datapath - data path
;          object   - object name
;          header   - FITS header
;          npix     - number of pixels
;          spec     - spectrum [erg/s/cm2/A]
;          sigspec  - error spectrum [erg/s/cm2/A]
;          sky      - sky spectrum [erg/s/cm2/A]
;          mask     - bad pixel mask
;          wave     - wavelength vector [A]
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 19, U of A
;       jm03apr17uofa - bug fix on reading the header - thanks to
;                       A. Marble; added WIDTH input
;       jm04jan04uofa - added FITS_INFO error checking
;       jm04jun16uofa - use double precision
;
; Copyright (C) 2001, 2003-2004, John Moustakas
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

function rd1dspec, speclist, width=width, datapath=datapath, $
  silent=silent, normalize=normalize

    if n_params() ne 1L then begin
       doc_library, 'rd1dspec'
       retall
    endif
    
    if (size(speclist[0],/type) ne 7L) or (strmatch(speclist[0],'*.fits*') eq 0B) then begin
       splog, 'SPECLIST list must be type string FITS files.'
       retall
    endif
    
    if not keyword_set(datapath) then datapath = './'

    if n_elements(width) eq 0L then width = 50.0 ; [Angstrom]
    
    speclist = strcompress(speclist,/remove)
    
; read the first extension of the first spectrum to initialize arrays

    if not keyword_set(silent) then splog, 'Reading '+speclist[0]+'.'
    spec = mrdfits(datapath+speclist[0],0,header,/silent)
    if (size(spec,/n_dimension) ne 1L) then begin
       splog, 'Spectrum '+strn(speclist[0])+' is not one-dimensional!'
       retall
    endif
    
    specsize = size(spec,/dimension)
    npix = specsize[0]
    
    nspec = n_elements(speclist)

    if nspec eq 1L then h = strarr(n_elements(header)) else h = ptr_new()
    speccube = {specname: '',           $ ; spectrum name
                datapath: '',           $ ; data path
                object:   '',           $ ; object name
                header:    h,           $ ; header
                npix:     0L,           $ ; number of pixels
                spec:     dblarr(npix), $ ; spectrum
                sigspec:  dblarr(npix), $ ; one-sigma spectrum
                sky:      dblarr(npix), $ ; sky spectrum
                mask:     intarr(npix), $ ; bad pixel mask
                wave:     dblarr(npix)}   ; wavelength vector
    if nspec gt 1L then speccube = replicate(speccube,nspec)

    fits_info, datapath+speclist[0], /silent, n_ext=n_ext

    speccube[0].specname = speclist[0]
    speccube[0].datapath = datapath
    speccube[0].object = strn(sxpar(header,'OBJECT'))
    speccube[0].npix = long(sxpar(header,'NAXIS1'))
    speccube[0].spec = spec
    if (n_ext ge 1L) then speccube[0].sigspec = mrdfits(datapath+speclist[0],1,/silent)
    if (n_ext ge 2L) then speccube[0].sky = mrdfits(datapath+speclist[0],2,/silent)
    if (n_ext ge 3L) then speccube[0].mask = mrdfits(datapath+speclist[0],3,/silent)
    speccube[0].wave = make_wave(header)

    if nspec eq 1L then speccube[0].header = header else speccube[0].header = ptr_new(header)

; read each spectrum

    for i = 1L, nspec-1L do begin
       
       if not keyword_set(silent) then splog, 'Reading ', speclist[i]+'.'
       spec = mrdfits(datapath+speclist[i],0,header,/silent)

       szt = size(spec,/dimension)
       if (size(spec,/n_dimension) ne 1L) then begin
          splog, 'Spectrum '+strn(speclist[i])+' is not one-dimensional!'
          heap_gc
          retall
       endif else if (specsize[0] ne szt[0]) then begin
          splog, 'Spectra are not the same dimension!'
          heap_gc
          retall
       endif
       
       speccube[i].specname = speclist[i]
       speccube[i].datapath = datapath
       speccube[i].object = strn(sxpar(header,'OBJECT'))
       speccube[i].header = ptr_new(header)
       speccube[i].npix = long(sxpar(header,'NAXIS1'))
       speccube[i].spec = spec
       if (n_ext ge 1L) then speccube[i].sigspec = mrdfits(datapath+speclist[i],1,/silent)
       if (n_ext ge 2L) then speccube[i].sky = mrdfits(datapath+speclist[i],2,/silent)
       if (n_ext ge 3L) then speccube[i].mask = mrdfits(datapath+speclist[i],3,/silent)
       speccube[i].wave = make_wave(header)

    endfor

    if keyword_set(normalize) then begin

        for i = 0L, nspec-1L do begin

           speccube[i].spec = im_normalize(speccube[i].spec,speccube[i].wave,$
             normwave=djs_mean(speccube[i].wave),binsize=2*width,const=norm,/mean)
           speccube[i].sigspec = speccube[i].sigspec/norm

        endfor

    endif

return, speccube
end
