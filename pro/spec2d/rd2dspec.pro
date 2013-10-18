;+
; NAME:
;       RD2DSPEC()
;
; PURPOSE:
;	Read a two dimensional FITS image into a structure cube. 
;
; CALLING SEQUENCE:
;	cube = rd2dspec(flist,datapath=,wset=,/silent)
;
; INPUTS:
;	flist - FITS file list
;
; OPTIONAL INPUTS:
;	datapath - path to the data to read in
;	
; KEYWORD PARAMETERS:
;	silent - suppress output to STDOUT
;
; OUTPUTS:
;	cube - data structure cube with the following fields:
;          fname     - file name [FLIST]
;          object    - object type (OBJECT)
;          type      - image type (IMAGETYP)
;          naxis1    - image length
;          naxis2    - image height
;          header    - FITS header (type pointer)
;          image     - two-dimensional FITS image
;          sigmamap  - corresponding two-dimensional error map
;          mask      - corresponding two-dimensional bad pixel map
;          sky       - corresponding two-dimensional sky image
;          wavemap   - corresponding two-dimensional wavelength map 
;          telluric_spec - telluric spectrum
;          telluric_wave - corresponding wavelength vector
;
; OPTIONAL OUTPUTS:
;          wset - trace set corresponding to WAVEMAP
;
; COMMENTS:
;       Use ICLEANUP() to free the pointer memory.  RETALL is called
;       if the images are not the same size. 
;
; EXAMPLE:
;
; PROCEDURES USED:
;       FITS_READ, SXPAR(), CWD(), FITS_INFO, MAKE_WAVE(),
;       STRUCT_ADDTAGS(), TRACESET2XY, DJS_LAXISGEN(), MRDFITS() 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2001 October 16, U of A
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm03dec08uofa - added support for reading in the
;                       two-dimensional sky image in the third FITS
;                       extension 
;       jm03dec24uofa - added support for reading in a telluric 
;                       spectrum in the fourth FITS extension
;       jm04jan14uofa - if the additional extensions do not exist then 
;                       do not try to read them in
;       jm05jun20uofa - simpler treatment of reading in the SKY and
;                       WAVEMAP extensions; experimented with using
;                       FITS_READ, but stayed with MRDFITS(); convert
;                       the 2D wavelength map trace set into a double
;                       precision 2D image; added WSET optional
;                       output; improved reading of the telluric 
;                       absorption spectrum, if present
;       jm06jan02uofa - bug on defining WSET1 reported by A. Marble 
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

function rd2dspec, flist, datapath=datapath, wset=wset, silent=silent

    if (n_elements(flist) eq 0L) then begin
       splog, 'Syntax - cube = rd2dspec(flist,datapath=,wset=,/silent)'
       retall
    endif
    
    if (size(flist[0],/type) ne 7L) or (strmatch(flist[0],'*.fits*') eq 0B) then begin
       splog, 'File list must be type string FITS files.'
       retall
    endif
    
    if (not keyword_set(datapath)) then datapath = cwd() ; current working directory

; read the first extension of the first image to get sizes

    if file_test(datapath+flist[0],/regular) eq 0L then begin
       splog, 'Image '+datapath+flist[0]+' not found.'
       retall
    endif

    if not keyword_set(silent) then splog, 'Reading ', flist[0]

;   fits_read, datapath+flist[0], image, header, exten_no=0L
    image = mrdfits(datapath+flist[0],0,header,/silent)

    if (size(image,/n_dimension) ne 2L) then begin
       splog, 'Image '+datapath+flist[0]+' is not two-dimensional.'
       retall
    endif
    
    imsz = size(image,/dimension)
    naxis1 = imsz[0]
    naxis2 = imsz[1]
    
    nflist = size(flist,/n_elements)
    template = {$
      fname:    '',                    $
      object:   '',                    $ ; object type (galaxy,standard,etc)
      type:     '',                    $ ; file type (zero,flat,comp,object)
      naxis1:   0L,                    $ ; x size
      naxis2:   0L,                    $ ; y size
      header:   ptr_new(),             $ ; header
      image:    fltarr(naxis1,naxis2), $ ; image
      sigmamap: fltarr(naxis1,naxis2), $ ; error map
      mask:     intarr(naxis1,naxis2), $ ; bad pixel mask (0 = good)
      sky:      fltarr(naxis1,naxis2), $ ; sky image
      wavemap:  dblarr(naxis1,naxis2)}   ; 2D wavelength map
      
    if (nflist gt 1L) then cube = replicate(template,nflist) else cube = template

    fits_info, datapath+flist[0], /silent, n_ext=n_ext

    cube[0].fname = flist[0]
    cube[0].object = string(sxpar(header,'OBJECT'))
    cube[0].type = string(sxpar(header,'IMAGETYP'))
    cube[0].naxis1 = long(sxpar(header,'NAXIS1'))
    cube[0].naxis2 = long(sxpar(header,'NAXIS2'))
    cube[0].image = image

    cube[0].header = ptr_new(header)

; error map    
    
    if (n_ext ge 1L) then begin
       cube[0].sigmamap = mrdfits(datapath+flist[0],1,/silent)
;      fits_read, datapath+flist[0], map, exten_no=1L
;      cube[0].sigmamap = temporary(map)
    endif

; bad pixel map    
    
    if (n_ext ge 2L) then begin
       cube[0].mask = mrdfits(datapath+flist[0],2,/silent)
;      fits_read, datapath+flist[0], map, exten_no=2L
;      cube[0].mask = temporary(map)
    endif

; sky spectrum    
    
    if (n_ext ge 3L) then begin
       cube[0].sky = mrdfits(datapath+flist[0],3,/silent)
;      fits_read, datapath+flist[0], map, exten_no=3L
;      cube[0].sky = temporary(map)
    endif

; wavelength map    
    
    if (n_ext ge 4L) then begin
       wset1 = mrdfits(datapath+flist[0],4,wh,/silent)
       if (size(wset1,/type) eq 8L) then begin
          nx = long(wset1.xmax-wset1.xmin+1L)
          ntrace = (size(wset1.coeff,/dimension))[1]
          xtemp = djs_laxisgen([nx,ntrace],iaxis=0L) + arm_double(wset1.xmin)
          traceset2xy, wset1, xtemp, map
          cube[0].wavemap = map
          wset = wset1
       endif
    endif

; telluric spectrum; MRDFITS() *must* be used over FITS_READ to
; properly handle the FITS extension header

    if (n_ext ge 5L) then begin
       telluric_spec = mrdfits(datapath+flist[0],5L,telluric_head,/silent)
;      fits_read, datapath+flist[0], telluric_spec, telluric_head, exten_no=5L
       telluric_wave = make_wave(telluric_head)
       ntellpix = n_elements(telluric_spec)
       cube = struct_addtags(cube,replicate({telluric_spec: fltarr(ntellpix), $
         telluric_wave: fltarr(ntellpix), telluric_head: $
         strarr(n_elements(telluric_head))},nflist))
       cube[0].telluric_spec = telluric_spec
       cube[0].telluric_wave = telluric_wave
       cube[0].telluric_head = telluric_head
    endif

; read each image

    for i = 1L, nflist-1L do begin
       
       if not keyword_set(silent) then splog, 'Reading ', flist[i]
       if file_test(datapath+flist[i],/regular) eq 0L then begin
          splog, 'Image '+datapath+flist[i]+' not found.'
          retall
       endif

;      fits_read, datapath+flist[i], image, header, exten_no=0L
       image = mrdfits(datapath+flist[i],0,header,/silent)

       szt = size(image,/dimension)
       if (size(image,/n_dimension) ne 2L) then begin
          splog, 'Image '+datapath+flist[i]+' is not two-dimensional.'
          heap_gc
          retall
       endif else if ((imsz[0] ne szt[0]) or (imsz[1] ne szt[1])) then begin
          splog, 'Images are not the same dimension!'
          heap_gc
          retall
       endif
       
       fits_info, datapath+flist[i], /silent, n_ext=n_ext

       cube[i].fname = flist[i]
       cube[i].object = string(sxpar(header,'OBJECT'))
       cube[i].type = string(sxpar(header,'IMAGETYP'))
       cube[i].naxis1 = long(sxpar(header,'NAXIS1'))
       cube[i].naxis2 = long(sxpar(header,'NAXIS2'))
       cube[i].image = image

       cube[i].header = ptr_new(header)

; error map       
       
       if (n_ext ge 1L) then begin
          cube[i].sigmamap = mrdfits(datapath+flist[i],1,/silent)
;         fits_read, datapath+flist[i], map, exten_no=1L
;         cube[i].sigmamap = temporary(map)
       endif

; bad pixel map

       if (n_ext ge 2L) then begin
          cube[i].mask = mrdfits(datapath+flist[i],2,/silent)
;         fits_read, datapath+flist[i], map, exten_no=2L
;         cube[i].mask = temporary(map)
       endif

; sky spectrum       

       if (n_ext ge 3L) then begin
          cube[i].sky = mrdfits(datapath+flist[i],3,/silent)
;         fits_read, datapath+flist[i], map, exten_no=3L
;         cube[i].sky = temporary(map)
       endif

; wavelength map    

       if (n_ext ge 4L) then begin

          wset1 = mrdfits(datapath+flist[i],4,/silent)
          traceset2xy, wset1, xtemp, map ; assumes XTEMP is the same as above!!
          cube[i].wavemap = temporary(map)

          wset = [ [wset], [wset1] ]

;         fits_read, datapath+flist[i], map, exten_no=4L
;         cube[i].wavemap = temporary(map)

       endif

; telluric spectrum; MRDFITS() *must* be used over FITS_READ to
; properly handle the FITS extension header

       if (n_ext ge 5L) then begin
          telluric_spec = mrdfits(datapath+flist[I],5L,telluric_head,/silent)
;         fits_read, datapath+flist[i], telluric_spec, telluric_head, exten_no=5L
          telluric_wave = make_wave(telluric_head)
          ntellpix = n_elements(telluric_spec)
          if (ntellpix ne n_elements(cube[0].telluric_spec)) then begin
             splog, 'Telluric spectra are not the same size!'
             heap_gc
             retall
          endif
          cube[i].telluric_spec = telluric_spec
          cube[i].telluric_wave = telluric_wave
          cube[i].telluric_head = telluric_head
       endif
    
    endfor
    
return, cube
end
