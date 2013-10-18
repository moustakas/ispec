;+
; NAME:
;       ISPEC_CREATE_PNG
;
; PURPOSE:
;       Generate PNG images of iSPEC2d FITS images for display with
;       ISPEC_WEBPAGE.  
;
; CALLING SEQUENCE:
;       ispec_create_png, fits_file, png_file=, png_thumbfile=, $
;          file_ext=, /mask
;
; INPUTS:
;       fits_file - scalar FITS file name
;
; OPTIONAL INPUTS:
;       png_file      - output PNG file name
;       png_thumbfile - output PNG file name of the thumbnail image 
;       file_ext      - FITS file extension number (default is 0, or 2
;                       if MASK=1)
;
; KEYWORD PARAMETERS:
;       mask - set this keyword when converting bad pixel masks rather
;              than the data spectra, because the byte scaling must be
;              treated differently
;
; OUTPUTS:
;       PNG files (see usage in ISPEC_WEBPAGE).
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR(), READFITS()
;
; COMMENTS:
;
; EXAMPLES:
;
; TODO:
;       Allow a vector of FITS files.
;
; MODIFICATION HISTORY:
;       J. Moustakas 2003 July 15, U of A, written based in *large*
;          part on Karl Gordon's KGCREATE_PNG 
;       jm04sep27uofa - use MRDFITS() rather than FITS_OPEN 
;       jm05jun16uofa - documented and cleaned up; better image 
;                       scaling; problems using FITS_READ, so use
;                       MRDFITS() 
;       jm05jun30uofa - problems with MRDFITS() reading unsigned
;                       integer FITS files, so use READFITS()
;
; Copyright (C) 2003-2005, John Moustakas
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

pro ispec_create_png, fits_file, png_file=png_file, png_thumbfile=png_thumbfile, $
  file_ext=file_ext, mask=mask

    nfits = n_elements(fits_file)
    if (nfits ne 1L) then begin
       print, 'Syntax - ispec_create_png, fits_file, png_file=, file_ext=, /mask'
       return
    endif

    if (n_elements(png_file) eq 0L) then png_file = repstr(repstr(repstr(fits_file,$
      '.fits','.png'),'.gz',''),'.fit','.png')
    if (n_elements(png_thumbfile) eq 0L) then png_thumbfile = repstr(png_file,'.png','_thumb.png')
    
    if (n_elements(file_ext) eq 0L) then file_ext = 0L

; read the image

    if (file_test(fits_file) eq 0L) then begin
       print, 'FITS File '+fits_file+' does not exist.'
       return
    endif

;   fits_read, fits_file, image, header, exten_no=file_ext ; DO NOT USE
;   image = mrdfits(fits_file,file_ext,header,/silent)     ; DO NOT USE
    image = float(readfits(fits_file,header,exten_no=file_ext,/silent))

; generate a full-scale PNG image and a thumb-sized PNG file
    
    disp_image = image

    size_image = size(image)
    if (size_image[1] gt size_image[2]) then begin
       ysize = 150L
       xsize = ysize*size_image[1]/size_image[2]
    endif else begin
       xsize = 150L
       ysize = xsize*size_image[2]/size_image[1]
    endelse

    disp_thumbimage = congrid(image,xsize,ysize)

; scale MASK images separately

    if keyword_set(mask) then begin

       disp_image = (disp_image gt 0.0)+1.0
       disp_thumbimage = (disp_thumbimage gt 0.0)+1.0

       scaled_image = bytscl(disp_image,min=1,max=2)
       scaled_thumbimage = bytscl(disp_image,min=1,max=2)

    endif else begin

       minvalue = -1.0
       maxvalue = +2.0 ; 2-sigma looks better than 1- or 3-sigma
       topvalue = !d.table_size-2L

       disp_image = (disp_image>0)^(1.0/2.0)
       stats = im_stats(disp_image)
       disp_image = (disp_image-stats.mean)/stats.sigma

       scaled_image = bytscl(disp_image,min=minvalue,max=maxvalue,top=topvalue)
       scaled_image = bytscl(topvalue-scaled_image,min=-10L,max=topvalue)

; thumbnail

       disp_thumbimage = (disp_thumbimage>0)^(1.0/2.0)
       stats = im_stats(disp_thumbimage)
       disp_thumbimage = (disp_thumbimage-stats.mean)/stats.sigma
              
       scaled_thumbimage = bytscl(disp_thumbimage,min=minvalue,max=maxvalue,top=topvalue)
       scaled_thumbimage = bytscl(topvalue-scaled_thumbimage,min=-10L,max=topvalue)
       
;      scaled_image = im_imgscl(disp_image,losig=losig,hisig=hisig,$
;        boxfrac=1.0,/sqrroot,log=0,/negative)
;      scaled_thumbimage = im_imgscl(disp_thumbimage,losig=losig,$
;        hisig=hisig,boxfrac=1.0,/sqrroot,log=0,/negative)

    endelse

; write out    
    
    loadct, 0, /silent
    tvlct, r, g, b, /get
    write_png, png_file, scaled_image, r, g, b
    write_png, png_thumbfile, scaled_thumbimage, r, g, b

return
end
