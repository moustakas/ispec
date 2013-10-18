;+
; NAME:
;	ITRIM
;
; PURPOSE:
;	Trim an image and update the header.
;
; CALLING SEQUENCE:
;       itrim, image, header=, trim=
;
; INPUTS:
;	image  - two dimensional image
;
; OPTIONAL INPUTS:
;	header - corresponding FITS header
;	trim   - [x1,x2,y1,y2] defining the good data region
;                (zero-indexed)
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	image  - (modified)
;	header - (if passed, then modified)
;
; COMMENTS:
;	If no header or trim region is passed then make the trim
;	region the whole image.
;
; EXAMPLE:
;	itrim, image, header=imheader, trim=[0,1200,0,125]
;
; PROCEDURES USED:
;	SXADDPAR, SXADDHIST
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 15, U of A
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm05jun17uofa - simplify the TRIM header keyword (no longer
;                       follows the IRAF convention); remove the
;                       CCDSEC header keyword
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

pro itrim, image, header=header, trim=trim

    imsize = size(image,/dimension)

;   on_error, 1
;   if (not keyword_set(header)) and (not keyword_set(trim)) then $
;     message, 'Need to pass a header or a trim region!'
;   if keyword_set(trim) then xy = long(trim) else $
;     xy = long(strsplit(sxpar(header,'TRIMSEC'),'[:,] ',/extract))-1L

    if n_elements(trim) eq 0L then xy = [0,imsize[0]-1,0,imsize[1]-1] else $
      xy = long(trim)
    
    xy[0] = xy[0] > 0L
    xy[1] = xy[1] < (imsize[0]-1L)
    xy[2] = xy[2] > 0L
    xy[3] = xy[3] < (imsize[1]-1L)

    ncols = xy[1]-xy[0]+1L
    nrows = xy[3]-xy[2]+1L    

    origimage = image ; original data

    image = origimage[xy[0]:xy[1],xy[2]:xy[3]] ; trim

; update the header.  IRAF indices start at one

    if keyword_set(header) then begin

       sxdelpar, header, 'DATASEC'
       sxdelpar, header, 'TRIMSEC'

       sxaddpar, header, 'NAXIS1', ncols
       sxaddpar, header, 'NAXIS2', nrows

       sxaddpar, header, 'TRIM', '['+strn(xy[0]+1)+':'+strn(xy[1]+1)+','+strn(xy[2]+1)+$
         ':'+strn(xy[3]+1)+']', ' trim region', before='HISTORY'

;      sxaddpar, header, 'TRIM', im_today()+' Trim section '+$
;        '['+strn(xy[0]+1)+':'+strn(xy[1]+1)+','+strn(xy[2]+1)+$
;        ':'+strn(xy[3]+1)+'].', before='HISTORY'
;      xy = xy+1L
;      sxaddpar, header, 'CCDSEC', '['+strn(xy[0])+':'+strn(xy[1])+','+$
;        strn(xy[2])+':'+strn(xy[3])+']'
       
    endif
    
return
end
