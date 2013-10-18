;+
; NAME:
;       REPAIR_LA_COSMIC_CRPIX
;
; PURPOSE:
;       Restore pixels improperly flagged as cosmic rays by
;       LA_COSMIC. 
;
; CALLING SEQUENCE:
;       repair_la_cosmic_crpix, repairfile, crpixfile, datapath=, /repair
;
; INPUTS:
;       repairfile - list of 2D spectra to repair [NREPAIR]
;       crpixfile  - pixel list to restore (see BADPIXFILE in
;                    ICCDPROC) [NREPAIR]
;
; OPTIONAL INPUTS:
;       datapath - path name to REPAIRFILE and CRPIXFILE
;
; KEYWORD PARAMETERS:
;       repair - modify the FITS file on disk with the good pixels
;                restored 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       MRDFITS(), SXPAR(), READCOL, MODFITS
;
; COMMENTS:
;       Occasionally, LA_COSMIC will flag pixels as cosmic rays even
;       though they are not.  This routine reads a list of pixels and
;       restores them to their pre-LA COSMIC values and updates the
;       FITS file on disk.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 24, U of A
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

pro repair_la_cosmic_crpix, repairfile, crpixfile, datapath=datapath, repair=repair

    nrepair = n_elements(repairfile)
    ncrpixfile = n_elements(crpixfile) 

    if (nrepair eq 0L) or (ncrpixfile eq 0L) then begin
       print, 'Syntax - repair_la_cosmic_crpix, repairfile, crpixfile, $'
       print, '   datapath=, /repair'
       return
    endif

    if (nrepair ne ncrpixfile) then begin
       splog, 'REPAIRFILE and CRPIXFILE must have the same number of elements.'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()

; call this routine recursively
    
    if (nrepair gt 1L) then begin
       for iobj = 0L, nrepair-1L do begin
          repair_la_cosmic_crpix, repairfile[iobj], crpixfile[iobj], $
            datapath=datapath, repair=repair
       endfor
       return
    endif

    if (file_test(datapath+crpixfile,/regular) eq 0L) then begin
       splog, 'CRPIXFILE '+datapath+crpixfile+' not found.'
       return
    endif

    if (file_test(datapath+repairfile,/regular) eq 0L) then begin
       splog, 'FITS file '+datapath+repairfile+' not found.'
       return
    endif

    splog, 'Reading '+datapath+repairfile+'.'
    repair_image = mrdfits(datapath+repairfile,0,header,/silent)
    repair_mask = mrdfits(datapath+repairfile,2,/silent)
    imsz = size(repair_image,/dimension)

    restore_file = sxpar(header,'ICOSMIC',count=ccount)
    if (ccount eq 0L) then begin
       splog, 'FITS File '+datapath+repairfile+' has not been cosmic-ray rejected.'
       return
    endif

    restore_image = mrdfits(datapath+restore_file,0,/silent)
    
; read the pixel file (zero-indexed)
    
    readcol, datapath+crpixfile, x1, x2, y1, y2, format='L,L,L,L', /silent, comment='#'
    xybad = long(transpose([[x1],[x2],[y1],[y2]]))

    ndim = size(xybad,/n_dimension)
    dims = size(xybad,/dimension)

    if (ndim eq 1L) then nlines = 1L else nlines = dims[1]
    splog, 'Reading '+string(nlines,format='(I0)')+' line(s) in pixel file '+crpixfile+'.'

    if (ndim[0] eq 1L) then jmax = 0L else jmax = dims[1]-1L
    
    for j = 0L, jmax do begin
       xy = xybad[*,j]
       repair_image[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)] = $
         restore_image[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)]
       repair_mask[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)] = 0
    endfor

    if keyword_set(repair) then begin
       splog, 'Repairing '+datapath+repairfile+'.'
       modfits, datapath+repairfile, repair_image, header, exten_no=0L
       modfits, datapath+repairfile, repair_mask, exten_no=2L
    endif

return
end
    
