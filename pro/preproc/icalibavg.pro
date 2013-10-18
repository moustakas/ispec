;+
; NAME:
;       ICALIBAVG
;
; PURPOSE:
;       Average calibration data with iterative sigma-clipping. 
;
; CALLING SEQUENCE:
;       icalibavg, caliblist, datapath=, sigrej=, outname=, $
;          calibavg=, /wfits, _extra=extra
;
; INPUTS:
;       caliblist - list of FITS files to average
;
; OPTIONAL INPUTS:
;       datapath  - I/O path
;       sigrej    - sigma rejection threshold (default 3.0)
;       outname   - output FITS file name of the average image 
;       extra     - extra parameters for DJS_AVSIGCLIP()
;
; KEYWORD PARAMETERS:
;       wfits     - write CALIBAVG to DATAPATH
;
; OUTPUTS:
;       See WFITS and CALIBAVG.
;
; OPTIONAL OUTPUTS:
;       calibavg  - average image
;
; COMMENTS:
;       This routine is basically a wrapper for DJS_AVSIGCLIP(), and
;       is meant to replace IRAF's IMCOMBINE task.
;
;       Note that all the images in CALIBLIST must be the same size. 
;
; PROCEDURES USED:
;       CWD(), IM_FITS_CUBE(), DJS_AVSIGCLIP(), SPLOG, WRITEFITS,
;       SXADDPAR 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 May 27, U of A, ISPEC v1.0.0 
;       jm05apr17uofa - documentation improved
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

pro icalibavg, caliblist, datapath=datapath, sigrej=sigrej, outname=outname, $
  calibavg=calibavg, wfits=wfits, _extra=extra

    ncalib = n_elements(caliblist)
    if ncalib eq 0L then begin
       print, 'Syntax - icalibavg, caliblist, datapath=, sigrej=, outname=, $'
       print, '   calibavg=, /wfits, _extra=extra'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(sigrej) eq 0L then sigrej = 3.0
    if n_elements(outname) eq 0L then outname = 'calibavg.fits'

    calibcube = im_fits_cube(caliblist,datapath=datapath)    
    calibavg = djs_avsigclip(calibcube.image,3,sigrej=sigrej,_extra=extra)

    header = calibcube[0].header
    sxaddpar, header, 'NCOMBINE', ncalib

    notblank = where(strcompress(header,/remove) ne '',nblank)
    if nblank ne 0L then header = header[notblank]
    
    if keyword_set(wfits) then begin
       splog, 'Writing '+datapath+outname+'.'
       writefits, datapath+outname, calibavg, header
    endif
    
return
end    
