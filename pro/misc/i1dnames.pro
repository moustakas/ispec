;+
; NAME:
;       I1DNAMES()
;
; PURPOSE:
;       Convert a 2D FITS file name to a 1D names, given the
;       extraction aperture.
;
; CALLING SEQUENCE:
;       newname = i1dnames(name,aperture=)
;
; INPUTS:
;       name - 2D FITS file name
;
; OPTIONAL INPUTS:
;       aperture - extraction aperture in arcsec
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       newname - iSPEC2d-formatted 1D FITS file name 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR()
;
; COMMENTS:
;       Converts wra.0121.fits --> wra.0121_030.ms.fits.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Jan 05, U of A - written
;       jm05jun17uofa - documented
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

function i1dnames, name, aperture=aperture

    nname = n_elements(name)
    if (nname eq 0L) then begin
       print, 'Syntax - newname = i1dnames(name,aperture=)'
       return, -1L
    endif

    if (n_elements(aperture) eq 0L) then aperture = 30.0
    
;   newname = repstr('w'+name,'.fits','_'+fmtaperture(aperture)+'.ms.fits')
    newname = repstr(name,'.fits','_'+fmtaperture(aperture)+'.ms.fits')
    
return, newname
end
