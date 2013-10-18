;+
; NAME:
;       FMTAPERTURE()
;
; PURPOSE:
;       Convert a floating-point extraction aperture into a string in 
;       iSPEC format.
;
; CALLING SEQUENCE:
;       strap = fmtaperture(aperture)
;
; INPUTS:
;       aperture - aperture [arcsec]
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       strap - string-formatted aperture
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2003 Jan 05, U of A
;
; Copyright (C) 2003, John Moustakas
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

function fmtaperture, aperture

    naperture = n_elements(aperture)
    if (naperture eq 0L) then begin
       print, 'Syntax - strap = fmtaperture(aperture)'
       return, -1L
    endif

    if (naperture gt 1L) then begin
       strap = strarr(naperture)
       for iap = 0L, naperture-1L do begin
          strap[iap] = fmtaperture(aperture[iap])
       endfor
       return, strap
    endif
    
    strap = string(aperture,format='(G0)')

    if (float(strap) mod fix(strap)) eq 0.0 then strap = string(strap,format='(I3.3)') else begin

       if (float(strap) lt 100.0) and (float(strap) ge 10.0) then strap = '0'+strap
       if float(strap) lt 10.0 then strap = '00'+strap

    endelse
    
return, strap
end
