;+
; NAME:
;       IBACKMODEL()
;
; PURPOSE:
;       Construct a continuum model for IBACKFIT().
;
; INPUTS:
;       x - see IBACKFIT()
;       p - see IBACKFIT()
;
; OPTIONAL INPUTS:
;       wave     - see IBACKFIT()
;       redcurve - see IBACKFIT()
;       nstar    - see IBACKFIT()
;       nback    - see IBACKFIT()
;       ndust    - see IBACKFIT()
;
; KEYWORD PARAMETERS:
;       doplot - generate a plot for debugging
;
; OUTPUTS:
;       model - continuum model [see IBACKFIT()]
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 15, U of A - written
;       jm04mar11uofa - added ALPHA parameter and documented 
;
; Copyright (C) 2002, 2004, John Moustakas
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

function ibackmodel, x, p, redcurve=redcurve

    xdim = size(x,/dim)
    rdim = size(redcurve,/dim)
    if (total(xdim eq rdim) ne 2.0) then begin
       splog, 'Dimensions of X and REDCURVE do not match'
       return, -1
    endif
    
    ntemplate = xdim[1]
    if (ntemplate+1 ne n_elements(p)) then begin
       splog, 'Dimensions of X+1 and P do not match'
       return, -1
    endif
    newflux = x*10.0D^(-0.4D*p[ntemplate]*redcurve)
    model = newflux # p[0L:ntemplate-1L]

return, model
end

