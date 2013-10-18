;+
; NAME:
;       BACKMODEL()
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

function backmodel, x, p, wave=wave, redcurve=redcurve, nstar=nstar, $
  nback=nback, ndust=ndust, normstar=normstar, doplot=doplot

    ntemplate = nstar+nback
    newflux = x
    nalpha = ndust

    if n_elements(normstar) eq 0L then normstar = dblarr(ntemplate)+1.0
    
; redden each stellar template when NDUST > 0

    for k = 0L, ndust-1L do begin
       alpha = p[ntemplate+ndust+k]
       ebv = p[ntemplate+k]
       newflux[*,k] = x[*,k] * 10.0^(-0.4D*alpha*ebv*redcurve)
    endfor

; make the model

    coeff = p[0L:ntemplate-1L] / normstar
    model = newflux # coeff

    if keyword_set(doplot) then begin
       plot, wave, newflux[*,0], xsty=3, ysty=3, ps=3
       for k = 1L, ntemplate-1L do oplot, wave, newflux[*,k], ps=3
       oplot, wave, model, ps=3, color=djs_icolor('yellow')
       cc = get_kbrd(1)
    endif
    
return, model
end

