;+
; NAME:
;       BALMER_MASK()
;
; PURPOSE:
;       Mask out Balmer absorption lines from a wavelength array. 
;
; CALLING SEQUENCE:
;       mask = balmer_mask(wave,good=,bad=,bdata=)
;
; INPUTS:
;       wave  - wavelength vector [Angstrom]
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       mask  - pixel mask for WAVE (1B is good and 0B is bad) 
;
; OPTIONAL OUTPUTS:
;       good  - indices in WAVE that have not been masked
;       bad   - indices in WAVE that have been masked
;       bdata - data structure from BALMER_DATA()
;
; COMMENTS:
;
; PROCEDURES USED:
;       BALMER_DATA()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 November 27, U of A
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

function balmer_mask, wave, good=good, bad=bad, bdata=bdata

    nwave = n_elements(wave)
    if nwave eq 0L then begin
       print, 'Syntax - mask = balmer_mask(wave,good=,bad=,bdata=)'
       return, -1
    endif

    mask = make_array(nwave,/byte,value=1)

    bdata = balmer_data(nbdata=nbalmer)

    for i = 0L, nbalmer-1L do $
      mask = mask and $
        ((wave lt bdata[i].lline) or (wave gt bdata[i].uline))
;       ((wave lt bdata[i].llimit+bdata[i].lwidth/2.0) or $
;        (wave gt bdata[i].ulimit-bdata[i].uwidth/2.0))
    
    good = where(mask eq 1B,ngood,comp=bad,ncomp=nbad)

return, mask
end
