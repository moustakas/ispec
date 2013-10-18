;+
; NAME:
;       IMASK_BITS()
;
; PURPOSE:
;       Return the bad pixel mask bit values assigned to each type of
;       bad pixel.
;
; CALLING SEQUENCE:
;       bit = imask_bits(maxpower=,/badpix,/saturation,$
;          /crsplits,/crpix,/skymask)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       maxpower - 
;
; KEYWORD PARAMETERS:
;       badpix     - dead pixel (=2) [ICCDPROC()]
;       saturation - saturated pixel (=2) [ICCDPROC()]
;       crsplits   - deviant pixel when co-adding cr-splits (=4) 
;                    [ICRCOMBINE] 
;       crpix      - cosmic ray pixel (=8) [ILA_COSMIC_WRAPPER] 
;       skymask    - sky subtraction residual pixel (=32) [ISKYMASK()]  
;
; OUTPUTS:
;       bit - bit number
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Dec 08, U of A
;       jm05apr19uofa - documented
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

function imask_bits, badpix=badpix, saturation=saturation, crsplits=crsplits, $
  crpix=crpix, skymask=skymask, maxpower=maxpower

    bit = 0
    maxpower = 5L ; this needs to be updated [see ICALIBRATE]

    if keyword_set(badpix)     then bit = fix(2^1) ; = 2 dead pixel [ICCDPROC()]
    if keyword_set(saturation) then bit = fix(2^1) ; = 2 saturated pixel [ICCDPROC()]
    if keyword_set(crsplits)   then bit = fix(2^2) ; = 4 deviant pixel when co-adding cr-splits [ICRCOMBINE]
    if keyword_set(crpix)      then bit = fix(2^3) ; = 8 cosmic ray pixel [ILA_COSMIC_WRAPPER]
                                                   ; = 16 extrapolated pixel in ICALIBRATE
    if keyword_set(skymask)    then bit = fix(2^5) ; = 32 sky subtraction residual pixels [ISKYMASK()] 

return, bit
end
