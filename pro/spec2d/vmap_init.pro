;+
; NAME:
;       VMAP_INIT()
;
; PURPOSE:
;       Create a two dimensional variance map from an image and the
;       noise properties of the CCD.
;
; CALLING SEQUENCE:
;       vmap = vmap_init(image,gain=,rdnoise=)
;
; INPUTS:
;       image - input image
;
; OPTIONAL INPUTS:
;       gain    - CCD gain [electron per ADU] (default 1.0)
;       rdnoise - CCD read noise [electron] (default 0.0) 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       vmap - variance image for IMAGE
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
;       J. Moustakas, 2001 July 10, U of A
;       jm03apr09uofa - checked out ISPEC v1.0.0
;
; Copyright (C) 2001, 2003, John Moustakas
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

function vmap_init, image, gain=gain, rdnoise=rdnoise

    if (n_elements(gain) eq 0) then gain = 1.0
    if (n_elements(rdnoise) eq 0) then rdnoise = 0.0
    
    vmap = image/gain + (rdnoise/gain)^2 ; [ADU^2]

return, vmap
end
