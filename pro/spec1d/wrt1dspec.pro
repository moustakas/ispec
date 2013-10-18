;+
; NAME:
;	WRT1DSPEC
;
; PURPOSE:
;	Write an iSPEC-format one-dimensional spectrum.
;
; INPUTS:
;	specname - output file name
;	spec     - spectrum
;	sigspec  - sigma spectrum
;	sky      - sky spectrum
;	mask     - bad pixel mask
;	header   - FITS header
;
; OPTIONAL INPUTS:
;	datapath - I/O data path
;	
; KEYWORD PARAMETERS:
;       nobadmask - do not write out the bad pixel mask spectrum 
;       gzip      - compress the output FITS file using SPAWN and GZIP 
;
; OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 21, U of A - written
;       jm04sep20uofa - NOBADMASK keyword added
;       jm08jan20nyu - enforce that all components of the spectra are
;                      written out in floating-point format; some
;                      documentation cleanup
;
; Copyright (C) 2002, 2004, 2008, John Moustakas
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

pro wrt1dspec, specname, spec, sigspec, sky, mask, header, $
  datapath=datapath, nobadmask=nobadmask, gzip=gzip

    if n_elements(specname) eq 0L then begin
       doc_library, 'wrt1dspec'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()

    mwrfits, float(spec), datapath+specname, header, /create
    mwrfits, float(sigspec), datapath+specname
    mwrfits, float(sky), datapath+specname
    if (not keyword_set(nobadmask)) then mwrfits, byte(mask), datapath+specname

    if keyword_set(gzip) then spawn, ['gzip -f '+datapath+specname], /sh

return
end
