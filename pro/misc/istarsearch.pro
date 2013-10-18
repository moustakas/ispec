;+
; NAME:
;       ISTARSEARCH()
;
; PURPOSE:
;       Find all objects that are in the standard-star database based
;       on their position.
;
; CALLING SEQUENCE:
;       starindx = istarsearch(ra,dec,epoch=,searchrad=)
;
; INPUTS:
;       ra  - right ascension (string, HMS) [NOBJ]
;       dec - declination (string, DMS) [NOBJ]
;
; OPTIONAL INPUTS:
;       epoch     - scalar coordinate epoch (default 2000.0)
;       searchrad - search radius (default 300.0 arcsec)
;       starinfo  - standard star database FITS table
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       starindx - indices of the standard stars in the coordinate
;                  list provided [NOBJ]
;
; OPTIONAL OUTPUTS:
;       objindx - index of RA and DEC corresponding to the matched
;                 stars 
;
; PROCEDURES USED:
;       IM_HMS2DEC(), PRECESS, MRDFITS(), IM_DJS_ANGLE_MATCH()
;
; COMMENTS:
;       If EPOCH is not 2000.0 then the coordinates are precessed. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 6, U of A
;       jm05apr14uofa - excised from ISPEC_MAKELISTS(), documented,
;                       and generalized
;       jm05jun21uofa - added STARINFO optional input; if STARINFO has
;                       been read in previously, it can be passed here
;                       to prevent it from being read in again
;       jm05jun28uofa - added OBJINDX optional output
;       jm05jul26uofa - bug fix when precessing coordinates (thanks to
;                       A. Marble)
;
; Copyright (C) 2002, 2005, John Moustakas
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

function istarsearch, ra, dec, epoch=epoch, searchrad=searchrad, $
  objindx=objindx, starinfo=starinfo

    nobj = n_elements(ra)
    if (nobj eq 0L) then begin
       print, 'Syntax - starindx = istarsearch(ra,dec,epoch=,$'
       print, '   searchrad=,objindx=,starinfo=)'
       return, -1L
    endif

    if (nobj ne n_elements(dec)) then begin
       print, 'RA and DEC must have the same number of elements.'
       return, -1L
    endif

    if (n_elements(epoch) eq 0L) then epoch = replicate(2000.0,nobj)

    if (nobj ne n_elements(epoch)) then begin
       print, 'RA and EPOCH must have the same number of elements.'
       return, -1L
    endif

    if (n_elements(searchrad) eq 0L) then searchrad = 300.0 ; [arcsec]
    
    raref = 15.0*im_hms2dec(ra) ; [degrees]
    decref = im_hms2dec(dec)    ; [degrees]

    doit = where(epoch ne 2000.0,ndoit)
    if (ndoit ne 0L) then for idoit = 0L, ndoit-1L do begin
       tra = raref[doit[idoit]] & tdec = decref[doit[idoit]]
       precess, tra, tdec, epoch[doit[idoit]], 2000.0
       raref[doit[idoit]] = tra & decref[doit[idoit]] = tdec
    endfor
    
    standardspath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/spectra')
    if (n_elements(starinfo) eq 0L) then $
      starinfo = mrdfits(standardspath+'ispec_standards.fits',1,/silent)

    ra = 15.0*im_hms2dec(starinfo.ra) ; [degrees]
    dec = im_hms2dec(starinfo.dec)    ; [degrees]

    ntot = im_djs_angle_match(raref,decref,ra,dec,mindx=mindx,$
      dtheta=searchrad/3600.0,units='degrees',mdist=mdist)
;   niceprint, lindgen(n_elements(raref)), raref, decref, mindx, mdist

    if (ntot eq 0L) then begin
       objindx = -1L
       return, -1L
    endif
    if (ntot ge 1L) then begin
       objindx = where(mindx ne -1L)
       starindx = mindx[objindx]
       if (ntot eq 1L) then starindx = starindx[0L]
    endif

return, starindx
end

