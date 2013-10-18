;+
; NAME:
;       ICHOOSE_WMAP()
;
; PURPOSE:
;       Select the wavelength solution nearest in (ra,dec) and time to
;       the given observation. 
;
; CALLING SEQUENCE:
;       wmapindx = ichoose_wmap(header,wheader,wmapname=,/silent)
;
; INPUTS:
;       header  - object FITS header
;       wheader - scalar or array of wavelength map header(s)
;                 generated by IRESTORE_WMAP()
;
; OPTIONAL INPUTS:
;       wmapname - wavelength map names corresponding to the trailing
;                  dimension in WHEADER 
;
; KEYWORD PARAMETERS:
;       silent - suppress messages to STDOUT
;
; OUTPUTS:
;       wmapindx - long integer index element of the appropriate
;                  wavelength map
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SXPAR(), IM_HMS2DEC(), DJS_DIFF_ANGLE()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Dec 08, U of A - excised from ICALIBRATE 
;       jm03dec16uofa - only choose wavelength maps whose arc lamps
;                       were observed the same day as the object
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

function ichoose_wmap, header, wheader, wmapname=wmapname, silent=silent

    ndim = size(wheader,/n_dimensions)
    if (ndim[0] eq 1L) then nwmap = 1L else begin
       dims = size(wheader,/dimensions)
       nwmap = dims[1]
    endelse

    if n_elements(wmapname) ne 0L then begin
       if (n_elements(wmapname) ne nwmap) then begin
          splog, 'WHEADER and WMAPNAME have incompatible dimensions.'
          return, -1L
       endif
    endif
    
    jdref = sxpar(header,'JD',count=jdcount)
    ratemp = sxpar(header,'RA',count=racount)
    detemp = sxpar(header,'DEC',count=decount)

    if (jdcount eq 0L) or (racount eq 0L) or (decount eq 0L) then begin
       splog, 'No JD, RA, or DEC header information in HEADER.'
       return, -1L
    endif
    
    raref = 15.0*im_hms2dec(ratemp)
    deref = im_hms2dec(detemp)

    diff = fltarr(nwmap)
    jddiff = fltarr(nwmap)
    
    for i = 0L, nwmap-1L do begin

       jd = sxpar(wheader[*,i],'JD',count=jdcount)
       ra = sxpar(wheader[*,i],'RA',count=racount)
       de = sxpar(wheader[*,i],'DEC',count=deccount)

       if (jdcount eq 0L) or (racount eq 0L) or (decount eq 0L) then begin
          if n_elements(wmapname) ne 0L then $
            splog, 'No JD, RA, or DEC header information for wavelength map '+wmapname[i]+'.' else $
            splog, 'No JD, RA, or DEC header information for wavelength map number '+string(i,format='(I0)')+'.'
          return, -1L
       endif

       jddiff[i] = jd-jdref
       diff[i] = djs_diff_angle(raref,deref,15.0*im_hms2dec(ra),im_hms2dec(de))

    endfor

    sameday = where((jddiff gt -0.5) and (jddiff lt +0.5),nsameday)
    if nsameday eq 0L then begin
       splog, 'WARNING: No wavelength maps were observed on the same day as this object!'
       sameday = lindgen(nwmap)
    endif

    diff = diff[sameday]
    mindiff = min(diff,wmapindx)
    tindx = sameday[wmapindx]
    
    if not keyword_set(silent) then begin
       
       if n_elements(wmapname) ne 0L then splog, 'Selecting wavelength map '+wmapname[tindx]+'.'
       splog, '   OBJECT, ARC (RA) : '+sxpar(header,'RA')+' '+sxpar(wheader[*,tindx],'RA')
       splog, '   OBJECT, ARC (DEC): '+sxpar(header,'DEC')+' '+sxpar(wheader[*,tindx],'DEC')

    endif
       
return, wmapindx
end