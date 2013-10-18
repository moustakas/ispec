;+
; NAME:
;       BALMERINDX()
;
; PURPOSE:
;       Return the indices of the unblended Balmer lines given a list
;       of wavelengths.  The indices of the complementary (assumed
;       forbidden) transitions are also returned.
;
; CALLING SEQUENCE:
;       balmerindx = balmerindx(linewave,blend,forbidden=,/allbalmer)
;
; INPUTS:
;       linewave - wavelengths in Angstroms [NLINE]
;
; OPTIONAL INPUTS:
;       blend - lines with the same BLEND string value are blended
;               [NLINE] 
;
; KEYWORD PARAMETERS:
;       allbalmer - return all the Balmer lines regardless of whether
;                   or not they are blended; in this case BLEND does
;                   not need to be passed
;
; OUTPUTS:
;       balmerindx - indices of the Balmer lines in LINEWAVE
;
; OPTIONAL OUTPUTS:
;       forbidden - indices of the non-Balmer lines in LINEWAVE 
;
; PROCEDURES USED:
;       REMOVE
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 09, U of A - excised from IFITSPEC 
;       jm05apr20uofa - documented and added error checking
;
; Copyright (C) 2004-2005, John Moustakas
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

function balmerindx, linewave, blend, forbidden=forbidden, allbalmer=allbalmer

    nline = n_elements(linewave)
    if (nline eq 0L) then begin
       print, 'Syntax - balmerindx = balmerindx(linewave,blend,forbidden=,/allbalmer)'
       return, -1L
    endif

    nblend = n_elements(blend)
    if (nblend eq 0L) and (not keyword_set(allbalmer)) then begin
       print, 'If ALLBALMER=0 then BLEND must be defined.'
       return, -1L
    endif

    forbidden = lindgen(nline)   ; default is to assume no Balmer lines 
    
; reject non-blended Balmer lines.  find the Balmer lines by matching
; the wavelength LINEWAVE to the nearest tenth of an Angstrom

    balmerwaves = [3889.049,3970.072,4101.734,4340.464,4861.325,6562.80] ; NOT GENERAL!
    nbalmer = n_elements(balmerwaves)

; BALMERINDX indexes all Balmer lines in LINEWAVE.  if no Balmer lines
; were fitted then return; if there are multiple Balmer lines at the
; same wavelength (e.g., fitting QSO, then return the first one)
    
    balmerindx = -1L
    for i = 0L, nbalmer-1L do balmerindx = [balmerindx,where(long(10*balmerwaves[i]) eq long(10*linewave))]
    balmerindx = balmerindx[1L:n_elements(balmerindx)-1L]
    
    pos = where(balmerindx ne -1L,nbalmerindx,comp=neg)
    if (nbalmerindx ne 0L) then balmerindx = balmerindx[pos] else return, -1L

; if ALLBALMER=1 then the user wants all the Balmer lines regardless
; of whether they are blended or not; also construct and return the
; FORBIDDEN line index    

    if keyword_set(allbalmer) then begin
       if (nbalmerindx eq nline) then forbidden = -1L else $
         remove, balmerindx, forbidden
       return, balmerindx
    endif

; now remove from BALMERINDX any blended Balmer lines (e.g., H-alpha) 
    
    rem = lonarr(nbalmerindx)-1L
    for i = 0L, nbalmerindx-1L do $
      if total(strmatch(blend,blend[balmerindx[i]])) gt float(1.0) then $
      rem[i] = balmerindx[i]

; if the only Balmer line in BALMERINDX is a blend then return

    match = where(rem ne -1L,nmatch)
    if (nmatch eq 1L) and (nbalmerindx eq 1L) then return, -1L
    if (nmatch ne 0L) then remove, match, balmerindx

    nbalmerindx = n_elements(balmerindx)
    
; finally, remove all the non-Balmer (forbidden) lines from the fully
; index linelist to generate the FORBIDDEN variable.  if NBALMERINDX
; equals NLINE that means that no forbidden lines were successfully
; fitted 

    if (nbalmerindx lt nline) then $
      remove, forbidden[balmerindx], forbidden else $
      forbidden = -1L

return, balmerindx
end

