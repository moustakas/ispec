;+
; NAME:
;       IREAD_LAMPLIST()
;
; PURPOSE:
;       Read the reference comparison line lists for various
;       elements. 
;
; CALLING SEQUENCE:
;       lamp = iread_lamplist(lampname,lamppath=,intensity_cut=,$
;          /keepblend,/silent)
;
; INPUTS:
;       lampname - comparison lamp name (e.g., 'HeAr') or a
;                  user-supplied lamp list (e.g., 'lamplist.dat')
;                  which is distinguished by the ".dat" extension
;
; OPTIONAL INPUTS:
;       lamppath      - data path to LAMPNAME
;       intensity_cut - only include lines brighter than this
;                       threshold (default 10.0)
;
; KEYWORD PARAMETERS:
;       keepblend - the default is to not use lines marked as blended;
;                   use this keyword to override this behavior
;       silent    - suppress messages to STDOUT
;
; OUTPUTS:
;       lamp - output data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       CWD(), READCOL, STRUCT_APPEND()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Jun 21, U of A - written, partially excised
;                                           from IARCFIT
;
; Copyright (C) 2004, John Moustakas
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

function iread_lamplist, lampname, lamppath=lamppath, intensity_cut=intensity_cut, $
  keepblend=keepblend, silent=silent

    if (n_elements(intensity_cut) eq 0L) then intensity_cut = 10.0
    if (n_elements(lamppath) eq 0L) then $
      lamppath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    
    if strmatch(lampname,'*.*') then begin

       crclist = lampname 

    endif else begin

       crc_ar = 'CRC04_Ar.dat'
       crc_cd = 'CRC04_Cd.dat'
       crc_he = 'CRC04_He.dat'
       crc_hg = 'CRC04_Hg.dat'
       crc_ne = 'CRC04_Ne.dat'

       indx_he = strmatch(lampname,'*He*',/fold)
       indx_ne = strmatch(lampname,'*Ne*',/fold)
       indx_ar = strmatch(lampname,'*Ar*',/fold)
       indx_cd = strmatch(lampname,'*Cd*',/fold)
       indx_hg = strmatch(lampname,'*Hg*',/fold)

       crclist = ''
       if indx_he then crclist = [crclist,crc_he]
       if indx_ne then crclist = [crclist,crc_ne]
       if indx_ar then crclist = [crclist,crc_ar]
       if indx_cd then crclist = [crclist,crc_cd]
       if indx_hg then crclist = [crclist,crc_hg]

       ncrclist = n_elements(crclist)
       if (ncrclist eq 1L) then begin
          splog, 'No matching line list for '+lampname+'.'
          return, -1L
       endif

       crclist = crclist[1L:ncrclist-1L]

    endelse

    ncrclist = n_elements(crclist)

    for ilist = 0L, ncrclist-1L do begin

       if (file_test(lamppath+crclist[ilist],/regular) eq 0L) then begin
          splog, 'Arc lamp file '+lamppath+crclist[ilist]+' not found.'
          return, -1L
       endif
       
       if not keyword_set(silent) then splog, 'Reading arc lamp file '+crclist[ilist]+'.'
       readcol, lamppath+crclist[ilist], wave, intensity, quality, element, $
         format='D,D,A,A', /silent, comment='#', delimiter='  '
       nline = n_elements(wave)

       lamp1 = replicate({lambda: 0.0D, intensity: 0.0D, $
         quality: '', element: '', good: 0B},nline)
       lamp1.element = element
       lamp1.lambda = wave
       lamp1.intensity = intensity
       lamp1.quality = quality

       if (ilist eq 0L) then lamp = lamp1 else lamp = struct_append(lamp,lamp1)
       
    endfor

    quality = strtrim(lamp.quality,2)
    intensity = lamp.intensity
    
    if keyword_set(keepblend) then $
      lamp.good = (quality ne 'BAD') else $
      lamp.good = (quality ne 'BAD') and (quality ne 'BLEND')

    lamp.good = lamp.good and (intensity gt intensity_cut)
    keep = where(lamp.good,nkeep)
    
    if (nkeep eq 0L) then begin
       splog, 'No good arc lines in the lamp list.'
       return, -1L
    endif else lamp = lamp[keep]
    
    srt = sort(lamp.lambda)
    lamp = lamp[srt]
    
return, lamp
end
    
