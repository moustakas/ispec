;+
; NAME:
;       INICE_LINENAME()
;
; PURPOSE:
;       Convert a spectral line name into a "nice" format.
;
; CALLING SEQUENCE:
;       nicename = inice_linename(name)
;
; INPUTS:
;       name - string line name [e.g., see usage in
;              IFIGURE_SPECLINEFIT()]. 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       nicename - nicely formatted line name
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Apr 22, U of A - written
;       jm05jun17uofa - documented
;       jm05jul22uofa - bug fix when converting [Ne III] 3869 
;
; Copyright (C) 2004, 2005, John Moustakas
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

function inice_linename, name

    nline = n_elements(name)
    nicename = strarr(nline)

    for iline = 0L, nline-1L do begin

       nm = strtrim(name[iline],2)
       
       if strmatch(nm,'H_*') then begin ; Balmer line

          split = strsplit(nm,'_',/extract)
          nsplit = n_elements(split)
          nicename[iline] = split[0]+'\'+strlowcase(split[1])
          if (nsplit gt 2L) then nicename[iline] = nicename[iline]+' '+strlowcase(strjoin(split[2L:nsplit-1L],' '))

       endif else begin ; forbidden line

          nicename[iline] = '['+repstr(strmid(nm,0,strpos(nm,'_'))+'] \lambda'+strmid(nm,strpos(nm,'_')+1,strlen(nm)),'_',' ')
          
;         ion = strpos(nm,'I')
;         element = strmid(nm,0,ion) ; strupcase(strmid(nm,0,ion))
;         therest = strmid(nm,ion,strlen(nm)-ion)
;         nicename[iline] = '['+element+' '+strtrim(repstr(therest,'_','] \lambda'),2)

       endelse

    endfor

    if (nline eq 1L) then nicename = nicename[0]
    
return, textoidl(nicename)
end
