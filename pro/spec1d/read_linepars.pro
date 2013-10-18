;+
; NAME:
;       READ_LINEPARS()
;
; PURPOSE:
;       Read a line list parameter file.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       linefile - name of the line list parameter file
;
; KEYWORD PARAMETERS:
;       strong - read ELINELIST_STRONG.DAT containing just the strong
;                nebular emission lines
;
; OUTPUTS:
;       linepars - array of structures with these fields [NLINE]
;          line   - line name
;          wave   - central wavelength
;          zindex - redshift code
;          windex - line width code
;          findex - flux ratio code
;          fvalue - flux ratio
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       See IFITSPEC() for more information.
;
; PROCEDURES USED:
;       READCOL, SPLOG
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 31, U of A, written
;       jm03nov04uofa - documented
;       jm04jan05uofa - improved error checking
;       jm04sep05uofa - added STRONG keyword
;       jm08aug11nyu - removed LINEPATH optional input; LINEFILE
;         should include the full path name
;
; Copyright (C) 2002-2004, John Moustakas
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

function read_linepars, linefile=linefile, strong=strong

    linepath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    if (n_elements(linefile) eq 0L) then begin
       if keyword_set(strong) then $
         linefile = linepath+'elinelist_strong.dat' else $
           linefile = linepath+'elinelist.dat'
    endif

    if (file_test(linefile,/regular) eq 0L) then begin
       splog, 'Line list file '+linefile+' not found.'
       return, -1L
    endif

    readcol, linefile, line, center, zindex, windex, findex, fvalue, $
      format = 'A,D,A,A,A,D', /silent, comment='#'

    linepars = {line: ' ', wave: 0.0D, zindex: '', $
      windex: '', findex: '', fvalue: 0.0D}
    nline = n_elements(line)
    linepars = replicate(linepars,nline)

    linepars.line = line
    linepars.wave = center
    linepars.zindex = zindex
    linepars.windex = windex
    linepars.findex = findex
    linepars.fvalue = fvalue

return, linepars
end    
