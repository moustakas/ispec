;+
; NAME:
;       IPARSE_QALOGFILE()
;
; PURPOSE:
;       Parse the QA photometry file LOGFILE written by ISENSFUNC(). 
;
; CALLING SEQUENCE:
;       photo = iparse_qalogfile(logfile,stats=)
;
; INPUTS:
;       logfile - string name
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       photo - output data structure
;
; OPTIONAL OUTPUTS:
;       stats - sensitivity function statistics data structure 
;
; PROCEDURES USED:
;       RSEX(), DJS_READILINES()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Apr 24, U of A - written
;       jm05jun23uofa - improved documentation and updated to
;                       read the ISENSFUNC log output using RSEX() 
;       jm05jun28uofa - vectorized
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

function iparse_qalogfile, logfile, stats=stats

    nlog = n_elements(logfile)
    if (nlog eq 0L) then begin
       print, 'Syntax - photo = iparse_qalogfile(logfile,stats=)'
       return, -1L
    endif

    if (nlog gt 1L) then begin
       for ilog = 0L, nlog-1L do begin
          photo1 = iparse_qalogfile(logfile[ilog],stats=stats1)
          if (ilog eq 0L) then begin
             photo = photo1
             stats = stats1
          endif else begin
             photo = struct_append(photo,photo1)
             stats = [ [stats], [stats1] ]
          endelse
       endfor
       return, photo
    endif
    
    if (file_test(logfile,/regular) eq 0L) then begin
       splog, 'Log file '+logfile+' not found.'
       retall
    endif

; read the log file

    photo = rsex(logfile)
    nstds = n_elements(photo)

; parse the sensitivity function statistics

    ostats = djs_readilines(logfile,indx=4)
    oresids = float(strsplit(strmid(ostats,strpos(ostats,':')+1,$
      strlen(ostats)),',',/extract))
    
    gstats = djs_readilines(logfile,indx=5)
    gresids = float(strsplit(strmid(gstats,strpos(gstats,':')+1,$
      strlen(gstats)),',',/extract))

    stats = {star: '', $
      omean: 0.0, omed: 0.0, omin: 0.0, omax: 0.0, osigma: 0.0, $ ; observed stats
      gmean: 0.0, gmed: 0.0, gmin: 0.0, gmax: 0.0, gsigma: 0.0}   ; grey-shift stats

    stats.omean = oresids[0]
    stats.omed = oresids[1]
    stats.omin = oresids[2]
    stats.omax = oresids[3]
    stats.osigma = oresids[4]
    
    stats.gmean = gresids[0]
    stats.gmed = gresids[1]
    stats.gmin = gresids[2]
    stats.gmax = gresids[3]
    stats.gsigma = gresids[4]
    
return, photo
end

