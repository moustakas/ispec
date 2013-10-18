;+
; NAME:
;   PARSE_ILINEFIT()
;
; PURPOSE:
;   Parse the output from ILINEFIT. 
;
; INPUTS:
;   linefit - output from ILINEFIT()
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   parsed  - parsed linefit structure
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine is really an internal support routine for
;   ISPECLINEFIT() and should not be called as a stand-alone program.
;   Essentially all this routine does is initialize the output
;   structure with the line fluxes and adds a galaxy, specfile, and a
;   mean redshift (and error) tags.
; 
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Mar 14, U of A - excised from ISPECLINEFIT()  
;   jm04may10uofa - compute the median, not the mean redshift 
;   jm08apr02nyu - occasionally, LINEDOF is 1.0, and so I was
;     dividing by zero when computing the reduced chi2; not sure
;     why/when this would happen, but LINEDOF is set to be greater
;     than equal to 1.0  
;
; Copyright (C) 2004, 2008, John Moustakas
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

function parse_ilinefit, linefit

    if (n_elements(linefit) eq 0L) then begin
       doc_library, 'parse_ilinefit'
       return, -1L
    endif

    nline = n_elements(linefit.linename)

    linename = strcompress(strupcase(linefit.linename),/remove)
    linewave = linefit.linewave

    for i = 0L, nline-1L do begin

; desired output tags

       tags = linename[i]+['','_BOX','_EW','_CHI2','_SIGMA','_SIGMA_INSTR',$
         '_SIGMA_TOTAL','_CONTINUUM','_WAVE','_LINEZ','_LIMIT','_EW_LIMIT']

; initialize the structure fields for each emission line fitted

       out1 = create_struct($
         tags[0], [linefit[i].linearea,linefit[i].linearea_err], $
         tags[1], [linefit[i].linebox,linefit[i].linebox_err], $
         tags[2], [linefit[i].lineew_area,linefit[i].lineew_area_err], $
         tags[3], linefit[i].linechi2,$
         tags[4], [linefit[i].linesigma,linefit[i].linesigma_err], $
         tags[5], linefit[i].linesigma_instr, $
         tags[6], linefit[i].linesigma_total, $
         tags[7], [linefit[i].linecontlevel,linefit[i].linecontlevel_err], $
         tags[8], linewave[i], $
         tags[9], [linefit[i].linez,linefit[i].linez_err], $
         tags[10],linefit[i].linelimit, $
         tags[11],linefit[i].lineew_limit)

; compute the reduced chi2
       
       if (linefit[i].linedof gt 0.0) and (linefit[i].linechi2 gt 0.0) then begin
          match = where(linename[i]+'_CHI2' eq tag_names(out1))
          out1.(match) = linefit[i].linechi2/((linefit[i].linedof-1.0)>1.0)
       endif

; grow the output structure       
       
       if (i eq 0L) then out = out1 else out = struct_addtags(out,out1)
       
    endfor

; add the mean redshift to the output structure

    good = where(linefit.linearea_err gt 0.0,ngood)
    if ngood ne 0L then begin
       z_line = djs_median(linefit[good].linez) ; median redshift
       z_line_err = djsig(linefit[good].linez)  ; standard error
    endif else begin
       z_line = 0.0
       z_line_err = -1.0
    endelse

    parsed = struct_addtags({linename: linename, $
      z_line: float(z_line), z_line_err: float(z_line_err)},out)

return, parsed
end
