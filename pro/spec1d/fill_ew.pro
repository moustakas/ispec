;+
; NAME:
;   FILL_EW()
;
; PURPOSE:
;   Compute the equivalent width from a measure of the continuum
;   and total line flux. 
;
; INPUTS:
;   linefit - output structure from ILINEFIT() or IABSLINEFIT() 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   linefit - (modified)
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;   IM_COMPUTE_ERROR()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2003 June 12, U of A - excised from IFITSPEC
;   jm06sep27nyu - compute the uncertainties in double precision 
;
; Copyright (C) 2003, 2006-2007, John Moustakas
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

function fill_ew, linefit

    nline = n_elements(linefit)
    if nline eq 0L then begin
       print, 'Syntax - linefit = fill_ew(linefit)'
       return, -1L
    endif

    for k = 0L, nline-1L do begin

       case linefit[k].linearea_err of

          -1.0: begin ; dropped from the fit by ILINEFIT()

             linefit[k].lineew_area_err = -1.0 
             linefit[k].lineew_box_err  = -1.0

             if tag_exist(linefit,'LINELIMIT') then linefit[k].lineew_limit = $
               linefit[k].linelimit / linefit[k].linecontlevel
             
          end

          -2.0: ; not measured

          -3.0: begin ; upper limit

             linefit[k].lineew_area = linefit[k].linearea / linefit[k].linecontlevel
             linefit[k].lineew_area_err = -3.0

             linefit[k].lineew_box = 0.0
             linefit[k].lineew_box_err = -3.0

          end

          else: begin
       
             linefit[k].lineew_area = linefit[k].linearea / linefit[k].linecontlevel
             linefit[k].lineew_area_err = im_compute_error(double(linefit[k].linearea),$
               double(linefit[k].linearea_err),double(linefit[k].linecontlevel),$
               double(linefit[k].linecontlevel_err),/quotient)

             linefit[k].lineew_box = linefit[k].linebox / linefit[k].linecontlevel
             linefit[k].lineew_box_err = im_compute_error(double(linefit[k].linebox),$
               double(linefit[k].linebox_err),double(linefit[k].linecontlevel),$
               double(linefit[k].linecontlevel_err),/quotient)

             if tag_exist(linefit,'LINELIMIT') then linefit[k].lineew_limit = $
               linefit[k].linelimit / linefit[k].linecontlevel
             
          end
             
       endcase

    endfor
    
return, linefit
end

