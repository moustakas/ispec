;+
; NAME:
;   ICREATE_LINEFIT()
;
; PURPOSE:
;   Initialize the iSPEC1d emission-line output structure and set
;   default return values.  
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   nline - optionally replicate LINEFIT NLINE times
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   linefit - see ILINEFIT() for details
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Apr 20, U of A - documented
;   jm05jul25uofa - output changed to floating-point precision 
;
; Copyright (C) 2005, 2008, John Moustakas
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

function icreate_linefit, nline

   linefit = create_struct($
    'linename'          ,   ' ', $
    'linewave'          ,   0.0, $
    'linez'             ,   0.0, $
    'linez_err'         ,  -1.0, $
    'linesigma'         ,   0.0, $ ; intrinsic line width [km/s]
    'linesigma_err'     ,  -2.0, $
    'linesigma_instr'   ,   0.0, $ ; instrumental line width [km/s]
    'linesigma_total'   ,   0.0, $ ; instrumental plus intrinsic line width [km/s]
    'linearea'          ,   0.0, $
    'linearea_err'      ,  -2.0, $
    'linelimit'         ,  -2.0, $
    'linebox'           ,   0.0, $
    'linebox_err'       ,  -2.0, $
    'linecontlevel'     ,   0.0, $
    'linecontlevel_err' ,  -2.0, $
    'lineew_area'       ,   0.0, $
    'lineew_area_err'   ,  -2.0, $
    'lineew_limit'      ,  -2.0, $
    'lineew_box'        ,   0.0, $
    'lineew_box_err'    ,  -2.0, $
    'linenpix'          ,   0,   $
    'linedof'           ,   0.0, $
    'linechi2'          ,  -2.0)
   if (n_elements(nline) ne 0L) then linefit = replicate(linefit,nline)

return, linefit
end
