;+
; NAME:
;       BALMER_DATA()
;
; PURPOSE:
;       Initialize parameters for the Balmer emission/absorption
;       lines.
;
; CALLING SEQUENCE:
;       bdata = balmer_data(nbdata=)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       bdata - data structure with all the relevant data on the first
;               four Balmer lines
;
; OPTIONAL OUTPUTS:
;       nbdata - number of Balmer lines (hard-wired to 4)
;
; PROCEDURES USED:
;
; COMMENTS:
;       Use the continuum windows for H-beta and H-delta defined by
;       Trager et al. (1998), modify the blue H-gamma window to avoid
;       the G-band, and define my own H-alpha windows.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Nov 26, U of A - written
;       jm05jul28uofa - new definition of the H-beta pseudo-continuum
;                       bands  
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

function balmer_data, nbdata=nbdata

    nbdata = 4L ; only 4 Balmer lines are supported
    
    bdata = {$
      line:  '', $
      label: '', $
      wave:   0.0, $
      llimit: 0.0, $
      lwidth: 0.0, $
      ulimit: 0.0, $
      uwidth: 0.0, $
      lline:  0.0, $
      uline:  0.0}

    bdata = replicate(bdata,nbdata)

    bdata.line   = ['H_delta', 'H_gamma', 'H_beta', 'H_alpha']
    bdata.label  = ['H\delta', 'H\gamma', 'H\beta', 'H\alpha'];+' Absorption'
    bdata.wave   = [ 4101.73 ,  4340.46 ,  4861.33,  6562.80 ]
    bdata.llimit = [ 4061.00 ,  4250.00 ,  4825.00,  6525.00 ]
    bdata.lwidth = [   38.00 ,    40.00 ,    20.00,    25.00 ]
    bdata.ulimit = [ 4145.00 ,  4394.00 ,  4900.00,  6630.00 ]
    bdata.uwidth = [   33.00 ,    53.00 ,    20.00,    60.00 ]
    bdata.lline  = [ 4083.00 ,  4320.00 ,  4848.00,  6533.00 ]
    bdata.uline  = [ 4122.00 ,  4363.00 ,  4877.00,  6593.00 ]

;   bdata.line   = ['H_delta', 'H_gamma', 'H_beta', 'H_alpha']
;   bdata.label  = ['H\delta', 'H\gamma', 'H\beta', 'H\alpha'];+' Absorption'
;   bdata.wave   = [ 4101.73 ,  4340.46 ,  4861.33,  6562.80 ]
;   bdata.llimit = [ 4061.00 ,  4250.00 ,  4838.00,  6525.00 ]
;   bdata.lwidth = [   38.00 ,    40.00 ,    20.00,    25.00 ]
;   bdata.ulimit = [ 4145.00 ,  4394.00 ,  4884.00,  6630.00 ]
;   bdata.uwidth = [   33.00 ,    53.00 ,    15.00,    60.00 ]
;   bdata.lline  = [ 4083.00 ,  4320.00 ,  4848.00,  6533.00 ]
;   bdata.uline  = [ 4122.00 ,  4363.00 ,  4877.00,  6593.00 ]

return, bdata
end
