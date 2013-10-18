;+
; NAME:
;       ISPEC_COUNT_CRPIX
;
; PURPOSE:
;       Count the number of pixels flagged as cosmic rays given a file
;       list. 
;
; CALLING SEQUENCE:
;       ispec_count_crpix, root=, file=
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       root - search for all file names that satisfy the condition
;              ROOT+'*.fits*' 
;       file - count the number of cosmic ray pixels in FITS file FILE 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       The number of bad pixels and the file list are printed to
;       STDOUT. 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       HEADFITS(), MRDFITS(), SXPAR()
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Dec 25, U of A
;
; Copyright (C) 2003, John Moustakas
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

pro ispec_count_crpix, root=root, file=file

    if n_elements(root) eq 0L then root = 'fw'
    
    if n_elements(file) eq 0L then $
      file = file_search(root+'*.fits*',count=count) else $
      count = n_elements(file)
    
    for i = 0L, count-1L do begin

       h = headfits(file[i])
       mask = mrdfits(file[i],2,/silent)
;      s = rd2dspec(file[i],/silent)
       crpix = long(total(mask / 8 mod 2))
       print, file[i], sxpar(h,'OBJECT'), crpix, format='(A25,2x,A20,2x,I5)'
       
    endfor
    
return
end    
