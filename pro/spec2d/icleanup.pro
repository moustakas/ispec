;+
; NAME:
;       ICLEANUP
;
; PURPOSE:
;       Free RD2DSPEC() pointers and array memory. 
;
; CALLING SEQUENCE:
;       icleanup, cube, /main
;
; INPUTS:
;       cube - data cube from RD2DSPEC()
;
; KEYWORD PARAMETERS:
;       main - call RETALL
;
; OUTPUTS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2001 August 23, U of A, written
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm03dec08uofa - added SKY structure entry
;       jm06jan26uofa - error checking for HEADER tag
;
; Copyright (C) 2001, 2003, John Moustakas
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

pro icleanup, cube, main=main

    ncube = size(cube,/n_elements)
    if ncube eq 0L then return
    
    for i = 0L, ncube-1L do begin

       if tag_exist(cube,'HEADER') then if (ptr_valid(cube[i].header))[0] then ptr_free, cube[i].header

       if tag_exist(cube[i],'image') then cube[i].image = 0.0
       if tag_exist(cube[i],'sky') then cube[i].sky = 0.0
       if tag_exist(cube[i],'vmap') then cube[i].vmap = 0.0
       if tag_exist(cube[i],'mask') then cube[i].mask = 0.0
       if tag_exist(cube[i],'spec') then cube[i].spec = 0.0
       if tag_exist(cube[i],'sigspec') then cube[i].sigspec = 0.0
       if tag_exist(cube[i],'wave') then cube[i].wave = 0.0

    endfor
    
    if keyword_set(main) then retall

return    
end

