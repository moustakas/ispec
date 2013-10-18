;+
; NAME:
;      ISPEC_VERSION()
;
; PURPOSE:
;      Return the current iSPEC version number.
;
; CALLING SEQUENCE:
;      ver = ispec_version()
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;      ver - string version number
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;      Based entirely on Dave Schlegel's IDLUTILS_VERSION() routine:
;
;      If this version is not tagged by CVS, then we return
;      'NOCVS:TOPLEVEL' where TOPLEVEL is the last directory in the
;      environment variable ${ISPEC_DIR}.  For example, if you are 
;      using a version of the code in the directory
;      '/u/schlegel/idlutils/v0_0', then this returns 'NOCVS:v0_0'. 
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;      J. Moustakas, 2003 April 7, U of A
;      jm03dec8uofa - upgraded to v1_1
;      jm05apr14uofa - upgraded to v2_0
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

function ispec_version

; the following expression in dollar signs is expanded by CVS and
; replaced by the tag name for this version

   name = '$Name: v2_1 $'

   words = str_sep(strcompress(name), ' ')

   if (words[0] EQ '$Name:' AND N_elements(words) EQ 3) then begin
      ver = words[1]
   endif else begin
      dirname = getenv('ISPEC_DIR')
      if (dirname NE '') then begin
         words = str_sep(dirname,'/')
         nword = N_elements(words)
         ver = 'NOCVS:' + words[nword-1]
      endif else begin
         ver = 'NOCVS:Unknown'
      endelse
   endelse

return, ver
end
