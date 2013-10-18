;+
; NAME:
;       IMAKE_KEYWORDFILE
;
; PURPOSE:
;       Generate a template header keywords file. 
;
; CALLING SEQUENCE:
;       imake_keywordfile, keywords, comments=, values=, $
;          datatype=, root=, datapath=, headfile=, imtype=
;
; INPUTS:
;       keywords - header keywords to add strarr(NKEYWORDS)
;
; OPTIONAL INPUTS:
;       comments - header comments corresponding to each keyword
;                  strarr(NKEYWORDS)
;       values   - values corresponding to each header keyword
;                  strarr(NKEYWORDS)
;       datatype - NOT SUPPORTED
;       root     - common prefix of all the FITS files to be searched 
;                  and written to HEADFILE
;       datapath - data path to the FITS files and to the directory 
;                  where HEADFILE should be written
;       headfile - output file name
;       imtype   - character string required of IMAGETYP
;                  (default='object') 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       A file DATAPATH+HEADFILE is written that looks like:
;
;          # PA, SCANLEN, APERTURE
;          # scan length [arcsec], slit position angle [degrees], slit aperture [arcsec]
;          n10069.fits,      90.0,      0.0,      3.0
;          n10070.fits,      90.0,      0.0,      3.0
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine generates a template header keywords file that
;       could then be modified by hand to change particular header
;       values.  Once the header file is satisfactory the FITS file
;       headers can be updated with IHEADER_KEYWORDS.
;       
;       All the FITS files in DATAPATH are searched and only files with
;       IMAGETYP OBJECT (or alternative supplied by keyword IMTYPE) are 
;       written to HEADFILE.
;       
;       If no values are given then they are assumed to be '0.0' for
;       all keywords.
;
; EXAMPLE:
;       We want to add the slit position angle (which we'll call 'PA')
;       to our FITS headers. The fits files have a common root of 'a.'.
;       Try 
;
;          IDL> imake_keywordfile, 'PA', comments='slit position angle', $
;          IDL> values='90.0', root='a.', datapath='mypath/', $
;          IDL> headfile='header_keywords.dat' 
;
;       The file 'header_keywords.dat' would then be given to
;       IHEADER_KEYWORDS for parsing and to update the FITS files with 
;       a header parameter of 'PA' and a value of 90 degrees.  If
;       different files have different position angles then
;       'header_keywords.dat' must be modified by hand.
;
; PROCEDURES USED:
;       CWD(), SXPAR(), HEADFITS(), SPLOG
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 8, U of A
;       A. Marble, 2002 December 11 - added IMTYPE keyword
;       jm03apr16uofa, checked out, ISPEC v1.0.0
;       jm05apr17uofa - replaced FINDFILE() with FILE_SEARCH() 
;
; Copyright (C) 2002-2003, 2005, John Moustakas
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

pro imake_keywordfile, keywords, comments=comments, values=values, $
  datatype=datatype, root=root, datapath=datapath, headfile=headfile, $
  imtype=imtype

    nkeys = n_elements(keywords)

    if (nkeys eq 0L) then begin
       print, 'Syntax - imake_keywordfile, keywords, comments=, values=, $'
       print, '   root=, datapath=, headfile=, imtype='
       return
    endif

    if n_elements(root) eq 0L then root = ''
    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(headfile) eq 0L then headfile = 'header_keywords.dat'
    if n_elements(imtype) eq 0L then imtype = 'object'

    ncomments = n_elements(comments)
    nvalues = n_elements(values)
;   ndatatype = n_elements(datatype)

    if ncomments ne 0L then if ncomments ne nkeys then begin
       splog, 'Incompatible dimensions in KEYWORDS and COMMENTS.'
       return
    endif
       
    if (nvalues eq 0L) then values = replicate('0.0',nkeys) else $
      if (nvalues ne nkeys) then begin
       splog, 'Incompatible dimensions in KEYWORDS and VALUES.'
       return
    endif

;   if ndatatype eq 0L then values = replicate('A15',nkeys) else if ndatatype ne nkeys then $
;     splog, 'Incompatible dimensions in KEYWORDS and DATATYPE.'

; generate the file list

    flist = file_search(root+'*.fits',count=fcount)
    if fcount eq 0L then begin
       splog, 'No FITS files found.'
       return
    endif
    mask = bytarr(fcount)

    object = strarr(fcount)
    for j = 0L, fcount-1L do begin
       h = headfits(flist[j])
       type = strlowcase(strcompress(sxpar(h,'IMAGETYP'),/remove))
       object[j] = sxpar(h,'OBJECT')
       if strmatch(type,'*'+imtype+'*') then mask[j] = 1B
    endfor
       
    good = where(mask,fcount)
    if (fcount ne 0L) then begin
       flist = flist[good] 
       object = object[good] 
    endif else begin
       splog, 'No files of image type '+strupcase(imtype)+'.'
       return
    endelse

    openw, lun1, datapath+headfile, /get_lun
    printf, lun1, '# '+strjoin(keywords,', ')
    if ncomments ne 0L then printf, lun1, '# '+strjoin(comments,', ') else printf, lun1, '# '
    printf, lun1, ' ' ; this needs to be blank for the datatype keyword
;   printf, lun1, '# '+strjoin(datatype,', ')
    for k = 0L, fcount-1L do printf, lun1, flist[k]+',   '+strjoin(string(values,format='(G0.0)'),',   ')+$
      '   # '+object[k]
    free_lun, lun1

return
end    
