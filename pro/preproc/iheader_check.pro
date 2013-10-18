;+
; NAME:
;       IHEADER_CHECK
;
; PURPOSE:
;       Verify FITS headers keywords before processing with ISPEC.
;
; CALLING SEQUENCE:
;       iheader_check, root=, datapath=
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       root     - search all FITS files that satisfy the string
;                  pattern ROOT+'*.fits' (default '')
;       datapath - I/O path
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       See COMMENTS.  
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       HEADFITS(), SXPAR(), CWD(), IM_HDR2STRUCT(),
;       STRUCT_TRIMTAGS(), SPLOG 
;
; COMMENTS:
;       Warning messages are printed to STDOUT.  Key keywords are,
;       OBJECT, IMAGETYP, DATE-OBS, EXPTIME, OBSERVAT, RA, DEC, EPOCH, 
;       UT.  DATE-OBS *must* be in FITS standard format of
;       YYYY-MM-DD (*not* YY/MM/DD).  Also, dark frames (identified by
;       having IMAGETYP = dark) must have a DARKTIME keyword.  
;
; EXAMPLE:
;       Check the headers for the night 1 data:
;
;   IDL> iheader_check, root='n1', datapath='./'
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 8, U of A
;       jm03apr16uofa - documented, ISPEC v1.0.0
;       jm05apr14uofa - replaced FINDFILE() with FILE_SEARCH();
;                       better error checking, more informative
;                       messages 
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

pro iheader_check, root=root, datapath=datapath

    if n_elements(root) eq 0L then root = ''
    if n_elements(datapath) eq 0L then datapath = cwd()

    flist = file_search(datapath+root+'*.fits*',count=fcount)
    if (fcount eq 0L) then begin
       print, 'No valid FITS files found.'
       return
    endif

    rootflist = file_basename(flist)
    
    keywords = [$
      'OBJECT',  $
      'IMAGETYP',$
      'DATE-OBS',$
      'EXPTIME', $
      'OBSERVAT',$
      'RA',      $
      'DEC',     $
      'EPOCH',   $
      'UT'       $
      ]
    nwords = n_elements(keywords)
    for iword = 0L, nwords-1L do keywords[iword] = idl_validname(keywords[iword],/convert_all)
    
    len = strn(max(strlen(keywords)))
    nwords = n_elements(keywords)
    
    for i = 0L, fcount-1L do begin

       head = headfits(flist[i])
       struct = im_hdr2struct(head)
       substruct = struct_trimtags(struct,select=keywords,format=replicate('A0',nwords))
       tags = tag_names(substruct)
       
;      struct_print, substruct
;      for i = 0L, n_elements(tags)-1L do print, tags[i], substruct.(i)

       missing = cmset_op(keywords,'AND',/not2,tag_names(substruct))
       if size(missing,/type) ne 3L then $
         splog, 'WARNING: keywords missing from '+rootflist[i]+': '+strjoin(strtrim(missing,2),', ')
;      notmissing = cmset_op(keywords,'AND',tag_names(substruct))

;      for j = 0L, nwords-1L do begin
;
;         value = sxpar(head,keywords[j],/silent,count=hcount)
;         keystr = string(keywords[j],format='(A'+len+')')
;         if hcount[0] eq 0L then message, 'Keyword '+keystr+' not found in '+rootflist[i]+'.', /info
;         if hcount[0] gt 1L then message, 'Multiple keyword entries for '+keystr+' found in '+rootflist[i]+'.', /info
;
;      endfor 

; check dark frames separately

;      imtype = sxpar(head,'IMAGETYP',count=good)
;      if good[0] eq 1L then begin
;         if strlowcase(strcompress(imtype,/remove)) eq 'dark' then begin
;            junk = sxpar(head,'DARKTIME',count=good)
;            if good[0] eq 0L then message, 'Keyword DARKTIME not found in dark frame '+rootflist[i]+'.', /info
;         endif
;      endif 
       
    endfor

return
end    
