;+
; NAME:
;       IHEADER_KEYWORDS
;
; PURPOSE:
;       Add header keywords to FITS headers.
;
; CALLING SEQUENCE:
;       iheader_keywords, keywordfile, datapath=, /update, /silent
;
; INPUTS:
;       keywordfile - input text file.  the first line should list the
;                     keywords to add, the next line should give the
;                     keyword comments.  the remaining file should
;                     list each file to be updated followed by the
;                     header values to add separated by commas.  for
;                     example (see IMAKE_KEYWORDFILE):
;
;    # SCANLEN, POSANGLE, APERTURE
;    # scan length [arcsec], slit position angle [degrees], slit aperture [arcsec]
;    
;    a.0101.fits,  0, 90, 2.5
;    a.0112.fits, 60, 90, 2.5
;
; OPTIONAL INPUTS:
;       datapath - path for I/O
;
; KEYWORD PARAMETERS:
;       update - set this keyword to actually overwrite the FITS
;                headers, otherwise the intended updates are printed
;                to the screen without modifying any files
;
; OUTPUTS:
;       The FITS headers are updated with the requested keywords.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       CWD(), SPLOG, DJS_READLINES(), HEADFITS(), DJS_MODFITS,
;       SXADDPAR 
;
; COMMENTS:
;       When updating the headers the data-type of the FITS files is
;       preserved [see DJS_MODFITS()].
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 7, U of A
;       jm03apr16uofa, documented, ISPEC v1.0.0
;
; Copyright (C) 2002-2003, John Moustakas
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

pro iheader_keywords, keywordfile, datapath=datapath, update=update, silent=silent
    
    if n_elements(keywordfile) eq 0L then begin
       print, 'Syntax - iheader_keywords, keywordfile, datapath=, /update, /silent'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    
    pushd, datapath

    if file_test(keywordfile,/regular) eq 0L then begin
       splog, 'Keyword file '+datapath+keywordfile+' not found.'
       return
    endif

    if not keyword_set(silent) then splog, 'Reading '+datapath+keywordfile+'.'
    data = djs_readlines(keywordfile,nhead=3,head=head)
    data = data[where(strcompress(data,/remove) ne '')]
    nobj = n_elements(data)

    keywords = strcompress(strsplit(strmid(head[0],strpos(head[0],'#')+1,$
      strlen(head[0])),',',/extract),/remove)
    nkeys = n_elements(keywords)

    comments = strsplit(strmid(head[1],strpos(head[1],'#')+1,strlen(head[1])),',',/extract)
    ncomments = n_elements(comments)

    if (nkeys ne ncomments) then message, 'Zero comments are not supported yet.'

    datatype = strsplit(strmid(head[2],strpos(head[2],'#')+1,strlen(head[2])),',',/extract)
    ndatatype = n_elements(datatype)

    if not keyword_set(silent) then splog, 'Found keywords '+strjoin(keywords,' ')+'.'
    if not keyword_set(silent) then splog, 'Found comments '+strjoin(comments,' ')+'.'

; crop trailing comments    
    
    junk = where(strmatch(data,'*#*') eq 1B,njunk) 
    if (njunk ne 0L) then for i = 0L, njunk-1L do data[junk[i]] = $
      strtrim(strmid(data[junk[i]],0,strpos(data[junk[i]],'#')))
    
    for i = 0L, nobj-1L do begin

       datasep = strsplit(data[i],',',/extract)
       
       fitsfile = datasep[0]
       h = headfits(fitsfile)
       h = h[where(strcompress(h,/remove) ne '')] ; clean up blank lines

;      for j = 0L, nkeys-1L do sxaddpar, h, keywords[j], $
;        strcompress(string(datasep[j+1],format='('+datatype[j]+')'),/remove), comments[j]

       for j = 0L, nkeys-1L do sxaddpar, h, keywords[j], $
         strcompress(datasep[j+1],/remove), comments[j], before='HISTORY'
;      splog, 'NOTE!!!!!!!!!!!!!!!!'
;      sxdelpar, h, 'PA'
       
       if keyword_set(update) then begin
          if not keyword_set(silent) then splog, 'Updating '+fitsfile+': '+$
            strjoin(strtrim(datasep[1:nkeys],2),', ')+'.'
          djs_modfits, fitsfile, 0, h
       endif else if not keyword_set(silent) then splog, 'To be updated: '+$
         fitsfile+': '+strjoin(strtrim(datasep[1:nkeys],2),', ')+'.'

    endfor
    
    popd
    
return
end
