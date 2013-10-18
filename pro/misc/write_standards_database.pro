;+
; NAME:
;       WRITE_STANDARDS_DATABASE
;
; PURPOSE:
;       Generate the master standard-star database.  
;
; CALLING SEQUENCE:
;       write_standards_database, /silent
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       See the code.
;
; PROCEDURES USED:
;       CWD()
;
; COMMENTS:
;       See the code comments for the priority list of standards.
;       This routine should be run whenever the standard-star
;       directory is updated.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 April 14, U of A
;       jm05jun01uofa - replace the CTIONEWCAL standards with those in
;                       BESSELL1999, which have the telluric absorption
;                       removed; added optional input FINALSTARINFO 
;
; Copyright (C) 2005, John Moustakas
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

pro write_standards_database, finalstarinfo, write=write

    rootdir = getenv('ISPEC_DIR')+'/'
    rootpath = filepath('',root_dir=rootdir,subdirectory='standards')
    datapath = rootpath+'data/'
    outpath = rootpath+'spectra/'
    cwd = cwd()

    if keyword_set(write) then begin
       
       if file_test(outpath,/directory) then begin
          splog, 'Press ENTER to remove all files in '+outpath+' or Q to quit.'
          cc = get_kbrd(1)
          if (strupcase(cc) eq 'Q') then return else spawn, ['/bin/rm -f '+outpath+'/*'], /sh
       endif

    endif

; recurvively identify all the STANDARDS_INFO.SEX files in DATAPATH 

    starfiles = file_search(datapath,'standards_info.sex',count=nstarfiles)

    for istar = 0L, nstarfiles-1L do begin

;      starpath = file_basename(file_dirname(starfiles[istar],/mark_directory))
;      starpath = strmid(starfiles[istar],0,strpos(starfiles[istar],'/',/reverse_search))
       
       starinfo1 = rsex(starfiles[istar])
       nstar = n_elements(starinfo1)
;      starinfo1 = struct_addtags(starinfo1,replicate({starpath: starpath},nstar))
       
       if (istar eq 0L) then starinfo = starinfo1 else $
         starinfo = struct_append(starinfo,starinfo1)

    endfor

; sort by RA; generate a uniq list of stars and correlate that list
; with the following priority list: calspec, bessell1999, spec50cal,
; oke1990, oke1983

    srt = sort(im_hms2dec(starinfo.ra))
    starinfo = starinfo[srt]
;   struct_print, starinfo
    
    ustarinfo = starinfo[uniq(strtrim(strupcase(starinfo.star),2))]
    nstar = n_elements(ustarinfo)

    for istar = 0L, nstar-1L do begin

       match = where(strtrim(ustarinfo[istar].star,2) eq strtrim(starinfo.star,2),nmatch)
       thisone = where(strmatch(starinfo[match].file,'*calspec*',/fold) eq 1B,nthisone)
       if (nthisone eq 0L) then thisone = where(strmatch(starinfo[match].file,'*bessell1999*',/fold) eq 1B,nthisone)
       if (nthisone eq 0L) then thisone = where(strmatch(starinfo[match].file,'*spec50cal*',/fold) eq 1B,nthisone)
       if (nthisone eq 0L) then thisone = where(strmatch(starinfo[match].file,'*oke1990*',/fold) eq 1B,nthisone)
       if (nthisone eq 0L) then thisone = where(strmatch(starinfo[match].file,'*oke1983*',/fold) eq 1B,nthisone)
       if (nthisone ne 1L) then message, 'We have a problem!'

       if (istar eq 0L) then finalstarinfo = starinfo[match[thisone]] else $
         finalstarinfo = struct_append(finalstarinfo,starinfo[match[thisone]])

    endfor
    struct_print, finalstarinfo

    nfinalstar = n_elements(finalstarinfo)

; write out    
    
    if keyword_set(write) then begin

       mwrfits, finalstarinfo, outpath+'ispec_standards.fits', /create
       for istar = 0L, nfinalstar-1L do $
         spawn, ['ln -fs ../data/'+finalstarinfo[istar].starpath+'/'+strtrim(finalstarinfo[istar].file,2)+' '+$
           outpath+strtrim(finalstarinfo[istar].file,2)], /sh
;        spawn, ['ln -fs '+datapath+finalstarinfo[istar].starpath+'/'+strtrim(finalstarinfo[istar].file,2)+' '+$
;          outpath+strtrim(finalstarinfo[istar].file,2)], /sh
       
    endif
    
return
end    
