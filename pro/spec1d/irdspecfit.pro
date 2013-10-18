;+
; NAME:
;       IRDSPECFIT()
;
; PURPOSE:
;       Read the one-dimensional spectral fitting results written out
;       by ISPECLINEFIT().
;
; CALLING SEQUENCE:
;       specfit = irdspecfit(object,objtagname=,root=,$
;          specfitpath=,mjdstr=,toc=,/silent)
;
; INPUTS:
;       object - read the spectral fitting results for this object
;       
;
; OPTIONAL INPUTS:
;       objtagname  - SPECFIT structure tag name in which to search for
;                     OBJECT (default 'GALAXY') 
;       root        - root suffix for SPECFIT (default
;                     MJDSTR+'_*'+ROOT+'*_specfit.fits.gz)
;       mjdstr      - see ROOT
;       specfitpath - full path name to SPECFIT
;
; KEYWORD PARAMETERS:
;       silent - suppress messages to STDOUT
;
; OUTPUTS:
;       toc - "table of contents" file for SPECFIT
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       MRDFITS(), DJS_FILEPATH(), SPLOG, MATCH_STRING(),
;       FILE_SEARCH() 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 15, U of A
;       jm05sep29uofa - added SPECDATA optional input
;-

function irdspecfit, object, specdata=specdata, objtagname=objtagname, $
  root=root, mjdstr=mjdstr, specfitpath=specfitpath, toc=toc, silent=silent

    nobject = n_elements(object)
    if nobject eq 0L then begin
       print, 'Syntax - specfit = irdspecfit()'
       return, -1L
    endif
    
    if (n_elements(objtagname) eq 0L) then objtagname = 'GALAXY'
    if (n_elements(root) eq 0L) then root = ''
    if (n_elements(specfitpath) eq 0L) then specfitpath = cwd()
    if (n_elements(mjdstr) eq 0L) then mjdstr = ''

; read the SPECFIT and corresponding SPECDATA FITS files 

    allfiles = file_search(djs_filepath(specfitpath+'*'+mjdstr+'_*'+root+'*_specfit.fits*'),count=fcount)
    specfitfile = allfiles[(reverse(sort(allfiles)))[0]]

    allfiles = file_search(djs_filepath(specfitpath+'*'+mjdstr+'_*'+root+'*_specdata.fits.gz'),count=fcount)
    specdatafile = allfiles[(reverse(sort(allfiles)))[0]]

    if (file_test(specdatafile,/regular) eq 0L) then begin
       splog, 'SPECDATAFILE '+specdatafile+' not found.'
       return, -1L
    endif
    
    if not keyword_set(silent) then splog, 'Reading '+specdatafile+'.'
    if (n_elements(specdata) eq 0L) then specdata = mrdfits(specdatafile,1,/silent)
    
    if arg_present(toc) then toc = mrdfits(specfitfile,1,/silent)
    
; match OBJECT in SPECDATA to find the appropriate FITS extension to
; read in SPECFIT

    tagindx = where(objtagname eq tag_names(specdata),ntagindx)
    if (ntagindx ne 1L) then begin
       message, 'Multiple or no tags found matching '+OBJTAGNAME+'.', /info
       return, -1L
    endif

    doit = match_string(object,specdata.(tagindx),index=indx,/exact,silent=silent)
    nindx = n_elements(indx)

    if (nobject eq 1L) then begin
       if (nindx gt nobject) then begin
          splog, 'WARNING! Multiple matches found for object '+strtrim(object[0],2)+'.'
          indx = indx[0]
          doit = doit[0]
       endif
       if (strcompress(doit,/remove) eq '') then return, -1L
    endif
    ext_no = indx+2

; loop on each object

    for k = 0L, nobject-1L do begin

       if strcompress(doit[k],/remove) ne '' then begin
          if not keyword_set(silent) then splog, 'Reading '+specfitfile+', extension '+$
            string(ext_no,format='(I0)')+'.'
          t0 = systime(1)

          specfit1 = mrdfits(specfitfile,ext_no[k],/silent)
;         print, 'Time to read extension '+string(ext_no[k],format='(I0)')+': '+strtrim(string(systime(1)-t0),2)+' seconds.'
          if k eq 0L then specfit = specfit1 else specfit = [ [specfit], [specfit1] ]
       endif

    endfor

    specfit = reform(specfit)

return, specfit
end
