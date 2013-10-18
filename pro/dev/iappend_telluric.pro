;+
; NAME:
;       IAPPEND_TELLURIC
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Dec 31, U of A
;-

pro iappend_telluric, filelist, sensname, datapath=datapath, $
  outpath=outpath, silent=silent, gzip=gzip, wfits=wfits
; given a sensitivity function and a file list append the telluric
; absorption spectrum to each file

; WARNING: if WFITS=1 then each file in FILELIST is over-written

    fcount = n_elements(filelist) ; number of images
    senscount = n_elements(sensname)
    if (fcount eq 0L) or (senscount ne 1L) then begin
       print, 'Syntax - iappend_telluric, filelist, sensname'
       return
    endif

    if not keyword_set(datapath) then datapath = cwd()
    if not keyword_set(outpath) then outpath = datapath

    outname = repstr(filelist,'.gz','') ; output file names

    sens = irdsensfunc(sensname,datapath=datapath,silent=silent)
    telluric_spec = sens.telluric_spec
    telluric_head = sens.telluric_head

    if (total(telluric_spec) eq 0.0) or (strcompress(telluric_head[0],/remove) eq '') then begin
       splog, 'Sensitivity function '+sensname+$
         ' does not have a valid telluric absorption spectrum.'
       return
    endif
    
    for i = 0L, fcount-1L do begin

       cube = rd2dspec(filelist[i],datapath=datapath,silent=silent)
       header = *cube.header

       splog, 'Successfully added telluric spectrum keyword TELLSPEC.'
       sxaddpar, header, 'TELLSPEC', 'T', ' appended telluric spectrum [T/F]', before='HISTORY'

       if keyword_set(wfits) then begin
          
          splog, 'Writing '+outpath+filelist[i]+'.'
          wrt2dspec, outname[i], cube.image, cube.sigmamap, cube.mask, header, $
            skyimage=cube.sky, telluric_spec=telluric_spec, telluric_head=telluric_head, $
            datapath=outpath, gzip=gzip
          
       endif
       
       icleanup, cube
    
    endfor

return
end    
