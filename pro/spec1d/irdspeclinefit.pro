;+
; NAME:
;       IRDSPECLINEFIT()
;
; PURPOSE:
;       Read the data structure output from PARSE_ISPECLINEFIT().  
;
; CALLING SEQUENCE:
;       speclinefit = irdspeclinefit(root=,$
;          speclinefile=,datapath=,/silent)
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       root     - search for the most recent file with the file name
;                  ROOT+'*_speclinefit.fits.gz' (default '') 
;       datapath - path name to SPECLINEFILE [default CWD()]
;
; KEYWORD PARAMETERS:
;       silent - suppress messages to STDOUT
;
; OUTPUTS:
;       speclinefit - output data structure structure written by
;                     PARSE_ISPECLINEFIT()
;
; OPTIONAL OUTPUTS:
;       speclinefile - file name read by this routine
;
; PROCEDURES USED:
;       CWD(), SPLOG, MRDFITS(), DJS_FILEPATH()
;
; COMMENTS:
;
; EXAMPLES:
;       IDL> speclinefit = irdspeclinefit(root='*integrated*',/silent)
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 14, U of A - written, based on an
;                                           earlier code
;-

function irdspeclinefit, root=root, speclinefile=speclinefile, $
  datapath=datapath, silent=silent
    
    if n_elements(root) eq 0L then root = ''
    if n_elements(datapath) eq 0L then datapath = cwd()

    pushd, datapath

    if n_elements(speclinefile) eq 0L then begin
       allfiles = findfile(djs_filepath(root+'*_speclinefit.fits.gz'),count=fcount)
       speclinefile = allfiles[(reverse(sort(allfiles)))[0]]
    endif

    if file_test(speclinefile,/regular) eq 0L then begin
       splog, 'File '+datapath+speclinefile+' not found.'
       return, -1L
    endif

    if not keyword_set(silent) then splog, 'Reading '+datapath+speclinefile+'.'
    speclinefit = mrdfits(datapath+speclinefile,1,/silent)

    popd
    
return, speclinefit
end
