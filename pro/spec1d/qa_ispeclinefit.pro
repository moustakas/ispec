;+
; NAME:
;       QA_ISPECLINEFIT
;
; PURPOSE:
;       Conduct some very simple tests of the quality assurance (QA)
;       of the fitting results from ISPECLINEFIT.
;
; CALLING SEQUENCE:
;       qa_ispeclinefit
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
; PROCEDURES USED:
;
;
; COMMENTS:
;       Check for: infinite fluxes; infinite EW's; 
;       
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 10, U of A - written based on an
;                                           earlier, less generalized
;                                           code 
;-

pro qa_ispeclinefit, nuclear=nuclear

    if keyword_set(nuclear) then root = 'nuclear_atlas'
    
    l = read_integrated(root=root,/silent)
    galaxy = strcompress(l.galaxy,/remove)
    ngalaxy = n_elements(galaxy)
    
    linename = strcompress(l[0].linename,/remove)
    nline = n_elements(linename)

    tags = tag_names(l)
    ntags = n_elements(tags)

    badindx = -1L
    
    for j = 0L, nline-1L do begin

       splog, 'Checking '+linename[j]+'.'

       keeptags = ['GALAXY','FIT_ID',linename[j]+['','_CONTINUUM','_EW','_SIGMA','_CHI2']]
       line = struct_trimtags(l,select=keeptags)
       line.fit_id = line.fit_id*3
       
       nomeasure = where((line.(2))[1,*] eq -2.0,nno)
       if (nno ne 0L) then begin
          splog, 'Line not measured.'
          struct_print, struct_trimtags(line[nomeasure],except=['*CHI2*,*CONTINUUM*'])
       endif

       upper = where((line.(2))[1,*] eq -3.0,nupper)
       if (nupper ne 0L) then begin
          splog, 'Upper limit.'
          struct_print, line[upper]
       endif

       cc = get_kbrd(1)
       print
       
    endfor

return
end    
