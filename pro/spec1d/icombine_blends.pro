;+
; NAME:
;       ICOMBINE_BLENDS()
;
; PURPOSE:
;       Combine the flux from unresolved blends, currently only
;       supporting [O II] 3726,3729.
;
; CALLING SEQUENCE:
;       linefit = icombine_blends(linefit)
;
; INPUTS:
;       linefit - linefit data structure from ILINEFIT  
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       outlinefit - updated LINEFIT structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       REMOVE, ICREATE_LINEFIT(), DJS_MEAN(), STRUCT_APPEND() 
;
; COMMENTS:
;       All doublets or duplicate lines should be of the form
;       OII_3727_1 and OII_3727_2 or H_BETA_1 and H_BETA_2, etc. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Jan 13, U of A - excised from IFITSPEC
;       jm04may09uofa - previously, if one line was dropped by MPFIT
;          (error code = -1.0) then the combined doublet was assigned
;          an error code of -1.0; but now we just take the flux of the
;          line that was *not* dropped
;       jm05apr20uofa - documented and added error checking
;
; Copyright (C) 2004-2005, 2007, John Moustakas
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

function icombine_blends, linefit

    nline = n_elements(linefit)
    if (nline eq 0L) then begin
       doc_library, 'icombine_blends'
       return, linefit
    endif

    keep = lindgen(nline)
    
    allnames = strcompress(linefit.linename,/remove)

    one = where((strmatch(allnames,'*_1') eq 1B),none)
    two = where((strmatch(allnames,'*_2') eq 1B),ntwo)
;   one = where((strmatch(allnames,'*_1') eq 1B) or (strmatch(allnames,'*_narrow',/fold) eq 1B),none)
;   two = where((strmatch(allnames,'*_2') eq 1B) or (strmatch(allnames,'*_broad',/fold) eq 1B),ntwo)

    if (none ne ntwo) or (none eq 0L) then begin
;      splog, ''
       return, linefit
    endif
    nmatch = none

    if ((none+ntwo) lt nline) then begin
       remove, [one,two], keep
       nkeep = n_elements(keep)
    endif else begin
       nkeep = 0L
    endelse
    
; alphabetize    
    
    lineone = allnames[one] & srtone = sort(lineone)
    linetwo = allnames[two] & srttwo = sort(linetwo)

    for k = 0L, nmatch-1L do begin

       line = [lineone[k],linetwo[k]]
;      newline = strmid(line[0],0,strlen(line[0])-2L) ; <-- NOT GENERAL!!
       newline = repstr(repstr(repstr(repstr(line[0],'_1',''),'_2',''),'_narrow',''),'_broad','')

       indx = [one[k],two[k]]
;      splog, 'Combining blends ['+strjoin(line,', ')+']'
       
       s = linefit[indx]
       s1 = s[0]
       s2 = s[1]

; initialize the output structure
       
       result = icreate_linefit(1)

       errflag = 0.0
       if (s1.linearea_err eq -1.0) or (s2.linearea_err eq -1.0) then errflag = -1.0  ; one line "dropped"
       if (s1.linearea_err eq -2.0) or (s2.linearea_err eq -2.0) then errflag = -2.0  ; "not measured"
       if (s1.linearea_err eq -1.0) and (s2.linearea_err eq -1.0) then errflag = -3.0 ; both lines "dropped"

       case errflag of

          0.0: begin            ; both lines are well-measured

; linearly add the redshift and add the error in quadrature
             result.linez = djs_mean(s.linez)
             result.linez_err = sqrt(s1.linez_err^2 + s2.linez_err^2)
             
; add the sigma line widths and errors in quadrature; compute the mean
; instrumental resolution

             dropped = where(s.linesigma eq -1.0,ndropped)
             if (ndropped eq 0L) then begin
                result.linesigma = sqrt(total(s.linesigma^2))
                result.linesigma_err = sqrt(total(s.linesigma_err^2))
             endif else begin
                result.linesigma = 0.0
                result.linesigma_err = -1.0
             endelse

             result.linesigma_instr = djs_mean(s.linesigma_instr)
             result.linesigma_total = sqrt(result.linesigma^2.0 + result.linesigma_instr^2.0)

; linearly add the Gaussian flux and add the error in quadrature
             result.linearea = s1.linearea + s2.linearea
             result.linearea_err = sqrt(s1.linearea_err^2 + s2.linearea_err^2)
             
; linearly add the box flux and add the error in quadrature
             result.linebox = s1.linebox + s2.linebox
             result.linebox_err = sqrt(s1.linebox_err^2 + s2.linebox_err^2)

; add the number of pixels linearly (not strictly correct)
             result.linenpix = total(s.linenpix)

; compute the mean number of degrees of freedom (not strictly correct)
             result.linedof = djs_mean(s.linedof)

; assign the chi2 to be the mean chi2
             result.linechi2 = djs_mean(s.linechi2)
          end

          -1.0: begin ; one line was dropped from the fit by MPFIT

             wgood = where(s.linearea_err gt 0.0)

; only use the redshift from the well-fitted line             
             result.linez = s[wgood].linez
             result.linez_err = s[wgood].linez_err

; add the sigma line widths in quadrature; compute the mean
; instrumental resolution
             result.linesigma = sqrt(total(s.linesigma^2))
             result.linesigma_err = s[wgood].linesigma ; <-- NOTE!
             result.linesigma_instr = djs_mean(s.linesigma_instr)
             result.linesigma_total = sqrt(result.linesigma^2.0 + result.linesigma_instr^2.0)

; linearly add the Gaussian flux and add the error in quadrature
             result.linearea = s[wgood].linearea         ; <-- NOTE
             result.linearea_err = s[wgood].linearea_err ; <-- NOTE
             
; linearly add the box flux and add the error in quadrature
             result.linebox = s[wgood].linebox         ; <-- NOTE
             result.linebox_err = s[wgood].linebox_err ; <-- NOTE

; add the number of pixels linearly (not strictly correct)
             result.linenpix = total(s.linenpix)

; compute the mean number of degrees of freedom (not strictly correct)
             result.linedof = djs_mean(s.linedof)

; assign the chi2 to be the mean chi2
             result.linechi2 = djs_mean(s.linechi2)

          end
          
;;;          -1.0: begin ; one line was dropped from the fit by MPFIT
;;;
;;;; linearly add the redshift and add the error in quadrature
;;;             result.linez = djs_mean(s.linez)
;;;             result.linez_err = sqrt(s1.linez_err^2 + s2.linez_err^2)
;;;
;;;; add the sigma line widths in quadrature; compute the mean
;;;; instrumental resolution
;;;             result.linesigma = sqrt(total(s.linesigma^2))
;;;             result.linesigma_err = -1.0
;;;             result.linesigma_instr = djs_mean(s.linesigma_instr)
;;;             result.linesigma_total = sqrt(result.linesigma^2.0 + result.linesigma_instr^2.0)
;;;
;;;; linearly add the Gaussian flux and add the error in quadrature
;;;             result.linearea = 0.0
;;;             result.linearea_err = -1.0
;;;             
;;;; linearly add the box flux and add the error in quadrature
;;;             result.linebox = 0.0
;;;             result.linebox_err = -1.0
;;;
;;;; add the number of pixels linearly (not strictly correct)
;;;             result.linenpix = total(s.linenpix)
;;;
;;;; compute the mean number of degrees of freedom (not strictly correct)
;;;             result.linedof = djs_mean(s.linedof)
;;;
;;;; assign the chi2 to be the mean chi2
;;;             result.linechi2 = djs_mean(s.linechi2)
;;;
;;;          end

          -2.0: begin           ; set the combined line as "not measured"

; simply copy the "unmeasured" line structure into RESULT          
             
             nm = where(s.linearea_err eq -2.0)
             result = s[nm[0]]
             
          end

          -3.0: begin ; both lines were dropped from the fit by MPFIT

; simply copy one of the "dropped" line structures into RESULT          
             
             result = s[0]

          end
          
       endcase

; fill the output structure with mean quantities 

       result.linename = newline
       result.linewave = djs_mean(s.linewave)
;      result.line_blend = 'good3727'

       if (k eq 0L) then result1 = result else result1 = [ [result1], [result] ]
       
    endfor

; concatenate RESULT to the input structure and return
    
    if (nkeep ne 0L) then $
      outlinefit = struct_append(reform(result1),linefit[keep]) else $
      outlinefit = reform(result1)

return, outlinefit
end

