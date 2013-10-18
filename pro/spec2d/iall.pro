;+
; NAME:
;	IALL
;
; PURPOSE:
;	Reduce an entire directory structure of data using IBATCH.
;
; CALLING SEQUENCE:
;       iall, paramfile, datapath=, procfile=, tracefile=, skyapfile=, $
;          crsplitfile=, crfile=, stdfile=, calibfile=, fluxcalibfile=, $
;          rednight=, senstitle=, sensname=, tellfits=, /noarcfit, $
;          sensinfo=, _extra=extra
;
; INPUTS:
;	paramfile   - vector or scalar input parameter file (see
;                     IBATCH)
;
; OPTIONAL INPUTS:
;	datapath    - path to the parameter file and the data
;	procfile    - text file with a list of FITS files (one per
;                     line) to process with ICCDPROC
;	tracefile   - text file with a list of FITS files (one per
;                     line) with which to trace the spatial distortion
;       skyapfile   - sky aperture file list generated by ISKYSELECT
;	crsplitfile - text file with a list of N cosmic ray split FITS 
;                     files (N files per line) to combine
;	crfile      - text file with a list of FITS files (one per
;                     line) to cosmic-ray reject with LA_COSMIC()
;	stdfile     - images with which to generate a sensitivity
;                     function
;	fluxcalibfile - text file with a list of FITS files (one per
;                     line) to flux calibrate with IFLUXCALIBRATE
;	calibfile   - text file with a list of FITS files (one per
;                     line) to calibrate with ICALIBRATE
;	rednight    - reduce a specific night of data (zero-indexed)
;       senstitle   - sensitivity function title array if
;                     n_elements(STDFILE) ne 0L, [NNIGHTS] string
;                     array 
;       sensname    - sensitivity function name array if
;                     n_elements(STDFILE) ne 0L, [NNIGHTS] string
;                     array 
;       tellfits    - telluric absorption FITS spectrum if
;                     n_elements(CALIBFILE) ne 0L, or
;                     n_elements(FLUXCALIBFILE) ne 0L [NNIGHTS] string
;                     array  
;	extra       - extra keywords for IBATCH
;	
; KEYWORD PARAMETERS:
;       noarcfit    - do not call IARCFIT in the initial reductions 
;
; OUTPUTS:
;       See IBATCH.
;
; OPTIONAL OUTPUTS:
;       sensinfo - structure array output from ISENSFUNC()
;
; COMMENTS:
;	This routine will do nothing unless at least one optional
;	input has been passed.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	READCOL, IBATCH, CWD()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 2 - written
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr23uofa - changed SENS keyword to MAKESENS; added
;                       SENSINFO optional output and SENSTITLE input 
;       jm03dec8uofa  - added support for 2D sky subtraction; removed
;                       the (redundant) keywords CCDPROC, SKYSUB,
;                       CRCLEAN, MAKESENS, and CALIBRATE
;       jm03dec31uofa - added TELLFILE optional input
;       jm04jun23uofa - suppress SKYFILE input; added NOARCFIT keyword
;       jm05jun20uofa - combining CR splits *and* cosmic-ray rejection
;                       now occurs before sky subtraction, because
;                       CR's can really mess up the sky fit
;       jm05jun23uofa - added FLUXCALIBFILE optional input; removed
;                       TELLFILE optional input; added SENSNAME and
;                       TELLFITS optional inputs
;       jm08jan18nyu  - if certain files are not found, try to keep
;                       going 
;
; Copyright (C) 2002-2005, 2008, John Moustakas
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

pro iall, paramfile, datapath=datapath, procfile=procfile, tracefile=tracefile, $
  skyapfile=skyapfile, crsplitfile=crsplitfile, crfile=crfile, stdfile=stdfile, $
  fluxcalibfile=fluxcalibfile, calibfile=calibfile, rednight=rednight, $
  senstitle=senstitle, sensname=sensname, tellfits=tellfits, noarcfit=noarcfit, $
  sensinfo=sensinfo, _extra=extra

    if not keyword_set(datapath) then datapath = cwd()

    if n_elements(paramfile) eq 0L then begin
       print, 'Syntax - iall, paramfile, datapath=, procfile=, tracefile=, $'
       print, '   skyapfile=, crsplitfile=, crfile=, stdfile=, calibfile=, $'
       print, '   fluxcalibfile=, rednight=, senstitle=, sensname=, tellfits=, $'
       print, '   /noarcfit, sensinfo='
       return
    endif
    
    nnights = n_elements(paramfile)

    if (n_elements(rednight) ne 0L) then begin
       rednight = (0L > long(rednight)) < nnights 
    endif else begin
       rednight = lindgen(nnights)
    endelse
    nreduce = n_elements(rednight)

    if (n_elements(noarcfit) eq 0L) then arcfit = 1L else arcfit = 0L
    
; ----------------------------------------------------------------------
; generate the flat field, ccd-process the data, optionally model the
; spatial distortion, and generate the wavelength map
    
    if (n_elements(procfile) ne 0L) then begin
       
       if (n_elements(procfile) ne nnights) then begin
          splog, 'PARAMFILE and PROCFILE have incompatible dimensions.'
          return
       endif

       if (n_elements(tracefile) eq 0L) then distortion = 0L else begin
          if (n_elements(tracefile) ne nnights) then begin
             splog, 'PARAMFILE and TRACEFILE have incompatible dimensions.'
             return
          endif
          distortion = 1L
       endelse 

       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Conducting initial reductions on night ",i0,"/",i0,".")', j+1, nreduce

          if file_test(datapath+procfile[i],/regular) eq 0L then begin
             splog, 'PROCFILE '+datapath+procfile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 
          endif else begin

             readcol, datapath+procfile[i], proclist, format='A', comment='#', /silent

             if distortion then begin
                if file_test(datapath+tracefile[i],/regular) then begin
                   readcol, datapath+tracefile[i], tracelist, format='A', comment='#', /silent
                endif else begin
                   splog, 'TRACEFILE '+datapath+tracefile[i]+' not found...continuing.'
;                  return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 
                endelse
             endif
             
             ibatch, paramfile[i], info, datapath=datapath, proclist=proclist, $
               tracelist=tracelist, /makeflat, /ccdproc, distortion=distortion, $
               arcfit=arcfit, _extra=extra
   
          endelse 

       endfor 

    endif 

; ----------------------------------------------------------------------
; model the spatial distortion; only do this if PROCFILE was not given 

    if (n_elements(tracefile) ne 0L) and (n_elements(procfile) eq 0L) then begin

       if (n_elements(tracefile) ne nnights) then begin
          splog, 'PARAMFILE and TRACEFILE have incompatible dimensions.'
          return
       endif

       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Modeling the spatial distortion for night ",i0,"/",i0,".")', j+1, nreduce

          if file_test(datapath+tracefile[i],/regular) eq 0L then begin
             splog, 'TRACEFILE '+datapath+tracefile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 
          endif else begin
             readcol, datapath+tracefile[i], tracelist, format='A', comment='#', /silent
             ibatch, paramfile[i], info, datapath=datapath, tracelist=tracelist, $
               /distortion, _extra=extra
          endelse 

       endfor 

    endif 

; ----------------------------------------------------------------------
; combine cosmic-ray split images, if they exist

    if (n_elements(crsplitfile) ne 0L) then begin
       
       if (n_elements(crsplitfile) ne nnights) then begin
          splog, 'PARAMFILE and CRSPLITFILE have incompatible dimensions.'
          return
       endif

       for j = 0L, nreduce-1L do begin
          
          i = rednight[j]
          splog, format='("Combining cosmic-ray split images for night ",i0,"/",i0,".")', j+1L, nreduce
          
          if file_test(datapath+crsplitfile[i],/regular) eq 0L then begin
             
             splog, 'No cosmic ray split file found...skipping.'
             
          endif else begin
             
             ibatch, paramfile[i], info, datapath=datapath, crsplitfile=crsplitfile[i], $
               _extra=extra, /crsplits

          endelse
               
       endfor 

    endif 

; ----------------------------------------------------------------------
; flag and remove cosmic rays

    if (n_elements(crfile) ne 0L) then begin

       if (n_elements(crfile) ne nnights) then begin
          splog, 'PARAMFILE and CRFILE have incompatible dimensions.'
          return
       endif

       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Flagging cosmic rays for night ",i0,"/",i0,".")', j+1L, nreduce

          if file_test(datapath+crfile[i],/regular) eq 0L then begin

             splog, 'CRFILE '+datapath+crfile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 

          endif else begin

             readcol, datapath+crfile[i], crlist, format='A', comment='#', /silent
             ibatch, paramfile[i], info, datapath=datapath, crlist=crlist, $
               /crclean, _extra=extra

          endelse 

       endfor

    endif
    
; ----------------------------------------------------------------------
; sky-subtract          
    
    if (n_elements(skyapfile) ne 0L) then begin

       if (n_elements(skyapfile) ne 0L) then if (n_elements(skyapfile) ne nnights) then begin
          splog, 'PARAMFILE and SKYAPFILE have incompatible dimensions.'
          return
       endif

       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Subtracting the sky for night ",i0,"/",i0,".")', j+1L, nreduce

          ibatch, paramfile[i], info, datapath=datapath, skyapfile=skyapfile[i], $
            /skysub, _extra=extra

       endfor 

    endif 
    
; ----------------------------------------------------------------------
; generate the sensitivity function

    if (n_elements(stdfile) ne 0L) then begin

       if (n_elements(stdfile) ne nnights) then begin
          splog, 'PARAMFILE and STDFILE have incompatible dimensions.'
          return
       endif

       if (n_elements(sensname) ne 0L) then begin
          if (n_elements(sensname) eq 1L) then begin
             sensname1 = replicate(sensname,nnights) 
          endif else begin
             if (n_elements(sensname) ne nnights) then begin
                splog, 'PARAMFILE and SENSNAME have incompatible dimensions.'
                return
             endif else sensname1 = sensname
          endelse
       endif
             
       if (n_elements(senstitle) ne nnights) then senstitle = replicate('',nnights)

       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          senstitle1 = senstitle[i]

          splog, format='("Generating the sensitivity function for night ",i0,"/",i0,".")', j+1L, nreduce

          if file_test(datapath+stdfile[i],/regular) eq 0L then begin

             splog, 'STDFILE '+datapath+stdfile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 

          endif else begin

             readcol, datapath+stdfile[i], stdlist, format='A', comment='#', /silent

             if (n_elements(sensname) ne 0L) then begin
                ibatch, paramfile[i], info, datapath=datapath, stdlist=stdlist, $
                  sensinfo=sensinfo1, senstitle=senstitle1, sensname=sensname1[i], $
                  /makesens, _extra=extra
             endif else begin
                ibatch, paramfile[i], info, datapath=datapath, stdlist=stdlist, $
                  sensinfo=sensinfo1, senstitle=senstitle1, /makesens, _extra=extra
             endelse
             
             if (j eq 0L) then sensinfo = sensinfo1 else sensinfo = [sensinfo,sensinfo1]

          endelse 

       endfor 
       
    endif
    
; ----------------------------------------------------------------------
; flux calibrate

    if (n_elements(fluxcalibfile) ne 0L) then begin

       if (n_elements(fluxcalibfile) ne nnights) then begin
          splog, 'PARAMFILE and FLUXCALIBFILE have incompatible dimensions.'
          return
       endif

       if (n_elements(tellfits) ne 0L) then begin
          if (n_elements(tellfits) ne nnights) then begin
             splog, 'PARAMFILE and TELLFITS have incompatible dimensions.'
             return
          endif
       endif else tellfits = replicate('',nnights)
             
       if (n_elements(sensname) ne 0L) then begin
          if (n_elements(sensname) eq 1L) then begin
             sensname1 = replicate(sensname,nnights) 
          endif else begin
             if (n_elements(sensname) ne nnights) then begin
                splog, 'PARAMFILE and SENSNAME have incompatible dimensions.'
                return
             endif else sensname1 = sensname
          endelse
       endif
             
       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Flux calibrating night ",i0,"/",i0,".")', j+1, nreduce

          if file_test(datapath+fluxcalibfile[i],/regular) eq 0L then begin

             splog, 'FLUXCALIBFILE '+datapath+fluxcalibfile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 

          endif else begin
             
             readcol, datapath+fluxcalibfile[i], caliblist, format='A', comment='#', /silent

             if (n_elements(sensname) ne 0L) then begin
                ibatch, paramfile[i], info, datapath=datapath, caliblist=caliblist, $
                  tellfits=tellfits[i], sensname=sensname1[i], _extra=extra, /fluxcalibrate
             endif else begin
                ibatch, paramfile[i], info, datapath=datapath, caliblist=caliblist, $
                  tellfits=tellfits[i], _extra=extra, /fluxcalibrate
             endelse

          endelse

       endfor 

    endif 

; ----------------------------------------------------------------------
; flux- and wavelength-calibrate, and distortion correct

    if (n_elements(calibfile) ne 0L) then begin

       if (n_elements(calibfile) ne nnights) then begin
          splog, 'PARAMFILE and CALIBFILE have incompatible dimensions.'
          return
       endif

       if (n_elements(tellfits) ne 0L) then begin
          if (n_elements(tellfits) ne nnights) then begin
             splog, 'PARAMFILE and TELLFITS have incompatible dimensions.'
             return
          endif
       endif else tellfits = replicate('',nnights)

       if (n_elements(sensname) ne 0L) then begin
          if (n_elements(sensname) eq 1L) then begin
             sensname1 = replicate(sensname,nnights) 
          endif else begin
             if (n_elements(sensname) ne nnights) then begin
                splog, 'PARAMFILE and SENSNAME have incompatible dimensions.'
                return
             endif else sensname1 = sensname
          endelse
       endif
             
       for j = 0L, nreduce-1L do begin

          i = rednight[j]
          splog, format='("Calibrating night ",i0,"/",i0,".")', j+1, nreduce

          if file_test(datapath+calibfile[i],/regular) eq 0L then begin

             splog, 'CALIBFILE '+datapath+calibfile[i]+' not found...continuing.'
;            return ; jm08jan16nyu - do not return, in case we're on, e.g., night 3/4 

          endif else begin
             
             readcol, datapath+calibfile[i], caliblist, format='A', comment='#', /silent

             if (n_elements(sensname) ne 0L) then begin
                ibatch, paramfile[i], info, datapath=datapath, caliblist=caliblist, $
                  tellfits=tellfits[i], sensname=sensname1[i], _extra=extra, /calibrate
             endif else begin
                ibatch, paramfile[i], info, datapath=datapath, caliblist=caliblist, $
                  tellfits=tellfits[i], _extra=extra, /calibrate
             endelse

          endelse

       endfor 

    endif 

; ----------------------------------------------------------------------    
; all done!  go extract some spectra!

return
end
