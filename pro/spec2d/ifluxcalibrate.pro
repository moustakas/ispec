;+
; NAME:
;       IFLUXCALIBRATE
;
; PURPOSE:
;       Flux-calibrate two-dimensional spectra using the
;       two-dimensional wavelength map.
;
; CALLING SEQUENCE:
;       ifluxcalibrate, caliblist, sensname=, extfile=, extpath=, $
;          tellfits=, datapath=, outpath=, /exclude_senserr, /gzip, $
;          /wfits
;
; INPUTS:
;	caliblist - file list of spectra to flux-calibrate 
;
; OPTIONAL INPUTS:
;	sensname  - FITS name of the sensitivity curve (ISENSFUNC)  
;       extfile   - name of the extinction file
;       extpath   - path to EXTFILE
;       tellfits  - name of the FITS file containing the telluric
;                   absorption spectrum (usually the output from
;                   ICONSTRUCT_TELLURIC) to append to the output
;	datapath  - path to the data (default CWD)
;	outpath   - output path (default CWD)
;
; KEYWORD PARAMETERS:
;       exclude_senserr - do not propagate the error in the
;                         sensitivity  function to the flux-calibrated
;                         error maps
;       gzip            - write GZIPPED FITS files
;	wfits           - write the flux-calibrated 2D spectra to
;                         OUTPATH 
;
; OUTPUTS:
;       Flux-calibrated spectra are written to OUTPATH as
;       'f'+CALIBLIST if WFITS=1.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       CWD(), REPSTR(), READCOL, SPLOG, IRDSENSFUNC(), SXPAR(),
;       RD2DSPEC(), TAG_EXIST(), FINDEX(), SXADDPAR, SXADDHIST,
;       WRT2DSPEC, ICLEANUP
;
; COMMENTS:
;       Extinction-correcting is optional (but not recommended!).
;       Only one TELLFITS file can be passed, and it is assumed that
;       the FITS file exists in DATAPATH.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jun 22, U of A, written, roughly based on
;          the the code in ICALIBRATE 
;       jm05jun23uofa - added TELLFITS optional input
;       jm05jun29uofa - added EXCLUDE_SENSERR keyword; removed IRAF
;                       EX-FLAG header keyword 
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

pro ifluxcalibrate, caliblist, sensname=sensname, extfile=extfile, $
  extpath=extpath, tellfits=tellfits, datapath=datapath, outpath=outpath, $
  exclude_senserr=exclude_senserr, gzip=gzip, wfits=wfits

    nspec = n_elements(caliblist)
    nsens = n_elements(sensname)
    
    if (nspec eq 0L) or (nsens eq 0L) then begin
       print, 'Syntax - ifluxcalibrate, caliblist, sensname=, extfile=, extpath=, $'
       print, '   tellfits=, datapath=, outpath=, /exclude_senserr, /gzip, /wfits'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = datapath

    outname = 'f'+repstr(caliblist,'.gz','') ; root output file names

; read in the extinction curve

    extcor = n_elements(extfile)
    if (extcor[0] eq 1L) then if (extfile[0] eq '') then extcor = 0L

    if (extcor ne 0L) then begin

       if (n_elements(extpath) eq 0L) then $
         extpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')

       if file_test(extpath+extfile,/regular) then begin
          
          splog, 'Reading the extinction file '+extpath+extfile+'.'
          readcol, extpath+extfile, extwave, extvals, format='F,F', /silent

       endif else begin

          splog, 'Extinction file '+extpath+extfile+' not found.'
          return

       endelse 

    endif else splog, 'WARNING: Spectra will not be extinction corrected.'
    
; read the telluric spectrum

    ntell = n_elements(tellfits)
    if (ntell eq 1L) then begin

       if file_test(datapath+tellfits,/regular) then begin
          
          splog, 'Reading the telluric spectrum '+datapath+tellfits+'.'
          tellspec = mrdfits(datapath+tellfits,0,tellhead,/silent)

       endif else begin

          splog, 'Telluric spectrum '+datapath+tellfits+' not found.'
          ntell = 0L

       endelse 

    endif
    
; restore the sensitivity function

    if file_test(datapath+sensname,/regular) then begin
       
       splog, 'Reading '+sensname+'.'

       senscube = irdsensfunc(sensname,datapath=datapath,/silent)
       senshead = senscube.header

; use the grey-shifted sensitivity function to flux-calibrate; the
; error in the grey-shifted sensitivity function (GREYSENSERR) is the
; uncertainty in the sensitivity function, which is propagated unless
; EXCLUDE_SENSERR=1; the scatter in the *observed* sensitivity
; function (SENSERR) is the absolute error in the spectrophotometry;
; typical values of both these uncertainties are written to the header 
       
       if (strmatch(sensname,'*grey*') eq 1B) then $
         tsens = senscube.greysens else $
         tsens = senscube.sens

       tsenserr = senscube.greysenserr ; senstivity function error [mag]
       senserr_abs = senscube.senserr  ; absolute error [mag]

; convert from magnitudes to counts
       
       sens = 10D^(0.4*tsens)
       if keyword_set(exclude_senserr) then senserr = sens*0.0 else $
         senserr = 0.4D * alog(10.0) * sens * tsenserr

       senswave = senscube.wave
       zptshift = sxpar(senshead,'ZPTSHIFT',count=zptcount) ; from ISENSFUNC()
       if (zptcount eq 0L) then zptshift = 0.0
       
    endif else begin
       
       splog, 'Unable to find sensitivity function '+sensname+'.'
       return

    endelse 

; loop on each 2D spectrum
    
    for i = 0L, nspec-1L do begin

       splog, 'Flux-calibrating '+caliblist[i]+'.'
       cube = rd2dspec(caliblist[i],datapath=datapath,wset=wset,/silent)

       image    = cube.image
       sigmap   = cube.sigmamap
       sky      = cube.sky
       header   = *cube.header
       ncols    = cube.naxis1
       nrows    = cube.naxis2

       if (tag_exist(cube,'WAVEMAP') eq 0L) then begin
          splog, 'Spectrum must be sky subtracted.'
          icleanup, cube
          continue
       endif

       wavemap = cube.wavemap
       dwavemap = wavemap - shift(wavemap,1) ; pixel size
       dwavemap[0,*] = dwavemap[1,*]
       
       exptime = sxpar(header,'EXPTIME',count=nexptime)
       airmass= sxpar(header,'AIRMASS',count=nairmass)

       if (nexptime eq 0L) or (nairmass eq 0L) then begin
          splog, 'Keywords EXPTIME and/or AIRMASS missing from header.'
          continue
       endif

; flux calibrate and extinction-correct each row according to its
; wavelength vector
       
       if extcor then $
         extinct = interpolate(extvals,findex(extwave,wavemap)) else $
         extinct = wavemap*0.0 + 1.0

       ext_image  = image  * 10.0^(0.4*airmass*extinct)
       ext_sky    = sky    * 10.0^(0.4*airmass*extinct)
       ext_sigmap = sigmap * 10.0^(0.4*airmass*extinct)

       sensfunc = interpolate(sens,findex(senswave,wavemap))
       sensfuncerr = sqrt(interpolate(senserr^2,findex(senswave,wavemap)))

       flux_image  = ext_image / sensfunc / exptime / dwavemap
       flux_sky    = ext_sky / sensfunc / exptime / dwavemap

       flux_sigmap = sqrt((ext_sigmap / sensfunc)^2 + $
         (ext_image * sensfuncerr / sensfunc^2)^2) / exptime / dwavemap

; update the header

       if extcor then begin
          sxaddpar, header, 'EXTCOR', 'T', ' extinction-corrected [T/F]', before='HISTORY'
;         sxaddpar, header, 'EX-FLAG', 0L, before='HISTORY'
          sxaddpar, header, 'IEXTNAME', extfile, ' extinction file name', before='HISTORY'
          sxaddhist, "'Extinction corrected to zero airmass "+im_today()+"'", header
       endif
       
       sxaddpar, header, 'FLUXCOR', 'T', ' flux-calibrated [T/F]', before='HISTORY'
       sxaddpar, header, 'ZUNITS', 'erg/s/cm2/A', ' flux units', before='HISTORY'

       sxaddpar, header, 'IFLUX', caliblist[i], ' input to IFLUXCALIBRATE', $
         before='HISTORY'
       sxaddpar, header, 'SENSNAME', sensname, ' sensitivity function name', before='HISTORY'
       sxaddpar, header, 'ZPTSHIFT', float(zptshift), format='(F12.4)', $
         ' sensitivity function zero point shift [mag]', before='HISTORY'

       sxaddpar, header, 'AMINERR', min(senserr_abs), format='(F12.4)', $
         ' minimum absolute spectrophotometric error [mag]', before='HISTORY'
       sxaddpar, header, 'AMAXERR', max(senserr_abs), format='(F12.4)', $
         ' maximum absolute spectrophotometric error [mag]', before='HISTORY'
       sxaddpar, header, 'AMEDERR', median(senserr_abs), format='(F12.4)', $
         ' median absolute spectrophotometric error [mag]', before='HISTORY'

       sxaddhist, "'Flux calibrated "+im_today()+"'", header
       if keyword_set(exclude_senserr) then sxaddhist, $
         "'Sensitivity function uncertainty not included.'", header
       
       if keyword_set(wfits) then begin

          splog, 'Writing '+datapath+outname[i]+'.'
          wrt2dspec, outname[i], float(flux_image), float(flux_sigmap), cube.mask, header, $
            skyimage=float(flux_sky), wset=wset, telluric_spec=tellspec, $
            telluric_head=tellhead, datapath=datapath, gzip=gzip
          
       endif

       icleanup, cube ; clean up memory
       
    endfor

return
end
