;+
; NAME:
;	IDOSPEC
;
; PURPOSE:
;	Extract one-dimensional spectra (wrapper for IDOEXTRACT). 
;
; CALLING SEQUENCE:
;       idospec, flist, specinfo, datapath=, outname=, outpath=, $
;          specnames=, traceinfo=, _extra=extra, /silent, /deredden, $
;          /telluric, /noplot, /gzip, /wfits
;
; INPUTS:
;	flist    - file list of wavelength-calibrated two dimensional
;                  spectra 
;    
; OPTIONAL INPUTS:
;	datapath  - path to the data (default PWD)
;	outname   - root of the output FITS file name (see
;                   MODIFICATION HISTORY)
;	outpath   - write the extracted spectra to this path
;	extra     - keywords to be passed to IEXTRACT or K_LAMBDA()
;	
; KEYWORD PARAMETERS:
;       silent   - do not print messages to the screen
;	deredden - de-redden the spectrum for reddening using the
;                  Schlegel, Finkbeiner, & Davis (1998) dust map and a
;                  reddening curve (the default is the Cardelli,
;                  Clayton, & Mathis (1989) reddening curve: see
;                  K_LAMBDA() for more choices)
;       telluric - correct the 1D spectrum for telluric absorption 
;       gzip     - compress the output spectra if WFITS=1
;	wfits    - write multiple extension FITS files to outpath
;
; OUTPUTS:
;	specinfo   - structure containing the following fields:
;	  spec     - one dimensional spectrum
;	  wave     - corresponding wavelength vector
;	  sigspec  - one dimensional variance spectrum
;	  header   - spectrum header
;	  skyinfo  - sky subtraction information structure 
;	  apinfo   - aperture extraction parameters
;
; OPTIONAL OUTPUTS:
;	specnames  - output file name(s) of the extracted spectra 
;	traceinfo  - structure containing the functional form of the
;                    trace (see IEXTRACT)
;
; COMMENTS:
;	The output spectrum name has an 'ms' appended to signify
;	"multiple spectra."  See the IRAF documentation in onedspec
;	for more information. 
;
; TODO: 
;       [1] Add an interactive extraction tool.
;       [2] Possibly take the average of the two reddening estimates. 
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD2DSPEC(), IDOEXTRACT, MWRFITS(), CWD(), WRT1DSPEC,
;	DUST_GETVAL(), GLACTC, IM_HMS2DEC(), SXADDHIST, SPLOG, 
;	ICLEANUP, FMTAPERTURE(), IM_TODAY(), K_LAMBDA(),
;	IDIVIDE_TELLURIC(), STRUCT_TRIMTAGS(), STRUCT_ADDTAGS(),
;	TAG_EXIST() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August/October 17-19, U of A
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr16uofa - added SKYSUB2D capability and SILENT keyword
;       jm03jul16uofa - developed SKYSUB2D by making the default sky
;                       aperture positions and sizes more intelligent 
;       jm03dec4uofa  - added GZIP keyword
;       jm03dec8uofa  - removed SKYSUB2D capability and updated with
;                       the new treatment of sky subtraction
;       jm03dec24uofa - added TELLURIC and NOPLOT keywords
;       jm04feb29uofa - bug fix!  if DEREDDEN=1 then the program was
;                       incorrectly not de-reddening the error
;                       spectrum 
;       jm05jun20uofa - name changed from ISPEC to IDOSPEC to support
;                       the new routine IDOEXTRACT
;       jm05jun23uofa - call IDIVIDE_TELLURIC() to remove the telluric
;                       absorption from the extracted 1D spectrum
;       jm05jul07uofa - OUTNAME is now just the root of the output 
;                       name: the extraction aperture and the suffix
;                       ('.ms.fits') is appended within this routine
;
; Copyright (C) 2001, 2003-2005, John Moustakas
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

pro idospec, flist, specinfo, datapath=datapath, outname=outname, outpath=outpath, $
  specnames=specnames, traceinfo=traceinfo, _extra=extra, silent=silent, $
  deredden=deredden, telluric=telluric, noplot=noplot, gzip=gzip, wfits=wfits, $
  galaxy=galaxy

    if n_elements(flist) eq 0L then begin
       print, 'Syntax - idospec, flist, specinfo, datapath=, outname=, outpath=, $'
       print, '   specnames=, traceinfo=, _extra=extra, /silent, /deredden, $'
       print, '   /telluric, /noplot, /gzip, /wfits'
       return
    endif

; read the data
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(outpath) eq 0L then outpath = datapath
    if n_elements(noplot) eq 0L then noplot = 0L

    nimage = n_elements(flist)

    noutname = n_elements(outname)
    if (noutname gt 0L) then if (noutname ne nimage) then begin
       splog, 'Dimensions of OUTNAME and FLIST do not agree.'
       return
    endif

; GALAXY is an internal keyword for my research

    ngalaxy = n_elements(galaxy)
    if (ngalaxy gt 0L) then if (ngalaxy ne nimage) then begin
       splog, 'Dimensions of GALAXY and FLIST do not agree.'
       return
    endif
    
; extract a one-dimensional spectrum

    specnames = strarr(nimage)
    for i = 0L, nimage-1L do begin

       cube = rd2dspec(flist[i],datapath=datapath,/silent)

       if (total(cube.sky) eq 0.0) then begin
          splog, 'Spectrum '+flist[i]+' must be sky subtracted.'
          continue
       endif

       if (total(cube.wavemap) eq 0.0) then begin
          splog, 'No wavelength map found for '+flist[i]+'.'
          continue
       endif

       if (not keyword_set(silent)) then splog, 'Extracting '+cube.fname+'.'
       idoextract, cube.image, cube.sigmamap, cube.mask, *cube.header, $
         specinfo, skyimage=cube.sky, wavemap=cube.wavemap, traceinfo=traceinfo, $
         noplot=noplot, _extra=extra
       header1d = specinfo.header
       
; add the name of the 2D image to the header

       sxaddpar, header1d, 'SPEC2D', flist[i], ' two-dimensional spectrum file name', $
         before='HISTORY'

       if keyword_set(deredden) then begin
          
          ra = 15.0*im_hms2dec(sxpar(header1d,'RA'))    ; [degree]
          de = im_hms2dec(sxpar(header1d,'DEC'))        ; [degree]
          glactc, ra, de, 2000.0, gl, gb, 1, /degree ; Galactic coordinates [degree]

; compute the average reddening from the Schlegel et al. & Burstein &
; Heiles (BH) maps; valid BH reddenings are always greater than -0.22
          
          ebv_sfd = dust_getval(gl,gb,/interp) ; E(B-V) color excess [Schlegel et al 1998]
          ebv_bh = dust_getval(gl,gb,map='BH') ; [Burstein & Heiles 1982]

          if (ebv_bh lt -0.22) then begin
             splog, 'WARNING: No valid Burstein & Heiles E(B-V) value!'
             ebv_bh = ebv_sfd 
          endif else ebv_bh = ebv_bh/4.0
          if (ebv_bh lt 1E-3) then ebv_bh = 0.0

          ebv_array = [ebv_sfd]
;         ebv_array = [ebv_sfd,ebv_bh]
          
          ebv = djs_mean(ebv_array)
          ebv_err = djsig(ebv_array)
          
          kl = k_lambda(specinfo.wave,/odonnell,R_V=3.1)

          specinfo.spec = specinfo.spec*10.0^(0.4*ebv*kl)
          specinfo.sigspec = sqrt((specinfo.sigspec*10.0^(0.4*ebv*kl))^2 + $
            (specinfo.spec*0.4*kl*alog(10)*ebv_err)^2)

          sxaddpar, header1d, 'EBV', float(ebv), $
            ' Galactic reddening [mag]', before='HISTORY'
          sxaddpar, header1d, 'EBV_ERR', float(ebv_err), $
            ' Galactic reddening error [mag]', before='HISTORY'
          sxaddhist, "'Spectrum corrected for Galactic extinction "+im_today()+"'", header1d
          if not keyword_set(silent) then splog, 'Correcting for Galactic E(B-V) = '+$
            string(ebv,format='(F6.4)')+' +/- '+string(ebv_err,format='(F6.4)')+'.'
          
       endif
       
       if keyword_set(telluric) then begin

          if ((tag_exist(telluric_flux,'CUBE')) eq 0B) then begin
             splog, 'No telluric spectrum found for '+flist[i]+'.'
             continue
          endif

stop          
          
          if (not keyword_set(silent)) then splog, 'Correcting for telluric absorption.'
          idivide_telluric, speccube=speccube, telluric_flux=cube[i].telluric_flux, $
            telluric_wave=cube[i].telluric_wave, vmaxshift=vmaxshift, $
            _extra=extra

          sxaddhist, "'Spectrum corrected for telluric absorption "+im_today()+"'", header1d
             
       endif

; extraction aperture (along the slit)

       strap = fmtaperture(specinfo.aperture)

       if (noutname ne 0L) then outname1 = outname[i] else begin
          outname1 = strmid(flist[i],0,strpos(flist[i],'.fits'))+'_'+strap+'.ms.fits'
          outname1 = repstr(outname1,'.gz','')
       endelse
       
       specnames[i] = outname1
       if keyword_set(gzip) then specnames[i] = specnames[i]+'.gz'
       
; update the header in SPECINFO

       specinfo = struct_addtags(struct_trimtags(specinfo,except='HEADER'),{header: header1d})
          
       if keyword_set(wfits) then begin

          if (ngalaxy ne 0L) then begin
             sxdelpar, header1d, 'GALAXY'
             sxaddpar, header1d, 'GALAXY', galaxy[i], after='OBJECT'
          endif

          print, 'Writing '+outpath+outname1+'.'
          wrt1dspec, outname1, specinfo.spec, specinfo.sigspec, specinfo.sky, $
            specinfo.mask, header1d, datapath=outpath, gzip=gzip, _extra=extra
          
       endif 
       
       icleanup, cube
       
    endfor

return
end
