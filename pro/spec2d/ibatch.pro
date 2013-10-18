;+
; NAME:
;	IBATCH
;
; PURPOSE:
;	Batch reduce one night of data using an input parameter file.  
;
; CALLING SEQUENCE:
;       ibatch, paramfile, info, datapath=, proclist=, tracelist=, $
;          skyapfile= crsplitfile=, crlist=, stdlist=, caliblist=, $
;          sensinfo=, _extra=extra, /ccdproc, /distortion, /arcfit, $
;          /skysub, /crsplits, /crclean, /makesens, /calibrate, $
;          /fluxcalibrate
;
; INPUTS:
;	paramfile   - input parameter file
;
; OPTIONAL INPUTS:
;	datapath    - path to the parameter file and the data
;	proclist    - string array list of the FITS files to process 
;                     with ICCDPROC
;	tracelist   - string array list of files with which to model
;                     the spatial distortion
;	skyapfile   - sky aperture file output from ISKYSELECT; two-
;                     dimensional sky subtraction
;	crsplitfile - text file which contains cosmic ray split images
;                     that must be combined; each line can have any
;                     number of two or more FITS names
;	crlist      - string array list of files from which to flag
;                     and clean cosmic rays
;	stdlist     - string array list of standard stars with which
;                     to construct a sensitivity function
;	caliblist   - string array list of processed date to
;                     wavelength and flux-calibrate
;	extra       - extra keywords for ICCDPROC, ICRCOMBINE, ISPEC
;                     LA_COSMIC_WRAPPER, IARCFIT, IFITDISTORTION,
;                     ISENSFUNC, ICALIBRATE, and IFLUXCALIBRATE;  can
;                     also be used to override any of the input
;                     parameters  
;	
; KEYWORD PARAMETERS:
;	makeflat      - create a flat field
;	ccdproc       - CCD processing (assumes the flat field exists) 
;       distortion    - determine the spatial distortion map
;	arcfit        - determine the wavelength solution map
;	skysub        - subtract the two-dimensional sky spectrum
;	crsplits      - combine cosmic-ray split images with rejection
;	crclean       - flag cosmic rays and interpolate over bad pixels
;	makesens      - generate a sensitivity function
;	calibrate     - wavelength and flux-calibrate
;	fluxcalibrate - flux-calibrate using IFLUXCALIBRATE
;
; OUTPUTS:
;	info - informational structure to verify that the correct
;              parameters were passed and read 
;
; OPTIONAL OUTPUTS:
;       sensinfo - structure output from ISENSFUNC(); for debugging 
;
; COMMENTS:
;       This routine is a wrapper to all the reduction utilities
;       available in ISPEC.  The keywords available are listed
;       above.  
;
;       In ${ISPEC_DIR}/etc/ibatch_template.txt there is an example
;       PARAMFILE.  
;
; EXAMPLES:
;	In all these examples we assume that a parameter file called
;	'ibatch.txt' exists which contains default parameters.  An
;	information structure called INFO is returned and contains the
;	values of the parameters used.
;
;	[1] Make a flat field with all defaults:
;
;		ibatch, 'ibatch.txt', info, /makeflat
;
;	[2] Given a flat field, process (bias, overscan, flatten) a
;	file list (proclist) of raw images:
;
;		ibatch, 'ibatch.txt', info, proclist=proclist, /ccdproc
;
;	[3] Derive a wavelength solution:
;
;		ibatch, 'ibatch.txt', info, /arcfit
;	
;	[4] To accomplish [1-3] try:
;
;		ibatch, 'ibatch.txt', info, proclist=proclist,
;		   /makeflat, /ccdproc, /arcfit
;
;	[5] To generate a sensitivity function we need to wavelength
;       calibrate the standard stars (given in stdlist), extract
;       one-dimensional spectra, and then fit the sensitivity curve.
;       Accomplish this procedure by typing:
;
;		ibatch, 'ibatch.txt', info, stdlist=stdlist, /makesens
;
; PROCEDURES USED:
;	ICCDPROC, IARCFIT, ISENSFUNC, DJS_READLINES(), CWD(),
;	IFITDISTORTION, ICALIBRATE, ICRCOMBINE, ILA_COSMIC_WRAPPER, 
;	IDOSPEC, ISENSFUNC, READCOL, SPLOG, ISKYSUBTRACT2D,
;	IFLUXCALIBRATE 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 December 14, U of A - written
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03jul14uofa - added SATURATION parameter to parameter file 
;       jm03dec07uofa - added GZIP support; added two dimensional sky
;                       subtraction ability
;       jm03dec31uofa - added TELLIST AND TELLURIC keywords
;       jm04jun23uofa - use SKYAPFILE rather than SKYLIST input to
;                       ISKYSUBTRACT2D 
;       jm04sep15uofa - IPARSE_PARAMFILE() updated to read LAMPNAME
;                       rather than LAMP
;       jm05apr15uofa - removed STDPATH and SKYSPECFILE parameter file
;                       keywords 
;       jm05jun17uofa - bug fix caused by removing the STDPATH
;                       keyword; combining CR splits *and* cosmic-ray
;                       rejection now occurs before sky subtraction,
;                       because CR's can really mess up the sky fit
;       jm05jun21uofa - altered MAKESENS procedure; now, IDOSPEC
;                       extracts wavelength-calibrated spectra using
;                       the 2D wavelength map, which are then passed
;                       to the new version of ISENSFUNC()
;       jm05jun22uofa - FLUXCALIBRATE keyword added; removed TELLLIST
;                       optional input and TELLURIC keyword
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

function iparse_paramfile, paramfile, info
; jm02sep27uofa
; parse the input parameter file
    
    paramtxt = djs_readlines(info.datapath+paramfile)
    paramtxt = paramtxt[where(strcompress(paramtxt,/remove) ne '')] ; remove blank lines

    nlines = n_elements(paramtxt)
    params = paramtxt

    comment = where(strmatch(paramtxt,'*#*') eq 1B,ncomment) ; crop trailing comments
    if ncomment ne 0L then for i = 0L, ncomment-1L do params[comment[i]] = $
      strtrim(strmid(paramtxt[comment[i]],0,strpos(paramtxt[comment[i]],'#'))) 

; parse the trim and overscan inputs
    
    indx = where(strmatch(params,'*TRIM*',/fold_case) eq 1B)
    trim = long((strsplit(params[indx],' ',/extract))[1:4])
    
    indx = where(strmatch(params,'*OVERSCAN*',/fold_case) eq 1B)
    overscan = long((strsplit(params[indx],' ',/extract))[1:4])

; parse the remaining parameters
    
    tagnames = tag_names(info)

    for i = 0L, nlines-1L do begin

       str = strcompress(strsplit(params[i],' ',/extract),/remove)
       str1 = str[0]
       w = where(strmatch(tagnames,'*'+str1+'*') eq 1B,nw)

       if (nw ne 0L) then begin
          
          if (strlowcase(str1) eq 'trim') or (strlowcase(str1) eq 'overscan') then begin

             str2 = strn(strmid(params[i],strlen(str1),strlen(params[i])-1L))
             info.(w) = long(strsplit(str2,' ',/extract))

          endif else begin

             str2 = strcompress(strmid(params[i],strlen(str1),strlen(params[i])-1L),/remove)
             if str2 ne '' then info.(w) = str2 ; non-blank value

          endelse

       endif

    endfor

return, info
end

pro ibatch, paramfile, info, datapath=datapath, proclist=proclist, $
  tracelist=tracelist, skyapfile=skyapfile, crsplitfile=crsplitfile, $
  crlist=crlist, stdlist=stdlist, caliblist=caliblist, sensinfo=sensinfo, $
  makeflat=makeflat, ccdproc=ccdproc, distortion=distortion, $
  arcfit=arcfit, skysub=skysub, crsplits=crsplits, crclean=crclean, $
  makesens=makesens, calibrate=calibrate, fluxcalibrate=fluxcalibrate, $
  _extra=extra

    nparam = n_elements(paramfile)
    
    if nparam eq 0L then begin
       print, 'Syntax - ibatch, paramfile, info, datapath=, proclist=, tracelist=, $'
       print, '   skyapfile= crsplitfile=, crlist=, stdlist=, caliblist=, $'
       print, '   sensinfo=, _extra=extra, /ccdproc, /distortion, /arcfit, $'
       print, '   /skysub, /crsplits, /crclean, /makesens, /calibrate, /fluxcalibrate'
       return
    endif

    if nparam gt 1L then begin
       splog, 'Multiple parameter files passed.'
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    
    if file_test(datapath+paramfile,/regular) eq 0L then begin
       splog, 'Unable to find '+datapath+paramfile+'.'
       return
    endif

; initialize the IBATCH information structure

    info = { $
      datapath:      datapath, $
      doplot:        1B, $
      gzip:          1B, $
      biasfile:      '', $
      saturation:    0.0,$
      domefile:      '', $
      domelist:      0B, $
      skyfile:       '', $
      darkfile:      '', $
      arcfile:       '', $
      arclist:       0B, $
      lampname:      '', $
      flatname:      '', $
      illumname:     '', $
      wmapname:      '', $
      tracename:     '', $
      sensname:      '', $
;     stdpath:       '', $
      extfile:       '', $
;     skyspecfile:   '', $
      badpixfile:    '', $
      gain:          0.0, $
      rdnoise:       0.0, $
      trim:          lonarr(4), $
      overscan:      lonarr(4), $
      pscale:        0.0, $
      bsorder_flat:  0L, $
      bsorder_sens:  0L, $
      minwave_guess: 0.0, $
      minwave_out:   0.0, $
      dwave:         0.0, $
;     aperture:      0.0, $
;     skyaperture:   0.0, $
      sigclip:       0.0, $
      objlim:        0.0, $
      sigfrac:       0.0, $
      makeflat:      keyword_set(makeflat), $
      ccdproc:       keyword_set(ccdproc),  $
      arcfit:        keyword_set(arcfit),   $
      makesens:      keyword_set(makesens), $
      crsplits:      keyword_set(crsplits), $
      crclean:       keyword_set(crclean),  $
      distortion:    keyword_set(distortion), $
      calibrate:     keyword_set(calibrate), $
      fluxcalibrate: keyword_set(fluxcalibrate), $
      proclist:      '', $
      objlist:       '', $
      stdlist:       '', $
      crlist:        '', $
      caliblist:     ''  $
      }

; parse the input parameter file

    info.datapath = datapath
    info = iparse_paramfile(paramfile,info)

; over-ride any of the input parameters using EXTRA

    if n_elements(extra) ne 0L then begin
       
       if tag_exist(extra,'DATAPATH') ne 0L then info.datapath = extra.datapath
       if tag_exist(extra,'DOPLOT') ne 0L then info.doplot = extra.doplot
       if tag_exist(extra,'GZIP') ne 0L then info.gzip = extra.gzip
       if tag_exist(extra,'BIASFILE') ne 0L then info.biasfile = extra.biasfile
       if tag_exist(extra,'SATURATION') ne 0L then info.saturation = extra.saturation
       if tag_exist(extra,'DOMEFILE') ne 0L then info.domefile = extra.domefile
       if tag_exist(extra,'DOMELIST') ne 0L then info.domelist = extra.domelist
       if tag_exist(extra,'SKYFILE') ne 0L then info.skyfile = extra.skyfile
       if tag_exist(extra,'DARKFILE') ne 0L then info.darkfile = extra.darkfile
       if tag_exist(extra,'ARCFILE') ne 0L then info.arcfile = extra.arcfile
       if tag_exist(extra,'ARCLIST') ne 0L then info.arclist = extra.arclist
       if tag_exist(extra,'LAMPNAME') ne 0L then info.lampname = extra.lampname
       if tag_exist(extra,'FLATNAME') ne 0L then info.flatname = extra.flatname
       if tag_exist(extra,'ILLUMNAME') ne 0L then info.illumname = extra.illumname
       if tag_exist(extra,'WMAPNAME') ne 0L then info.wmapname = extra.wmapname
       if tag_exist(extra,'TRACENAME') ne 0L then info.tracename = extra.tracename
       if tag_exist(extra,'SENSNAME') ne 0L then info.sensname = extra.sensname
;      if tag_exist(extra,'STDPATH') ne 0L then info.stdpath = extra.stdpath
       if tag_exist(extra,'EXTFILE') ne 0L then info.extfile = extra.extfile
;      if tag_exist(extra,'SKYSPECFILE') ne 0L then info.skyspecfile = extra.skyspecfile
       if tag_exist(extra,'BADPIXFILE') ne 0L then info.badpixfile = extra.badpixfile
       if tag_exist(extra,'GAIN') ne 0L then info.gain = extra.gain
       if tag_exist(extra,'RDNOISE') ne 0L then info.rdnoise = extra.rdnoise
       if tag_exist(extra,'TRIM') ne 0L then info.trim = extra.trim
       if tag_exist(extra,'OVERSCAN') ne 0L then info.overscan = extra.overscan
       if tag_exist(extra,'PSCALE') ne 0L then info.pscale = extra.pscale
       if tag_exist(extra,'BSORDER_FLAT') ne 0L then info.bsorder_flat = extra.bsorder_flat
       if tag_exist(extra,'BSORDER_SENS') ne 0L then info.bsorder_sens = extra.bsorder_sens
       if tag_exist(extra,'MINWAVE_GUESS') ne 0L then info.minwave_guess = extra.minwave_guess
       if tag_exist(extra,'MINWAVE_OUT') ne 0L then info.minwave_out = extra.minwave_out
       if tag_exist(extra,'DWAVE') ne 0L then info.dwave = extra.dwave
;      if tag_exist(extra,'APERTURE') ne 0L then info.aperture = extra.aperture
;      if tag_exist(extra,'SKYAPERTURE') ne 0L then info.skyaperture = extra.skyaperture
       if tag_exist(extra,'SKYAPFILE') ne 0L then skyapfile = extra.skyapfile
       if tag_exist(extra,'SIGCLIP') ne 0L then info.sigclip = extra.sigclip
       if tag_exist(extra,'OBJLIM') ne 0L then info.objlim = extra.objlim
       if tag_exist(extra,'SIGFRAC') ne 0L then info.sigfrac = extra.sigfrac

       if tag_exist(extra,'GREY') eq 0L then grey = 1L

    endif

; ----------------------------------------------------------------------
    if keyword_set(makeflat) then begin ; generate the master flat

       splog, 'MAKEFLAT: generating the flat field and illumination pattern.'

; DOMEFILE is a file name with a list of dome flats
       
       if info.domelist then begin 

          if file_test(datapath+info.domefile,/regular) then $
            readcol, datapath+info.domefile, domefile, format='A', /silent, comment='#' else begin
             splog, 'Dome file '+datapath+domefile+' not found!'
             return
          endelse
          if file_test(datapath+info.flatname,/regular) then $
            readcol, datapath+info.flatname, flatname, format='A', /silent, comment='#' else begin
             splog, 'Flat name list file '+datapath+flatname+' not found!'
             return
          endelse

       endif else begin

          domefile = info.domefile
          flatname = info.flatname

       endelse
       
       iccdproc, biasfile=info.biasfile, domefile=domefile, skyfile=info.skyfile, $
         darkfile=info.darkfile, datapath=datapath, trim=info.trim, overscan=info.overscan, $
         bsorder=info.bsorder_flat, flatname=flatname, illumname=info.illumname, $
         badpixfile=info.badpixfile, doplot=info.doplot, gzip=info.gzip, pscale=info.pscale, $
         gain=info.gain, rdnoise=info.rdnoise, saturation=info.saturation, /makeflat, _extra=extra

    endif

; ----------------------------------------------------------------------
    if keyword_set(ccdproc) then begin ; process all the objects

       if n_elements(proclist) eq 0L then begin
          splog, 'No list of images to process specified.'
          return
       endif

       splog, 'CCDPROC: CCD-processing the file list.'

; DOMEFILE is a file name with a list of dome flats
       
       if info.domelist then begin 

          if file_test(datapath+info.flatname,/regular) then $
            readcol, datapath+info.flatname, flatname, format='A', /silent, comment='#' else begin
             splog, 'Flat name list file '+datapath+flatname+' not found!'
             return
          endelse

       endif else flatname = info.flatname
       
       iccdproc, proclist, biasfile=info.biasfile, darkfile=info.darkfile, datapath=datapath, $
         trim=info.trim, overscan=info.overscan, flatname=flatname, $
         illumname=info.illumname, badpixfile=info.badpixfile, objlist=objlist, $
         pscale=info.pscale, gain=info.gain, rdnoise=info.rdnoise, _extra=extra, $
         saturation=info.saturation, doplot=info.doplot, gzip=info.gzip, /clean, /wfits

; update the information structure

       info.proclist = strjoin(proclist,' ')
       if n_elements(objlist) ne 0L then info.objlist = strjoin(objlist,' ')
       
    endif 
    
; ----------------------------------------------------------------------
    if keyword_set(distortion) then begin ; model the trace

       if n_elements(tracelist) eq 0L then begin
          splog, 'No list of trace images specified.'
          return
       endif

       splog, 'DISTORTION: modeling the spatial distortion.'
       
       tracename = info.tracename
       ifitdistortion, tracelist, datapath=datapath, tracename=tracename, $
         doplot=info.doplot, tinfo=tinfo, finaltrace=finaltrace, _extra=extra
       
       info.tracename = tracename ; update
       
    endif

; ----------------------------------------------------------------------
    if keyword_set(arcfit) then begin ; derive the wavelength solution

       splog, 'ARCFIT: determining the wavelength solution.'
       
; ARCFILE is a file name with a list of arc lamps
       
       if info.arclist then begin 

          if file_test(datapath+info.arcfile,/regular) then $
            readcol, datapath+info.arcfile, arcfile, format='A', /silent, comment='#' else begin
             splog, 'Arc lamp name list file '+datapath+info.arcfile+' not found!'
             return
          endelse

          if file_test(datapath+info.wmapname,/regular) then $
            readcol, datapath+info.wmapname, wmapname, format='A', /silent, comment='#' else begin
             splog, 'Wavelength map list file '+datapath+info.wmapname+' not found!'
             return
          endelse

       endif else begin

          arcfile = info.arcfile ; arc lamp image
          wmapname = info.wmapname

       endelse

       tracename = info.tracename ; distortion/trace structure

       narc = n_elements(arcfile)
       for iarc = 0L, narc-1L do begin

; if a distortion map exists then un-distort the arc lamp but do not
; crop any extrapolated rows
          
; jm06jul31uofa - removed; it is not required if 2D sky subtraction is
;                 used 
          
;         if file_test(datapath+tracename,/regular) ne 0L then begin
;            icalibrate, arcfile[iarc], tracename=tracename, crop_frac=1.0, $
;              gzip=info.gzip, /wfits
;            arcfile[iarc] = 'd'+arcfile[iarc]
;         endif

          iarcfit, arcfile[iarc], wmapname=wmapname[iarc], datapath=datapath, lamp=info.lampname, $
            xcoeff=[info.minwave_guess,info.dwave], doplot=info.doplot, _extra=extra, /write

       endfor
          
       info.wmapname = strjoin(wmapname,' ') ; update
       
    endif

; ----------------------------------------------------------------------
    if keyword_set(crsplits) then begin ; combine cr-split images

       if file_test(datapath+crsplitfile,/regular) then begin

; read the cr-split file.  remove blank and commented lines
          
          crfiles = djs_readlines(datapath+crsplitfile)
          crfiles = crfiles[where(strcompress(crfiles,/remove) ne '')]

          good = where(strmatch(crfiles,'#*') ne 1B,ngood)
          if ngood ne 0L then crfiles = crfiles[good] else begin
             splog, 'No valid CRSPLIT file names read in '+datapath+crsplitfile+'.'
             return
          endelse
          
          nsplits = n_elements(crfiles)
;         crlist = strarr(nsplits) ; output files
          
          for i = 0L, nsplits-1L do begin

             crinput = strsplit(crfiles[i],' ',/extract) ; cr-split image
             icrcombine, crinput, gzip=info.gzip, /wfits, _extra=extra
;            crlist[i] = crimname

          endfor

       endif else begin

          splog, 'CRSPLITS: cosmic ray split file not found...returning.'
          return

       endelse

;      info.crlist = strjoin(crlist,' ')
       
    endif

; ----------------------------------------------------------------------
    if keyword_set(crclean) then begin ; flag cosmic rays and write out the cleaned images

       if n_elements(crlist) eq 0L then begin
          splog, 'No list of cosmic ray images specified.'
          return
       endif

       splog, 'CRCLEAN: cleaning cosmic rays.'
       
       ila_cosmic_wrapper, crlist, datapath=datapath, sigclip=info.sigclip, $
         objlim=info.objlim, sigfrac=info.sigfrac, gain=info.gain, $
         readn=info.rdnoise, gzip=info.gzip, /wfits, _extra=extra

    endif
    
; ----------------------------------------------------------------------
    if keyword_set(skysub) then begin ; flag cosmic rays and write out the cleaned images

       if (n_elements(skyapfile) eq 0L) then begin
          splog, 'No sky aperture file given.'
          return
       endif

       splog, 'SKYSUB: sky subtracting the two-dimensional sky spectrum.'

; ARCFILE is a file name with a list of arc lamps

       if info.arclist then begin 

          if file_test(datapath+info.arcfile,/regular) then $
            readcol, datapath+info.arcfile, arcfile, format='A', /silent, comment='#' else begin
             splog, 'Arc lamp name list file '+datapath+info.arcfile+' not found!'
             return
          endelse

          if file_test(datapath+info.wmapname,/regular) then $
            readcol, datapath+info.wmapname, wmapname, format='A', /silent, comment='#' else begin
             splog, 'Wavelength map list file '+datapath+info.wmapname+' not found!'
             return
          endelse

       endif else begin

          arcfile = info.arcfile   ; arc lamp image
          wmapname = info.wmapname

       endelse

       iskysubtract2d, skyapfile, datapath=datapath, wmapname=wmapname, $
         tracename=info.tracename, gzip=info.gzip, doplot=info.doplot, $
         /wfits, _extra=extra

    endif
    
; ----------------------------------------------------------------------
    if keyword_set(makesens) then begin ; generate the sensitivity curve

       if n_elements(stdlist) eq 0L then begin
          splog, 'No list of standard star images specified.'
          return
       endif

       splog, 'MAKESENS: generating the sensitivity function.'

; wavelength-calibrate       

; WMAPNAME is a file name with a list of wavelength maps

       if (info.arclist) and (strcompress(info.wmapname,/remove) ne '') then begin 

          if file_test(datapath+info.wmapname,/regular) then $
            readcol, datapath+info.wmapname, wmapname, format='A', /silent, comment='#' else begin
             splog, 'Wavelength map list file '+datapath+info.wmapname+' not found!'
             return
          endelse

       endif else wmapname = info.wmapname

; jm08jan17nyu: do not pass _EXTRA because if SENSNAME is defined then
; it will overwrite the call SENSNAME=''
       
       icalibrate, stdlist, outcube, datapath=datapath, extfile='', $
         tracename='', wmapname=wmapname, sensname='', gzip=info.gzip, $
         minwave=info.minwave_out, dwave=info.dwave, /wfits;, _extra=extra

; extract one-dimensional counts spectra; jm08jan17nyu: some issues
; with optimal extraction, so do not do it (e.g., 05apr stars)

       newstdlist = strtrim(outcube.outname)
       ispec, newstdlist, datapath=datapath, minwave=info.minwave_out, $
         dwave=info.dwave, specnames=specnames, noplot=abs(info.doplot-1), $
         _extra=extra, /tracespec, /seeing, gzip=info.gzip, optimal=0, /wfits
;      idospec, stdlist, datapath=datapath, minwave=info.minwave_out, $
;        dwave=info.dwave, specnames=specnames, noplot=abs(info.doplot-1), $
;        _extra=extra, /tracespec, /seeing, gzip=info.gzip, /wfits

; derive the sensitivity function

       sensname = info.sensname
       sensinfo = isensfunc(specnames,datapath=datapath,sensname=sensname,$
         extfile=info.extfile,bsorder=info.bsorder_sens,grey=grey,$
         doplot=info.doplot,/write,_extra=extra)
       
; update the information structure

       info.sensname = sensname 
       info.stdlist = strjoin(stdlist,' ')
       
    endif 

; ----------------------------------------------------------------------
    if keyword_set(calibrate) then begin ; flux- and wavelength-calibrate the data

       if n_elements(caliblist) eq 0L then begin
          splog, 'No list of files specified to calibrate.'
          return
       endif

       splog, 'CALIBRATE: undistorting, wavelength-calibrating, and flux-calibrating the data.'
       
; WMAPNAME is a file name with a list of wavelength maps

       if (info.arclist) and (strcompress(info.wmapname,/remove) ne '') then begin 

          if file_test(datapath+info.wmapname,/regular) then $
            readcol, datapath+info.wmapname, wmapname, format='A', /silent, comment='#' else begin
             splog, 'Wavelength map list file '+datapath+info.wmapname+' not found!'
             return
          endelse

       endif else wmapname = info.wmapname

       icalibrate, caliblist, outcube, datapath=datapath, extfile=info.extfile, $
         tracename=info.tracename, wmapname=wmapname, sensname=info.sensname, $
         minwave=info.minwave_out, dwave=info.dwave, gzip=info.gzip, $
         /wfits, _extra=extra

    endif
    
; ----------------------------------------------------------------------
    if keyword_set(fluxcalibrate) then begin ; flux-calibrate the data

       if n_elements(caliblist) eq 0L then begin
          splog, 'No list of files specified to calibrate.'
          return
       endif

       splog, 'FLUXCALIBRATE: flux-calibrating.'

       ifluxcalibrate, caliblist, datapath=datapath, extfile=info.extfile, $
         sensname=info.sensname, gzip=info.gzip, /wfits, _extra=extra
       
    endif
    
return
end
