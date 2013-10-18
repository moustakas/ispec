;+
; NAME:
;	ICALIBRATE
;
; PURPOSE:
;	Rectify a two-dimensional spectrum both spatially and in the
;	wavelength dimension, and flux-calibrate.
;
; CALLING SEQUENCE:
;       icalibrate, caliblist, outcube, tracename=, wmapname=, $
;          sensname=, extfile=, extpath=, datapath=, minwave=, $
;          maxwave=, dwave=, crop_frac=, outpath=, tellfits=, $
;          /exclude_senserr, /gzip, /wfits
;
; INPUTS:
;	caliblist - file list of data to calibrate
;    
; OPTIONAL INPUTS:
;	tracename - name of the trace/distortion structure
;                   (IFITDISTORTION) 
;	wmapname  - name of the pixel-to-wavelength map (IARCFIT)
;	sensname  - FITS name of the sensitivity curve (ISENSFUNC) 
;       extfile   - name of the extinction file
;       extpath   - path to EXTFILE
;	datapath  - path to the data (default PWD)
;	minwave   - starting wavelength for each row of the output
;                   spectrum [Angstrom]
;	maxwave   - ending wavelength for each row of the output
;                   spectrum [Angstrom] 
;	dwave     - dispersion for each row of the output spectrum
;                   [Angstrom/pixel] 
;       crop_frac - if more than CROP_FRAC percent of the columns of
;                   any row within NEDGE rows of the edge of the CCD
;                   were extrapolated, then crop those rows (defaul
;                   0.05). for no cropping set CROP_FRAC = 1.0
;	outpath   - write the calibrated data to this path (default
;                   datapath) 
;       tellfits  - name of the FITS file containing the telluric
;                   absorption spectrum (usually the output from
;                   ICONSTRUCT_TELLURIC) to append to the output
;       extra     - keywords for PAD_EXTRAPOLATE()
;	
; KEYWORD PARAMETERS:
;       exclude_senserr - do not propagate the error in the
;                         sensitivity  function to the flux-calibrated
;                         error maps
;       gzip            - write GZIPPED FITS files
;	wfits           - write out the calibrated data 
;
; OUTPUTS:
;	outcube - output data cube
;
; COMMENTS:
;	This routine will do one or all of distortion-correct,
;	wavelength calibrate, flux calibrate, and extinction correct.
;	The spatial distortion correction is based on the results of
;	IFITDISTORTION and passed through the keyword TRACENAME.
;	Wavelength calibration can be done with one or more arc lamp.
;	If multiple arcs are passed then the (RA,DEC) of the science
;	object and the arcs is minimized.  Finally the data can be
;	flux calibrated and/or extinction corrected using a
;	user-specified sensitivity curve and extinction file.
;
;       The input spectra are written out according to the following
;       prefixes:  'd' for distortion, 'w' for wavelength calibration,
;       'e' for extinction correction, and 'f' for flux calibration.
;       If both extinction correction *and* flux calibration are being
;       applied then the 'f' prefix superseeds the 'e' prefix.  For
;       example, 'ra.0001.fits' would be written out as
;       'fwdra.0001.fits'.  
;
;       Interpolated pixels are assigned a bad pixel mask bit number
;       according to IMASK_BITS().  
;
;       Only one TELLFITS can be passed, and it is assumed that the
;       FITS file exists in DATAPATH.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD2DSPEC(), ICLEANUP, CMRESTORE, TRACESET2XY, MRDFITS(),
;	SXADDPAR, SXADDHIST, SXPAR(), MWRFITS(), MAKE_WAVE(), CWD(),
;	WRT2DSPEC, SPLOG, HMS2DEC(), DJS_DIFF_ANGLE(),
;	PAD_EXTRAPOLATE(), REMOVE, IRDSENSFUNC(), IMASK_BITS(),
;	IRESTORE_WMAP(), REPSTR(), ICHOOSE_WMAP(), XY2TRACESET
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 17, U of A - written
;	jm02feb20uofa - distortion and wavelength-calibration combined
;	   into one procedure using triangulation
;       jm02oct23uofa - identified a bug in TRIANGULATE and/or TRIGRID  
;       jm02nov07uofa - added support for multiple wavelength maps
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr14uofa - added CROP_FRAC parameter
;       jm03apr29uofa - added correct propagation of sensitivity
;                       function error but see NOSENSERR
;       jm03dec07uofa - added GZIP keyword; made IRESTORE_WMAP its own
;                       routine
;       jm03dec08uofa - added support for 2D sky subtraction
;       jm03dec30uofa - added TELLURIC keyword
;       jm05jun23uofa - changed NOSENSERR keyword to EXCLUDE_SENSERR;
;                       the default is to propagate the error in the
;                       sensitivity function; removed TELLURIC keyword
;                       in favor of the TELLFITS optional input (see
;                       implementation in IFLUXCALIBRATE) 
;       jm05jun29uofa - removed IRAF EX-FLAG header keyword
;       jm06jan24uofa - replace the IMAGE output with OUTCUBE
;
; Copyright (C) 2001-2003, 2005, John Moustakas
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

pro icalibrate, caliblist, outcube, tracename=tracename, wmapname=wmapname, $
  sensname=sensname, extfile=extfile, extpath=extpath, datapath=datapath, $
  minwave=minwave, maxwave=maxwave, dwave=dwave, crop_frac=crop_frac, $
  outpath=outpath, tellfits=tellfits, exclude_senserr=exclude_senserr, $
  gzip=gzip, wfits=wfits, _extra=extra

    if (n_elements(caliblist) eq 0L) then begin
       print, 'Syntax - icalibrate, caliblist, image, tracename=, wmapname=, $'
       print, '   sensname=, extfile=, extpath=, datapath=, minwave=, maxwave=, $'
       print, '   dwave=, crop_frac=, outpath=, tellfits=, /exclude_senserr, /gzip, $'
       print, '   /wfits'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = datapath

    if n_elements(extpath) eq 0L then $
      extpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
       
    tcalib = n_elements(tracename)
    wcalib = n_elements(wmapname)
    fcalib = n_elements(sensname)
    extcor = n_elements(extfile)

    if (tcalib[0] eq 1L) then if (tracename[0] eq '') then tcalib = 0L
    if (wcalib[0] eq 1L) then if (wmapname[0] eq '') then wcalib = 0L
    if (fcalib[0] eq 1L) then if (sensname[0] eq '') then fcalib = 0L
    if (extcor[0] eq 1L) then if (extfile[0] eq '') then extcor = 0L
    
    if (tcalib eq 0L) and (wcalib eq 0L) and (fcalib eq 0L) and (extcor eq 0L) then begin
       splog, 'No calibration data passed [TRACENAME, WMAPNAME, SENSNAME, EXTFILE].'
       return
    endif

    nimage = n_elements(caliblist)       ; number of images
    outname = repstr(caliblist,'.gz','') ; root output file names

; extrapolation parameters
    
    nedge = 20L
    if n_elements(crop_frac) eq 0L then crop_frac = 0.05
    
; restore the trace/distortion structure

    if (tcalib ne 0L) then begin

       if file_test(datapath+tracename,/regular) then begin

          splog, 'Restoring '+tracename+'.'
          cmrestore, datapath+tracename, tracefit, /quiet
          traceinfo = tracefit

          npad = traceinfo.npad
          dmap = traceinfo.dmap ; distortion map padded by 2*NPAD pixels spatially 
          
       endif else begin

          splog, 'Unable to find trace structure '+tracename+'.'
          return

       endelse

    endif

; initialize the wavelength parameters and restore the wavelength
; map(s) 

    if (wcalib ne 0L) then begin
    
       if n_elements(minwave) ne 0L then minwave1 = minwave
       if n_elements(maxwave) ne 0L then maxwave1 = maxwave

       awavemap = irestore_wmap(wmapname,datapath=datapath,minwave_in=minwave1,$
         maxwave_in=maxwave1,dwave=adwave,minwave_out=aminwave,maxwave_out=$
         amaxwave,wheader=awheader)
       if awavemap[0] eq -1L then return

       if n_elements(dwave) ne 0L then adwave = adwave*0.0+dwave

    endif
       
; read in the extinction curve

    if (extcor ne 0L) then begin

       if file_test(extpath+extfile,/regular) then begin
          
          splog, 'Reading the extinction file '+extpath+extfile+'.'
          readcol, extpath+extfile, extwave, extvals, format='F,F', /silent

       endif else begin

          splog, 'Extinction file '+extpath+extfile+' not found.'
          return

       endelse 

    endif

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

    if (fcalib ne 0L) then begin
    
       if file_test(datapath+sensname,/regular) then begin
          
          splog, 'Reading '+sensname+'.'

          senscube = irdsensfunc(sensname,datapath=datapath,/silent)
          senshead = senscube.header
          
; this may be a little confusing.  use either the observed or the
; grey-shifted sensitivity function to flux-calibrate, but always use
; the error in the *grey-shifted* function as the error in the
; sensitivity function.  basically this is because this routine only
; propagates *relative* spectrophotometric errors (from the
; grey-shifted function) and not absolute errors
          
          if strmatch(sensname,'*grey*') eq 1B then $
            tsens = senscube.greysens else $
            tsens = senscube.sens

          tsenserr = senscube.greysenserr ; senstivity function error [mag]
          senserr_abs = senscube.senserr  ; absolute error [mag]

;         relerr = 10.0^(0.4*tsenserr)-1.0     ; [fraction]
;         relmaxerr = 100.0*max(relerr)        ; maximum relative error [percent]
;         relmeanerr =  100.0*djs_mean(relerr) ; mean relative error [percent]

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

    endif else begin

       telluric_spec = 0.0
       telluric_head = ''

    endelse

; loop on each image

    for i = 0L, nimage-1L do begin

       cube = rd2dspec(caliblist[i],datapath=datapath,wset=wset)

       ncols = cube.naxis1
       nrows = cube.naxis2
       colaxis = findgen(ncols)
       rowaxis = findgen(nrows)
    
       image = cube.image
       sigmap = cube.sigmamap
       skyimage = cube.sky
       mask = cube.mask
       header = *cube.header

       oimage = image           ; original data
       oskyimage = skyimage
       osigmap = sigmap
       omask = mask

; if not distortion-correcting then initialize a simple distortion map
; with no padding (extrapolation)

; ---------------------------------------------------------------------------       
; ---------------------------------------------------------------------------       
       if (tcalib eq 0L) then begin

          npad = 0L
          dmap = (fltarr(ncols)+1) # findgen(nrows)

       endif
; ---------------------------------------------------------------------------       
; ---------------------------------------------------------------------------       

; if either wavelength calibrating or distortion correcting then
; initialize the inputs to TRIGRID()
       
       if (tcalib ne 0L) or (wcalib ne 0L) then begin
       
; if wavelength calibrating then the horizontal spacing (GS) of the
; output image is DWAVE over the limits [MINWAVE,MAXWAVE], otherwise
; the spacing is constant in pixels.  the vertical spacing limits are
; dictated by the distortion map.  NEWXMAP is the desired output
; wavelength map.
       
          if (wcalib ne 0L) then begin ; wavelength-calibrating

; minimize the difference in RA, DEC to select the appropriate
; wavelength map

             indx = ichoose_wmap(header,awheader,wmapname=wmapname)
             wavemap = reform(awavemap[*,*,indx])
             wmapname_out = wmapname[indx]
             dwave = (adwave[indx])[0]
             minwave = (aminwave[indx])[0]
             maxwave = (amaxwave[indx])[0]

             splog, format='("Wavelength specifications: [",F7.2,", ",F8.2,", ",F4.2,", ",I4,"].")', $
               minwave, maxwave, dwave, ncols

             xmap = wavemap
             newxmap = (findgen((maxwave-minwave)/dwave+1)*dwave+minwave) # $
               (fltarr(nrows+2*npad)+1)
             
             gs = [dwave,1.0]                              ; [horizontal,vertical] spacing
             limits = [minwave,-npad,maxwave,npad+nrows-1] ; [horizontal,vertical] limits

; pad and extrapolate the wavelength map
             
; ---------------------------------------------------------------------------             
; ---------------------------------------------------------------------------             
             xmap = pad_extrapolate(wavemap,npad,ncoeff=5,_extra=extra)
; ---------------------------------------------------------------------------             
; ---------------------------------------------------------------------------             

             outname[i] = 'w'+outname[i] ; output image name
       
          endif else begin ; no wavelength map

             xmap = findgen(ncols) # (fltarr(nrows+2*npad)+1)
             newxmap = xmap

             gs = [1.0,1.0]
             limits = [0,-npad,ncols-1,npad+nrows-1]
          
          endelse

          if (tcalib ne 0L) then begin ; distortion-correcting
             
; ---------------------------------------------------------------------------             
; ---------------------------------------------------------------------------             
; pad and extrapolate the image and the sigma map spatially; pad the
; bad pixel mask and make every extrapolated pixel bad
             
             image = pad_extrapolate(oimage,npad,ncoeff=2,_extra=extra)
             skyimage = pad_extrapolate(oskyimage,npad,ncoeff=2,_extra=extra)
             sigmap = pad_extrapolate(osigmap,npad,ncoeff=2,_extra=extra) ; <-- wrong?!?

             mask = make_array(ncols,nrows+2*npad,value=-1.0,/int)
             mask[*,npad:npad+nrows-1] = omask ; 0 is good
; ---------------------------------------------------------------------------             
; ---------------------------------------------------------------------------             

             outname[i] = 'd'+outname[i] ; output image name

          endif

; undistort and wavelength calibrate.  the bit of code using
; INTERPOL() is effectively the same as the TRIGRID() code if not
; distortion correcting.  QUINTIC interpolation takes substantially
; longer but I am not sure which is more correct
          
;;        splog, 'Interpolating '+cube.fname+'.'
;;
;;        interpimage = float(newxmap*0.0)
;;        interpskyimage = float(newxmap*0.0)
;;        interpsigmap = float(newxmap*0.0)
;;        interpmask = fix(newxmap*0.0)
;;
;;        for k = 0L, nrows+2*npad-1L do begin
;;           interpimage[*,k] = interpol(image[*,k],xmap[*,k],newxmap[*,k])
;;           interpskyimage[*,k] = interpol(skyimage[*,k],xmap[*,k],newxmap[*,k])
;;           interpsigmap[*,k] = interpol(sigmap[*,k],xmap[*,k],newxmap[*,k])
;;           interpmask[*,k] = interpol(mask[*,k],xmap[*,k],newxmap[*,k])
;;        endfor

          splog, 'TRIANGULATING '+cube.fname+'.'

          triangulate, xmap, dmap, tr, b
          interpimage = trigrid(xmap,dmap,image,tr,gs,limits,quintic=quintic)
          interpskyimage = trigrid(xmap,dmap,skyimage,tr,gs,limits,quintic=quintic)
          interpsigmap = trigrid(xmap,dmap,sigmap,tr,gs,limits,quintic=quintic)
;         interpmask = trigrid(xmap,dmap,mask,tr,gs,limits,missing=-1.0)

;        djs_plot, findgen(n_elements(interpimage[*,0])), interpimage[*,111], $
;          ps=10, ysty=3, xsty=3
;        djs_oplot, colaxis, image[*,111], ps=10, color='green'
          
; in order to properly rectify the bad pixel mask we have to generate
; a bad pixel mask that is comprised of only one bad pixel value [see
; IMASK_BITS()], then run TRIGRID() through that mask and then co-add
; all the masks.

          outmask = fix(interpimage*0)
          junk = imask_bits(maxpower=maxpower)
          for j = 1L, long(maxpower) do begin

             bitmask = mask*0
             w = where(mask eq fix(2.0^j),nw)
             if (nw ne 0L) then begin
                bitmask[w] = 1.0
                ibitmask = trigrid(xmap,dmap,bitmask,tr,gs,limits,missing=-1.0)
                neww = where(ibitmask gt 0.0,nneww)
                if (nneww ne 0L) then outmask[neww] = outmask[neww] + fix(2.0^j)
             endif

          endfor

; crop the padded rows
          
          image = float(interpimage[*,npad:npad+nrows-1])
          skyimage = float(interpskyimage[*,npad:npad+nrows-1])
          sigmap = float(interpsigmap[*,npad:npad+nrows-1])
          mask = fix(outmask[*,npad:npad+nrows-1])

; if more than CROP_FRAC percent of the columns of any row within NEDGE
; rows of the edge of the CCD were extrapolated, then crop those rows
; as well  

          if (tcalib ne 0L) then begin
          
             newrows = long(rowaxis)
             
             junk = total((mask gt 0.0),1)
             extrap = where(junk gt crop_frac*ncols,nextrap)
             if nextrap ne 0L then keep = where((extrap lt nedge) or (extrap gt nrows-nedge),nextrap)

             if nextrap ne 0L then begin

                extrap = extrap[keep]
                remove, extrap, newrows

                nrows = n_elements(newrows)
                image = image[*,newrows]
                skyimage = skyimage[*,newrows]
                sigmap = sigmap[*,newrows]
                mask = mask[*,newrows]
                
                sxaddpar, header, 'NAXIS2', nrows ; update the header          
                sxaddpar, header, 'CRPIX2', nrows/2L+1L

             endif

; ---------------------------------------------------------------------------
; assign any extrapolated pixels an appropriate bit number
; ---------------------------------------------------------------------------
;            w = where((mask ne fix(0.0)) and (mask ne fix(2.0^1)) and $
;              (mask ne fix(2.0^2)) and (mask ne fix(2.0^3)) and $
;              (mask ne fix(6)) and (mask ne fix(10)) and (mask ne fix(12)) and $
;              (mask ne fix(14)),nw)
;            if nw ne 0L then begin
;               message, 'This needs to be generalized!'
;               mask[w] = fix(2.0^4)
;            endif
; ---------------------------------------------------------------------------

          endif
             
          if (wcalib ne 0L) then begin

; wavelength calibrate the wavelength map for testing; do not crop
; padded rows because we use the padded INTERPMAP below to construct
; the output WSET

             interpmap = trigrid(xmap,dmap,xmap,tr,gs,limits,quintic=quintic)
;            interpmap = float(interpmap[*,npad:npad+nrows-1])
             nx = (size(interpmap,/dimension))[0]

;; ----------------------------------------------------------------------          
;; test plots
;; ----------------------------------------------------------------------          
;
;; this plot shows the sense of the shift applied
;          
;             djs_plot, xmap[*,nrows/2L]-interpmap[*,nrows/2L], xsty=3, ysty=3, ps=10
;
;; these plots are the killers! clearly the correct interpolation is
;; being applied here!
;             
;            djs_plot, xmap[*,nrows/2L], image[*,nrows/2L], xsty=3, ysty=3, ps=10, xrange=[4400,4510]
;            djs_oplot, interpmap[*,nrows/2L], interpimage[*,nrows/2L], ps=10, color='green'
;
;            djs_plot, xmap[ncols/2L,*], xsty=3, ysty=3, ps=10
;            djs_oplot, interpmap[ncols/2L,*], xsty=3, ysty=3, ps=10, color='green'
;
;; ----------------------------------------------------------------------          

; check for extrapolated columns

             good = lindgen(nx)
             for ix = 0L, nx-1L do begin
                zero = where((interpmap[ix,*] eq 0.0),nzero)
                if (nzero gt crop_frac*nrows) then good[ix] = -1L
             endfor 

             extrap = where((good eq -1L),nextrap,comp=keep,ncomp=ngood)
             if (ngood eq 0L) then message, 'This should never happen!'
             good = good[keep]
stop             
             if (nextrap ne 0L) then begin
                splog, 'Cropping '+strn(nextrap)+' extrapolated column(s) in '+outname[i]+'.'
                image = image[good,*]
                skyimage = skyimage[good,*]
                sigmap = sigmap[good,*]
                mask = mask[good,*]
                interpmap = interpmap[good,*]
                nx = ngood
             endif

; old, less reliable code (jm08jan17nyu)
;
;            im1d = total(interpmap,2)
;            imzero = where(im1d eq float(0),nzero,comp=good,ncomp=ngood)
;            if nzero ne 0L then begin
;               splog, 'Cropping '+strn(nzero)+' extrapolated column(s) in '+outname[i]+'.'
;               image = image[good,*]
;               skyimage = skyimage[good,*]
;               sigmap = sigmap[good,*]
;               mask = mask[good,*]
;               interpmap = interpmap[good,*]
;               nx = ngood
;            endif
      
             w = where(interpmap gt 0.0)
             minwave = min(interpmap[w])
             maxwave = max(interpmap[w])

; compute the new dispersion - thanks to A. Marble 03may18
             
;            dwave = interpmap[nx/2+1,nrows/2]-interpmap[nx/2,nrows/2]
             dwave = (interpmap[nx-1,nrows/2] - interpmap[0,nrows/2]) / (nx-1.0)
             if (dwave lt 0.0) then message, 'DWAVE should never be zero!'

             splog, 'Starting wavelength '+strn(minwave,length=7)+'.'
             splog, 'Ending   wavelength '+strn(maxwave,length=7)+'.'
             splog, 'Spectral dispersion '+strn(dwave,length=7)+'.'
             splog, 'Number of columns '+string(nx,format='(I4)')+'.'

; generate the output wavelength map trace set
             
             xxmap = rebin(findgen(nx),nx,nrows)
             xy2traceset, xxmap, interpmap[*,npad:nrows+2*npad-npad-1L], wset, ncoeff=2L, /silent
;            xy2traceset, xxmap, interpmap, wset, ncoeff=2L, /silent ; linear solution; jm05jul07uofa 
             
; update the headers, including particular IRAF keywords

             sxaddpar, header, 'NAXIS1', nx

;            sxaddpar, header, 'WAT0_001', 'system=world', before='HISTORY' 
;            sxaddpar, header, 'WAT1_001', 'wtype=linear', before='HISTORY' 
;            sxaddpar, header, 'WAT2_001', 'wtype=linear label=Wavelength units=Angstroms', before='HISTORY' 
;            sxaddpar, header, 'WCSDIM', 2, before='HISTORY' 
;            sxaddpar, header, 'DC-FLAG', 0, ' [0 or 1] linear dispersion sampling', before='HISTORY' 
             
             sxaddpar, header, 'CRVAL1', float(minwave), ' wavelength at CRPIX1', before='HISTORY'
             sxaddpar, header, 'CRPIX1', float(1.0), ' reference pixel number', before='HISTORY'
             sxaddpar, header, 'CD1_1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
             sxaddpar, header, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
             sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type', before='HISTORY'

             sxaddpar, header, 'IWAVE', caliblist[i], ' input to ICALIBRATE [wavelength calibration]', $
               before='HISTORY'
             sxaddpar, header, 'WMAPNAME', wmapname_out, ' wavelength map name', before='HISTORY'

             sxaddhist, "'Linearized wavelength scale applied "+im_today()+"'", header

          endif
             
          if (tcalib ne 0L) then begin
          
             sxaddpar, header, 'ITRACE', caliblist[i], ' input to ICALIBRATE [distortion correction]', $
               before='HISTORY'
             sxaddpar, header, 'TRACENAM', tracename, ' distortion map name', before='HISTORY'

             sxaddhist, "'Distortion correction applied "+im_today()+"'", header

          endif
          
       endif ; close the condition of TCALIB or WCALIB

; check interpolated images
       
;      djs_plot, image[*,63], xsty=3, ysty=3, xr=[500,600], ps=10
;      djs_oplot, oimage[*,63], ps=10, color='green'
       
; verify that the image has been wavelength calibrated then extinction
; correct 
       
       if (fcalib ne 0L) or (extcor ne 0L) then begin

          wave = make_wave(header,cd1_1=dwave) ; wavelength vector
          if wave[0] eq -1L then begin
             splog, 'Image '+cube.fname+' (#'+strn(i)+') is not wavelength calibrated.'
             icleanup, cube, /main
          endif
          
       endif

       if extcor then begin ; extinction correct

          junk = sxpar(header,'EXTCOR',count=extcount) ; has this spectrum been extinction corrected already?
          if extcount ne 1L then begin

             airmass= sxpar(header,'AIRMASS')
             extinct = interpol(extvals,extwave,wave) # (fltarr(nrows)+1.0) ; [NCOLS,NROWS]

             splog, 'Extinction correcting '+cube.fname+' to zero airmass.'
             image = image * 10D^(0.4*airmass*extinct)
             skyimage = skyimage * 10D^(0.4*airmass*extinct)
             sigmap = sigmap * 10D^(0.4*airmass*extinct)

; update the header

             sxaddpar, header, 'EXTCOR', 'T', ' extinction-corrected [T/F]', before='HISTORY'
;            sxaddpar, header, 'EX-FLAG', long(0), before='HISTORY' ; extinction-corrected
             sxaddpar, header, 'IEXTNAME', extfile, ' extinction file name', before='HISTORY'
             sxaddhist, "'Extinction corrected to zero airmass "+im_today()+"'", header

             if (fcalib eq 0L) then outname[i] = 'e'+outname[i]

          endif else begin

             splog, 'Image '+cube.fname+', has already been extinction corrected.'

          endelse
          
       endif

; flux calibrate
       
       if fcalib ne 0L then begin

          sensfunc = interpol(sens,senswave,wave) # (fltarr(nrows)+1)       ; [NCOLS,NROWS]
          sensfuncerr = interpol(senserr,senswave,wave) # (fltarr(nrows)+1) ; [NCOLS,NROWS]

          exptime = sxpar(header,'EXPTIME')

          junk = sxpar(header,'EXTCOR',count=extcount)
          if extcount ne 1L then splog, 'WARNING: Image '+cube.fname+' has not been extinction corrected.'

          splog, 'Flux calibrating '+cube.fname+'.'

          sigmap = sqrt((sigmap / sensfunc)^2 + $
            (image * sensfuncerr / sensfunc^2)^2) / exptime / dwave

          image = image / sensfunc / exptime / dwave
          skyimage = skyimage / sensfunc / exptime / dwave

; update the header

          sxaddpar, header, 'FLUXCOR', 'T', ' flux-calibrated [T/F]', before='HISTORY'
          sxaddpar, header, 'ZUNITS', 'erg/s/cm2/A', ' flux units', before='HISTORY'

          sxaddpar, header, 'IFLUX', caliblist[i], ' input to ICALIBRATE [flux calibration]', $
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
          
          outname[i] = 'f'+outname[i]

       endif 

; return an output data cube       
       
       if arg_present(outcube) then begin
          outcube1 = {fname: cube.fname, outname: outname[i], object: cube.object, $
            type: cube.type, naxis1: ncols, naxis2: nrows, $
            image: float(image), sigmamap: float(sigmap), $
            mask: mask, sky: float(skyimage), wavemap: wavemap, $
            wsettags: tag_names(wset)}
          outcube1 = struct_addtags(outcube1,wset)
          if (i eq 0L) then outcube = outcube1 else outcube = struct_append(outcube,outcube1)
       endif
    
; write out the calibrated images

       if keyword_set(wfits) then begin
          
          splog, 'Writing '+outpath+outname[i]+'.'
          wrt2dspec, outname[i], float(image), float(sigmap), mask, header, $
            skyimage=float(skyimage), wset=wset, telluric_spec=tellspec, $
            telluric_head=tellhead, datapath=outpath, gzip=gzip
          
       endif

       icleanup, cube ; clean up memory
    
    endfor

    if fcalib ne 0L then icleanup, senscube
    
return
end
