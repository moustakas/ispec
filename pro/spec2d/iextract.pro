;+
; NAME:
;	IEXTRACT
;       
; PURPOSE:
;	Extract a one-dimensional spectrum.
;
; CALLING SEQUENCE:
;       iextract, image, sigmamap, mask, header, specinfo, $
;          skyimage=, aperture=, loaperture=, upaperture, $
;          skyaperture=, skylower=, skyupper=, flanking=, $
;          inputspec=, refrow=, refwave=, fluxrange=, $
;          /tracespec, /noskysub, /noskyshift, /noskymask, /seeing, $
;          /noplot, /meanprofile, /optimal, /starfit, $
;          /silent, _extra=extra
;
; INPUTS:
;	image       - two-dimensional wavelength calibrated image frame
;	sigmamap    - corresponding two-dimensional error map
;	mask        - corresponding two-dimensional bad pixel mask
;	header      - corresponding FITS header for the image frame
;
; OPTIONAL INPUTS:
;	skyimage    - two-dimensional sky image frame
;	aperture    - extraction aperture diameter [arcsec]
;	loaperture  - lower extraction aperture radius [arcsec]
;	upaperture  - upper extraction aperture radius [arcsec]
;	skyaperture - two-element sky aperture diameter corresponding 
;                     to [SKYLOWER,SKYUPPER] [arcsec]; (the default is
;                     20% of NROWS-1); if SKYAPERTURE is a scalar then
;                     use the same aperture for both sky windows 
;       skylower    - distance from row 0 to start the lower sky
;                     aperture of diameter SKYAPERTURE [arcsec] (the
;                     default is 0.0)
;       skyupper    - distance from row NROWS-1 to end the upper sky
;                     aperture of diameter SKYAPERTURE [arcsec] (the
;                     default is 0.0)
;       flanking    - use FLANKING sky apertures (SKYLOWER and
;                     SKYUPPER are not used) [arcsec]; place the sky
;                     aperture a distance FLANKING+SKYAPERTURE away
;                     from APERTURE
;	inputspec   - input spectrum to subtract from the extracted
;                     one dimensional spectrum.  a structure with the
;                     following fields:
;	  wave           - one dimensional wavelength vector
;	  spec           - corresponding one dimensional spectrum
;         invvar         - corresponding inverse variance spectrum 
;       refrow      - extract the spectrum centered on REFROW.  if
;                     both REFROW and REFWAVE are set then REFROW
;                     takes precedence (the default is the row
;                     corresponding to the flux-weighted center of the
;                     mean spatial profile)
;	refwave     - reference wavelength.  center the extraction
;                     aperture at the column corresponding to REFWAVE
;                     (the default is the middle wavelength)
;	fluxrange   - vertical plot range for the spatial profile
;                     window [MIN,MAX]
;	extra       - keywords for ITRACESPEC, ISKYSUBTRACT, and 
;                     ISKYSHIFT
;	
; KEYWORD PARAMETERS:
;	tracespec   - trace the spectrum as a function of wavelength 
;	noskysub    - do not sky subtract
;	noskyshift  - do not use the night sky lines to improve the
;                     wavelength calibration
;	noskymask   - do not mask sky subtraction residuals
;       seeing      - estimate the seeing by fitting a Gaussian
;                     function to the spatial PSF (useful for point
;                     sources) 
;	noplot      - do not generate extraction windows
;       meanprofile - set this to examine the mean spatial profile
;                     while still being able to set REFROW; this
;                     keyword is stronger than REFWAVE
;       optimal     - optimal extraction (weight by variance)
;       starfit     - model and subtract a point source (e.g., a star)
;                     from the 2D spectrum before extraction [see
;                     I2DSTARFIT() for details and relevant 
;                     parameters] 
;       silent      - do not print messages to STDOUT
;
; OUTPUTS:
;	specinfo   - structure containing the following fields:
;	  spec     - one dimensional spectrum
;	  wave     - corresponding wavelength vector
;	  vspec    - one dimensional variance spectrum
;	  header   - spectrum header
;	  skyinfo  - sky subtraction information structure (see
;                    ISKYSUBTRACT) 
;	  apinfo   - aperture extraction parameters
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;	The spectra must be wavelength calibrated.  Partial pixels
;	corresponding to small extraction apertures are treated
;	properly.  
;
;       If LOAPERTURE, UPAPERTURE, and APERTURE are specified (silly!)
;       then APERTURE is overridden.
;
; EXAMPLE:
;	
; TODO:
;       [1] Set the extraction apertures interactively. 
;       [2] Find, trace, and extract multiple apertures on a single
;           image (for multi-object spectroscopy).
;       [3] Optionally allow an input trace.
;
; PROCEDURES USED:
;	SXPAR(), MAKE_WAVE(), MAKE_ARCSECAXIS(), FLAM_UNITS(),
;	GET_ELEMENT, ITRACESPEC, LEGEND, DJS_PLOTLIMITBOX,
;	ISKYSUBTRACT, ISKYSHIFT, COMBINE1FIBER, READCOL, SPLOG,
;	STRN(), SXADDHIST, ISKYMASK(), SXADDPAR, SXDELPAR, DJS_PLOT,
;	IM_HMS2DEC(), DELVARX, IEXTRACT_OPTIMAL, IM_WINDOW 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August & October 8, U of A
;	jm01dec27uofa - developed tracing algorithm
;       jm02oct16uofa - documented and optimized
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr16uofa - added SILENT keyword
;       jm03jul16uofa - if REFROW is not set then centroid on the
;                       mean spatial profile
;       jm03sep12uofa - added MEANPROFILE keyword
;       jm03dec8uofa  - added support for 2D sky subtraction and the
;                       NOSKYMASK keyword
;       jm04sep01uofa - added OPTIMAL keyword and IEXTRACT_OPTIMAL,
;                       developed by Andy Marble
;       jm05jan19uofa - previously, if either NOSKYSHIFT or NOSKYSUB
;                       were set, then no sky-line shift was applied;
;                       now, this only occurs if NOSKYSHIFT by itself
;                       is set 
;       jm05jul07uofa - HELIOCOR keyword removed; SEEING keyword
;                       added; LOAPERTURE and UPAPERTURE optional
;                       inputs added; generally cleaned up 
;       jm05jul20uofa - added STARFIT keyword
;       jm07jun20nyu  - fixed SEEING value written in header (the
;                       Gaussian sigma was being written, rather than
;                       the FWHM)
;
; Copyright (C) 2001-2005, 2007, John Moustakas
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

pro iextract, image, sigmamap, mask, header, specinfo, skyimage=skyimage, $
  aperture=aperture, loaperture=loaperture, upaperture=upaperture, $
  skyaperture=skyaperture, skylower=skylower, skyupper=skyupper, $
  flanking=flanking, inputspec=inputspec, refrow=refrow, refwave=refwave, $
  fluxrange=fluxrange, tracespec=tracespec, noskysub=noskysub, $
  noskyshift=noskyshift, noskymask=noskymask, seeing=seeing, $
  noplot=noplot, meanprofile=meanprofile, optimal=optimal, starfit=starfit, $
  silent=silent, _extra=extra

    if (n_elements(image) eq 0L) or (n_elements(sigmamap) eq 0L) or $
      (n_elements(mask) eq 0L) or (n_elements(header) eq 0L) then begin
       print, 'Syntax - iextract, image, sigmamap, mask, header, specinfo, $'
       print, '   skyimage=, aperture=, loaperture=, upaperture=, skyaperture=, $'
       print, '   skylower=, skyupper=, flanking=, inputspec=, refrow=, refwave=, $'
       print, '   fluxrange=, /tracespec, /noskysub, /noskyshift, /noskymask, $'
       print, '   /seeing, /noplot, /meanprofile, /optimal, /starfit, $'
       print, '   /silent, _extra=extra'
       return
    endif

    if (n_elements(skyimage) ne 0L) then if (total(skyimage) eq 0.0) then delvarx, skyimage
    
    light = 2.99792458D5 ; speed of light [km/s]

    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    midrow = nrows/2L
    midcol = ncols/2L
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)
    
    header1d = header
    objname = sxpar(header1d,'OBJECT',count=count_objname)
    if (count_objname eq 0L) then objname = ''

    wave = make_wave(header1d,crval1=crval1) ; wavelength vector
    dwave = sxpar(header1d,'CD1_1',count=count)
    if (count eq 0L) then $
      message, 'The image of '+objname+' is not wavelength calibrated.' else $

    arcsecaxis = make_arcsecaxis(header1d,pscale=pscale,refpix=refpix)
    if (long(arcsecaxis[0]) eq -1L) then begin
       specinfo = -1L
       return
    endif
    
    fluxcor = sxpar(header1d,'FLUXCOR',count=nfluxcor)
    if (nfluxcor eq 1L) then begin
       ytitle = 'f_{\lambda} (10^{-17} '+flam_units()+')'
       scale = 1E15
    endif else begin 
       ytitle = 'Counts'
       scale = 1.0
    endelse

; determine the reference column based on the reference wavelength

    minlam = min(wave)
    maxlam = max(wave)
    
    if n_elements(refwave) eq 0L then $
      refwave1 = (maxlam-minlam)/2.0+minlam else $
      refwave1 = (minlam > float(refwave)) < maxlam

    get_element, wave, refwave1, refcol

    if (not keyword_set(noplot)) then doplot = 1

; subtract one or more point sources (e.g., stars) from the 2D image  

    if keyword_set(starfit) then begin    

       starimage = i2dstarfit(image,sigmamap,header1d,_extra=extra,doplot=doplot)
       image = image - starimage
       
    endif
    
; trace    
    
    if keyword_set(tracespec) then begin 

       itracespec, image, header, refrow=refrow, refcol=refcol, $
         traceinfo=traceinfo, tracewindow=tracewindow, doplot=doplot, $
         _extra=extra

       center = traceinfo.trace
       refrow = center[midcol]

    endif else begin ; do not trace

       if n_elements(refrow) eq 0L then begin

          if n_elements(refwave) eq 0L then begin

             sprofile = total(image,1)
             frac = fix(0.05*nrows)
             mmax = max(sprofile[frac:nrows-frac]) ; initial guess
             cen = where(sprofile eq mmax)

          endif else begin
             
;            flag = 1L
             mmax = max(image[refcol,*],cen)

          endelse
             
          for i = 0L, 2L do $      
            fmax = trace_fweight(image[refcol,*],float(cen),0)

          center = replicate(cen,ncols)
          refrow = long(center[refcol])
       
       endif else center = replicate(float(refrow),ncols)

    endelse

; default aperture extraction parameters

    if (n_elements(aperture) eq 0L) then begin
       aperture = 0.1*nrows*pscale                   ; 10% aperture [arcsec]
       aperture = float(round(100.0*aperture)/100.0) ; round to 2 decimal places
    endif

    if ((n_elements(loaperture) ne 0L) and (n_elements(upaperture) eq 0L)) or $
      ((n_elements(loaperture) eq 0L) and (n_elements(upaperture) ne 0L)) then begin
       splog, 'Both LOAPERTURE and UPAPERTURE must be given.'
       delvarx, loaperture, upaperture
    endif

    if (n_elements(loaperture) eq 0L) and (n_elements(upaperture) eq 0L) then begin
       loaperture = aperture/2.0
       upaperture = aperture/2.0
    endif else aperture = loaperture + upaperture

    lower = (center-loaperture/pscale) > 0L         ; start aperture window [pixel]
    upper = (center+upaperture/pscale) < (nrows-1L) ; end aperture window [pixel]
       
; overplot the aperture limits on the trace

    if keyword_set(tracespec) and not keyword_set(noplot) then begin
       wset, tracewindow
       oplot, wave, lower, line=2, thick=2.0
       oplot, wave, upper, line=2, thick=2.0
    endif

; define the sky windows

    if n_elements(skylower) eq 0L then skylower = 0.0*nrows*pscale ; [arcsec]
    if n_elements(skyupper) eq 0L then skyupper = skylower

    nskyap = n_elements(skyaperture)
    if (nskyap eq 0L) then skyaperture = 0.20*nrows*pscale ; 20% aperture [arcsec]
    if (nskyap eq 0L) or (nskyap eq 1L) then skyaperture = replicate(skyaperture,2)

    skylower = round(10.0*skylower)/10.0
    skyupper = round(10.0*skyupper)/10.0
    skyaperture = round(10D0*skyaperture)/10.0

    if nskyap gt 2L then message, 'Unsupported number of sky apertures!'

    skyap = long(skyaperture/pscale)>1                       ; [pixel]
    skylow = replicate(long(skylower/pscale),ncols)          ; [pixel]
    skyup = replicate(long(nrows-1.0-skyupper/pscale),ncols) ; [pixel]

    nflank = n_elements(flanking) ; over-ride SKYLOW and SKYUP
    if nflank ne 0L then begin

       if nflank gt 2L then message, 'Unsupported number of dimensions in FLANKING!'
       
       if (nflank eq 1L) then flanking = replicate(flanking,2)
       flanking = float(flanking)

       skylow = long(center - (loaperture+flanking[0]+skyaperture[0])/pscale)
       skyup = long(center + (upaperture+flanking[1]+skyaperture[1])/pscale)

    endif else flanking = [0.0,0.0]
    
; pack the aperture extraction parameters into a structure
    
    apinfo = {aperture: float(aperture), loaperture: loaperture, $
      upaperture: upaperture, center: center, lower: lower, $
      upper: upper, skyap: skyap, skylow: long(skylow), $
      skyup: long(skyup), flanking: flanking}

; give a warning message if the number of sky rows is much less than
; the number of aperture rows (see Robertson 1986)
    
; ultimately make the extraction interactive at this point
    
; ----------------------------------------------------------------------
; spatial profile
; ----------------------------------------------------------------------
;   nsum = 10
;   sprofile = total(image[midcol-nsum/2:midcol+nsum/2,*],1) ; sum columns to increase the S/N

    if (n_elements(flag) eq 1L) or keyword_set(meanprofile) then begin
       sprofile = scale*total(image,1)/ncols
       refrowtext = 'Mean Profile/'+strn(refrow)
    endif else begin
       sprofile = scale*image[refcol,*]
       refrowtext = strn(refrow)+'/'+string(refwave1,format='(I0)');+' \AA'
    endelse

    if n_elements(fluxrange) eq 0L then begin
       fluxrange = [min(sprofile)-(max(sprofile)-min(sprofile))*0.1,$
         max(sprofile)+(max(sprofile)-min(sprofile))/2.0]
    endif else fluxrange = scale*fluxrange

    if not keyword_set(noplot) then begin

       airmass = sxpar(header1d,'AIRMASS',count=nairmass)
       if (nairmass eq 0L) then airmass = 'Unknown' else $
         airmass = strtrim(string(airmass,format='(F12.3)'),2)
       
       if (!d.window ne 0L) then im_window, 0, xratio=0.49, /square
       skyaptxt = strtrim(string(skyaperture[0],format='(F12.2)'),2)+', '+$
         strtrim(string(skyaperture[1],format='(F12.2)'),2)
       legtxt = ['Aperture '+strtrim(string(apinfo.aperture,format='(F12.2)'),2),$
         'Sky '+skyaptxt,refrowtext,'Airmass = '+strn(airmass,length=6)]

       if keyword_set(tracespec) then begin

; plot SPORFILE in the upper panel; plot the average spatial profile
; of the lower and upper 10% of the columns in the lower two panels;
; below we overplot the traced apertures
          
          nfrac = fix(0.01*ncols)
          sprofilelo = scale*total(image[0L:nfrac-1L,*],1)/nfrac
          sprofilehi = scale*total(image[ncols-nfrac:ncols-1L,*],1)/nfrac

          meanlo = fix(total(colaxis[0L:nfrac-1L])/nfrac)
          meanhi = fix(total(colaxis[ncols-nfrac:ncols-1L])/nfrac)
          meanrowlo = fix(total(apinfo.center[0L:nfrac-1L])/nfrac)
          meanrowhi = fix(total(apinfo.center[ncols-nfrac:ncols-1L])/nfrac)
          meanlowave = wave[meanlo]
          meanhiwave = wave[meanhi]

          legtxtlo = ['Lower Mean Profile/'+string(meanrowlo,format='(I0)'),$
            string(meanlo,format='(I0)')+'/'+string(meanlowave,format='(I0)')];+' \AA']
          legtxthi = ['Upper Mean Profile/'+string(meanrowhi,format='(I0)'),$
            string(meanhi,format='(I0)')+'/'+string(meanhiwave,format='(I0)')];+' \AA']

          djs_plot, rowaxis, sprofile, xsty=11, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.5, xtitle='', ytitle='', $
            yrange=fluxrange, xthick=2.0, ythick=2.0, $
            position=[0.15,0.5,0.95,0.9] ; xmargin=[8,3], ymargin=[4,3] ; position=[0.25,0.1,0.95,0.9]
          axis, xaxis=1, xthick=2.0, charthick=2.0, xsty=1, $
            charsize=1.5, xtitle='Position (arcsec)', xrange=minmax(interpol(arcsecaxis,rowaxis,!x.crange))
          xyouts, 0.08, 0.5, textoidl(ytitle), align=0.5, orientation=90, charsize=1.5, $
            charthick=2.0, /normal
          legend, objname, /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, legtxt, /right, /top, box=0, charsize=1.5, charthick=2.0
; sky windows
          if (not keyword_set(noskysub)) and (n_elements(skyimage) eq 0L) then begin
             ymin = !y.crange[0] & ymax = (!y.crange[1]-ymin)/1.75+ymin
             djs_plotlimitbox, [apinfo.skylow[midcol],apinfo.skylow[midcol]+apinfo.skyap[0]], $
               [ymin,ymax], line=2, thick=3.0, color=djs_icolor('blue') ; sky box
             djs_plotlimitbox, [apinfo.skyup[midcol],apinfo.skyup[midcol]-apinfo.skyap[1]], $
               [ymin,ymax], line=2, thick=3.0, color=djs_icolor('blue')
          endif

; mean aperture windows
          ymin = !y.crange[0] & ymax = (!y.crange[1]-ymin)/1.25+ymin
          if keyword_set(tracespec) then cindx = midcol else cindx = refcol
          djs_plotlimitbox, [apinfo.lower[cindx],apinfo.upper[cindx]], [ymin,ymax], $
            line=0, thick=3.0, color=djs_icolor('green') ; extraction box
          oplot, [apinfo.center[cindx],apinfo.center[cindx]], $ ; vertical line
            ymax*[0.95,1.05], line=0, thick=2.0, color=djs_icolor('green')
          
          djs_plot, rowaxis, sprofilelo, xsty=1, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.3, xtitle='Row Number', ytitle='', $
            yrange=fluxrange, xthick=2.0, ythick=2.0, position=[0.15,0.1,0.55,0.44], $
            /noerase
; lower aperture windows
          djs_plotlimitbox, [apinfo.lower[0L],apinfo.upper[0L]], [ymin,ymax], $
            line=2, thick=2.0, color=djs_icolor('green')
          legend, legtxtlo, /right, /top, box=0, charsize=1.3, charthick=2.0

          djs_plot, rowaxis, sprofilehi, xsty=1, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.3, xtitle='Row Number', ytitle='', $
            yrange=fluxrange, xthick=2.0, ythick=2.0, position=[0.55,0.1,0.95,0.44], $
            ytickname=replicate(' ',10), /noerase
; upper aperture windows
          djs_plotlimitbox, [apinfo.lower[ncols-1L],apinfo.upper[ncols-1L]], [ymin,ymax], $
            line=2, thick=2.0, color=djs_icolor('green')
          legend, legtxthi, /right, /top, box=0, charsize=1.3, charthick=2.0

       endif else begin
          
          djs_plot, rowaxis, sprofile, xsty=11, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.5, xtitle='Row Number', ytitle=ytitle, $
            yrange=fluxrange, xthick=2.0, ythick=2.0, xmargin=[8,3], ymargin=[4,3]
          axis, xaxis=1, xthick=2.0, charthick=2.0, xsty=1, $
            charsize=1.5, xtitle='Position (arcsec)', xrange=minmax(interpol(arcsecaxis,rowaxis,!x.crange))
          legend, objname, /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, legtxt, /right, /top, box=0, charsize=1.5, charthick=2.0
; sky windows
          if (not keyword_set(noskysub)) and (n_elements(skyimage) eq 0L) then begin
             ymin = !y.crange[0] & ymax = (!y.crange[1]-ymin)/1.75+ymin
             djs_plotlimitbox, [apinfo.skylow[midcol],apinfo.skylow[midcol]+apinfo.skyap[0]], $
               [ymin,ymax], line=2, thick=3.0, color=djs_icolor('blue') ; sky box
             djs_plotlimitbox, [apinfo.skyup[midcol],apinfo.skyup[midcol]-apinfo.skyap[1]], $
               [ymin,ymax], line=2, thick=3.0, color=djs_icolor('blue')
          endif
; aperture windows
          ymin = !y.crange[0] & ymax = (!y.crange[1]-ymin)/1.25+ymin
          if keyword_set(tracespec) then cindx = midcol else cindx = refcol
          djs_plotlimitbox, [apinfo.lower[cindx],apinfo.upper[cindx]], [ymin,ymax], $
            line=0, thick=3.0, color=djs_icolor('green') ; extraction box
          oplot, [apinfo.center[cindx],apinfo.center[cindx]], $ ; vertical line
            ymax*[0.95,1.05], line=0, thick=2.0, color=djs_icolor('green')

       endelse
          
    endif

; ----------------------------------------------------------------------
; sky subtraction
; ----------------------------------------------------------------------

    if (n_elements(skyimage) eq 0L) then begin

       iskysubtract, image, sigmamap, header, apinfo, skyinfo, noskysub=noskysub, _extra=extra
       imnosky = skyinfo.imnosky

    endif else begin

       iskysubtract, image, sigmamap, header, apinfo, skyinfo, /noskysub, _extra=extra
       imnosky = image
       skyinfo.skyfit = skyimage

    endelse

; ----------------------------------------------------------------------
; estimate the seeing
; ----------------------------------------------------------------------

    if keyword_set(seeing) then begin

       subimage = image[*,floor(min(apinfo.lower))>0L:ceil(max(apinfo.upper))<(nrows-1L)]
       
; fit a Gaussian function at NSEEING positions along columns, median
; averaging NMED rows to increase the S/N at to prevent spurious fits

       nmed = 5L
       nseeing = 10L
       seeing = fltarr(nseeing)

       colspace = 0.9*ncols/nseeing
       colpos = fix(findgen(nseeing)*colspace+colspace/2.0)

       for isee = 0L, nseeing-1L do begin

          sprofile = djs_median(subimage[(colpos[isee]-nmed)>0L:(colpos[isee]+nmed)<(ncols-1L),*],1)
          gaussaxis = findgen(n_elements(sprofile))
          
          gfit = mpfitpeak(gaussaxis,sprofile,a,nterm=4L,/gaussian,/positive)
          seeing[isee] = 2D0*sqrt(2.0*alog(2.0))*a[2]*pscale

;         plot, gaussaxis, sprofile, ps=4, xsty=3, ysty=3, thick=2.0, xthick=2.0, $
;           ythick=2.0, charsize=1.8, charthick=2.0, xtitle='Row', ytitle='Flux/Counts', $
;           title='Column '+string(colpos[isee],format='(I0)')
;         djs_oplot, gaussaxis, gfit, color='green', ps=10, thick=3.0
;         cc = get_kbrd(1)

       endfor

; update the header       

       sxaddpar, header1d, 'SEEING', float(median(seeing)), $
         ' median FWHM seeing [arcsec]', before='HISTORY'

    endif

; ----------------------------------------------------------------------
; spectrum extraction
; ----------------------------------------------------------------------

    if keyword_set(optimal) then begin

       iextract_optimal, imnosky, sigmamap, mask, apinfo, skyimage=skyimage, $ 
         spec1d=spec1d, sig1d=sig1d, mask1d=mask1d, sky1d=sky1d, _extra=extra

    endif else begin

; sum along the spatial dimension to extract the one dimensional
; spectrum
       
       spec1d = fltarr(ncols)
       sig1d = spec1d*0.0
       sky1d = spec1d*0.0
       mask1d = intarr(ncols)

; sum the spectrum using partial pixels (conserving flux)
       
;      npix = fix(aperture/pscale)+1L
       npix = fix(apinfo.upper)-fix(apinfo.lower)+1L

       up = floor(apinfo.upper)
       lo = ceil(apinfo.lower)

       for k = 0L, ncols-1L do begin

          pixels = lindgen(npix[k])+fix(apinfo.lower[k])

          weights = replicate(1.0D,npix[k])
          weights[0L] = lo[k]-apinfo.lower[k]
          weights[npix[k]-1L] = apinfo.upper[k]-up[k]
          if (k eq 0L) then outweights = weights else outweights = [outweights,weights]
          
          maskcut = reform(mask[k,pixels])

          spec1d[k] = total(weights*imnosky[k,pixels])
          sig1d[k] = sqrt(total(weights*sigmamap[k,pixels]^2.0))
          sky1d[k] = total(weights*skyinfo.skyfit[k,pixels])

          indx = iindx_badpix(maskcut,bit=bit)
          mask1d[k] = bit
          
       endfor

    endelse

; ----------------------------------------------------------------------    
; update the output header
; ----------------------------------------------------------------------    
    
    sxaddpar, header1d, 'NAXIS', 1L
    sxaddpar, header1d, 'NAXIS1', ncols
    sxdelpar, header1d, 'NAXIS2'

    sxaddpar, header1d, 'APERLO', apinfo.loaperture, format='(F12.2)', $
      ' lower extraction aperture radius [arcsec]', before='HISTORY'
    sxaddpar, header1d, 'APERUP', apinfo.upaperture, format='(F12.2)', $
      ' upper extraction aperture radius [arcsec]', before='HISTORY'
    sxaddpar, header1d, 'APERWID', apinfo.aperture, format='(F12.2)', $
      ' extraction aperture diameter [arcsec]', before='HISTORY'
    sxaddpar, header1d, 'APERCEN', round(apinfo.center[midcol]), $
      ' extraction aperture center [pixel]', before='HISTORY'
    sxaddpar, header1d, 'MEDSNR', median(spec1d/sig1d), $
      ' median signal-to-noise per pixel', before='HISTORY'
;   sxaddpar, header1d, 'EXTRACT', apinfo.aperture, format='(F12.2)', $
;     ' extraction aperture diameter [arcsec]', before='HISTORY'
;   sxaddpar, header1d, 'EXSHIFT', (refrow-refpix)*pscale, $
;     ' extraction center [arcsec]', before='HISTORY'
    
    sxaddhist, "'One-dimensional spectrum extracted "+im_today()+"'", header1d

; ----------------------------------------------------------------------    
; improve the wavelength solution using the sky lines
; ----------------------------------------------------------------------    

    if keyword_set(noskyshift) then begin
;   if keyword_set(noskyshift) or keyword_set(noskysub) then begin
       if not keyword_set(silent) then splog, 'No sky line shift applied.' 
    endif else begin
       skyshift = iskyshift(sky1d,header1d,_extra=extra)
       wave = wave + skyshift*dwave
    endelse

; ----------------------------------------------------------------------    
; repair sky-subtraction residuals
; ----------------------------------------------------------------------    

    if keyword_set(noskymask) or keyword_set(noskysub) then begin
       if not keyword_set(silent) then splog, 'Sky subtraction residuals not repaired.' 
    endif else begin
       spec1d = iskymask(spec1d,sig1d,wave,mask=mask1d,doplot=(keyword_set(noplot)-1))
    endelse

; ----------------------------------------------------------------------    
; subtract an input spectrum from the extracted spectrum (useful to
; remove stellar contamination or to pass a sky-spectrum)
; ----------------------------------------------------------------------    
    
    if n_elements(inputspec) ne 0L then begin

       if not keyword_set(silent) then splog, 'Subtracting an input spectrum.'

; need to check for wavelength overlap here!

;      splog, format='("Subtracting an input spectrum . . .",$)'
       combine1fiber, alog10(inputspec.wave), inputspec.spec, inputspec.invvar, $
         newloglam=alog10(wave), newflux=newflux
       spec1d = spec1d - newflux

    endif

; pack the extraction into a structure

    specinfo = {$
      spec:    spec1d,      $
      wave:    wave,        $
      sigspec: sig1d,       $
      sky:     sky1d,       $
      mask:    fix(mask1d), $
      header:  header1d,    $
      refrow:  refrow,      $
      refcol:  refcol,      $
      refwave: refwave1}
    specinfo = create_struct(create_struct(specinfo,apinfo),skyinfo)

    if not keyword_set(noplot) then begin

       if !d.window ne 2L then im_window, 2, xratio=0.49, /square
;      plot, specinfo.wave, specinfo.spec/sqrt(specinfo.vspec), xsty=3, ysty=3, ps=10, $
       djs_plot, specinfo.wave, scale*specinfo.spec, xsty=3, ysty=3, ps=10, $
         xtitle='Wavelength (\AA)', ytitle=ytitle, charthick=2.0, charsize=1.5, $
         xmargin=[10,3], yrange=minmax(scale*specinfo.spec)*[1.0,1.1]
;      oplot, specinfo.wave, specinfo.sigspec, line=2

       snr = median(specinfo.spec/specinfo.sigspec) ; median S/N ratio

       scanlen = sxpar(header1d,'SCANLEN',count=count_scanlen) 
       if (count_scanlen ne 0L) then begin
          legtext = [sxpar(header1d,'OBJECT'),strn(aperture,format='(F6.2)')+'" Aperture',$
            strn(scanlen,format='(F6.2)')+'" Scan','S/N = '+strtrim(string(snr,format='(F10.1)'),2)] 
       endif else begin
          legtext = [sxpar(header1d,'OBJECT'),strn(aperture,format='(F6.2)')+'" Aperture',$
            'S/N = '+strtrim(string(snr,format='(F10.1)'),2)]
       endelse

       legend, strtrim(legtext,2), /right, /top, box=0, charthick=2.0, charsize=1.3

    endif
    
return
end
