;+
; NAME:
;	IDOEXTRACT
;       
; PURPOSE:
;	Extract a one-dimensional spectrum.
;
; CALLING SEQUENCE:
;       idoextract, image, sigmamap, mask, header, specinfo, $
;          skyimage=, aperture=, loaperture=, upaperture=, $
;          inputspec=, refrow=, fluxrange=, /tracespec, $
;          /noskyshift, /noskymask, /noplot, /meanprofile, $
;          /seeing, /optimal, /starfit, /silent, _extra=extra
;
; INPUTS:
;	image       - 2D sky-subtracted spectrum
;	sigmamap    - corresponding 2D error map
;	mask        - corresponding 2D bad pixel mask
;	header      - corresponding FITS header
;
; OPTIONAL INPUTS:
;	skyimage    - 2D sky spectrum
;       wavemap     - 2D wavelength map
;	aperture    - extraction aperture diameter [arcsec]
;	loaperture  - lower extraction aperture radius [arcsec]
;	upaperture  - upper extraction aperture radius [arcsec]
;	inputspec   - data structure containing an input spectrum to
;                     subtract from the extracted 1D spectrum:
;	  wave           - 1D wavelength vector
;	  spec           - corresponding 1D spectrum
;         invvar         - corresponding inverse variance spectrum 
;       refrow      - extract the spectrum centered on REFROW.  (the
;                     default is the row corresponding to the
;                     flux-weighted center of the mean spatial
;                     profile) 
;	fluxrange   - vertical plot range for the spatial profile
;                     window [MIN,MAX]
;	minwave     - starting wavelength for the output 1D spectrum
;                     [Angstrom] 
;	maxwave     - ending wavelength for the output 1D spectrum
;                     [Angstrom] 
;	dwave       - wavelength dispersion for the output 1D spectrum
;                     [Angstrom/pixel]
;	extra       - keywords for ITRACESPEC, IEXTRACT_OPTIMAL, and 
;                     ISKYSHIFT
;	
; KEYWORD PARAMETERS:
;	tracespec   - trace the spectrum as a function of wavelength 
;	noskyshift  - do not use the night sky lines to improve the
;                     wavelength calibration
;	noskymask   - do not mask sky subtraction residuals
;	noplot      - do not generate extraction windows
;       meanprofile - set this to examine the mean spatial profile
;                     while still being able to set REFROW; this
;                     keyword is stronger than REFWAVE
;       seeing      - estimate the seeing by fitting a Gaussian
;                     function to the spatial PSF (useful for point
;                     sources) 
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
;	  sigspec  - one dimensional error spectrum
;	  sky      - corresponding sky spectrum
;	  header   - spectrum header
;	  apinfo   - aperture extraction parameters
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       IDOEXTRACT expects a ccd-processed and sky-subtracted
;       spectrum. 
;
;       If LOAPERTURE, UPAPERTURE, and APERTURE are specified then
;       APERTURE is overridden.
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
;	ISKYSHIFT, COMBINE1FIBER, READCOL, SPLOG, SXADDHIST,
;	ISKYMASK(), SXADDPAR, SXDELPAR, DJS_PLOT, IM_HMS2DEC(),
;	OBSERVATORY, DELVARX, IEXTRACT_OPTIMAL, IM_WINDOW, REM_DUP(),
;	DJS_MEDIAN(), MPFITPEAK(), IM_TODAY()  
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
;       jm05jun20uofa - routine exported from IEXTRACT to the new name
;                       IDOEXTRACT due to all the major changes,
;                       including: 
;          * Input 2D spectrum must be sky-subtracted and the
;              appropriate wavelength map must be given (as of v2.0
;              this is handled by ISKYSUBTRACT2D)
;          * IM_WINDOW is used to make resolut-independent windows 
;          * SKYAPERTURE, SKYLOWER, SKYUPPER, FLANKING, REFWAVE, and
;              HELIOCOR keywords removed
;          * WAVEMAP, MINWAVE, MAXWAVE, and DWAVE optional inputs
;              added 
;          * SEEING keyword added
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

pro idoextract, image, sigmamap, mask, header, specinfo, skyimage=skyimage, $
  wavemap=wavemap, aperture=aperture, loaperture=loaperture, upaperture=upaperture, $
  inputspec=inputspec, refrow=refrow, fluxrange=fluxrange, minwave=minwave1, $
  maxwave=maxwave1, dwave=dwave1, tracespec=tracespec, noskyshift=noskyshift, $
  noskymask=noskymask, noplot=noplot, meanprofile=meanprofile, seeing=seeing, $
  optimal=optimal, starfit=starfit, silent=silent, _extra=extra

    if (n_elements(image) eq 0L) or (n_elements(sigmamap) eq 0L) or $
      (n_elements(mask) eq 0L) or (n_elements(header) eq 0L) then begin
       print, 'Syntax - idoextract, image, sigmamap, mask, header, specinfo, $'
       print, '   skyimage=, aperture=, loaperture=, upaperture=, inputspec=, $'
       print, '   refrow=, fluxrange=, /tracespec, /noskyshift, /noskymask, $'
       print, '   /noplot, /meanprofile, /seeing, /optimal, /starfit, $'
       print, '   /silent, _extra=extra'
       return
    endif

; preliminaries    
    
    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    midrow = nrows/2L
    midcol = ncols/2L
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)
    refcol = midcol
    
    header1d = header
    objname = sxpar(header1d,'OBJECT',count=count_objname)
    if (count_objname eq 0L) then objname = ''

    arcsecaxis = make_arcsecaxis(header1d,pscale=pscale,refpix=refpix)
    if (long(arcsecaxis[0]) eq -1L) then begin
       specinfo = -1L
       return
    endif

    fluxcor = sxpar(header1d,'FLUXCOR',count=nfluxcor)
    if (nfluxcor eq 1L) then begin
       ytitle = 'f_{\lambda} (10^{-15} '+flam_units()+')'
       scale = 1E15
    endif else begin 
       ytitle = 'Counts'
       scale = 1.0
    endelse
    
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

          sprofile = total(image,1)
          frac = fix(0.05*nrows)
          mmax = max(sprofile[frac:nrows-frac]) ; initial guess
          cen = where(sprofile eq mmax)

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
       oplot, colaxis, lower, line=2, thick=2.0
       oplot, colaxis, upper, line=2, thick=2.0
    endif

; pack the aperture extraction parameters into a structure
    
    apinfo = {aperture: float(aperture), loaperture: loaperture, $
      upaperture: upaperture, center: center, lower: lower, $
      upper: upper}

; ultimately make the extraction interactive at this point
    
; ----------------------------------------------------------------------
; spatial profile
; ----------------------------------------------------------------------

    if keyword_set(meanprofile) then begin
       sprofile = scale*total(image,1)/ncols
       refrowtext = 'Mean Profile/'+string(refrow,format='(I0)')
    endif else begin
       sprofile = scale*image[refcol,*]
       refrowtext = string(refrow,format='(I0)')+'/'+string(refcol,format='(I0)')
    endelse

    if n_elements(fluxrange) eq 0L then begin
       fluxrange = [min(sprofile)-(max(sprofile)-min(sprofile))*0.1,$
         max(sprofile)+(max(sprofile)-min(sprofile))/2.0]
    endif else fluxrange = scale*fluxrange

    if (not keyword_set(noplot)) then begin

       airmass = sxpar(header1d,'AIRMASS',count=nairmass)
       if (nairmass eq 0L) then airmass = 'Unknown' else $
         airmass = strtrim(string(airmass,format='(F12.3)'),2)
       
       if (!d.window ne 0L) then im_window, 0, xratio=0.49, /square
       legtxt = 'Aperture = '+strtrim(string(apinfo.aperture,format='(F12.2)'),2)+'"'

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

          legtxtlo = 'Column '+string(meanlo,format='(I0)')+$
            ' Mean Profile/'+string(meanrowlo,format='(I0)')
          legtxthi = 'Column '+string(meanhi,format='(I0)')+$
            ' Mean Profile/'+string(meanrowhi,format='(I0)')

          djs_plot, rowaxis, sprofile, xsty=11, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.5, xtitle='', ytitle='', $
            yrange=fluxrange, xthick=2.0, ythick=2.0, $
            position=[0.17,0.5,0.95,0.9] ; xmargin=[8,3], ymargin=[4,3] ; position=[0.25,0.1,0.95,0.9]
          axis, xaxis=1, xthick=2.0, charthick=2.0, xsty=1, $
            charsize=1.5, xtitle='Position (arcsec)', xrange=minmax(interpol(arcsecaxis,rowaxis,!x.crange))
          xyouts, 0.05, 0.5, textoidl(ytitle), align=0.5, orientation=90, charsize=1.5, $
            charthick=2.0, /normal
          legend, objname, /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, legtxt, /right, /top, box=0, charsize=1.5, charthick=2.0

; mean aperture windows
          ymin = !y.crange[0] & ymax = (!y.crange[1]-ymin)/1.25+ymin
          if keyword_set(tracespec) then cindx = midcol else cindx = refcol
          djs_plotlimitbox, [apinfo.lower[cindx],apinfo.upper[cindx]], [ymin,ymax], $
            line=0, thick=3.0, color=djs_icolor('green') ; extraction box
          oplot, [apinfo.center[cindx],apinfo.center[cindx]], $ ; vertical line
            ymax*[0.95,1.05], line=0, thick=2.0, color=djs_icolor('green')
          
          djs_plot, rowaxis, sprofilelo, xsty=1, ysty=1, thick=2.0, ps=10, $
            charthick=2.0, charsize=1.3, xtitle='Row Number', ytitle='', $
            yrange=fluxrange, xthick=2.0, ythick=2.0, position=[0.17,0.1,0.55,0.44], $
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
            yrange=fluxrange, xthick=2.0, ythick=2.0, position=[0.17,0.1,0.95,0.9]
          axis, xaxis=1, xthick=2.0, charthick=2.0, xsty=1, $
            charsize=1.5, xtitle='Position (arcsec)', xrange=minmax(interpol(arcsecaxis,rowaxis,!x.crange))
          legend, objname, /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, legtxt, /right, /top, box=0, charsize=1.5, charthick=2.0
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

       message, 'Optimal extraction not yet supported.'
       specinfo = -1L
       return

;      iextract_optimal, imnosky, sigmamap, mask, apinfo, skyimage=skyimage, $ 
;        spec1d=spec1d, sig1d=sig1d, mask1d=mask1d, sky1d=sky1d, _extra=extra

    endif else begin

; extraction aperture weights
       
       apweights = image*0D0

       npix = fix(apinfo.upper)-fix(apinfo.lower)+1L
       up = floor(apinfo.upper)
       lo = ceil(apinfo.lower)

       for k = 0L, ncols-1L do begin

          appixels = lindgen(npix[k])+fix(apinfo.lower[k])

          apweights[k,appixels] = 1.0
          apweights[k,appixels[0]] = lo[k]-apinfo.lower[k]
          apweights[k,appixels[npix[k]-1L]] = apinfo.upper[k]-up[k]
          
       endfor

; determine the default output wavelength vector       

       if (n_elements(dwave1) eq 0L) then $
         dwave = fix(100.0*(max(wavemap[*,nrows/2L]) - $
           min(wavemap[*,nrows/2L]))/(ncols-1))/100D else $
         dwave = dwave1
       dwave = arm_double(dwave)

       if (n_elements(minwave1) eq 0L) then $
         minwave = float(ceil(max(wavemap[0,*]))) else $
         minwave = minwave1 > ceil(max(wavemap[0,*]))

       if (n_elements(maxwave1) eq 0L) then maxwave1 = minwave + dwave*(ncols-1)
       maxwave = maxwave1 < float(floor(min(wavemap[ncols-1,*])))

       npix1d = (maxwave-minwave)/dwave + 1L
       wave = minwave + findgen(npix1d)*dwave

       spec1d = wave*0.0
       sig1d = spec1d*0.0
       sky1d = spec1d*0.0
       mask1d = bytarr(npix1d) ; not supported

; wavelength weights

       thesepixels = where(apweights gt 0.0,nthesepixels)
       if (nthesepixels eq 0L) then message, 'This is not good.'

       srted = sort(wavemap[thesepixels])

       twodpixels = im_convert_index(thesepixels[srted],apweights)

       bigflux      = image[thesepixels[srted]]
       bigsky       = skyimage[thesepixels[srted]]
       bigferr      = sigmamap[thesepixels[srted]]
       bigapweights = apweights[thesepixels[srted]]

       bigwave          = wavemap[thesepixels[srted]]
       dwavematrix      = wavemap-shift(wavemap,1)
       dwavematrix[0,*] = dwavematrix[1,*]
       bigdwave         = dwavematrix[thesepixels[srted]]

;      plot, bigwave, bigflux, ps=4, xsty=3, ysty=3, xrange=[4650,4750] ; [5320.0,5380] ; [4150,4200]
;      plot, bigwave, bigferr, ps=4, xsty=3, ysty=3, xrange=[4000,4500]

       for k = 0L, npix1d-1L do begin

; the left side (low-wavelength end) of UPPIXELS is inside the desired
; output pixel, but the right side (high-wavelength end) is beyond
; WAVE[K] + DWAVE; LOPIXELS enter *into* the output pixel from the
; "left" (low wavelength side); however these "where" statements will
; double-count pixels that are smaller than the output dispersion, and
; that are roughly centered on the output wavelength; so flag and
; remove duplicates; also, we have to account for pixels that are
; roughly centered on the output pixel, but that is larger than the
; output pixel; compute the pixel weights
          
          uppixels = where((bigwave ge wave[k]) and (bigwave le wave[k]+dwave),nuppixels)
          if (nuppixels ne 0L) then begin
             upweights = ((wave[k] + dwave - bigwave[uppixels]) / bigdwave[uppixels]) < 1.0
             if (n_elements(pixels) eq 0L) then begin
                pixels = uppixels
                weights = upweights
             endif else begin
                pixels = [pixels,uppixels]
                weights = [weights,upweights]
             endelse
          endif

          lopixels = where((bigwave+bigdwave ge wave[k]) and (bigwave+bigdwave le wave[k]+dwave),nlopixels)
          if (nlopixels ne 0L) then begin
             loweights = ((bigwave[lopixels]+bigdwave[lopixels] - wave[k]) / bigdwave[lopixels]) < 1.0
             if (n_elements(pixels) eq 0L) then begin
                pixels = lopixels
                weights = loweights
             endif else begin
                pixels = [pixels,lopixels]
                weights = [weights,loweights]
             endelse
          endif

          morepixels = where((bigwave le wave[k]) and (bigwave+bigdwave ge wave[k]+dwave),nmorepixels)
          if (nmorepixels ne 0L) then begin
             moreweights = (bigdwave[morepixels]/dwave) < 1.0
             if (n_elements(pixels) eq 0L) then begin
                pixels = morepixels
                weights = moreweights
             endif else begin
                pixels = [pixels,morepixels]
                weights = [weights,moreweights]
             endelse
          endif
          
          uniqindx = rem_dup(pixels)
          pixels = pixels[uniqindx]
          weights = weights[uniqindx]

          final_weights = weights*bigapweights[pixels]

; construct the output spectrum

          spec1d[k] = total(final_weights*bigflux[pixels])
;         sig1d[k] = sqrt(total((bigferr[pixels])^2))
;         sig1d[k] = sqrt(total((final_weights*bigferr[pixels])^2))
          sig1d[k] = sqrt(total(final_weights*(bigferr[pixels])^2))
          sky1d[k] = total(final_weights*bigsky[pixels])

; delete variables!

          delvarx, pixels, weights
          
; debugging code for a.0108 in the 98mar data - jm05jun22uofa          
;         
;         if (k eq 206L) or (k eq 300L) then begin
;;        if (k eq 628L) or (k eq 632L) then begin
;;        if (k eq 206L) or (k eq 200L) then begin
;
;            plotsym, 0, 1.5, /fill
;            djs_oplot, bigwave[lopixels], bigferr[lopixels], ps=8, color='red'
;            djs_oplot, bigwave[lopixels]+bigdwave[lopixels], bigferr[lopixels], ps=8, color='red'
;            noweight = where(loweights eq 0.0,nno)
;            if (nno ne 0L) then begin
;               tvcircle, 1.0, bigwave[lopixels[noweight]], bigflux[lopixels[noweight]], /data
;               tvcircle, 1.0, bigwave[lopixels[noweight]]+bigdwave[lopixels[noweight]], $
;                 bigflux[lopixels[noweight]], /data
;            endif
;
;            plotsym, 8, 1.5;, /fill
;            djs_oplot, bigwave[uppixels], bigferr[uppixels], ps=8, color='green'
;            djs_oplot, bigwave[uppixels]+bigdwave[uppixels], bigferr[uppixels], ps=8, color='green'
;            noweight = where(upweights eq 0.0,nno)
;            if (nno ne 0L) then begin
;               tvcircle, 1.0, bigwave[uppixels[noweight]], bigferr[uppixels[noweight]], /data
;               tvcircle, 1.0, bigwave[uppixels[noweight]]+bigdwave[uppixels[noweight]], $
;                 bigferr[uppixels[noweight]], /data
;            endif
;            
;            djs_oplot, wave[k]*[1,1], !y.crange, line=2, thick=2
;            djs_oplot, (wave[k]+dwave)*[1,1], !y.crange, line=2, thick=2
;            
;            niceprint, bigferr[lopixels], loweights, bigferr[lopixels]*loweights, $
;              bigwave[lopixels], bigwave[lopixels]+bigdwave[lopixels], $
;              reform(twodpixels[0,lopixels]), reform(twodpixels[1,lopixels])
;            print
;
;            niceprint, bigferr[uppixels], upweights, bigferr[uppixels]*upweights, $
;              bigwave[uppixels], bigwave[uppixels]+bigdwave[uppixels], $
;              reform(twodpixels[0,uppixels]), reform(twodpixels[1,uppixels])
;
;            print, wave[k], wave[k]+dwave, spec1d[k], sig1d[k], spec1d[k]/sig1d[k], $
;              sqrt(total((loweights*bigferr[lopixels])^2)), $
;              sqrt(total((upweights*bigferr[uppixels])^2)), $
;              n_elements(lopixels), n_elements(uppixels)
;
;            stop
;
;         endif
         
; debugging code for a.0121 in the 98mar data - jm05jun22uofa          
;         
;         if (k eq 385L) or (k eq 399L) then begin
;;        if (k eq 628L) or (k eq 632L) then begin
;;        if (k eq 206L) or (k eq 200L) then begin
;
;            plotsym, 0, 1.5, /fill
;            djs_oplot, bigwave[lopixels], bigflux[lopixels], ps=8, color='red'
;            djs_oplot, bigwave[lopixels]+bigdwave[lopixels], bigflux[lopixels], ps=8, color='red'
;            noweight = where(loweights eq 0.0,nno)
;            if (nno ne 0L) then begin
;               tvcircle, 1.0, bigwave[lopixels[noweight]], bigflux[lopixels[noweight]], /data
;               tvcircle, 1.0, bigwave[lopixels[noweight]]+bigdwave[lopixels[noweight]], $
;                 bigflux[lopixels[noweight]], /data
;            endif
;
;            plotsym, 8, 1.5;, /fill
;            djs_oplot, bigwave[uppixels], bigflux[uppixels], ps=8, color='green'
;            djs_oplot, bigwave[uppixels]+bigdwave[uppixels], bigflux[uppixels], ps=8, color='green'
;            noweight = where(upweights eq 0.0,nno)
;            if (nno ne 0L) then begin
;               tvcircle, 1.0, bigwave[uppixels[noweight]], bigflux[uppixels[noweight]], /data
;               tvcircle, 1.0, bigwave[uppixels[noweight]]+bigdwave[uppixels[noweight]], $
;                 bigflux[uppixels[noweight]], /data
;            endif
;            
;            djs_oplot, wave[k]*[1,1], !y.crange, line=2, thick=2
;            djs_oplot, (wave[k]+dwave)*[1,1], !y.crange, line=2, thick=2
;            
;            niceprint, bigflux[lopixels], loweights, bigflux[lopixels]*loweights, $
;              bigwave[lopixels], bigwave[lopixels]+bigdwave[lopixels], $
;              reform(twodpixels[0,lopixels]), reform(twodpixels[1,lopixels])
;            print
;
;            niceprint, bigflux[uppixels], upweights, bigflux[uppixels]*upweights, $
;              bigwave[uppixels], bigwave[uppixels]+bigdwave[uppixels], $
;              reform(twodpixels[0,uppixels]), reform(twodpixels[1,uppixels])
;
;            print, wave[k], wave[k]+dwave, total(spec1d[k]), total(loweights*bigflux[lopixels]), $
;              total(upweights*bigflux[uppixels]), n_elements(lopixels), n_elements(uppixels)
;
;            stop
;
;         endif

       endfor
       
    endelse 

; ----------------------------------------------------------------------    
; update the output header
; ----------------------------------------------------------------------    

    sxaddpar, header1d, 'NAXIS', 1L
    sxaddpar, header1d, 'NAXIS1', npix1d
    sxdelpar, header1d, 'NAXIS2'

    sxaddpar, header1d, 'CRVAL1', float(minwave), ' wavelength at CRPIX1', before='HISTORY'
    sxaddpar, header1d, 'CRPIX1', float(1.0), ' reference pixel number', before='HISTORY'
    sxaddpar, header1d, 'CD1_1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
    sxaddpar, header1d, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
    sxaddpar, header1d, 'CTYPE1', 'LINEAR', ' projection type', before='HISTORY'
;   sxaddpar, header1d, 'DC-FLAG', 0, ' log-linear flag'

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
       if (not keyword_set(silent)) then splog, 'No sky line shift applied.' 
    endif else begin
       skyshift = iskyshift(sky1d,header1d,_extra=extra)
       wave = wave + skyshift*dwave
    endelse

; ----------------------------------------------------------------------    
; repair sky-subtraction residuals
; ----------------------------------------------------------------------    

    if keyword_set(noskymask) then begin
       if (not keyword_set(silent)) then splog, 'Sky subtraction residuals not repaired.' 
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

; pack the extraction information into a structure

    specinfo = {$
      spec:    spec1d,      $
      wave:    wave,        $
      sigspec: sig1d,       $
      sky:     sky1d,       $
      mask:    fix(mask1d), $
      header:  header1d,    $
      refrow:  refrow,      $
      refcol:  refcol}
    specinfo = create_struct(specinfo, apinfo)
      
    if (not keyword_set(noplot)) then begin

       if !d.window ne 2L then im_window, 2, xratio=0.49, /square
;      plot, specinfo.wave, specinfo.spec/sqrt(specinfo.vspec), xsty=3, ysty=3, ps=10, $
       djs_plot, specinfo.wave, scale*specinfo.spec, xsty=3, ysty=3, ps=10, $
         xtitle='Wavelength (\AA)', ytitle=ytitle, charthick=2.0, charsize=1.5, $
         xmargin=[10,3], yrange=minmax(scale*specinfo.spec)*[1.0,1.1]
;      oplot, specinfo.wave, specinfo.sigspec, line=2

       snr = median(specinfo.spec/specinfo.sigspec) ; median S/N ratio

       scanlen = sxpar(header1d,'SCANLEN',count=count_scanlen) 
       if (scanlen gt 0.0) then begin
          legtext = [sxpar(header1d,'OBJECT'),strtrim(string(aperture,format='(F12.2)'),2)+'" Aperture',$
            strtrim(string(scanlen,format='(F12.2)'),2)+'" Scan','S/N = '+strtrim(string(snr,format='(F10.1)'),2)]
       endif else begin
          legtext = [sxpar(header1d,'OBJECT'),strtrim(string(aperture,format='(F12.2)'),2)+'" Aperture',$
            'S/N = '+strtrim(string(snr,format='(F10.1)'),2)]
       endelse

       legend, strtrim(legtext,2), /right, /top, box=0, charthick=2.0, charsize=1.3

    endif
    
return
end
