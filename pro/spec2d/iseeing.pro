;+
; NAME:
;	ISEEING
;
; PURPOSE:
;       Read a list of two dimensional standard star spectra and
;       compute the seeing.
;
; CALLING SEQUENCE:
;       iseeing, specfile, seeing=, energy=, datapath=, $
;          title=, psname=, refwave=, slit_width=, /silent, $
;          /postscript, /doplot
;
; INPUTS:
;	specfile - column-format text file input list of standard
;                  stars, one per line, or a SPECLIST (identified by
;                  string searching for *.FITS) in SPECFILE
;
; OPTIONAL INPUTS:
;	datapath - path to SPECFILE
;	title    - plot title
;       psname   - postscript file name (default 'iseeing.ps')
;       refwave  - measure the seeing at these wavelengths (can be a
;                  scalar or a vector) (see COMMENTS)
;       slit_width - compute the amount of energy enclosed within
;                    SLIT_WIDTH based on the fit to the Gaussian
;                    profile (default 4.5")
;
; KEYWORD PARAMETERS:
;       silent     - do not print the seeing measurements to STDOUT
;	postscript - generate postscript output
;	doplot     - generate a seeing plot
;
; OUTPUTS:
;       seeing     - seeing FWHM at each REFWAVE (arcsec)
;       energy     - encircled energy within SLIT_WIDTH
;
; OPTIONAL OUTPUTS:
;	Postscript output of the seeing as a function of time if 
;	POSTSCRIPT is set.
;
; COMMENTS:
;       The seeing and enclosed energy are measured at three equally
;       spaced wavelengths based on the input (standard star) spectra
;       unless these wavelengths are specified via REFWAVE.  
;
; EXAMPLE:
;
; INTERNAL ROUTINES:
;       IFIT_SEEING(), SEEING_PLOT
;
; PROCEDURES USED:
;	CWD(), RD2DSPEC(), GET_ELEMENT, MPFITPEAK(), MAKE_WAVE(),
;	IFORAGE(), DJS_PLOT, DJS_OPLOT, LEGEND, DFPSPLOT, DFPSCLOSE,
;	SPLOG, READCOL, ICLEANUP, NICEPRINT
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 October 8, U of A - based on
;          SEEING_HISTOGRAM, an older code 
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr24uofa - calculate encircled energy; if only one star
;                       is passed then plot the Gaussian fit
;       jm03jul15uofa - 2D spectra do not have to be the same size 
;
; Copyright (C) 2002-2003, John Moustakas
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

function ifit_seeing, image, wave, pscale, refwave=refwave, $
  slit_width=slit_width, energy=energy, profplot=profplot

    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)
    
    if n_elements(refwave) eq 0L then refwave = [4000.0,5500.0,6500.0]
    nref = n_elements(refwave)

    box = 15.0                  ; [arcsec]
    sigma = fltarr(nref)
    seeing = fltarr(nref)
    energy = fltarr(nref)

    get_element, wave, refwave, refcol

    if keyword_set(profplot) then begin

       djs_plot, [0], [0], /nodata, xrange=[-20,20*(nref+1)], yrange=[0,1.2], $
         xsty=3, ysty=3, xthick=2.0, ythick=2.0, thick=2.0, charsize=2.0, $
         charthick=2.0, xtitle='Arbitrary Axis', ytitle='Normalized Flux'
       
    endif

    colors = ['blue','green','red']
    
    for i = 0L, nref-1L do begin
       
       sprofile = image[refcol[i],*]
       smax = max(sprofile,refrow)
       arcaxis = (rowaxis-refrow)*pscale ; [arcsec]
       
; fit a Gaussian to the spatial profile

       local = where((arcaxis gt -box) and (arcaxis lt box))
       profit = mpfitpeak(arcaxis[local],sprofile[local],a,nterm=4,/gaussian,/positive,yfit=yfit)

       sigma[i] = a[2]                          ; [arcsec]
       seeing[i] = 2.0*sqrt(2.0*alog(2.0))*a[2] ; FWHM [arcsec]
       energy[i] = errorf(slit_width/2.0/a[2]/sqrt(2.0))
       
       if keyword_set(profplot) then begin

          shift = -a[1]+20*(i+1)
          norm = max(sprofile-a[3])
          djs_oplot, arcaxis+shift, (sprofile-a[3])/norm, ps=10, thick=3.0
          djs_oplot, arcaxis[local]+shift, profit/norm, ps=10, color=colors[i<2L], thick=2.0
          
       endif

    endfor

    if keyword_set(profplot) then begin

       legstr = ''
       for j = 0L, nref-1L do $
         legstr = [legstr,textoidl([$
         'Wavelength = '+string(refwave[j],format='(G0.0)')+' \AA',$
         'Seeing (1\sigma) = '+string(sigma[j],format='(F5.3)')+'"',$
         'Seeing (FWHM) = '+string(seeing[j],format='(F5.3)')+'"',$
         'Energy ('+string(slit_width,format='(G0.0)')+'" Slit) = '+string(100*energy[j],format='(F6.2)')+'%',''])]
       legstr = legstr[1:n_elements(legstr)-1L]
       
       legend, legstr, /left, /top, box=0, charsize=1.3, charthick=2.0

    endif
    
return, seeing
end

pro seeing_plot, julday, seeing, refwave, title=title, postscript=postscript, psname=psname

    if keyword_set(postscript) then begin
       !p.font = 0
       dfpsplot, psname, /square, /color, /isolatin1
    endif else window, 0, xs=500, ys=500

    dummy = label_date(date_format=['%N-%D'])
    djs_plot, julday, seeing[*,0], ps=-2, xsty=3, ysty=3, yrange=minmax(seeing)*[0.95,1.1], $
      charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, line=0, xtitle='UT Date', $
      ytitle='FWHM Seeing (arcsec)', color=fsc_color('blue',100), thick=2.0, title=title, $
      xtickformat='LABEL_DATE'
    djs_oplot, julday, seeing[*,1], ps=-4, line=2, color=fsc_color('green',101), thick=2.0
    djs_oplot, julday, seeing[*,2], ps=-6, line=3, color=fsc_color('red',102), thick=2.0

    legend, string(refwave,format='(I5)')+' \AA ', psym=[-2,-4,-6], $
      linestyle=[0,2,3], thick=2.0, /right, /top, $
      box=0, charthick=2.0, charsize=1.5, color=[fsc_color('blue',100),$
      fsc_color('green',101),fsc_color('red',102)]

    if keyword_set(postscript) then begin
       dfpsclose
       !p.font = -1
    endif
    
return
end    

pro iseeing, specfile, seeing=seeing, energy=energy, datapath=datapath, $
  title=title, psname=psname, refwave=refwave, slit_width=slit_width, $
  silent=silent, postscript=postscript, doplot=doplot

    if n_elements(specfile) eq 0L then begin
       print, 'Syntax - iseeing, specfile, seeing=, energy=, datapath=, $'
       print, '   title=, psname=, refwave=, slit_width=, /silent, $'
       print, '   /postscript, /doplot'
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(psname) eq 0L then psname = 'iseeing.ps'
    if n_elements(doplot) eq 0L then doplot = 1L

    if n_elements(slit_width) eq 0L then slit_width=4.5 ; [arcsec]
    
    pushd, datapath

    if (strmatch(specfile,'*fits*'))[0] eq 1B then speclist = specfile else begin
    
       if file_test(specfile,/regular) eq 0L then begin
          splog, 'File '+datapath+specfile+' not found.'
          return
       endif
    
       readcol, specfile, speclist, format='A', comment='#', /silent

    endelse
    
; verify that the data have been wavelength-calibrated, then read them 

    forage = iforage(speclist)
    if size(forage,/type) ne 8L then begin
       splog, 'There was an error in IFORAGE().'
       return
    endif
    nspec = n_elements(speclist)

    now = where(forage.crval1 eq float(0),nnow)
    if nnow ne 0L then begin
       splog, 'The following spectra have not been wavelength calibrated:'
       niceprint, forage[now].file
       return
    endif

    julday = forage.jd

; reference wavelengths [Angstrom].  select NSAMP wavelengths equally
; spaced over the average spectral range, unless REFWAVE is specified

    if n_elements(refwave) eq 0L then begin

       nref = 3L
       refwave = fltarr(nref)
       rcd1_1 = djs_mean(forage.cd1_1)
       rnaxis1 = djs_mean(forage.naxis1)
       rcrval1 = djs_mean(forage.crval1)
       
       refwave[0] = rnaxis1/(2.5*nref)*rcd1_1+rcrval1           ; first wavelength
       refwave[1] = rnaxis1/2.0*rcd1_1+rcrval1                  ; middle of the wavelength range
       refwave[2] = rnaxis1*rcd1_1+rcrval1-(refwave[0]-rcrval1) ; ending wavelength

       refwave = float(round(refwave))
       
    endif else nref = n_elements(refwave)
    
    seeing = fltarr(nspec,nref) ; [arcsec]
    energy = fltarr(nspec,nref) ; fraction

    if keyword_set(doplot) and (nspec eq 1L) then profplot = 1L    
    
    for j = 0L, nspec-1L do begin
    
       cube = rd2dspec(speclist[j],datapath=datapath,silent=silent)
       image = cube.image
       header = *cube.header
       wave = make_wave(header)
       
       seeing[j,*] = ifit_seeing(image,wave,forage[j].cd2_2,refwave=refwave,$
         slit_width=slit_width,profplot=profplot,energy=energy1)

       energy[j,*] = energy1

       icleanup, cube
       
    endfor

    if keyword_set(doplot) and (nspec gt 1L) then begin
    
; sort by Julian date

       srt = sort(julday)
       julday = julday[srt]
       seeing = seeing[srt,*]
       seeing_plot, julday, seeing, refwave, title=title

    endif
       
    if keyword_set(postscript) and (nspec gt 1L) then seeing_plot, julday, $
      seeing, refwave, psname=psname, title=title, /postscript
    
    popd
    
return
end
