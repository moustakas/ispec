;+
; NAME:
;	ISENSFUNC()
;
; PURPOSE:
;	Use standard star observations to derive a sensitivity
;	function.  
;
; CALLING SEQUENCE:
;       info = isensfunc(stdlist,datapath=,sensname=,extfile=,$
;          bsorder_sens=,nord_sens=,zptshift=,senstitle=,$
;          grey=,slit_width=,/doplot,/write,_extra=extra)
;
; INPUTS:
;	stdlist   - file list of calibrated standard stars
;
; OPTIONAL INPUTS:
;	datapath     - file path to the standard star observations
;	sensname     - root name of the sensitivity curve
;       extfile      - name of the extinction file to use (default 
;                      'KPNOEXTINCT.DAT')
;	bsorder_sens - order of the b-spline fit to the sensitivity 
;                      function (default 50)
;	nord_sens    - order of the spline at every b-spline
;                      breakpoint (default 3) 
;       zptshift     - scalar sensitivity function zero point shift in
;                      magnitudes; shifts the sensitivity function by
;                      a constant amount (relative to the mean
;                      sensitivity at the mean wavelength of the
;                      data); useful for absolute spectrophotometry 
;       senstitle    - plot title (default none)
;	grey         - apply a grey shift either to the most sensitive
;                      star (GREY=1) or the mean sensitivity (GREY=2)
;                      at the average wavelength of the standards
;                      (default GREY=1)
;       slit_width   - slit width used to compute the encircled energy
;                      [arcsec] 
;	extra        - extra keywords for BSPLINE_ITERFIT() and 
;                      ISTARSEARCH() 
;
; KEYWORD PARAMETERS:
;	doplot    - generate diagnostic and quality assurance plots 
;	write     - generate postscript output of all the diagnostic
;                   plots, write the final sensitivity function as a
;                   FITS file, and create a log file
;
; OUTPUTS:
;       info - output data structure with the following fields
;          specname  - input FITS spectrum name
;          starname  - star name from the header
;          tablename - name of the standard in the STARINFO table 
;          date      - date of observation
;          airmass   - airmass of observation
;          posangle  - slit position angle
;          parangle  - parallactic angle
;          greyshift - either maximum (if GREY=1) or mean (if GREY=2)
;                      grey shift
;          greywave  - wavelength of the grey shift
;          senszero  - median sensitivity at the mean wavelength of
;                      the sensitivity function in units of 
;                      [counts/s/A]/[erg/s/cm2/A]
;          psym      - plot symbol
;          color     - star color
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The plot legend shows the five most deviant stars with the
;       date of observation, the airmass, and the grey shift.
;
;       The error spectra and bad pixel masks of the standards are not
;       incorporated into generating the sensitivity function.
;
;       Assuming that the shifts in the sensitivity functions are due
;       to clouds (and not slit loss), then the absolute
;       spectrophotometric error is given by the scatter of the
;       standards with respect to the sensitivity function fit with no
;       grey shift.  If you are not concerned with absolute
;       spectrophotometry then shift each sensitivity function to the
;       maximum sensitivity (GREY=1).  Alternatively, grey shifting to
;       the mean sensitivity (GREY=2) and flux-calibrating gives the
;       relative spectrophotometric error as the scatter about that
;       mean curve.
;
;       If ZPTSHIFT is non-zero then the SENSZERO field in INFO (see
;       OUTPUTS) is *not* modified, but the output sensitivity
;       functions and plots have the zero point shift applied.
;
; EXAMPLE:
;
; INTERNAL ROUTINES:
;       CLEANSENS, SENSFUNCPLOT, SENSFUNCFIT, WRITE_SENSFUNCFIT 
;
; PROCEDURES USED:
;	BSPLINE_ITERFIT(), BSPLINE_VALU(), READCOL, MAKE_WAVE(),
;	SXPAR(), MKHDR, CWD(), DFPSPLOT, DFPSCLOSE, DJS_ITERSTAT,
;	DJS_PLOT, DJS_OPLOT, IM_SYMBOLS, ANGSTROM(), IM_LEGEND,
;	MRDFITS(), JD2FITSDATE(), DJS_MEAN(), FSC_COLOR(),
;	GET_ELEMENT, SXADDPAR, SXADDHIST, MWRFITS, SPLOG, MINMAX(),
;	STRUCT_PRINT, STRUCT_TRIMTAGS(), DJS_MEDIAN(), REPSTR(),
;	FLAM_UNITS(), TELLURIC_MASK(), IM_WINDOW, IM_BEST_SYMBOLS(),
;	ISTARSEARCH(), DJSIG()
;
; DATA FILES:
;	${ISPEC_DIR}/etc/kpnoextinct.dat
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August/October, U of A - written 
;       jm02jun12uofa - verify that the headers have all the required
;                       fields 
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr21uofa - major updates to grey-shifting capabilities;
;                       expanded documentation
;       jm03apr22uofa - more major additions.  added WRITE_SENSFUCFIT
;                       routine.  calculate the error spectrum in the
;                       sensitivity function.  make grey-shifting
;                       mandatory and write out separate sensitivity
;                       functions 
;       jm03apr23uofa - several bug fixes in STARSEARCH()
;       jm03dec30uofa - added TELLURIC keyword
;       A.R.Marble, 2004 January 6, U of A - plot color bugs fixed
;       jm04jul11uofa - use a black plot background for contrast 
;       jm05jan14uofa - reject bandpasses from the fitting that are
;                       affected by telluric absorption
;       jm05jan19uofa - when computing the grey-shift wavelength,
;                       exclude masked bandpasses
;       jm05jun21uofa - major rewrite:
;          * obsolete STDPATH keyword removed
;          * obsolete HEADCHECK() and STARSEARCH() functions removed;
;               new ISTARSEARCH() routine now used to locate standards  
;          * SEARCHRAD optional input now superseeded by ISTARSEARCH()  
;          * the call to ISEEING (and any subsequent references to the
;               seeing or encircled energy) has been removed since no
;               wavelength-calibrated 2D images are created as of
;               iSPEC2d v2.0 
;          * the TELLURIC keyword and all associated code were removed 
;          * use IM_WINDOW to spawn monitor resolution independent
;               windows 
;          * SLIT_WIDTH optional input added to compute the encircled
;               energy if an estimate of the seeing appears in the
;               header 
;          * drastic improvement in the computation of the sensitivity
;               function error; the error in the mean (fit minus all
;               the standards) across the wavelength range is computed
;               and a low-order b-spline is fitted; this is done
;               separately for the observed and grey-shifted
;               sensitivity functions; set DOPLOT=1 for a QA plot
;          * improved log file output
;          * default b-spline order increased to 50 due to the higher
;               sampling of the flux standards
;          * better contrast for the sensitivity function fit (thanks
;               to Andy Marble) and thicker plot symbols
;          * improved selection of break points and better accounting
;               of gaps and exclusion of telluric regions in the
;               standard star spectra 
;          * improved memory management and speed
;
; Copyright (C) 2001-2005, John Moustakas
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

pro cleansens, sensinfo
; clean up memory

    nsens = n_elements(sensinfo)
    if (nsens eq 0L) then return

    for k = 0L, nsens-1L do begin

       if ptr_valid(sensinfo[k].bandwave) then ptr_free, sensinfo[k].bandwave
       if ptr_valid(sensinfo[k].sens) then ptr_free, sensinfo[k].sens
       if ptr_valid(sensinfo[k].bandmask) then ptr_free, sensinfo[k].bandmask
       if ptr_valid(sensinfo[k].binsize) then ptr_free, sensinfo[k].binsize

    endfor

return
end    
    
pro sensfuncfit, sensinfo, bsorder_sens, nord_sens, sensres, $
  fitsens, fitwave, fitmask, coarsefit, senswave, sensfit, $
  sensfiterr, sset, doplot=doplot, greyshifted=greyshifted, $
  _extra=extra
; fit a b-spline to the sensitivity function in logarithm space
; because it is smoother

    nstds = n_elements(sensinfo)
    
    fitwave = fltarr(total(sensinfo.nbins))
    fitsens = fitwave*0.0
    fitmask = byte(fitwave*0.0)
    fitbinsize = fitwave*0.0
    
    for k = 0L, nstds-1L do begin

       if (k eq 0L) then pbins = 0L else pbins = long(total(sensinfo[0:k-1].nbins))
       nbins = sensinfo[k].nbins

       fitwave[pbins:nbins+pbins-1] = *sensinfo[k].bandwave
       fitsens[pbins:nbins+pbins-1] = *sensinfo[k].sens
       fitmask[pbins:nbins+pbins-1] = *sensinfo[k].bandmask
       fitbinsize[pbins:nbins+pbins-1] = *sensinfo[k].binsize

    endfor

; check for NaNs or masked points; reject break points that are inside
; masked bins
    
    flag = where((finite(fitsens) eq 0B) or (fitsens lt -900.0),$
      nflag,comp=good,ncomp=ngood)
    if (nflag ne 0L) then fitmask[flag] = 0B

    fitvar = variance(fitsens[good],/double)
    if (fitvar eq 0) then fitvar = 1.0
    fitinvvar = 0.0 * fitsens + 1.0/fitvar

; construct break points that are equally spaced among the non-masked
; wavelengths 
    
    bkptres = (max(fitwave)-min(fitwave))/bsorder_sens ; break point resolution [Angstrom]
    bkpt = findgen(bsorder_sens)*bkptres+min(fitwave)

    nbkpt = n_elements(bkpt)
    keep = lonarr(nbkpt)+1L

    masked = where(fitmask eq 0B,nmasked)
    if (nmasked ne 0L) then begin
       
       srtwave = sort(fitwave[masked])
       uindx = uniq(fitwave[masked],srtwave)
       ufitwave = fitwave[masked[uindx]]
       ubinsize = fitbinsize[masked[uindx]]

       for ibkpt = 0L, nbkpt-1L do $
         for iwave = 0L, n_elements(ufitwave)-1L do $
         if (bkpt[ibkpt] ge ufitwave[iwave]-ubinsize[iwave]/2.0) and $
           (bkpt[ibkpt] le ufitwave[iwave]+ubinsize[iwave]/2.0) then $
         keep[ibkpt] = 0L

    endif

    bkpt = bkpt[where(keep)] ; this should never crash

    sset = bspline_iterfit(fitwave,fitsens,invvar=fitinvvar*fitmask,bkpt=bkpt,$
      nord=nord_sens,outmask=sensmask,upper=3.0,lower=3.0,/silent,_extra=extra)
    coarsefit = bspline_valu(fitwave,sset,lower=3.0,upper=3.0,mask=sensmask)
    
    senswave = findgen((max(fitwave)-min(fitwave))/sensres+1.0)*sensres+min(fitwave) 
    sensfit = bspline_valu(senswave,sset)

; determine the error in the sensitivity function fit by computing the
; standard error over a wavelength interval that is much larger than
; the binsize; then fit a low-order b-spline to the error function;
; this procedure ensures that the uncertainty varies smoothly, and not
; susceptible to spurious bins

    errorsampling = 0.01
    nerrorbins = 1.0/errorsampling
    errorbinsize = errorsampling*(max(fitwave)-min(fitwave)) ; [Angstroms]
    errorbins = findgen(nerrorbins+1L)*errorbinsize+min(fitwave)

    sensbinerror = fltarr(nerrorbins)-999.0
    
    for ibin = 0L, nerrorbins-2L do begin

       inrange = where((fitwave gt errorbins[ibin]) and (fitwave le errorbins[ibin+1L]),ninrange)
       if (ninrange gt 2L) then begin
          diff = coarsefit[inrange]-fitsens[inrange]
          if keyword_set(greyshifted) then $
            sensbinerror[ibin] = djsig(diff,sigrej=3.0)/sqrt(ninrange) else $
            sensbinerror[ibin] = djsig(diff,sigrej=3.0)
       endif
       
    endfor

; fit a smooth function    
    
    goodbins = where(sensbinerror gt -900.0,ngoodbins)
    if (ngoodbins gt 10L) then begin

       bsorder_senserr = 3.0
       bkptres = (max(errorbins[goodbins])-min(errorbins[goodbins]))/bsorder_senserr
       bkpt = findgen(bsorder_senserr)*bkptres+min(errorbins[goodbins])

       errorsset = bspline_iterfit(errorbins[goodbins],sensbinerror[goodbins],bkpt=bkpt,/silent)
       sensfiterr = bspline_valu(senswave,errorsset)

;      poly_senserr = 3L
;      im_poly_iter, errorbins[goodbins], sensbinerror[goodbins], poly_senserr, 3.0, coeff=coeff
;      sensfiterr = poly(senswave,coeff)

    endif else sensfiterr = sensfit*0.0

    if keyword_set(doplot) then begin

       if keyword_set(greyshifted) then begin
          title = 'Grey-Shifted Sensitivity Function' 
          ytitle = 'Sensitivity Function Error (mag)'
          im_window, 3, xratio=0.35, /square
       endif else begin
          title = 'Observed Sensitivity Function'
          ytitle = 'Absolute Error (mag)'
          im_window, 5, xratio=0.35, /square ; 5 won't interfere with ITRACESPEC
       endelse
       
       plot, errorbins[goodbins], sensbinerror[goodbins], ps=4, xthick=2.0, $
         ythick=2.0, charsize=1.2, charthick=2.0, xtitle=textoidl('Wavelength (\AA)'), $
         ytitle=ytitle, xsty=3, ysty=3, title=title
       oplot, senswave, sensfiterr, line=0, thick=2.0

    endif

return
end

pro write_sensfuncfit, info, senswave, sensfit, sensfiterr, sensres, $
  sensname, zptshift, histnote=histnote
; jm03apr22uofa

    nstds = n_elements(info)

    mkhdr, senshead, sensfit, /extend
    sxaddpar, senshead, 'OBJECT', 'Sensitivity Function'
    sxaddpar, senshead, 'CRVAL1', min(senswave), ' wavelength at the reference pixel'
    sxaddpar, senshead, 'CRPIX1', 1.0, ' reference pixel number'
    sxaddpar, senshead, 'CD1_1', sensres, ' dispersion in Angstroms per pixel'
    sxaddpar, senshead, 'CDELT1', sensres, ' dispersion at the reference pixel'
    sxaddpar, senshead, 'CTYPE1', 'LINEAR', ' projection type'
    sxaddpar, senshead, 'INSTAR', nstds, ' number of standard stars'
    sxaddpar, senshead, 'ZPTSHIFT', zptshift, ' zero point shift [mag]'
    if (n_elements(histnote) ne 0L) then sxaddhist, histnote, senshead

; log of the individual stars
    
    for i = 1L, nstds do sxaddpar, senshead, 'ISTAR'+string(i,format='(I2.2)'), $
      strcompress(info[i-1L].tablename,/remove), before='HISTORY'

; add the seeing measurements to the header
    
    for i = 1L, nstds do sxaddpar, senshead, 'SEEING'+string(i,format='(I2.2)'), $
      info[i-1L].seeing, format='(F5.3)', ' FHWM [arcsec]', before='HISTORY'
    
    splog, 'Writing '+sensname+'.'
    mwrfits, sensfit+zptshift, sensname, senshead, /create
    mwrfits, sensfiterr, sensname ; extension 1
    
return
end

pro sensfuncplot, sensinfo, info, fitsens1, fitwave1, fitmask, coarsefit1, $
  senswave, sensfit, sensfiterr, bsorder_sens, nord_sens, sset, $
  zptshift, greyleg=greyleg, senstitle=senstitle, residuals=residuals, $
  thick=thick

; this hack is needed to get the colors on the sensitivity function
; fit right both for display purposes *and* postscript output
    
    nstds = n_elements(info)

; use a black background for maximum contrast

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
    
; residuals plot

    good = where(fitmask eq 1B)
    fitsens = fitsens1[good]
    fitwave = fitwave1[good]
    coarsefit = coarsefit1[good]
    
    residuals = fitsens - coarsefit
    percent_residuals = 100.0*(10.0^(0.4*residuals)-1.0)
;   percent_residuals = 100.0*((10^(0.4*fitsens)-10^(0.4*coarsefit))/10^(0.4*fitsens))
    maxdev = fltarr(nstds)

    djs_iterstat, abs(residuals), sigrej=5.0, sigma=rsig, mean=rmean

    plot, [0], [0], xsty=3, ysty=11, xtitle=textoidl('Wavelength (\AA)'), $
      xrange=minmax(fitwave), ytitle='Residuals (mag)', charsize=1.1, charthick=thick, $
      xthick=thick, ythick=thick, position=[0.13,0.12,0.88,0.35], yminor=3, $
      yrange=(rmean+5.0*rsig)*[-1,1], /noerase, color=fsc_color('white',!d.table_size-55)
      
    axis, /yaxis, yrange=interpol(percent_residuals,residuals,!y.crange), ysty=1, $
      charsize=1.1, charthick=thick, ythick=thick, ytitle='Residuals (%)', $
      yminor=3, color=fsc_color('white',!d.table_size-55)
    for k = 0L, nstds-1L do begin
       diff = *sensinfo[k].sens - bspline_valu(*sensinfo[k].bandwave,sset)
       maxdev[k] = max(abs(diff))
       im_symbols, info[k].psym, color=fsc_color(info[k].color,!d.table_size-k-2), $ ;fill=1, $
         thick=(thick-2L)>2L, fill=info[k].psymfill, psize=info[k].psymsize
       good = where(*sensinfo[k].bandmask eq 1B,ngood,comp=bad,ncomp=nbad)
       if (ngood ne 0L) then oplot, (*sensinfo[k].bandwave)[good], diff[good], $
         ps=8, thick=(thick-2L)>2L, symsize=0.8
;      if (nbad ne 0L) then oplot, (*sensinfo[k].bandwave)[bad], diff[bad], ps=7, syms=2.0
    endfor
    oplot, [!x.crange[0],!x.crange[1]], [0,0], line=0, thick=thick, color=fsc_color('white',!d.table_size-55)

; sensitivity function plot.  apply the zero point shift in the plot
; by just shifting YRANGE

    mnmx = max(fitsens-zptshift)-min(fitsens-zptshift)
    yrange = [min(fitsens-zptshift)-mnmx*0.30,max(fitsens-zptshift)+mnmx*0.30]
;   yrange = minmax(fitsens-zptshift)*[1.0/1.001,1.018]

    plot, [0], [0], xsty=3, ysty=11, ytitle=textoidl('2.5 log [Counts s^{-1} '+$
      '\AA^{-1}] / ['+flam_units()+']'), charsize=1.1, charthick=thick, $
      yrange=yrange, xthick=thick, ythick=thick, position=[0.13,0.35,0.88,0.95], $
      xtickname=replicate(' ',10), xrange=minmax(fitwave), /noerase, $
      title=senstitle, color=fsc_color('white',!d.table_size-55)
    axis, yaxis=1, ysty=3, ythick=thick, charthick=thick, charsize=1.1, $
      yrange=yrange, ytitle=textoidl('2.5 log [Counts s^{-1} \AA^{-1}] / ['+flam_units()+']'), $
      color=fsc_color('white',!d.table_size-55)
    for k = 0L, nstds-1L do begin
       im_symbols, info[k].psym, color=fsc_color(info[k].color,!d.table_size-k-2), $ ; fill=1
         thick=(thick-2L)>2L, fill=info[k].psymfill, psize=info[k].psymsize
       good = where(*sensinfo[k].bandmask eq 1B,ngood,comp=bad,ncomp=nbad)
       if (ngood ne 0L) then oplot, (*sensinfo[k].bandwave)[good], (*sensinfo[k].sens)[good], $
         ps=8, thick=(thick-2L)>2L, symsize=0.8
;      if (nbad ne 0L) then oplot, (*sensinfo[k].bandwave)[bad], (*sensinfo[k].sens)[bad], ps=7, syms=2.0
    endfor
    oplot, senswave, sensfit, line=0, thick=thick+2, color=fsc_color('black',!d.table_size-56)
    oplot, senswave, sensfit, line=0, thick=thick, color=fsc_color('white',!d.table_size-55)

; overplot the b-spline breakpoints

    bkplot = where(sset.bkmask eq 1L,comp=nobkplot)
    djs_oplot, sset.fullbkpt[bkplot], bkplot*0.0+!y.crange[0]+0.025*(!y.crange[1]-!y.crange[0]), $
      psym=7, symsize=1.0
    if (nobkplot[0] ne -1L) then djs_oplot, sset.fullbkpt[nobkplot], $
      bkplot*0.0+!y.crange[0]+0.025*(!y.crange[1]-!y.crange[0]), psym=4

    im_legend, ['b-spline '+string(long(bsorder_sens),format='(I0)'),'order '+$
      string(nord_sens,format='(I0)'),string(nstds,format='(I0)')+' stars',$
      'ZP Shift = '+string(zptshift,format='(G0.0)')+' mag', $
      greyleg], /left, /top, box=0, charsize=1.3, charthick=thick, $
      textcolor=fsc_color('white',!d.table_size-55)

    im_legend, [textoidl('\Delta')+' = '+string(mean(residuals),format='(F6.3)')+' (mag)',$
      textoidl('\sigma_{\Delta}')+' = '+string(stddev(residuals),format='(F6.3)')+' (mag)',$
      'Mean Scatter = '+string(stddev(percent_residuals),format='(F5.1)')+'%'], $
;     'Max  Scatter = '+string(max(abs(percent_residuals)),format='(F4.1)')+'%'], $
;     'Scatter = '+string((10.0^stddev(residuals)-1)*100.0,format='(F4.1)')+'%'], $
      /right, /top, box=0, charsize=1.3, charthick=thick, textcolor=fsc_color('white',!d.table_size-55)

; choose the most deviant stars and put them on the legend in reverse
; deviance order!

    indx = (reverse(sort(maxdev)))[0:4<(nstds-1L)]
    
    c = lonarr(n_elements(indx))
    for k = 0L, n_elements(indx)-1L do c[k] = fsc_color(info[indx[k]].color,!d.table_size-k-2)
    
    starlabel = string(info[indx].id,format='(I3)')+' '+$
      info[indx].starname+' '+info[indx].date+' '+$
      string(info[indx].airmass,format='(F7.4)')+' '+$
      string(info[indx].greyshift,format='(F7.4)')+' '+$
      string(info[indx].seeing,format='(F5.3)')+' '+$
      string(100*info[indx].energy,format='(I3)')+'%'
    
    im_legend, [starlabel], /right, /bottom, box=0, charsize=1.0, charthick=(thick-2L)>2L, $
      textcolor=c, color=c, psym=info[indx].psym, thick=(thick-2L)>2L, $
      fill=info[indx].psymfill
;     textcolor=[fsc_color('white',200),c], color=[fsc_color('white',200),c], psym=[-1,info[indx].psym]

return
end

function isensfunc, stdlist, datapath=datapath, sensname=sensname, extfile=extfile, $
  bsorder_sens=bsorder_sens, nord_sens=nord_sens, zptshift=zptshift, senstitle=senstitle, $
  blockwave=blockwave, blocktrans=blocktrans, blockfactor=blockfactor, grey=grey, $
  slit_width=slit_width, doplot=doplot, write=write, _extra=extra

    nstds = n_elements(stdlist)
    if nstds eq 0L then begin
       print, 'Syntax - info = isensfunc(stdlist,datapath=,sensname=,extfile=,$'
       print, '   bsorder_sens=,nord_sens=,zptshift=,senstitle=,grey=,slit_width=,$'
       print, '   /doplot,/write,_extra=extra)'
       return, -1L
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(nord_sens) eq 0L) then nord_sens = 3L
    if (n_elements(bsorder_sens) eq 0L) then bsorder_sens = 50.0 else $
      bsorder_sens = float(bsorder_sens)
    if (bsorder_sens le 0.0) then begin
       splog, 'BSORDER_SENS must be positive.'
       return, -1L
    endif

    if (n_elements(zptshift) eq 0L) then zptshift = 0.0 ; [mag]
    
    if (n_elements(senstitle) eq 0L) then senstitle = ''
    if (n_elements(senstitle) ne 1L) then senstitle = senstitle[0] ; this error arose jm04sep10uofa

    if (n_elements(grey) eq 0L) then grey = 1L else grey = long(grey)
    if (grey ne 1L) and (grey ne 2L) then grey = 1L

    if (n_elements(sensname) eq 0L) then sensname = 'sensfunc.fits'
    greysensname = repstr(sensname,'.fits','_grey.fits')
    
    sensroot = strmid(sensname,0,strpos(sensname,'.fits'))
    psname = datapath+'qaplot_'+sensroot+'.ps'
    greypsname = repstr(psname,'.ps','_grey.ps')
    
; read the extinction file
    
    if (n_elements(extfile) eq 0L) then extfile = 'kpnoextinct.dat'
    etcpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')

    if (file_test(etcpath+extfile,/regular) eq 0L) then begin
       splog, 'Extinction file '+etcpath+extfile+' not found.'
       return, -1L
    endif
    
    if keyword_set(write) then openw, lun, datapath+'qalog_'+sensroot+'.datalog', /get_lun
    
    info1 = {$
      id:          0L,$
      specname:    '',$
      starname:    '',$
      tablename:   '',$
      date:        '',$
      jd:          0D,$ ; julian date
      ra:          '',$ ; HMS
      dec:         '',$ ; DMS
      epoch:      0.0,$
      minwave:    0.0,$
      maxwave:    0.0,$
      dwave:      0.0,$
      exptime:    0.0,$
      airmass:    0.0,$
      posangle:   0.0,$
      parangle:   0.0,$
      seeing:     0.0,$ ; arcsec
      energy:     0.0,$ ; encircled energy
      greyshift:  0.0,$
      greywave:   0.0,$
      senszero:   0.0,$
      psym:        0L,$
      psymfill:    0L,$
      psymsize:   0.0,$
      color:       ''}

    sensinfo1 = {$
      bandwave: ptr_new(), $
      sens:     ptr_new(), $
      bandmask: ptr_new(), $
      binsize:  ptr_new(), $
      nbins:    0L}

; read the iSPEC standards FITS table and the extinction file 
    
    rootpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/spectra')
    splog, 'Reading '+rootpath+'ispec_standards.fits'
    starinfo = mrdfits(rootpath+'ispec_standards.fits',1,/silent)

    splog, 'Reading the extinction file '+etcpath+extfile+'.'
    readcol, etcpath+extfile, extwave, extvals, format='F,F', /silent

; oversampling variables    
    
    finerpix = 50.0 ; oversampling factor
    sensres = 5.0   ; resolution of the output sensitivity function [Angstrom]
    
; loop on each standard star observation

    for i = 0L, nstds-1L do begin

       if (file_test(datapath+stdlist[i],/regular) eq 0L) then begin
          splog, 'File '+datapath+stdlist[i]+' not found.'
          continue
       endif

       info1.id = i
       info1.specname = stdlist[i]
       
; read the spectrum
       
       print, format='("Reading ",A0," . . . ",$)', strtrim(info1.specname,2)
;      splog, 'Reading '+info1.specname+'.'
       counts = mrdfits(datapath+info1.specname,0,header,/silent)
       
       wave = make_wave(header)
       if (wave[0] eq -1L) then continue ; insufficient wavelength information
       npix = n_elements(wave)

; fill the info structure with header information

       jd = sxpar(header,'JD')
       date = jd2fitsdate(jd)
       
       info1.starname = strcompress(strupcase(sxpar(header,'OBJECT')),/remove)
       info1.date     = date
       info1.jd       = jd
       info1.ra       = sxpar(header,'RA') ; 15.0*im_hms2dec(sxpar(header,'RA'))
       info1.dec      = sxpar(header,'DEC') ; im_hms2dec(sxpar(header,'DEC'))
       info1.epoch    = sxpar(header,'EPOCH')
       info1.minwave  = float(sxpar(header,'CRVAL1'))
       info1.dwave    = float(sxpar(header,'CD1_1'))
       info1.maxwave  = float(info1.minwave + info1.dwave * (sxpar(header,'naxis1')-1.0))
       info1.exptime  = float(sxpar(header,'EXPTIME'))
       info1.airmass  = float(sxpar(header,'AIRMASS'))
       info1.parangle = float(sxpar(header,'PARANGLE'))
       posangle = sxpar(header,'POSANGLE',count=pacount)
       if (pacount eq 1L) then info1.posangle = float(posangle)

; retrieve the seeing and compute the encircled energy       
       
       seeing = sxpar(header,'SEEING',count=nseeing)
       if (nseeing ne 0L) then begin
          if (n_elements(slit_width) eq 0L) then slit_width = 0.0
          info1.seeing = seeing
          info1.energy = errorf(slit_width/2.0/seeing/sqrt(2.0))
       endif

       if keyword_set(write) then printf, lun, info1.starname, info1.specname
       print, format='(A0,$)', strtrim(info1.starname,2)

; extinction correct to zero airmass 
       
       extinct = interpol(extvals,extwave,wave,/quadratic)
       truecounts = counts * 10.0D^(0.4D*info1.airmass*extinct)

       starindx = istarsearch(info1.ra,info1.dec,epoch=info1.epoch,$
         starinfo=starinfo,_extra=extra)
       if (starindx[0] eq -1L) then begin
          print, ' . . . no standard star found.'
          continue
       endif else print, ' . . . found '+strtrim(starinfo[starindx].star,2)+$
         ' ('+strtrim(starinfo[starindx].starpath,2)+').'

       stdfile = strtrim(starinfo[starindx].file,2)
       info1.tablename = strtrim(starinfo[starindx].star,2)
       info1.color = strtrim(starinfo[starindx].starcolor,2)
       
; all the standard star data files have the same format
; (UPDATE_STANDARDS_DATABASE)

       readcol, rootpath+stdfile, bandwave, stdmag, binsize, $
         format='D,D', comment='#', /silent
       nbands = n_elements(bandwave)

       bandstart = bandwave - binsize/2.0
       bandend = bandwave + binsize/2.0

; include an optional correction for second-order contamination

       stdflux = 10D^(-0.4*(stdmag+48.6D)) * 2.99793D18 / bandwave / bandwave ; [erg/s/cm2/A]

       if (n_elements(blockwave) ne 0L) and (n_elements(blocktrans) ne 0L) then begin

          if (n_elements(blockwave) ne n_elements(blocktrans)) then $
            splog, 'BLOCKWAVE and BLOCKTRANS must have the same number of elements!' else begin

             if (n_elements(blockfactor) eq 0L) then blockfactor = 0.1
             
             linterp, blockwave, blocktrans, bandwave, blocktrans_band, /nointerp 
             blocktrans_band = blocktrans_band > 0          
             order2 = blocktrans_band * stdflux * float(blockfactor)
             linterp, bandwave*2, order2, bandwave, extraflux

;            djs_plot, bandwave, stdflux, xsty=3, ysty=3, xr=[3000,9000]
;            djs_oplot, bandwave, stdflux+extraflux, color='red'

             stdflux = stdflux + extraflux

          endelse
       endif 

; reject bandpasses that are outside of the spectral range

       good = where((bandstart gt min(wave)) and (bandend lt max(wave)),nbins)

       stdflux = stdflux[good]
       bandwave = bandwave[good]
       bandstart = bandstart[good]
       stdmag = stdmag[good]
       bandend = bandend[good]
       binsize = binsize[good]

       extband = interpol(extvals,extwave,bandwave,/quadratic)

       sensinfo1.nbins = nbins
       sensinfo1.bandwave = ptr_new(bandwave)
       sensinfo1.binsize = ptr_new(binsize)

; sum the counts in each bin using partial pixels and construct the
; sensitivity curve; reject bandpasses with a sensitivity equal to
; -999.0; the latter will occur if there are gaps in the standard-star
; observations

       sens = fltarr(nbins)-999.0
       nanmask = bytarr(nbins)+1B ; 1B is good
       
       for j = 0L, nbins-1L do begin

          wavearray = bandstart[j] + findgen(binsize[j]*finerpix)/finerpix + 1.0/finerpix/2.0
          
          bincounts = interpol(truecounts,wave,wavearray)
          starcounts = total(bincounts)/binsize[j]/finerpix
;         if (j gt 10L) and (j lt 20L) then starcounts = -100.0
          
          value = starcounts / info1.exptime / info1.dwave / stdflux[j]
          if (value gt 0.0) then sens[j] = 2.5*alog10(value) else nanmask[j] = 0B

;         if (bandwave[j] gt 4845.0) and (bandwave[j] lt 4865.0) then $
;           print, bandwave[j], stdflux[j], starcounts, sens[j]
          
          if keyword_set(write) then printf, lun, bandwave[j], stdflux[j], starcounts

       endfor

; reject bandpasses affected by telluric absorption

       tmask = telluric_mask(wave,good=tfree,bad=tbad)

       telluric_bandmask = bytarr(nbins)+1B ; 1B is good
       for iband = 0L, nbins-1L do begin
          inband = where((wave gt bandstart[iband]) and (wave lt bandend[iband]),ninband)
          if (ninband eq 0L) then telluric_bandmask[iband] = 0B else $
            if total((tmask[inband] eq 0B)) ne 0.0 then telluric_bandmask[iband] = 0B
       endfor

; reject bandpasses that contain the Balmer absorption lines and the
; 5577 night sky line

       rejectwaves = [4101.0,4340.0,4861.0,5577.0,6563.0]
       rejectdwave = 10.0            
       rejectmask = bytarr(nbins)+1B ; 1B is good

       for ireject = 0L, n_elements(rejectwaves)-1L do begin
          rejbands = where((bandstart ge rejectwaves[ireject]-rejectdwave) and $
            (bandend le rejectwaves[ireject]+rejectdwave),nreject)
          if (nreject ne 0L) then rejectmask[rejbands] = 0B
       endfor

       bandmask = (telluric_bandmask + nanmask + rejectmask) eq 3B
;      bandmask = (telluric_bandmask + nanmask) eq 2B

       sensinfo1.sens = ptr_new(sens)
       sensinfo1.bandmask = ptr_new(bandmask)
;      niceprint, bandwave, sens, long(telluric_bandmask), long(nanmask), long(bandmask)

; grow the information structure

       if (n_elements(info) eq 0L) then info = info1 else info = [ [info], [info1] ]
       if (n_elements(sensinfo) eq 0L) then begin
          sensinfo = sensinfo1
       endif else begin
          sensinfo = [ [sensinfo], [sensinfo1] ]
       endelse

    endfor 

    nstds = n_elements(info)
    if (nstds eq 0L) then begin
       splog, 'No good standard stars found.'
       cleansens, sensinfo
       return, -1L
    endif
    
    info = reform(info)
    sensinfo = reform(sensinfo)

    if keyword_set(write) then free_lun, lun ; close the data log file

; assign a plot symbol to each unique standard star

    bigpsym = im_best_symbols(fill=bigfill,psize=bigsize)
    
    srt = sort(info.tablename)
    unique = uniq(info.tablename,srt)
    nuniq = n_elements(unique)

    if (nuniq gt n_elements(bigpsym)) then begin
       splog, 'More symbols must be defined in IM_BEST_SYMBOLS() before proceeding.'
       cleansens, sensinfo & cleansens, sensinfo1
       return, -1L
    endif
    
    psymvec = bigpsym[0L:nuniq-1L] ; assign plot symbols
    fillvec = bigfill[0L:nuniq-1L]
    sizevec = bigsize[0L:nuniq-1L]
    
    for k = 0L, nuniq-1L do begin
       w = where(strmatch(strtrim(info.tablename,2),strtrim(info[unique[k]].tablename,2),/fold_case) eq 1B,nw)
       info[w].psym = psymvec[k]
       info[w].psymfill = fillvec[k]
       info[w].psymsize = sizevec[k]
    endfor

; fit, write out, and plot the observed sensitivity function 

    splog, 'Fitting a b-spline to the observed sensitivity curve.'
    splog, 'Zero-point shift = '+string(zptshift,format='(G0.0)')+' mag.'

    sensfuncfit, sensinfo, bsorder_sens, nord_sens, sensres, $
      _extra=extra, fitsens, fitwave, fitmask, coarsefit, $
      senswave, sensfit, sensfiterr, sset, doplot=doplot, $
      greyshifted=0L

    if keyword_set(doplot) then begin
       im_window, 0, xratio=0.49, /square
       sensfuncplot, sensinfo, info, fitsens, fitwave, fitmask, coarsefit, $
         senswave, sensfit, sensfiterr, bsorder_sens, nord_sens, sset, $
         zptshift, greyleg='No Grey Shift', senstitle=senstitle, thick=2.0
    endif

; generate postscript output of the observed sensitivity function 
    
    if keyword_set(write) then begin

       splog, 'Writing '+psname+'.'
       dfpsplot, psname, /square, /color, /isolatin1
       sensfuncplot, sensinfo, info, fitsens, fitwave, fitmask, coarsefit, $
         senswave, sensfit, sensfiterr, bsorder_sens, nord_sens, sset, $
         zptshift, greyleg='No Grey Shift', senstitle=senstitle, residuals=resid, $
         thick=8.0
       dfpsclose

       write_sensfuncfit, info, senswave, sensfit, sensfiterr, sensres, $
         datapath+sensname, zptshift, histnote="'Observed Sensitivity Function'"

    endif
    
; compute some numbers for grey-shifting.  calculate the median
; sensitivity centered on the bin nearest to the average wavelength of
; all the standards, plus or minus one bin.  this procedure should
; make the shift independent of the wavelength range of the data, and
; impervious to bad bins from cosmic rays; account for masked
; bandpasses 

    gwave = (djs_mean(info.maxwave)-djs_mean(info.minwave))/2.0+djs_mean(info.minwave)
    info.greywave = gwave

    for i = 0L, nstds-1L do begin
       good = where(*sensinfo[i].bandmask eq 1L)
       get_element, (*sensinfo[i].bandwave)[good], gwave, xx
       info[i].senszero = median(((*sensinfo[i].sens)[good])[xx-1L:xx+1L])
;      info[i].senszero = interpol(*sensinfo[i].sens,*sensinfo[i].bandwave,5556.0)
    endfor

    maxsens = max(info.senszero,maxindx)  ; maximum sensitivity 
    meansens = djs_mean(info.senszero) ; mean sensitivity

    print, format='("Grey-shifting ",A0, " standards centered on ", A0, " Angstrom to the",$)', $
      string(nstds,format='(I0)'), string(gwave,format='(I4)')
    
    if (grey eq 1L) then begin  ; grey-shift to the maximum sensitivity
       print, ' maximum sensitivity.'
       info.greyshift = info[maxindx].senszero-info.senszero
       greyleg = 'Maximum Grey Shift' & ghistnote = 'Maximum'
    endif
    
    if (grey eq 2L) then begin  ; grey-shift to the mean sensitivity
       print, ' mean sensitivity.'
       info.greyshift = meansens-info.senszero
       greyleg = 'Mean Grey Shift' & ghistnote = 'Mean'
    endif

; do the grey-shifting now
    
    for i = 0L, nstds-1L do begin

       greysens = *sensinfo[i].sens + info[i].greyshift
       ptr_free, sensinfo[i].sens & sensinfo[i].sens = ptr_new(greysens)

    endfor

; fit, write out, and plot the grey-shifted sensitivity function 

    splog, 'Fitting a b-spline to the grey-shifted sensitivity curve.'
    splog, 'Zero-point shift = '+string(zptshift,format='(G0.0)')+' mag.'

    sensfuncfit, sensinfo, bsorder_sens, nord_sens, sensres, $
      _extra=extra, fitsens, fitwave, fitmask, coarsefit, $
      senswave, sensfit, sensfiterr, sset, doplot=doplot, $
      greyshifted=1L

    if keyword_set(doplot) then begin
       im_window, 2, xratio=0.49, /square
       sensfuncplot, sensinfo, info, fitsens, fitwave, fitmask, coarsefit, $
         senswave, sensfit, sensfiterr, bsorder_sens, nord_sens, sset, $
         zptshift, greyleg=greyleg, senstitle=senstitle, thick=2.0
    endif
    
; generate postscript output of the grey-shifted sensitivity function 
    
    if keyword_set(write) then begin

       splog, 'Writing '+greypsname+'.'

       dfpsplot, greypsname, /square, /color, /isolatin1
       sensfuncplot, sensinfo, info, fitsens, fitwave, fitmask, coarsefit, $
         senswave, sensfit, sensfiterr, bsorder_sens, nord_sens, sset, $
         zptshift, greyleg=greyleg, senstitle=senstitle, residuals=greyresid, $
         thick=8.0
       dfpsclose

       write_sensfuncfit, info, senswave, sensfit, sensfiterr, sensres, $
         datapath+greysensname, zptshift, histnote="'"+ghistnote+" Grey-Shifted Sensitivity Function'"

    endif
    
; crop the output information structure to only useful fields

    outinfo = struct_trimtags(info,select=['SPECNAME','STARNAME','TABLENAME','DATE',$
      'AIRMASS','PARANGLE','SEEING','ENERGY','GREYSHIFT','SENSZERO','PSYM','COLOR'])
    
    if keyword_set(write) then begin 

; compute residuals of the fit
       
       residstr = string(djs_mean(resid),format='(F10.7)')+', '+$
         string(djs_median(resid),format='(F10.7)')+', '+$
         string(min(resid),format='(F10.7)')+', '+$
         string(max(resid),format='(F10.7)')+', '+$
         string(stddev(resid),format='(F10.7)')

       if n_elements(grey) ne 0L then $
         greyresidstr = string(djs_mean(greyresid),format='(F10.7)')+', '+$
         string(djs_median(greyresid),format='(F10.7)')+', '+$
         string(min(greyresid),format='(F10.7)')+', '+$
         string(max(greyresid),format='(F10.7)')+', '+$
         string(stddev(greyresid),format='(F10.7)') else $
         greyresidstr = string(0.0,format='(F10.7)')+', '+$
         string(0.0,format='(F10.7)')+', '+string(0.0,format='(F10.7)')+', '+$
         string(0.0,format='(F10.7)')+', '+string(0.0,format='(F10.7)')

; write the log file

       openw, lun, datapath+'qalog_'+sensroot+'.log', /get_lun
       printf, lun, '## ISENSFUNC log file '+im_today()
       printf, lun, '## Starting and ending wavelength: '+string(minmax(senswave),$
         format='(F7.2,2x,F7.2)')
       printf, lun, '## b-spline fit (order,nbkpoints): '+string(nord_sens,format='(I0)')+', '+$
         string(n_elements(sset.fullbkpt),format='(I0)')
       printf, lun, '## '+greyleg
       printf, lun, '## Observed Residuals   (Mean, Median, Min, Max, Stdv): '+residstr
       printf, lun, '## Grey-Shift Residuals (Mean, Median, Min, Max, Stdv): '+greyresidstr
       printf, lun, '#  1 SPECNAME'
       printf, lun, '#  2 STARNAME'
       printf, lun, '#  3 TABLENAME'
       printf, lun, '#  4 DATE'
       printf, lun, '#  5 AIRMASS'
       printf, lun, '#  6 PARANGLE'
       printf, lun, '#  7 SEEING'
       printf, lun, '#  8 ENERGY'
       printf, lun, '#  9 GREYSHIFT'
       printf, lun, '# 10 SENSZERO'
       printf, lun, '# 11 PSYM'
       printf, lun, '# 12 COLOR'
       struct_print, outinfo, lun=lun, /no_head
       free_lun, lun    
       
    endif

    if keyword_set(doplot) then struct_print, struct_trimtags(info,select=$
      ['ID','SPECNAME','STARNAME','TABLENAME','DATE','AIRMASS','POSANGLE','PARANGLE',$
      'SEEING','ENERGY','GREYSHIFT','SENSZERO'])

    cleansens, sensinfo ; cleanup memory
    cleansens, sensinfo1
    
return, info
end
