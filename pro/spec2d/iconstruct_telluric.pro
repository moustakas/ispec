;+
; NAME:
;       ICONSTRUCT_TELLURIC
;
; PURPOSE:
;       Construct a telluric absorption-line spectrum from a list of
;       one-dimensional (telluric) standard star spectra.
;
; CALLING SEQUENCE:
;       iconstruct_telluric, tlist, datapath=, tellmethod=, $
;          contmethod=, tcorr_bsorder=, tcorr_nord=, tcorr_sigclip=, $
;          tcorr_ncoeff=, tcorr_pixspace=, plottitle=, psname=, $
;          tellfits=, tellwave=, tellcorr=, /debug, /doplot, /write, $
;          _extra=extra
;
; INPUTS:
;       tlist - list of one-dimensional telluric spectra [NTELL] 
;
; OPTIONAL INPUTS:
;       datapath   - I/O path
;
;       tellmethod - telluric correction method (see COMMENTS) 
;          1: global telluric correction
;          2: local telluric correction (default)
;
;       contmethod - continuum normalization method (see COMMENTS) 
;          1: fit the full continuum in the 6750-8250 A wavelength
;             range using any of the fitting techniques implemented in
;             IM_FITCONTINUUM(); CONTMETHOD=1 is *required* by
;             TELLMETHOD=1, and is the default for TELLMETHOD=2 
;          2: local Legendre fit to the continuum (only available for
;             TELLMETHOD=2) 
;          3: local b-spline fit to the continuum (only available for 
;             TELLMETHOD=2) 
;
;       tcorr_bsorder  - order of the b-spline fit either to the full
;                        super-sampled telluric spectrum if
;                        TELLMETHOD=1, or to each individual telluric
;                        absorption band if TELLMETHOD=2 (default 10) 
;       tcorr_nord     - order of the fit between breakpoints (default
;                        4) 
;       tcorr_sigclip  - sigma-clipping threshold (default 3.0)
;       tcorr_ncoeff   - number of Legendre coefficients (default 4)  
;       tcorr_pixspace - pixel spacing of the b-spline knots (default 5) 
;       plottitle      - QA plot title 
;       psname         - if DOPLOT=1 or WRITE=1 then generate a
;                        postscript QA plot file with this file name
;                        (default 'qaplot_telluric_spectrum.ps')
;       tellfits       - if WRITE=1 then write out the average
;                        telluric correction spectrum (default
;                        'telluric_spectrum.fits')
;       extra          - extra keywords for RD1DSPEC(),
;                        BSPLINE_ITERFIT(), IM_FITCONTINUUM() 
;
; KEYWORD PARAMETERS:
;       debug    - generate a plot useful for debugging then wait for
;                  a keystroke
;       doplot   - generate a QAPLOT
;       write    - write out the mean telluric correction spectrum 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       tellwave - telluric correction spectrum wavelength vector
;                  [Angstroms] 
;       tellcorr - telluric correction spectrum [unitless]
;
; PROCEDURES USED:
;       RD1DSPEC(), ICLEANUP, SPLOG, BSPLINE_ITERFIT(),
;       BSPLINE_VALU(), FUNC_FIT(), FLEGENDRE(), DJS_REJECT(),
;       DJS_PLOT, DJS_OPLOT, MWRFITS, DFPSPLOT, DFPSCLOSE, ANGSTROM(),
;       REPSTR(), IM_FITCONTINUUM(), PLOTSYM, SXADDPAR, SXADDHIST 
;
; COMMENTS:
;       Two algorithms have been implemented for generating the mean
;       telluric correction spectrum.  The first algorithm,
;       TELLMETHOD=1, assumes that the only absorption features in the
;       spectra provided in TLIST are telluric features.
;       Consequently, after the low-frequency continuum has been
;       divided out (CONTMETHOD=1), the residual spectrum *is* the
;       telluric correction spectrum.  This method is the "global"
;       telluric correction, and depends sensitively on having already
;       masked non-telluric absorption-line features from the input
;       spectra.
;
;       The second algorithm is the "local" telluric correction method
;       (TELLMETHOD=2).  This technique uses pre-defined bandpasses to
;       find the prominent telluric absorption features.
;       Consequently, almost any stellar spectrum can be used to model
;       the telluric absorption-line spectrum.  In this method you can
;       use one of three methods to normalize the continuum.  The
;       first method (CONTMETHOD=1), and the default, fits the
;       large-scale continuum in the same way as in TELLMETHOD=1.  The
;       other two methods fit the continuum locally around each
;       telluric band individually, using either a Legendre
;       polynomial (CONTMETHOD=2) or a b-spline (CONTMETHOD=3). 
;
;       Once the telluric correction spectrum has been generated for
;       each object in TLIST, then a high-order b-spline is fitted to
;       all the spectra simultaneously to generate the mean telluric
;       correction spectrum, which is subsequently written out and a
;       QA plot is generated.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 12, U of A, written, based in large
;          part on David Schlegel's TELLURIC_CORR SDSS routine 
;       jm04sep10uofa - significant improvements; adopt two distinct
;                       algorithms to construct the mean telluric
;                       correction spectrum; added TCORR_SIGCLIP optional
;                       input 
;       jm04nov07uofa - substantial revision and updates, and improved
;                       documentation 
;       jm05jun24uofa - TELLFILE optional input changed to TELLFITS
;                       for clarity
;
; Copyright (C) 2004-2005, John Moustakas
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

pro iconstruct_telluric, tlist, datapath=datapath, tellmethod=tellmethod, $
  contmethod=contmethod, tcorr_bsorder=tcorr_bsorder, tcorr_nord=tcorr_nord, $
  tcorr_sigclip=tcorr_sigclip, tcorr_ncoeff=tcorr_ncoeff, tcorr_pixspace=tcorr_pixspace, $
  plottitle=plottitle, psname=psname, tellfits=tellfits, tellwave=tellwave, $
  tellcorr=tellcorr, debug=debug, doplot=doplot, write=write, _extra=extra

    ntell = n_elements(tlist)
    if (ntell eq 0L) then begin
       print, 'Syntax - iconstruct_telluric, tlist, datapath=, tellmethod=, $'
       print, '   contmethod=, tcorr_bsorder=, tcorr_nord=, tcorr_sigclip=, $'
       print, '   tcorr_ncoeff=, tcorr_pixspace=, plottitle=, psname=, $'
       print, '   tellfits=, tellwave=, tellcorr=, /debug, /doplot, /write, $'
       print, '   _extra=extra'
       return
    endif

    splog, 'Using '+string(ntell,format='(I0)')+' telluric standards.'
    
    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(tcorr_bsorder) eq 0L) then tcorr_bsorder = 10.0
    if (n_elements(tcorr_nord) eq 0L) then tcorr_nord = 4
    if (n_elements(tcorr_sigclip) eq 0L) then tcorr_sigclip = 3.0
    if (n_elements(tcorr_ncoeff) eq 0L) then tcorr_ncoeff = 4
    if (n_elements(tcorr_pixspace) eq 0) then tcorr_pixspace = 5.0
    if (n_elements(plottitle) eq 0L) then plottitle = ''
    if (n_elements(psname) eq 0L) then psname = 'qaplot_telluric_spectrum.ps'
    if (n_elements(tellfits) eq 0L) then tellfits = 'telluric_spectrum.fits'

    if (n_elements(tellmethod) eq 0L) then tellmethod = 2L
    if (tellmethod lt 1L) or (tellmethod gt 2L) then tellmethod = 2L
    if (tellmethod eq 1L) then contmethod = 1L

    if (n_elements(contmethod) eq 0L) then contmethod = 1L
    if (contmethod lt 1L) or (contmethod gt 3L) then contmethod = 1L

    if keyword_set(write) then begin
       doplot = 0L
       debug = 0L
       postthick = 5.0
    endif else postthick = 2.0

; construct a structure describing which pixels to use for the local
; continuum-fitting, and which to use to construct the
; telluric-correction; these values are primarily used for
; TELLMETHOD=2, but they are also used to generate the zoomed-in QA
; plot 

    tellbands1 = { TELLBAND, $
      twave1: 6850., $
      twave2: 6960., $
      cwave1: [6750.0,6960.0,0.0], $ ; [6600., 6950., 0], $
      cwave2: [6850.0,7060.0,0.0]}   ; [6860., 7200., 0] }
    tellbands2 = { TELLBAND, $
      twave1: 7150., $
      twave2: 7350., $
      cwave1: [7050., 7115., 7340.], $
      cwave2: [7160., 7130., 7440.] }
    tellbands3 = { TELLBAND, $
      twave1: 7560., $
      twave2: 7720., $
      cwave1: [7400., 7700., 0], $
      cwave2: [7580., 8000., 0] }
    tellbands4 = { TELLBAND, $
      twave1: 8105., $
      twave2: 8240., $
      cwave1: [8000., 8225., 0], $
      cwave2: [8105., 8325., 0] }
;   tellbands5 = { TELLBAND, $
;     twave1: 8530., $
;     twave2: 8865., $
;     cwave1: [8490., 8865., 0], $
;     cwave2: [8530., 8905., 0] }
;   tellbands6 = { TELLBAND, $
;     twave1: 8644., $
;     twave2: 8697., $
;     cwave1: [8604., 8697., 0], $
;     cwave2: [8644., 8737., 0] }

    tellbands = [tellbands1, tellbands2, tellbands3, tellbands4]

; the following parameters describe the output telluric correction
; spectrum
    
    telldwave = 1.0 ; [Angstrom/pixel]

    mintellwave = 6750.0
    maxtellwave = 8250.0
;   mintellwave = min((tellbands.cwave1)[where(tellbands.cwave1 gt 0.0)])
;   maxtellwave = max(tellbands.cwave2)

    tellwave = findgen((maxtellwave-mintellwave)/telldwave+1)*telldwave+mintellwave
    
; read all the telluric standard star spectra, interpolating them all
; onto a common wavelength vector (for convenience, the reference
; spectrum is taken as the zeroth spectrum in TLIST)

    tspec1 = rd1dspec(tlist[0],datapath=datapath,_extra=extra)
    if (size(tspec1,/type) ne 8L) then return

    nbigpix = tspec1.npix
    bigwave = tspec1.wave

    inrange = where((bigwave gt mintellwave) and (bigwave lt maxtellwave),npix)
    if (npix eq 0L) then begin
       splog, 'There are no telluric features in the spectral range.'
       return
    endif
    
    bigflux = dblarr(npix,ntell)
    bigferr = dblarr(npix,ntell)

    bigwave = bigwave[inrange]
    bigflux[*,0] = tspec1.spec[inrange]
    bigferr[*,0] = tspec1.sigspec[inrange]
    
    icleanup, tspec1

    for itell = 1L, ntell-1L do begin

       tspec = rd1dspec(tlist[itell],datapath=datapath,_extra=extra)
       if (size(tspec,/type) ne 8L) then return

       bigflux[*,itell] = interpol(tspec.spec,tspec.wave,bigwave)
       bigferr[*,itell] = sqrt(interpol(tspec.sigspec^2,tspec.wave,bigwave))
       
       icleanup, tspec

    endfor
    
    bigivar = 1.0/bigferr^2

; initialize some b-spline parameters    
    
    tcorr_nbkpts = long(tcorr_bsorder) + 2L ; fix(npix/tcorr_bsorder) + 2L
    tcorr_bkpt = findgen(tcorr_nbkpts) * (max(bigwave) - min(bigwave)) / (tcorr_nbkpts-1) + min(bigwave)
;   print, tcorr_bkpt
    
; ---------------------------------------------------------------------------
; CONTMETHOD=1: continuum normalization in [MINTELLWAVE-MAXTELLWAVE] 
; ---------------------------------------------------------------------------

    if (contmethod eq 1L) then begin
    
; loop on each object in TLIST and divide out the continuum

       fitwave = bigwave # replicate(1.0,ntell) ; [NPIX,NTELL]
       fitflux = bigflux*0.0
       fitferr = bigferr*0.0
       
       for itell = 0L, ntell-1L do begin

          continuum = im_fitcontinuum(bigwave,bigflux[*,itell],$
            doplot=debug,/tellmask,/silent,_extra=extra)

          fitflux[*,itell] = bigflux[*,itell] / continuum
          fitferr[*,itell] = bigferr[*,itell] / continuum
          
          if keyword_set(debug) and (ntell ne 1L) then begin
             splog, 'Press any key to continuum.'
             cc = get_kbrd(1)
          endif

       endfor

       fitivar = 1.0 / fitferr^2
       
    endif

; ---------------------------------------------------------------------------
; TELLMETHOD=1: global telluric correction method
; ---------------------------------------------------------------------------

    if (tellmethod eq 1L) then begin

       splog, 'TELLMETHOD = '+string(tellmethod,format='(I0)')+' - Global telluric correction.'
       
; in this implementation the normalized telluric standards are stacked
; and a b-spline is fitted to the super-sampled spectrum

       tellset = bspline_iterfit(fitwave,fitflux,lower=tcorr_sigclip,$
         upper=tcorr_sigclip,invvar=fitivar,yfit=yfit,bkpt=tcorr_bkpt,$
         nord=tcorr_nord,/silent,_extra=extra)
       tellcorr = bspline_valu(tellwave,tellset)

    endif

; ---------------------------------------------------------------------------
; TELLMETHOD=2: local telluric correction method     
; ---------------------------------------------------------------------------

    if (tellmethod eq 2L) then begin

       splog, 'TELLMETHOD = '+string(tellmethod,format='(I0)')+' - Local telluric correction.'

       tellcorr = tellwave*0.0D; + 1.0
       
; if CONTMETHOD=2 or CONTMETHOD=3 then fit the continuum locally
; around each telluric feature

       if (contmethod eq 2L) or (contmethod eq 3L) then begin
          
          fitwave = bigwave # replicate(1.0,ntell)
          fitflux = dblarr(npix,ntell) + 1.0
          fitferr = dblarr(npix,ntell) + 1.0
          fitivar = dblarr(npix,ntell) + 1.0

       endif
          
; loop on each telluric band

       for iband = 0L, n_elements(tellbands)-1L do begin

; set masks to 1 where continuum and telluric bands are defined, then
; loop on each telluric standard

          for itell = 0L, ntell-1L do begin

             cmask = bytarr(npix)
             tmask = bytarr(npix)

             for i=0L, n_elements(tellbands[iband].cwave1)-1L do $
               cmask = cmask or (bigwave ge tellbands[iband].cwave1[i] and $
                 bigwave le tellbands[iband].cwave2[i])

             for i=0L, n_elements(tellbands[iband].twave1)-1L do $
               tmask = tmask or (bigwave ge tellbands[iband].twave1[i] and $
                 bigwave le tellbands[iband].twave2[i])

             indc = where(cmask ne 0,nindc)
             indt = where(tmask ne 0,nindt)

             if (nindc ne 0L) and (nindt ne 0L) then begin

; ---------------------------------------------------------------------------
; legendre fit to continuum             
; ---------------------------------------------------------------------------

                if (contmethod eq 2L) then begin
                   
; rescale X axis to be between 0 and 1 in the fitting region
                   
                   xmin = min(bigwave[indc],max=xmax)
                   xfit = (bigwave - xmin) / (xmax - xmin)

; iterate the fit using rejection at the 6-sigma level. at most,
; reject 25% of the points
                   
                   iiter = 0
                   maxiter = 5
                   qdone = 0
                   inmask = bigivar[indc,itell] NE 0
                   outmask = 0

                   while (not keyword_set(qdone) and iiter le maxiter) do begin

                      res = func_fit(xfit[indc],bigflux[indc,itell],tcorr_ncoeff,$
                        invvar=bigivar[indc,itell],function_name='flegendre')
                      continuum = flegendre(xfit,tcorr_ncoeff) # res

                      qdone = djs_reject(continuum[indc]-bigflux[indc,itell],bigivar[indc,itell],$
                        lower=tcorr_sigclip,upper=tcorr_siglip,maxrej=ceil(0.25*n_elements(indc)),$
                        inmask=inmask,outmask=outmask)

                      iiter = iiter + 1

                   endwhile 

                endif 

; ---------------------------------------------------------------------------                   
; spline fit to continuum
; ---------------------------------------------------------------------------                   

                if (contmethod eq 3L) then begin 

                   nbkpts = fix(npix/tcorr_pixspace) + 2L
                   bkpt = findgen(nbkpts) * (max(bigwave) - min(bigwave)) / (nbkpts-1) + min(bigwave)

                   istart = (where(bkpt gt min(bigwave[indc])))[0]
                   istart = (istart - 1) > 0
                   iend = (where(bkpt gt max(bigwave[indc])))[0]
                   if (iend eq -1L) then iend = nbkpts-1L

                   sset = bspline_iterfit(bigwave[indc],bigflux[indc,itell],$
                     invvar=bigivar[indc,itell],maxiter=10,upper=tcorr_sigclip,$
                     lower=tcorr_sigclip,nord=nord,bkpt=bkpt[istart:iend],$
                     maxrej=maxrej)

                   continuum = bspline_valu(bigwave,sset)

                endif

                if (contmethod eq 2L) or (contmethod eq 3L) then begin

                   fitflux[indt,itell] = bigflux[indt,itell] / continuum[indt]
                   fitferr[indt,itell] = bigferr[indt,itell] / continuum[indt]
                   fitivar[indt,itell] = bigivar[indt,itell] * (continuum[indt])^2

                   if keyword_set(debug) then begin
                      
                      djs_plot, bigwave, bigflux[*,itell], xsty=3, ysty=3, ps=10, $
                        xrange=[tellbands[iband].cwave1[0],max(tellbands[iband].cwave2)], $
                        charsize=1.3, charthick=postthick, xthick=postthick, ythick=postthick, $
                        title=repstr(repstr(repstr(tlist[itell],'_','-'),'.gz'),'.fits'), $
                        xtitle='Wavelength [\AA]', ytitle='Intensity'
                      djs_oplot, bigwave[indt], bigflux[indt,itell], color='cyan', ps=10
                      djs_oplot, bigwave[indc], bigflux[indc,itell], color='red', ps=10
                      djs_oplot, bigwave, continuum, color='green', thick=postthick
                      cc = get_kbrd(1)

                   endif        ; close debug IF statement

                endif           ; close CONTMETHOD IF statement 
                   
             endif              ; close TELLBAND "in range" IF statement 
             
          endfor                ; close object loop

; fit a b-spline to each absorption band, using all the spectra
; simultaneously 

          if (nindt ne 0L) then begin

             tellset = bspline_iterfit(fitwave[indt,*],fitflux[indt,*],$
               lower=tcorr_sigclip,upper=tcorr_sigclip,invvar=fitivar[indt,*],$
               yfit=yfit,everyn=1,/silent)

;            istart = (where(tcorr_bkpt gt min(bigwave[indt])))[0]
;            istart = (istart - 1L) > 0L
;
;            iend = (where(tcorr_bkpt gt max(bigwave[indt])))[0]
;            if (iend eq -1L) then iend = tcorr_nbkpts-1L
;
;            tellset = bspline_iterfit(fitwave[indt,*],fitflux[indt,*],$
;              lower=tcorr_sigclip,upper=tcorr_sigclip,invvar=fitivar[indt,*],$
;              yfit=yfit,bkpt=tcorr_bkpt[istart:iend],/silent)

             srt = sort(fitwave[indt,*])
             lintwave = (fitwave[indt,*])[srt]

             linterp, lintwave, yfit[srt], tellwave, tellcorr1, missing = 0.0
             tellcorr = tellcorr + tellcorr1

          endif
             
       endfor                   ; close TELLBAND loop

       zero = where(tellcorr eq 0.0,nzero)
       if (nzero ne 0L) then tellcorr[zero] = 1.0
          
    endif                       ; close TELMETHOD=2 loop
    
    tellcorr = tellcorr < 1.0
    
; generate the median telluric spectrum
    
;      if (ntell gt 1L) then $
;        tellcorr_median = djs_median(fitflux,2) < 1.0 else $
;        tellcorr_median = fitflux
;
;      linterp, fitwave, tellcorr_median, tellwave, tellcorr, missing=1.0

    if keyword_set(doplot) or keyword_set(write) then begin
       
       if keyword_set(write) then begin
          splog, 'Generating postscript output '+datapath+psname+'.'
          dfpsplot, datapath+psname, /color, /landscape
          colors = ['navy','red']
       endif else colors = ['cyan','red']

       plotsym, 0, 0.2, /fill
       djs_plot, fitwave, fitflux, ps=8, xsty=3, ysty=3, charsize=1.3, $
         charthick=postthick, xthick=postthick, ythick=postthick, $
         xtitle='Wavelength [\AA]', ytitle='Correction factor', $
         color=colors[0], title=plottitle
       djs_oplot, tellwave, tellcorr, color=colors[1], thick=postthick, ps=10

       tellcolors = ['purple','blue','dark green','pink']
       if (tellmethod eq 2L) then begin

          for iband = 0L, n_elements(tellbands)-1L do begin

             djs_oplot, tellbands[iband].twave1*[1,1], !y.crange, line=2, $
               thick=postthick, color=tellcolors[iband]
             djs_oplot, tellbands[iband].twave2*[1,1], !y.crange, line=2, $
               thick=postthick, color=tellcolors[iband]
             
          endfor
          
       endif
       
       if keyword_set(write) then dfpsclose

    endif 

; write out the telluric correction spectrum
    
    if keyword_set(write) then begin

       mkhdr, tellhead, tellwave
       sxdelpar, tellhead, 'COMMENT' & sxdelpar, tellhead, 'DATE'
       sxaddpar, tellhead, 'CRVAL1', tellwave[0], ' central wavelength of first pixel'
       sxaddpar, tellhead, 'CRPIX1', 1, ' starting pixel (1-indexed)'
       sxaddpar, tellhead, 'CD1_1', tellwave[1]-tellwave[0], ' dispersion [Angstrom/pixel]'
       sxaddpar, tellhead, 'CDELT1', tellwave[1]-tellwave[0], ' dispersion [Angstrom/pixel]'
       sxaddpar, tellhead, 'CTYPE1', 'LINEAR'
;      sxaddpar, tellhead, 'DC-FLAG', 0, ' log-linear flag'
       sxaddhist, "'Telluric correction spectrum generated "+im_today()+"'", tellhead

       splog, 'Writing '+datapath+tellfits+'.'
       mwrfits, float(tellcorr), datapath+tellfits, tellhead, /create

    endif

return
end
    
