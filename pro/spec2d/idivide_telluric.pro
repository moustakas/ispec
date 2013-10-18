;+
; NAME:
;       IDIVIDE_TELLURIC
;
; PURPOSE:
;       Divide a spectrum by an appropriately scaled and shifted
;       telluric absorption-line spectrum.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;       specfile - 
;       tellfile -
;
; OPTIONAL INPUTS:
;       datapath - path name to SPECFILE
;       tellpath - path name to TELLFILE
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 12, U of A - based on an earlier code 
;-

function scale_function, x, p, xwave=xwave, yspec=yspec
; x - telluric spectrum
; xwave - wavelength of the telluric spectrum
; p[0] - wavelength shift
; p[1] - amplitude scale factor

    z = p[0]
    scale = p[1]

    shiftwave = xwave*(1.0+z)
    scaleflux = scale * x
    
    shiftmodel = interpol(scaleflux,xwave,shiftwave)

    model = yspec / shiftmodel - 1.0
    
return, model
end

pro idivide_telluric, specfile, tellfile, datapath=datapath, tellpath=tellpath, $
  vmaxshift=vmaxshift, method=method, suffix=suffix, postscript=postscript, $
  write=write

    light = 2.99792458D5        ; speed of light [km/s]
    
    nspec = n_elements(specfile)
    if (nspec eq 0L) then begin
       print, 'Syntax - idivide_telluric, specfile, tellfile, datapath=, $'
       print, '   tellpath='
       return
    endif

    ntell = n_elements(tellfile)
    if (ntell ne 1L) then begin
       splog, 'TELLFILE must be a scalar file name.'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(tellpath) eq 0L) then tellpath = cwd()
    if (n_elements(vmaxshift) eq 0L) then vmaxshift = 1500.0 ; [km/s]
    if (n_elements(method) eq 0L) then method = 2L

    mch = machar()
    eps = mch.eps

    tellbands1 = { TELLBAND, $
      twave1: 6850., $
      twave2: 6960., $
      cwave1: [6600., 6950., 0], $
      cwave2: [6860., 7200., 0] }
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
    tellbands = [tellbands1, tellbands2, tellbands3, tellbands4]
    nbands = n_elements(tellbands)
    
; read the telluric absorption spectrum

    if (file_test(tellpath+tellfile,/regular) eq 0L) then begin
       splog, 'Telluric spectrum '+tellpath+tellfile+' not found.'
       return
    endif

    telluric_flux = mrdfits(tellpath+tellfile,0,telluric_header,/silent)
    telluric_wave = make_wave(telluric_header)

    if (n_elements(suffix) eq 0L) then suffix = '' else suffix = '_'+suffix
    psname = 'qaplot_idivide_telluric'+suffix+'.ps'
    
    if keyword_set(write) then postscript = 1L
    
    if keyword_set(postscript) or keyword_set(write) then begin
       splog, 'Generating postscript output '+datapath+psname+'.'
       dfpsplot, datapath+psname, /color, /landscape
       postthick = 5.0
    endif else postthick = 2.0

; loop on each spectrum

    for ispec = 0L, nspec-1L do begin

;      psname = repstr(repstr(specfile[ispec],'.gz',''),'.fits','.ps')
       
       scube = rd1dspec(specfile[ispec],datapath=datapath)
       wave = scube.wave
       flux = scube.spec
       ferr = scube.sigspec
       ivar = 1.0/ferr^2.0
       npix = scube.npix
       
       corrflux = flux          ; corrected spectrum
       tellfit = flux*0.0+1.0   ; best-fitting telluric spectrum

       file = repstr(repstr(specfile[ispec],'.gz',''),'.fits','')
       plottitle = scube.object+' ['+file+']'

; interpolate the telluric spectrum onto the wavelength spacing of the
; object spectrum, replacing missing values with 1.0

       linterp, telluric_wave, telluric_flux, wave, tflux, missing=1.0

       tmask = bytarr(npix)
       for iband = 0L, nbands-1L do for i = 0L, n_elements(tellbands[iband].twave1)-1L do $
         tmask = tmask or (wave ge tellbands[iband].twave1[i] and $
           wave le tellbands[iband].twave2[i])
       bigindx = where(tmask ne 0L,nbigindx)

       if (nbigindx eq 0L) then begin
          splog, 'Something is wrong with the telluric spectrum.'
          stop
       endif

       trange = where((wave gt min(telluric_wave)) and (wave lt max(telluric_wave)),nrange)

; ---------------------------------------------------------------------------       
; fit all bands simultaneously       
; ---------------------------------------------------------------------------       
       
       if (method eq 1L) then begin

          parinfo = {$
            value:   1.0D, $
            limits:  [0.5D,1.5D], $
            limited: [1L,1L]}
          parinfo = replicate(parinfo,2L) ; shift + scale factor

; starting wavelength shift          
          
          parinfo[0].value = 0.0
          parinfo[0].limits = [-vmaxshift,+vmaxshift]/light
          parinfo[0].limited = [1L,1L]
          
          plottitle = plottitle+' - METHOD 1'

          functargs = {xwave: wave[bigindx], yspec: flux[bigindx]}

          result = mpfitfun('scale_function',tflux[bigindx],flux[bigindx],$
            weights=weights,parinfo=parinfo,functargs=functargs,status=status,$
            covar=covar,perror=perror,yfit=yfit,nfev=nfev,niter=niter,/quiet)
          splog, 'MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
          if (status EQ 5) then splog, 'Warning: Maximum number of iterations reached: ', niter

          corrflux[bigindx] = yfit

          djs_plot, wave, flux, ps=10, xsty=3, ysty=3, xrange=[6500,8400], color='red'
          djs_oplot, wave, corrflux, ps=10, color='green', thick=2.0
          djs_oplot, wave, flux/tflux, ps=10, color='cyan', thick=1.0

stop          
          
       endif
       
       if (method eq 2L) then begin ; fit each band individually

          parinfo = {$
            value:   0.0D, $
            limits:  [0.0D,0.0D], $
            limited: [0L,0L]}
          parinfo = replicate(parinfo,4L)

          parinfo[0].value = 1E-12 ; wavelength shift
          parinfo[0].limits = [-vmaxshift,+vmaxshift]/light ; starting wavelength shift
          parinfo[0].limited = [1,1]
          
          parinfo[1].value = 1.0 ; scale factor
          parinfo[1].limits = [0.5,1.5]
          parinfo[1].limited = [0,0]
          
          parinfo[2].value = 1.0 ; intercept

          parinfo[3].value = 0.1 ; slope

          plottitle = plottitle+' - METHOD 2'

          result = dblarr(2*nbands)
          
          for iband = 0L, nbands-1L do begin
          
             tmask = bytarr(npix)

             for i = 0L, n_elements(tellbands[iband].twave1)-1L do $
               tmask = tmask or (wave ge tellbands[iband].twave1[i] and $
                 wave le tellbands[iband].twave2[i])

             indt = where(tmask ne 0,nindt)
             
             functargs = {xwave: wave[indt], yspec: flux[indt]}

             bandresult = mpfitfun('scale_function',tflux[indt],flux[indt],$
               weights=weights,parinfo=parinfo,functargs=functargs,status=status,$
               covar=covar,perror=perror,yfit=yfit,nfev=nfev,niter=niter,/quiet)
             splog, 'MPFIT nfev=', nfev, ' niter=', niter, ' status=', status
             if (status EQ 5) then splog, 'Warning: Maximum number of iterations reached: ', niter

             print, bandresult
             
             result[2*iband] = bandresult[0]
             result[2*iband+1] = bandresult[1]

;            yfit = yfit / djs_mean(yfit)                       

             tellfit[indt] = yfit < 1.0
             corrflux[indt] = flux[indt]/yfit

;            plot, wave[indt], flux[indt]/median(flux[indt]), ysty=3, ps=10
;            djs_oplot, wave[indt], yfit, ps=10, color='green'
;            cc = get_kbrd(1)
             
          endfor
             
       endif

;      niceprint, result
       
; generate the QA plots

       !p.multi = [0,2,2]

       xmargin = 75.0           ; extra plot range in Angstroms

       for iband = 0L, nbands-1L do begin

          xlo = min(tellbands[iband].twave1) - xmargin
          xhi = max(tellbands[iband].twave2) + xmargin
          
          plotindx = where(wave ge xlo and wave LE xhi)
          norm = max(flux[plotindx])>max(corrflux[plotindx])

          yrange = [min(flux[plotindx])<min(corrflux[plotindx]),$
            max(flux[plotindx])>max(corrflux[plotindx])]*[0.95,1.05]

          if (plotindx[0] NE -1) then begin

             djs_plot, wave[plotindx], flux[plotindx], ps=10, xsty=3, ysty=3, $
               xrange=[xlo,xhi], ymargin=[4,4], xtitle='\lambda [\AA]', $
               ytitle='Flux', xthick=postthick, color='grey', yrange=yrange, $
               ythick=postthick, charsize=1.2, charthick=postthick, $
               thick=postthick+2
             djs_oplot, wave[plotindx], corrflux[plotindx], ps=10, color='red', thick=postthick
             djs_oplot, wave[plotindx], tellfit[plotindx], ps=10, color='dark green', thick=postthick

             djs_oplot, tellbands[iband].twave1*[1,1], !y.crange, line=2, thick=postthick
             djs_oplot, tellbands[iband].twave2*[1,1], !y.crange, line=2, thick=postthick
             
             if (iband eq 0) then xyouts, 0.5, 0.95, plottitle, $
               charsize=1.5, align=0.5, charthick=postthick, /normal

          endif

       endfor

       !p.multi = 0

       if (not keyword_set(postscript)) then begin
          splog, 'Press any key to continue.'
          cc = get_kbrd(1)
       endif

    endfor

    if keyword_set(postscript) or keyword_set(write) then begin
       dfpsclose
;      spawn, ['gzip -f '+datapath+psname], /sh
    endif
       
return
end    
