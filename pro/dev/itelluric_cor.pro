;+
; NAME:
;       ITELLURIC_COR()
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
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
; COMMON BLOCKS:
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
;       J. Moustakas, 2003 Dec 30, U of A
;-

function itelluric_cor, telluric_wave, telluric_spec, wave, spec, doplot=doplot
    
; interpolate the telluric spectrum onto the wavelength spacing of the
; object spectrum, replacing missing values with 1.0    

    linterp, telluric_wave, telluric_spec, wave, tspec, missing=1.0

; index the telluric absorption wavelengths    
    
    tellmask = telluric_mask(wave,bband=bband,bfreeband=bfreeband)
    const = djs_median(spec[bband]/tspec[bband])
    normspec = spec/const
    
;   plot, wave, normspec, xrange=[6760,6900], xsty=3, ysty=3, ps=10
;   djs_oplot, wave[bband], tspec[bband], ps=10, color='yellow'
    
; initialize the grid of constants; allow the depth of the telluric
; absorption to vary by +/- 10%

    cmin = 0.90
    cmax = 1.1
    dc = 0.02
    ncgrid = round((cmax-cmin)/dc)+1.0
    cgrid = findgen(ncgrid)*dc+cmin

    chi2array = fltarr(ncgrid)
    zshift = fltarr(ncgrid)
    
    for i = 0L, ncgrid-1L do begin

; amplitude and velocity shift
       
       scale_tspec = tspec
       scale_tspec[bband] = cgrid[i]*tspec[bband]
       
       zans = im_ztweak(scale_tspec[bband],wave[bband],normspec[bband],wave[bband],$
         vmin=-1500.0,vmax=1500.0,minwave=6800.0,maxwave=6950.0,doplot=0,/silent)
;      cc = get_kbrd(1)
       
       zshift[i] = zans.zshift
       newwave = wave*(1.0+zshift[i])
       linterp, newwave, scale_tspec, wave, shift_tspec, missing=1.0

       chi2array[i] = total((normspec[bband]-shift_tspec[bband])^2.0)

;      wset, 0
;      plot, wave, normspec, xrange=[6760,6900], xsty=3, ysty=3, ps=10
;      djs_oplot, wave, normspec/shift_tspec, ps=10, color='cyan'
;      djs_oplot, wave, scale_tspec, ps=10, color='yellow'
;      djs_oplot, wave, shift_tspec, ps=10, color='red'
;      cc = get_kbrd(1)

    endfor

    chi2array = chi2array/max(chi2array)
    findchi2min, cgrid, chi2array, minchi2, cgrid_best;, /doplot
    zshift_best = interpol(zshift,cgrid,cgrid_best)

    scale_tspec = tspec
    scale_tspec[bband] = cgrid_best*tspec[bband]

    linterp, wave*(1.0+zshift_best), scale_tspec, wave, best_tspec, missing=1.0
    corspec = spec/best_tspec

    if keyword_set(doplot) then begin

       window, 5, xsize=300, ysize=300
       plot, wave, spec/const, xrange=[6800,6950<max(wave)], xsty=3, ysty=3, ps=10, $
         xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $
         xtitle='Wavelength [\AA]', ytitle='Normalized Flux'
;      djs_oplot, wave[bband], spec[bband]/const, ps=10, color='red'
       djs_oplot, wave, best_tspec, ps=10, color='red'
;      djs_oplot, wave[bband], tspec[bband]*const, ps=10, color='green'
       djs_oplot, wave, corspec/const, ps=10, color='yellow'

    endif
    
return, corspec
end    
