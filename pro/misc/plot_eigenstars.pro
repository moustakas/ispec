pro plot_eigenstars, ps=ps
;jm02jan25uofa
    
    eigen = mrdfits('spEigenStar-51984.fits',0,h)
    nstar = (size(eigen,/dimension))[1]
    
    coeff1 = sxpar(h,'coeff1')
    coeff0 = sxpar(h,'coeff0')
    naxis = sxpar(h,'naxis1')
    wave = 10.0^(coeff0+coeff1*findgen(naxis))
    minwave = min(wave)
    maxwave = max(wave)

; first panel

    if keyword_set(ps) then begin
       ps_open, 'sdss_eigenstars_panel1', /ps_fonts, /portrait
       device, /inches, /times
    endif else window, 0, xs=550, ys=550

; left side
    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3], xsty=3, ysty=3, $
      /nodata, position=[0.05,0.7,0.5,1.0], xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,0]
    xyouts, 7500, 2.5, sxpar(h,'NAME0'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3], xsty=3, ysty=3, $
      /nodata, /noerase, position=[0.05,0.4,0.5,0.7], xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Relative Flux'
    oplot, wave, eigen[*,1]
    xyouts, 7500, 2.5, sxpar(h,'NAME1'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3.5], xsty=3, ysty=3, $
      /nodata, /noerase, position=[0.05,0.1,0.5,0.4], $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Wavelength (\AA)'
    oplot, wave, eigen[*,2]
    xyouts, 7500, 2.5, sxpar(h,'NAME2'), /data, align=0.5, charsize=1.8, charthick=2.0

; right side
    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,2.5], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.7,0.95,1.0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10),       charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    axis, yaxis=1, ysty=3, yrange=[0,2.5], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,3]
    xyouts, 7500, 2, sxpar(h,'NAME3'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.4,0.95,0.7], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10),       charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    axis, yaxis=1, ysty=3, yrange=[0,3], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,4]
    xyouts, 7500, 2.5, sxpar(h,'NAME4'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,2.0], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.1,0.95,0.4], ytickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Wavelength (\AA)'
    axis, yaxis=1, ysty=3, yrange=[0,2.0], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,5]
    xyouts, 7500, 1.5, sxpar(h,'NAME5'), /data, align=0.5, charsize=1.8, charthick=2.0

    if keyword_set(ps) then ps_close
    
; second panel
    
    if keyword_set(ps) then begin
       ps_open, 'sdss_eigenstars_panel2', /ps_fonts, /portrait
       device, /inches, /times
    endif else window, 2, xs=550, ys=550

; left side
    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,2.0], xsty=3, ysty=3, $
      /nodata, position=[0.05,0.7,0.5,1.0], xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,6]
    xyouts, 7500, 1.5, sxpar(h,'NAME6'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,8], xsty=3, ysty=3, $
      /nodata, /noerase, position=[0.05,0.4,0.5,0.7], xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Relative Flux'
    oplot, wave, eigen[*,7]
    xyouts, 7500, 7, sxpar(h,'NAME7'), /data, align=0.5, charsize=1.8, charthick=2.0

    djs_plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,2.5], xsty=3, ysty=3, $
      /nodata, /noerase, position=[0.05,0.1,0.5,0.4], $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Wavelength (\AA)'
    oplot, wave, eigen[*,8]
    xyouts, 7500, 2.0, sxpar(h,'NAME8'), /data, align=0.5, charsize=1.8, charthick=2.0

; right side
    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.7,0.95,1.0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10),       charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    axis, yaxis=1, ysty=3, yrange=[0,3], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,9]
    xyouts, 7500, 2.5, sxpar(h,'NAME9'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,3], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.4,0.95,0.7], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10),       charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    axis, yaxis=1, ysty=3, yrange=[0,3], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,10]
    xyouts, 7500, 2.5, sxpar(h,'NAME10'), /data, align=0.5, charsize=1.8, charthick=2.0

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,2.5], xsty=3, ysty=11, $
      /nodata, /noerase, position=[0.5,0.1,0.95,0.4], ytickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Wavelength (\AA)'
    axis, yaxis=1, ysty=3, yrange=[0,2.5], charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0
    oplot, wave, eigen[*,11]
    xyouts, 7500, 2., sxpar(h,'NAME11'), /data, align=0.5, charsize=1.8, charthick=2.0

    if keyword_set(ps) then ps_close

stop

return
end
