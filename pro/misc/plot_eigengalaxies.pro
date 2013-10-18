pro plot_eigengalaxies, ps=ps
;jm02jan25uofa
    
    eigen = mrdfits('spEigenGal-51984.fits',0,h)
    ngal = (size(eigen,/dimension))[1]
    
    coeff1 = sxpar(h,'coeff1')
    coeff0 = sxpar(h,'coeff0')
    naxis = sxpar(h,'naxis1')
    wave = 10.0^(coeff0+coeff1*findgen(naxis))
    minwave = min(wave)
    maxwave = max(wave)

; first panel

    if keyword_set(ps) then begin
       ps_open, 'sdss_eigengalaxies', /ps_fonts;, /portrait
       device, /inches, /times
    endif else window, 0, xs=550, ys=550

; left side
    plot, [0], [0], xrange=[minwave,maxwave], yrange=[0,1E4], xsty=3, ysty=3, $
      /nodata, position=[0.1,0.55,0.5,1.0], xtickname=replicate(' ',10), $
      charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Relative Flux'
    oplot, wave, eigen[*,0]

    plot, [0], [0], xrange=[minwave,maxwave], yrange=[-4000,0], xsty=3, ysty=3, $
      /nodata, /noerase, position=[0.1,0.1,0.5,0.55],  xtitle='Wavelength (\AA)', $
      charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Relative Flux'
    oplot, wave, eigen[*,1]

; right side
    plot, [0], [0], xrange=[minwave,maxwave], xsty=3, ysty=11, yrange=[-3000,1000], $
      /nodata, /noerase, position=[0.5,0.55,0.9,1.0], xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10),       charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0
    axis, yaxis=1, ysty=3, yrange=[-3000,1000], charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, $
      ytitle='Relative Flux'
    oplot, wave, eigen[*,2]

    plot, [0], [0], xrange=[minwave,maxwave], xsty=3, ysty=11, $
      /noerase, /nodata, position=[0.5,0.1,0.9,0.55], yrange=[-1200,200], $
      ytickname=replicate(' ',10),       charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, $
      xtitle='Wavelength (\AA)'
    axis, yaxis=1, ysty=3, yrange=[-1200,200], charsize=1.3, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Relative Flux'
    oplot, wave, eigen[*,3]

    if keyword_set(ps) then ps_close

stop

return
end
