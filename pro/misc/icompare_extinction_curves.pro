pro icompare_extinction_curves, ps=ps
; jm07may10nyu - compare the various extinction curves

    path = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    file = ['kpnoextinct.dat','ctioextinct.dat','lasillaextinct.dat']

    if keyword_set(ps) then begin
       dfpsplot, 'icompare_extinction_curves.ps', /square, /color
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.4, /square
       postthick = 2.0
    endelse

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal
    
    xrange = [4000,9200] ; [3200,10400]
    yrange = [0.005,1.0]

    xtitle = 'Wavelength (\AA)'
    ytitle = 'Extinction (mag/AM)'

    xaxis = findgen((xrange[1]-xrange[0])/0.05+1)*0.05+xrange[0]
    
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), position=pos, /ylog

; overplot the telluric bands

    polyfill, [6850,6850,6960,6960], [yrange[0],yrange[1],yrange[1],yrange[0]], $
      color=fsc_color('light grey',10), noclip=0
    polyfill, [7150,7150,7350,7350], [yrange[0],yrange[1],yrange[1],yrange[0]], $
      color=fsc_color('light grey',10), noclip=0
    polyfill, [7560,7560,7720,7720], [yrange[0],yrange[1],yrange[1],yrange[0]], $
      color=fsc_color('light grey',10), noclip=0
    polyfill, [8105,8105,8240,8240], [yrange[0],yrange[1],yrange[1],yrange[0]], $
      color=fsc_color('light grey',10), noclip=0

    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xthick=postthick, ythick=postthick, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xrange=xrange, yrange=yrange, position=pos, /ylog
    
; KPNO
    im_symbols, 108, psize=0.9, fill=1, thick=postthick, color=djs_icolor('red')
    readcol, path+file[0], w, e, /silent, format='F,F'
    djs_oplot, w, e, ps=-8, color='red', thick=postthick, line=0
; CTIO
    im_symbols, 106, psize=1.2, fill=1, thick=postthick, color=djs_icolor('dark green')
    readcol, path+file[1], w, e, /silent, format='F,F'
    djs_oplot, w, e, ps=-8, color='dark green', thick=postthick, line=2
; LA SILLA
    im_symbols, 115, psize=1.5, fill=1, thick=postthick, color=djs_icolor('blue')
    readcol, path+file[2], w, e, /silent, format='F,F'
    djs_oplot, w, e, ps=-8, color='blue', thick=postthick, line=1

; legend

    im_legend, ['KPNO','CTIO','La Silla'], /left, /bottom, psym=[108,106,115], /fill, $
      color=djs_icolor(['red','dark green','blue']), charsize=2.0, charthick=postthick, $
      box=0, symsize=1.5, line=[0,2,1]

    if keyword_set(ps) then dfpsclose

return
end
    
