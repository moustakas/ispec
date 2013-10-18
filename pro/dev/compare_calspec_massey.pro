pro compare_calspec_massey
; jm04mar21uofa 
; compare the Massey et al. (1988) and CALSPEC absolute calibrations

    light = 2.99793D18 ; [A/s]
    
    cal = read_standards_table(datapath=calpath,/calspec)
    mas = read_standards_table(datapath=maspath,/spec50cal)

; find all stars in common

    doit = match_string(cal.star,mas.star,/exact,index=indx,/silent)
    good = where(indx ne -1L,ngood)
    niceprint, cal[good].star, mas[indx[good]].star

    res = {mean: 0.0, median: 0.0, sigma: 0.0}
    res = replicate(res,ngood)
    
    for i = 0L, ngood-1L do begin

       calfile = strtrim(cal[good[i]].file,2)
       masfile = strtrim(mas[indx[good[i]]].file,2)
       
       readcol, calpath+calfile, calwave, calmag, /silent
       calflux = 10^(-0.4*(calmag+48.59)) * light / calwave / calwave ; [erg/s/cm2/A]

       readcol, maspath+masfile, maswave, masmag, /silent
       masflux = 10^(-0.4*(masmag+48.59)) * light / maswave / maswave ; [erg/s/cm2/A]

       imasflux = interpol(masflux,maswave,calwave)
       
       resid = (calflux-imasflux)/calflux

       xrange = [3500,7500]
       get_element, calwave, xrange, ww
       resid = resid[ww[0]:ww[1]]
       wave = calwave[ww[0]:ww[1]]
       
       stats = im_stats(resid)
       yrange = (abs(stats.min)>stats.max)*[-1,1]
       title = strupcase(cal[good[i]].star)
       strstats = strtrim(string(stats.mean,format='(F12.5)'),2)+' +/- '+$
         strtrim(string(stats.sigma,format='(F12.5)'),2)

       res[i].mean = stats.mean
       res[i].sigma = stats.sigma
       res[i].median = stats.median
       
       djs_plot, wave, resid, xsty=3, ysty=3, xrange=xrange, ps=10, $
         charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, yrange=yrange
       djs_oplot, !x.crange, [0,0], line=0, thick=2.0
       legend, [title,strstats], /left, /top, box=0, charsize=1.5, charthick=2.0
       
       cc = get_kbrd(1)
       
    endfor

    splog, 'Mean offset : '+strtrim(string(djs_mean(res.mean),format='(F12.5)'),2)
    splog, 'Mean scatter: '+strtrim(string(djs_mean(res.sigma),format='(F12.5)'),2)
    
return
end
    
