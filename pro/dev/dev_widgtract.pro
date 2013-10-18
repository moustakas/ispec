pro extract_event, event

    widget_control, event.top, get_uvalue = apinfo


stop
    
    case event_name of

       'FSC_FIELD': begin
;          'center': apinfo.center = *event.value
;          'aperture':
;          'lower':
;          'upper':
;          'skyaperture':
;          'skylower':
;          'skyupper':
       end

       'WIDGET_BUTTON': begin
;         'extract':
;         'done': ; widget_control, event.top, /destroy
       end

    endcase

return
end

pro widgtract, image, vmap, apinfo
; jm01oct4uofa
; interactive widget-based aperture extraction routine
; user-interactive widget-based spectrum extraction routine

;    imsize = size(image,/dimension)
;    ncols = imsize[0]
;    nrows = imsize[1]
;
;; compute the mean spatial profile of the image
;    
;    mprof = total(image,1)/ncols
;    rowaxis = findgen(nrows)
;    
;; spatial profile widget
;
;    base_profile = widget_base(title='Spatial Profile',/column)
;    draw_profile = widget_draw(base_profile,scr_xsize=550L,scr_ysize=550L)
;    widget_control, base_profile, /realize
;    widget_control, draw_profile, get_value=windx_profile
;
;    wset, windx_profile
;; ----------------------------------------------------------------------
;; spatial profile
;; ----------------------------------------------------------------------
;    plot, rowaxis, mprof, xsty=3, ysty=3, thick=2.0, ps=10, $
;      charthick=2.0, charsize=1.5, xtitle='Row', ytitle='Counts', $
;      yrange=[min(mprof),max(mprof)*1.2], xrange=[1,nrows], title=title
;; ----------------------------------------------------------------------
;; aperture window
;; ----------------------------------------------------------------------
;    ymax = max(mprof)
;    oplot, [apinfo.lower,apinfo.center], [ymax,ymax]*1.1, line=0, thick=2.0
;    oplot, [apinfo.center,apinfo.upper], [ymax,ymax]*1.1, line=0, thick=2.0
;    oplot, [apinfo.center,apinfo.center], [ymax*1.05,ymax*1.15], line=0, thick=2.0
;; ----------------------------------------------------------------------
;; sky windows
;; ----------------------------------------------------------------------
;    oplot, [apinfo.skylower-apinfo.skyaperture/2.0,apinfo.skylower+apinfo.skyaperture/2.0], $
;      [ymax,ymax]*1.1, line=3, thick=2.0
;    oplot, [apinfo.skyupper-apinfo.skyaperture/2.0,apinfo.skyupper+apinfo.skyaperture/2.0], $
;      [ymax,ymax]*1.1, line=3, thick=2.0
;
;; extraction parameter widget
;    
;    base_params = widget_base(title='Extraction',/column);,scr_xsize=300L)
;    f1 = fsc_inputfield(base_params,title='Center:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.center,decimal=2.0,/float,uvalue='center',event_pro='extract_event')
;    f2 = fsc_inputfield(base_params,title='Aperture:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.aperture,decimal=2.0,/float,uvalue='aperture',event_pro='extract_event')
;    f3 = fsc_inputfield(base_params,title='Lower:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.lower,decimal=2.0,/float,uvalue='lower',event_pro='extract_event')
;    f4 = fsc_inputfield(base_params,title='Upper:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.upper,decimal=2.0,/float,uvalue='upper',event_pro='extract_event')
;    f5 = fsc_inputfield(base_params,title='Sky aperture:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.skyaperture,decimal=2.0,/float,uvalue='skyaperture',event_pro='extract_event')
;    f6 = fsc_inputfield(base_params,title='Sky lower:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.skylower,decimal=2.0,/float,uvalue='skylower',event_pro='extract_event')
;    f7 = fsc_inputfield(base_params,title='Sky upper:',labelsize=80,/positive,/cr_only,xsize=8,$
;                        value=apinfo.skyupper,decimal=2.0,/float,uvalue='skyupper',event_pro='extract_event')
;
;    extract_button = widget_button(base_params,value='Extract',uvalue='extract',xsize=8,units=0)
;    done_button = widget_button(base_params,value='Done',uvalue='done',xsize=8,units=0)
;    widget_control, base_params, /realize, set_uvalue=apinfo
;
;    xmanager, 'rk_extract_user', base_params, event_handler='extract_event';, /no_block
    
; plot the extracted spectrum.  plot the average spatial profile

;    base_spec = widget_base(title='Extracted Spectrum',/column)
;    draw_spec = widget_draw(base_spec,scr_xsize=700,scr_ysize=400)
;    widget_control, base_spec, /realize
;    widget_control, draw_spec, get_value=windx_spec
;
;    wset, windx_spec
;    plot, wave[*,midrow], image[*,midrow], xsty=3, ysty=3, thick=2.0, $
;      charthick=2.0, charsize=1.5, ytitle='Counts', title=title, $
;      xtickname=replicate(' ',10), position=[0.1,0.52,0.95,0.92], ps=10
;;   legend, ['Spectrum'], /left, /top, box=0, charsize=1.5, charthick=2.0
;    plot, wave[*,midrow], sqrt(vmap[*,midrow]), xsty=3, ysty=3, thick=2.0, $
;      charthick=2.0, charsize=1.5, ytitle='Counts', position=[0.1,0.12,0.95,0.52], $
;      xtitle='Wavelength', /noerase, ps=10
;;   legend, ['Uncertainty'], /left, /top, box=0, charsize=1.5, charthick=2.0

; once "extract" is selected print out a useful message giving the
; extraction parameters and write to the log if specified    

return
end
