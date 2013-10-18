pro rkexspec, image, vmap, mask, header, wave

    imsize = size(image,/dimension)
    xsize = imsize[0]
    ysize = imsize[1]
    iminfo = sstretch(image,npix=1000,/silent)
    
    state = {base_id: 0L,                $ ; rkexspec base id
             base_min_size: [600L,400L], $
             draw_window_id: 0L,         $
             draw_widget_id: 0L,         $
             location_bar_id: 0L,        $
             done_id: 0L,                $ ; done button
             xrange: [0.0,1.0],          $
             yrange: [0.0,1.0],          $
             position: [0.15,0.15,0.95,0.95], $
             draw_window_size: [600L,512L], $
             mouse: [0L,0L], $
             mphys: [0L,0L], $
             base_pad: [0L,0L]}

; define the widgets

    base = widget_base(title='rkexspec',/row,/base_align_right,uvalue='rkexspec_base')
;                      scr_xsize=800L,scr_ysize=500L)
    label_base = widget_base(base,/column,/base_align_left,uvalue='label_base',group_leader=base)
;   done_id = widget_button(label_base,value='Done',uvalue='done')
    draw_base = widget_base(base,/column,/base_align_right,uvalue='image_base',group_leader=base)
    imdraw = widget_draw(draw_base,uvalue='image_window')
;   spec = widget_draw(draw_base,uvalue='spectrum_window');,xsize=xsize+50L,ysize=ysize)
;   vspec = widget_draw(draw_base,uvalue='vspectrum_window',xsize=xsize+50L,ysize=ysize+50L)

    f1 = fsc_inputfield(label_base,title='Aperture',/doublevalue,event_pro='rkextract_event', $
                        uvalue='aperture',value=20.0,/cr_only,/positive,decimal=2,$
                        labelsize=120,xsize=5,/frame)
    f2 = fsc_inputfield(label_base,title='Lower',/doublevalue,event_pro='rkextract_event', $
                        uvalue='aperture',value=20.0,/cr_only,/negative,decimal=2,$
                        labelsize=120,xsize=5,/frame)
    f3 = fsc_inputfield(label_base,title='Upper',/doublevalue,event_pro='rkextract_event', $
                        uvalue='aperture',value=20.0,/cr_only,/positive,decimal=2,$
                        labelsize=120,xsize=5,/frame)

    widget_control, base, /realize
;   xmanager, 'rkexspec', base, event_handler='rkexspec_event', group_leader=base

    widget_control, imdraw, get_value = wsetmp & wset, wsetmp
    plot, image[595,*], xsty=3, ysty=3

    
    state.base_id = base
;   state.done_id = done_id
    state.draw_window_id = draw_base

; display the two-dimensional spectrum in the top window
    
    widget_control, imdraw, get_value = wsetmp
    wset, wsetmp
    tv, bytscl(image,min=iminfo.min,max=iminfo.max)
    
; display the extracted spectrum and variance spectrum in the next two
; windows

    widget_control, spec, get_value = wsetmp
    wset, wsetmp
    plot, wave[*,63], image[*,63], position=[0.08,0.05,1,0.95], xsty=8, ysty=11, $
      xthick=2.0, ythick=2.0, yminor=3, charthick=2.0, ytitle='Counts', $
      xtickname=replicate(' ',10), xrange=minmax(wave)
    
;    widget_control, vspec, get_value = wsetmp
;    wset, wsetmp
;    plot, wave[*,63], sqrt(vmap[*,63]), position=[0.08,0.2,1,0.95], xsty=8, ysty=11, $
;      xthick=2.0, ythick=2.0, yminor=3, charthick=2.0, ytitle='Counts'

stop


return
end
