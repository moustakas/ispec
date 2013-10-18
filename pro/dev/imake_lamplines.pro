;+
; NAME:
;       IMAKE_LAMPLINES
;
; PURPOSE:
;       Generate a arc lamp line list file.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 July 11, U of A, written based on Andy
;          Marble's IARCFIT_WRAPPER
;-

pro imake_lamplines, arcfile, input_lampfile, result, datapath=datapath, $
  output_lampfile=output_lampfile, lambda0=lambda0, dispersion=dispersion, $
  nolines=nolines, _extra=extra

    narcfile = n_elements(arcfile)
    if (narcfile ne 1L) then begin
       splog, 'ARCFILE must be a scalar.'
       return
    endif
    
    ninput_lampfile = n_elements(input_lampfile)
    if (ninput_lampfile ne 1L) then begin
       splog, 'INPUT_LAMPFILE must be a scalar.'
       return
    endif

    if (n_elements(lambda0) eq 0L) then begin
       splog, 'LAMBDA0 must be defined.'
       return
    endif

    if (n_elements(dispersion) eq 0L) then begin
       splog, 'DISPERSION must be defined.'
       return
    endif

    if (n_elements(output_lampfile) eq 0L) then output_lampfile = 'lamplines.dat'
    if (n_elements(datapath) eq 0L) then datapath = cwd()
    
; read lamp line information

    if (file_test(datapath+input_lampfile,/regular) eq 0L) then begin
       splog, 'Input lampfile '+datapath+input_lampfile+' not found.'
       return
    endif
    
    readcol, datapath+input_lampfile, format='A,A,A', arcwavetxt, quality, $
      element, comment='#', /silent
    arcwave = double(arcwavetxt)
    nwave = n_elements(arcwave)
    
; read in first lamp image
    
    arclamp = rd2dspec(arcfile,datapath=datapath)
    arc = arclamp.image
    ncol = arclamp.naxis1
    colaxis = lindgen(ncol)
    flux = total(arc,2)
    flux = 1000.0*flux/max(flux)
    
    result1 = {$
      linepix:   0.0D, $
      intensity: 0.0D, $
      linewave:  0.0D, $
      element:    '',  $
      quality:    '',  $
      lineindx:    0L}

    disp0 = dispersion
    cntr = 0L

    window, 0, xs=800, ys=500
    
; loop over all supplied lines
    
    splog, 'Select negative intensity to skip a line.'
    
    for i=0L, nwave-1L do begin

       if (arcwave[i] ge lambda0) and (arcwave[i] le lambda0+dispersion*ncol) then begin

          splog, 'Please identify '+element[i]+' at '+string(arcwave[i],format='(F7.2)')+' A.'
          
          clr = replicate('green',nwave)
          clr[i] = 'red'

          user_xrange = (arcwave[i]-lambda0)/dispersion+50.0*[-1.0,+1.0]
          user_xrange[0] = user_xrange[0]>(-10)
          user_xrange[1] = user_xrange[1]<(ncol-1+10)
          
          if keyword_set(nolines) then begin
             arm_zoomplot, colaxis, flux, xsty=3, ysty=3, /single, coordinates=coord, $
               oplotline_x=(arcwave[i]-lambda0)/dispersion, oplotline_color=clr[i], $
               oplotline_linestyle=0, oplotline_thick=2, psym=10, /silent, $
               xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, user_xrange=user_xrange, $
               xtitle='Column Number [Pixel]', ytitle='Normalized Intensity'
          endif else begin
             arm_zoomplot, colaxis, flux, xsty=3, ysty=3, /single, coordinates=coord, $
               oplotline_x=(arcwave-lambda0)/dispersion, oplotline_color=clr, $
               oplotline_linestyle=0, oplotline_thick=2, psym=10, /silent, $
               xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, user_xrange=user_xrange, $
               xtitle='Column Number [Pixel]', ytitle='Normalized Intensity'
          endelse
          
          if (coord[1] gt 0.0) then begin

             result1.linepix = coord[0]
             result1.intensity = coord[1]
             result1.linewave = arcwave[i]
             result1.element = element[i]
             result1.quality = quality[i]
             result1.lineindx = i

; if there are 3 or more measurements then update the starting
; wavelength and dispersion; compute the local dispersion and starting
; wavelength 

             if (cntr ge 3L) then begin

;               coeff = robust_linefit(result[(cntr-3L)>0L:(cntr-1L)].linepix,result[(cntr-3L)>0L:(cntr-1L)].linewave)
                coeff = robust_linefit(result.linepix,result.linewave)

                lambda0 = coeff[0]
                dispersion = coeff[1]

                splog, 'New starting wavelength: '+string(lambda0,format='(G0.0)')
                splog, 'New dispersion         : '+string(dispersion,format='(G0.0)')

             endif else begin

                lambda0 = arcwave[i]-result1.linepix*dispersion

             endelse
             
             if (cntr eq 0L) then result = result1 else $
               result = struct_append(result,result1)
             
             cntr = cntr + 1L
             
          endif else begin

             splog, 'Skipping '+element[i]+' at '+string(arcwave[i],format='(F7.2)')+' A.'

          endelse
          
       endif 

    endfor

; open output file

    if file_test(datapath+output_lampfile,/regular) then begin
       splog, 'Overwrite existing file '+datapath+output_lampfile+' [Y/N]? ', format='(A,$ )'
       cc = get_kbrd(1)
       if strupcase(cc) ne 'Y' then return
    endif

    openw, lun, datapath+output_lampfile, /get_lun
    arm_filestamp, lun, notes='Wavelength, Intensity, GOOD/BLEND/BAD, Element, Pixel'
    struct_print, struct_trimtags(result,select=['LINEWAVE','INTENSITY','QUALITY','ELEMENT','LINEPIX'],$
      format=['F11.4','I6','A6','A6','F11.2']), lun=lun, /no_head
    free_lun, lun
    
return
end

