;+
; NAME: 
;       IARCFIT_MAKEFILE
;       
; PURPOSE:
;       Generate a file containing line positions for use in obtaining
;       wavelength solutions for arc calibration lamp images.
;
; CALLING SEQUENCE: 
;       iarcfit_makefile, filelist, outfile, linelist=, lambda0=, $
;          dispersion=, [/nolines]
;
; INPUTS:
;       filelist   - array of 1 or more filenames, each a text file
;                    listing arc calibration lamp filenames
;       outfile    - name of file to be generated
;       linelist   - text file containing two columns: line
;                    wavelength, line identification (eg, 'HeI')
;       lambda0    - initial guess at starting wavelength of arc
;                    calibration lamp files
;       dispersion - dispersion of arc calibration lamp files
;       
; OPTIONAL INPUTS:
;
; KEYWORDS:
;       nolines    - suppress plotting of lines not currently being
;                    identified
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;       ARM_ZOOMPLOT, RD2DSPEC(), ARM_FILESTAMP, READCOL
;
; COMMENTS: 
;       The output file is generated in the following way.  Lines are
;       plotted on top of the first lamp file spectrum at the
;       wavelengths specified in LINELIST.  The use can zoom in on the
;       plot as desired by selecting opposing corners of the desired
;       region with the left mouse button.  The right mouse button
;       pans back out to show the original plot.  The user indicates
;       the center/peak of the emission feature corresponding to the
;       red line (how well they line up depends on the accuracy of
;       LAMBDA0 and DISPERSION) by clicking on it with the middle
;       mouse button.  This continues until the position of each line
;       has been determined.  Then the relative offset for each
;       additional lamp file is similarly determined.  The strongest
;       line is indicated by a red line.  The user uses the middle
;       mouse button to indicate the center of the line.  The
;       difference between that position and that of the red line is
;       the offset for that file.
; 
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., 2005 Feb 25
;
; Copyright (C) 2005, Andy R. Marble
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

pro iarcfit_makefile, filelist, outfile, linelist=linelist, $
  lambda0=lambda0, dispersion=dispersion, nolines=nolines

; error checking

    if n_elements(outfile) eq 0L then message,'OUTFILE not specified.'
    if n_elements(filelist) eq 0L then message, 'FILELIST not specified.'
    if n_elements(lambda0) eq 0L then message,'LAMBDA0 not defined.'
    if n_elements(dispersion) eq 0L then message,'DISPERSION not defined.'
    if n_elements(linelist) eq 0L then message,'LINELIST not specified.'
       
; read lamp line information
    
    readcol, linelist, f='A,A,A', arcwavetxt, arcclass, arctype, comment='#'
    arcwave = double(arcwavetxt)
    nwave = n_elements(arcwave)
    
; open output file
    
    if file_search(outfile) ne '' then $
      spawn, 'cp '+outfile+' '+outfile+'.old'
    openw, lun, outfile, /get_lun
    arm_filestamp, lun, notes='Wavelength, Intensity, GOOD/BLEND/BAD, Pixel'
    
; read in first lamp image
    
    readcol, filelist[0], f='A', comment='#', list
    firstlamp = (rd2dspec(list[0])).image
    ncol = n_elements(firstlamp[*,0])
    spec = total(firstlamp,2)
    spec = spec/max(spec)*1000.
    
    position = fltarr(2,nwave)
    
    cntr = 0L
;;    prevwave = lambda0
    ;;   prevpxl = 0
    
; loop over all supplied lines
    
    disp0=dispersion
    
    for i=0,nwave-1 do begin
       
;;       if i gt 0 then begin
;;          prevwave = arcwave[i-1]
;;          prevpxl = position[0,i-1]
;;       endif
       
       if arcwave[i] ge lambda0 and arcwave[i] le lambda0+dispersion*ncol then begin
          
          cntr = cntr + 1L
          
          clr = replicate('green',nwave)
          clr[i] = 'red'
          
          if keyword_set(nolines) then begin
             arm_zoomplot, spec, xst=3, yst=3, /single, coordinates=coord, $
               oplotline_x=(arcwave[i]-lambda0)/dispersion, oplotline_color=clr[i], $
               oplotline_linestyle=0,oplotline_thick=2
          endif else begin
             arm_zoomplot, spec, xst=3, yst=3, /single, coordinates=coord, $
               oplotline_x=(arcwave-lambda0)/dispersion, oplotline_color=clr, $
               oplotline_linestyle=0,oplotline_thick=2
          endelse
          position[*,i] = coord
          lambda0 = arcwave[i]-coord[0]*dispersion
          
;;          if cntr ge 2 then begin
;;             dispersion = (arcwave[i]-arcwave[i-1])/(position[0,i]-position[0,i-1])
          if cntr eq 1L then wave0 = arcwave[i]-coord[0]*disp0
          ;;         endif
          printf, lun, f='(A,2x,i4,2x,A,2x,A,1x,f9.3)', $
            arcwavetxt[i],round(coord[1]),arcclass[i], arctype[i],coord[0]
          
       endif
       
    endfor
    
    printf, lun
    printf, lun, 'The intial wavlength and dispersion values are: '+ $
      string(wave0,f='(f9.3)')+' '+string(disp0,f='(f9.3)')
    
    
    printf, lun
    printf, lun, '# Comparison Lamp Frame, Measured Lamp Frame, Relative Pixel Offset'
    printf, lun
    
    wh = (where(position[1,*] eq max(position[1,*])))[0]
    clr = replicate('green',nwave)
    clr[wh] = 'red'
    
    for i=0,n_elements(filelist)-1 do begin
       
       readcol, filelist[i], f='A', comment='#', list
       
       for j=0,n_elements(list)-1 do begin
          
; read in first lamp image
          
          lamp = (rd2dspec(list[j])).image
          spec = total(lamp,2)
          spec = spec/max(spec)*1000.
          
          arm_zoomplot, spec, xst=3, yst=3, /single, coordinates=coord, $
            oplotline_x=(arcwave-arcwave[wh])/dispersion+position[0,wh], oplotline_color=clr, $
            oplotline_linestyle=0,oplotline_thick=2
          if i+j eq 0 then begin
             first = coord[0]
             firstfile = list[0]             
          endif
          printf, lun, f='(A,2x,A,1x,f9.3)', list[j], firstfile, first - coord[0]
          
       endfor
       
    endfor
    
    close, lun
    free_lun, lun
    
 end
