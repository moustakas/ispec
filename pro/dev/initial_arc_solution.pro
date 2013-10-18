;+
; NAME:
;       INITIAL_ARC_SOLUTION
;
; PURPOSE:
;       Derive the first-order "guess" at the wavelength solution by
;       marking and labeling comparison lamp lines.  
;
; CALLING SEQUENCE:
;       initial_arc_solution, arcfile, datapath=, row=, nmed_arc=, $
;          outfile=
;
; INPUTS:
;       arcfile - two-dimensional arc lamp (post-ICCDPROC)
;
; OPTIONAL INPUTS:
;       datapath - I/O path
;       row      - median average NMED_ARC rows centered on ROW
;       nmed_arc - see ROW
;       outfile  - output file name (default initial_arc_solution.dat) 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       A text file is written out with the following information in ?
;       columns: 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       CWD(), RD2DSPEC(), DJS_MEDIAN(), DJS_PLOT, SPLOG, MPFITPEAK(),
;       MPFITPEAK_GAUSS(), ROBUST_LINEFIT(), STRUCT_APPEND(), ICLEANUP 
;
; TODO:
;       [1] Handle a list of arc lamps.
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Jun 23, U of A - written
;       jm04jul07uofa - documented and generalized
;-

pro initial_arc_solution, arcfile, datapath=datapath, row=row, $
  nmed_arc=nmed_arc, outfile=outfile

    narcfile = n_elements(arcfile)
    if (narcfile eq 0L) then begin
       print, 'Syntax - initial_arc_solution, arcfile, datapath=, row=, $'
       print, '   nmed_arc=, outfile='
       return
    endif
    
    if (n_elements(outfile) eq 0L) then outfile = 'initial_arc_solution.dat'

; read ARCFILE    
    
    if (n_elements(datapath) eq 0L) then datapath = cwd()

    cube = rd2dspec(arcfile,datapath=datapath)

    arc = cube.image
    arcinvvar = 1.0/cube.sigmamap^2.0D ; inverse variance
    archead = *cube.header
    ncols = cube.naxis1
    nrows = cube.naxis2

    colaxis = findgen(ncols)
    
; median filter NMED_ARC rows centered on ROW
    
    if n_elements(nmed_arc) eq 0L then nmed_arc = 5L
    if (n_elements(row) eq 0L) then row = nrows/2L else row = long(row)

    flux = djs_median(arc[*,row-nmed_arc/2L:row+nmed_arc/2L],2)
    mask = byte(flux*0.0)

; initialize the data structure

    result1 = {$
      linewave:  0.0, $
      intensity: 0.0, $
      linepix:   0.0}
    
    device, get_screen_size=scrnsize
    window, 0, xsize=0.95*scrnsize[0], ysize=0.5*scrnsize[1]
    
    bigfit = flux*0.0
    
    done = 0L
    count = 0L
    
    while (done eq 0L) do begin

       if (count eq 0L) then begin
          splog, '   left button: select opposite corners of zoom region'
          splog, ' middle button: identify the line center and peak'
          splog, '  right button: restore original plotting range'
       endif

       if (count eq 0L) then begin
          arm_zoomplot, colaxis, flux, xsty=3, ysty=3, ps=10, charsize=1.5, $
            charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Column', ytitle='Flux', $
            /single, coordinates=coords, /silent, thick=2.0, oplot_psym=10
       endif else begin
          arm_zoomplot, colaxis, flux, xsty=3, ysty=3, ps=10, charsize=1.5, $
            charthick=2.0, xthick=2.0, ythick=2.0, xtitle='Column', ytitle='Flux', $
            /single, coordinates=coords, /silent, oplot_x=colaxis, oplot_y=bigfit, $
            oplot_color='orange', oplot_psym=10
       endelse       

; fit a Gaussian function       
       
       gflux = flux

       bad = where(mask eq 1B,nbad)
       if (nbad ne 0L) then gflux[bad] = 0.0

       x1 = fix(coords[0]-10)
       x2 = fix(coords[0]+10)       

       gfit = mpfitpeak(colaxis[x1:x2],gflux[x1:x2],a,nterms=5,/positive,/gauss)

       xrange = fix(a[1]+[-5.0,5.0]*a[2])
       
       subaxis = colaxis[xrange[0]:xrange[1]]
       subflux = flux[xrange[0]:xrange[1]]

; compute the flux-weighted centroid then re-evaluate the Gaussian       

       center1 = total(subflux*indgen(n_elements(subaxis))) / total(subflux) + xrange[0]
       result1.pixel = center1
       result1.intensity = a[0]

       a[1] = center1
       gfit = mpfitpeak_gauss(subaxis,a)
       bigfit[xrange[0]:xrange[1]] = gfit

;      xrange = fix(a[1]+[-5.0,5.0]*a[2]) ; re-compute the range
       mask[xrange[0]:xrange[1]] = (mask[xrange[0]:xrange[1]] + 1B) ge 1B

; plot the results       
       
       djs_oplot, subaxis, gfit, ps=10, color='orange'
       djs_oplot, center1*[1,1], !y.crange, color='blue', thick=2.0, line=2

       linewave1 = 0.0
       read, prompt='Line wavelength (negative to quit)? ', linewave1

       if (linewave1 gt 0.0) then begin

          result1.linewave = linewave1

          if (count eq 0L) then result = result1 else $
            result = struct_append(result,result1)

       endif else begin

          done = 1L

       endelse
       
; print the "robust" starting wavelength and dispersion
       
       if (count ge 3L) then begin

          coeff = robust_linefit(result.pixel,result.linewave)

          print, 'Starting wavelength: '+string(coeff[0],format='(F9.4)')
          print, 'Dispersion         : '+string(coeff[1],format='(F7.5)')
          
       endif
          
       count = count + 1L

    endwhile

    icleanup, cube

; initialize the output file and verify that it does not exist 
    
    if file_test(datapath+outfile,/regular) then begin
       splog, 'Overwrite existing file '+datapath+outfile+' [Y/N]? ', format='(A,$ )'
       cc = get_kbrd(1)
       if strupcase(cc) ne 'Y' then return
    endif

    openw, lun, outfile, /get_lun
    printf, lun, '# File generated by INITIAL_ARC_SOLUTION '+im_today()
    struct_print, struct_trimtags(result,select='*',format=['F7.2','F12.2','F7.2']), $
      lun=lun, /no_head
    free_lun, lun
    
stop

return
end
