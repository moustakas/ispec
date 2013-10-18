;+
; NAME:
;	IARCFIT
;
; PURPOSE:
;	Generate a two dimensional map between pixel position and
;	wavelength from an arc calibration lamp.
;
; CALLING SEQUENCE:
;       iarcfit, arcfile, wset, xcen, ycen, datapath=, lampname=, $
;          wmapname=, row=, nmed_arc=, xcoeff=, dxcoeff=, func_arc=, $
;          norder_arc=, ncoeff_arc=, mintol=, arc_pixel=, arclambda=, $
;          /write, /debug, /doplot, _extra=extra
;
; INPUTS:
;	arcfile  - arc lamp FITS file name
;
; OPTIONAL INPUTS:
;	datapath   - I/O path name
;	lampname   - comparison lamp name (e.g., 'HeAr') or a
;                    user-supplied lamp list (e.g., 'lamplist.dat')
;                    which is distinguished by the ".dat" extension 
;       wmapname   - output wavelength map name
;       row        - row to use for initial wavelength solution
;                    (default to the middle row) 
;	nmed_arc   - number of rows to median filter around ROW for
;                    initial wavelength solution (default 5)
;	xcoeff     - two-element array giving the initial
;                    guesses [starting wavelength (Angstrom),
;                    dispersion (Angstrom/pixel)] (default
;                    [3600.0,2.75])
;	dxcoeff    - deviations from XCOEFF (default [50.0,0.1])
;	func_arc   - name of fitting function (default 'legendre')
;	norder_arc - order of the fit of column number versus
;                    wavelength (default 4)
;       ncoeff_arc - number of coefficients for tracing lamp lines
;                    (default 7) 
;	mintol     - minimum tolerance for acceptable arc lines
;                    (default 0.5 Angstrom)
;       arc_pixel  - pixel positions of all lamp lines (useful for
;                    forcing iarcfit to use weak lines in
;                    user-supplied line list)
;       extra      - extra keywords
;
; KEYWORD PARAMETERS:
;	write    - write the wavelength map and the QA plots to disk
;       debug    - plot the functional fits to the arc line traces 
;       doplot   - generate QA plots on the screen (will be generated
;                  anyway if WRITE=1)
;
; OUTPUTS:
;	wset    - 2D wavelength map (pixel -> lambda)
;	xcen    - column pixel positions of arc lines [nrows,nlines] 
;	ycen    - row positions of arc lines [nrows,nlines]
;
; OPTIONAL OUTPUTS:
;	arclambda - wavelengths of good lamp lines (Angstrom)
;
; COMMENTS:
;       An NCOEFF_ARC Legendre polynomial is fitted to the row-dependent
;       arc line position, while an NORDER_ARC Legendre polynomial is
;       fitted in the wavelength dimension.
;
;	Also see FITARCIMAGE in Dave Schlegel's IDLSPEC2D package.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD2DSPEC(), READCOL, IARCFIT_GUESS(), DJS_MEDIAN(),
;	TRACE_CRUDE(), FIND_NPEAKS(),  TRACESET2XY(), TRACESET2PIX(),
;	XY2TRACESET, TRACE_FWEIGHT(), TRACE_GWEIGHT(), ICLEANUP,
;	DFPSPLOT, DFPSCLOSE, CWD(), GET_ELEMENT, CLEANPLOT,
;	LINEID_PLOT, DJS_ITERSTAT, ERRPLOT, LEGEND, CMSAVE,
;	IREAD_LAMPLIST(), IM_WINDOW
;
; DATA FILES:
;	${ISPEC_DIR}/etc/CRC04_*.dat
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August, U of A - written
;       jm02nov15uofa - added a quality assurance plot
;       jm03jan12uofa - added support for a HeNeAr lamp
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm03dec07uofa - added DEBUG keyword
;       jm04jun21uofa - call IREAD_LAMPLIST(); changed LAMP keyword to
;                       LAMPNAME; relegated XDIFF and GAUSS keywords;
;                       added NCOEFF_ARC, ARC_PIXEL, and DEBUG
;                       keywords; implemented the changes by A. Marble
;                       from 2004-Apr-8
;       jm05jun17uofa - make sure NMED_ARC does not exceed the size of
;                       the image (thanks to A. Marble)
;       jm05jun21uofa - use IM_WINDOW to spawn monitor resolution 
;                       independent windows
;
; Copyright (C) 2001-2005, John Moustakas
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

pro iarcfit, arcfile, wset, xcen, ycen, datapath=datapath, lampname=lampname, $
  wmapname=wmapname, row=row, nmed_arc=nmed_arc, xcoeff=xcoeff, dxcoeff=dxcoeff, $
  func_arc=func_arc, norder_arc=norder_arc, ncoeff_arc=ncoeff_arc, mintol=mintol, $
  arc_pixel=arc_pixel, arclambda=arclambda, write=write, debug=debug, doplot=doplot, $
  _extra=extra

    narc = n_elements(arcfile)
    if narc eq 0L then begin
       print, 'Syntax - iarcfit, arcfile, wset, xcen, ycen, datapath=, lampname=, $'
       print, '   wmapname=, row=, nmed_arc=, xcoeff=, dxcoeff=, func_arc=, $'
       print, '   norder_arc=, ncoeff_arc=, mintol=, arc_pixel=, arclambda=, $'
       print, '   /write, /debug, /doplot, _extra=extra'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(lampname) eq 0L then lampname = 'HeAr'
    if n_elements(mintol) eq 0L then mintol = 0.5
    if n_elements(xcoeff) eq 0L then xcoeff = [3600.0,2.75]
    if n_elements(dxcoeff) eq 0L then dxcoeff = [50.0,0.1]
    if n_elements(nmed_arc) eq 0L then nmed_arc = 5L
    if n_elements(norder_arc) eq 0L then norder_arc = 4L
    if n_elements(ncoeff_arc) eq 0L then ncoeff_arc = 7L
    if n_elements(func_arc) eq 0L then func_arc = 'legendre' else func_arc = strlowcase(func_arc)

    if (narc gt 1L) then begin

       if n_elements(wmapname) eq 0L then wmapname = replicate('',narc) else $
         if n_elements(wmapname) ne narc then $
         message, 'Incompatible dimensions in ARCFILE and WMAPNAME.'
       
       for k = 0L, narc-1L do begin
          
          iarcfit, arcfile[k], wset, xcen, ycen, datapath=datapath, lampname=lampname, $
            wmapname=wmapname[k], row=row, nmed_arc=nmed_arc, xcoeff=xcoeff, dxcoeff=dxcoeff, $
            func_arc=func_arc, norder_arc=norder_arc, ncoeff_arc=ncoeff_arc, mintol=mintol, $
            arc_pixel=arc_pixel, arclambda=arclambda, write=write, debug=debug, doplot=doplot, $
            _extra=extra

       endfor 
       return

    endif

    cube = rd2dspec(arcfile,datapath=datapath,_extra=extra)

    arc = cube.image
    arcinvvar = 1.0/cube.sigmamap^2.0D ; inverse variance
    archead = *cube.header
    ncols = cube.naxis1
    nrows = cube.naxis2

    if (n_elements(row) eq 0L) then row = nrows/2L else row = long(row)

; read in the lamplist for wavelength calibration into a structure

    lamp = iread_lamplist(lampname,_extra=extra)
    if (size(lamp,/type) ne 8L) then return
    
; median filter NMED_ARC rows centered on ROW; this median spectrum
; determines the initial wavelength solution
    
    spec = djs_median(arc[*,(row-nmed_arc/2L)>0L:(row+nmed_arc/2L)<(nrows-1L)],2)
    spec = 1E4*spec/max(spec)

; initialize the input wavelength solution based on the dispersion and
; starting wavelength.  the linear transformations below are being
; computed to be compatible with the way the tracesets are stored

    xrange = ncols-1.0
    xmid = 0.5*xrange

    acoeff = fltarr(4)
    acoeff[0] = xcoeff[0] + xmid * xcoeff[1]
    acoeff[1] = xcoeff[1] * xrange / 2.0
    acoeff[2:3] = [5.0,0.0] ; [5.0,-0.8]

    dcoeff = fltarr(4)
    dcoeff[0] = dxcoeff[0] + xmid * dxcoeff[1]
    dcoeff[1] = dxcoeff[1] * xrange / 2.0
    dcoeff[2:3] = [10.0,2.0]
    
    splog, 'Searching for the initial wavelength solution.'
;   nsteps = [1, 80, 40, 1]
    aset = iarcfit_guess(spec,lamp.lambda,100*lamp.intensity,acoeff=acoeff,$
      dcoeff=dcoeff,nsteps=nsteps,bestcor=bestcor,model=model)
    initcoeff = [aset.coeff[0]-xmid*xcoeff[1],aset.coeff[1]*2.0/xrange,aset.coeff[2:3]]
    splog, 'Initial wavelength fit = ', initcoeff, format='(a,99f12.5)'

;   traceset2xy, aset, xdum, rowlambda
;   splog, 'Wavelength range in ROW ', minmax(rowlambda)
;   splog, 'Dispersion in ROW ', (max(rowlambda)-min(rowlambda))/ncols

; trim the lamplist to only those wavelengths within the wavelength
; range.  predict the pixel positions of the lines in our lamp list
; using the initial wavelength solution and attempt to trace those
; lines on the arc lamp image

    if n_elements(arc_pixel) eq 0L then $
      xstart = traceset2pix(aset,lamp.lambda,/silent) else xstart = arc_pixel

    qtrim = (xstart gt 1L) and (xstart lt ncols-2L) and lamp.good
    itrim = where(qtrim,nlamp)
    if nlamp eq 0L then message, 'No arc lines in wavelength range.'
    xstart = xstart[itrim]
    lamp = lamp[itrim]
    
; trace arc lines as a function of row number on the 2D arc image.
; allow for a shift of up to 2 pixels in the initial centers, but only
; 0.5 pixels while tracing, starting at the middle row

    splog, 'Tracing ', nlamp, ' arc lines.'
    xcen = trace_crude(arc,arcinvvar,xstart=xstart,ystart=row,radius=radius,$
      nave=1.0,nmed=1.0,maxerr=maxerr,maxshift0=1.0D,maxshifte=0.2D,$
      yset=ycen)
;   print, total(xcen,1)/nrows

; trace from the bottom of the CCD using the previous trace as the
; initial line centers

    xcen = trace_crude(arc,arcinvvar,xstart=xcen[0,*],ystart=0L,radius=radius,$
      nave=1.0,nmed=1.0,maxerr=maxerr,maxshifte=0.2D,maxshift0=1.0D,xerr=xerr,$
      yset=ycen)
;   print, total(xcen,1)/nrows

    if keyword_set(debug) then begin
    
       djs_plot, findgen(ncols), arc[*,row], xsty=3, ysty=3, ps=10, $
         charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, ytitle='Counts', $
         xtitle='Column Number';, xrange=[100,200]
       djs_oplot, findgen(ncols), arc[*,20], color='yellow', ps=10
       djs_oplot, findgen(ncols), arc[*,nrows-20], color='red', ps=10
       for i = 0L, nlamp-1L do djs_oplot, [xcen[row,i],xcen[row,i]], $
         !y.crange, line=1, thick=0.5, color='cyan'
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
       
    endif

; fix bad traces (optional)

;   xcen = trace_fix(xcen,minsep=minsep,ngrow=ngrow,ycen=ycen,xerr=xerr)

; for each line, fit a legendre function to the trace to make the
; center a smoothly varying function and to interpolate over bad rows

    xy2traceset, ycen, xcen, crudeset, ncoeff=ncoeff_arc, yfit=xcrudefit, /silent

    if keyword_set(debug) then begin

       print
       print, 'left  mouse button deselects nearest value and refits.'
       print, 'right mouse button advances.' 
       print
    
       xranges = fltarr(nlamp)
       for i = 0L, nlamp-1L do xranges[i] = max(xcen[*,i])-min(xcen[*,i])

       for i = 0L, nlamp-1L do begin

          use = lonarr(nrows) + 1L
          inmask = xcen * 0L + 1L
          prev = nrows + 1L

          while total(use) lt prev do begin

             djs_plot, findgen(nrows), xcen[*,i], ps=4, xsty=3, ysty=3, $
               charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
               yr=total(minmax(xcen[*,i]))/2.+[-1,1]*max(xranges)/2., $
               ytitle='Column Position (pixel)', xtitle='Row Number'
             if total(use) lt nrows then djs_oplot, [where(use eq 0L),where(use eq 0L)], $
               [xcen[where(use eq 0L),i], xcen[where(use eq 0L),i]], ps=7,th=2
             djs_oplot, findgen(nrows), xcrudefit[*,i], line=0, thick=2.0, color='red'
             legend, 'Line '+string(i+1,format='(I3.3)')+'/'+string(nlamp,format='(I3)'), $
               /left, /top, box=0, charsize=1.5, charthick=2.0

             prev = total(use)
             cursor, xmouse, ymouse, 3, /device
             if !mouse.button eq 1L then begin

                coords = convert_coord(findgen(nrows), xcen[*,i], /data, /to_device)
                dist = sqrt((coords[0,*]-xmouse)^2 + (coords[1,*]-ymouse)^2)

                wh = where(dist eq min(dist[where(use)]))
                use[wh] = 0L
                
                inmask[*,i] = use
                
                xy2traceset, ycen, xcen, crudeset, ncoeff=ncoeff_arc, $
                  yfit=xcrudefit, /silent, inmask = inmask

             endif
             
          endwhile
          
       endfor

    endif
    
; iterate the flux-weighted centers.  in the last iteration, use the
; formal errors in the arc image

    xnew = trace_fweight(arc,xcrudefit,ycen)
    xnew = trace_fweight(arc,xnew,ycen)

    if (keyword_set(gauss)) then $
      xcen = trace_gweight(arc,xnew,ycen,sigma=1.0,invvar=arcinvvar,xerr=xerr) else $
      xcen = trace_fweight(arc,xnew,ycen,radius=2.0,invvar=arcinvvar,xerr=xerr)

; to obtain the final two-dimensional wavelength map fit column
; position (pixels) as a function of line list wavelength (Angstrom).
; use the formal errors in the centroids and do an NORDER_ARC order
; fit to the function between wavelength and column number at each row 

    if (nlamp eq 1L) then $
      ytmp = transpose(lamp.lambda*(dblarr(nrows)+1)) else $
      ytmp = lamp.lambda # (dblarr(nrows)+1)

    xy2traceset, transpose(double(xcen)), ytmp, wset, func=func_arc, $
      ncoeff=(norder_arc+1), maxiter=nlamp, maxrej=1.0, /sticky, xmin=0, $
      xmax=ncols-1, yfit=yfit, /silent

;   plot, xcen[63,*], ytmp[*,63]-yfit[*,63], ps=4

; check for bad or mismatched lines and reject them

    resid2d = ytmp - yfit
    resid1d = total(resid2d,2)/nrows ; mean residuals

    repeat begin
      
       dev = max(abs(resid1d),indx)
;      print, 'Bad line - off by '+strn(dev,length=5)+' Angstrom.'

       lmask = bytarr(nlamp)+1B
       lmask[indx] = 0B
       w = where(lmask,nlamp)
       lamp = lamp[w]
       xcen = xcen[*,w]
       
       if (nlamp eq 1L) then $
         ytmp = transpose(lamp.lambda*(dblarr(nrows)+1)) else $
         ytmp = lamp.lambda # (dblarr(nrows)+1)

       xy2traceset, transpose(double(xcen)), ytmp, wset, func=func_arc, ncoeff=(norder_arc+1), $
         maxiter=nlamp, maxrej=1.0, /sticky, xmin=0, xmax=ncols-1, yfit=yfit, /silent
       
       resid2d = ytmp - yfit
       resid1d = total(resid2d,2)/nrows ; mean residuals
       
    endrep until total(abs(resid1d) gt mintol) eq 0

    splog, 'Final arcfit complete based on '+string(nlamp,format='(I0)')+' lines.'

;   niceprint, lamp.lambda, ytmp[*,nrows/2L], yfit[*,nrows/2L], ytmp[*,nrows/2L]-yfit[*,nrows/2L]
    
    arclambda = lamp.lambda                                ; wavelengths of good arc lines
    xarc = transpose(traceset2pix(wset,arclambda,/silent)) ; arc pixel positions [nrows,nlamp]
    xdiff = xcen-xarc ; residuals between the traced and the predicted positions

    traceset2xy, wset, xmap, ymap
    xwave = reform(ymap[*,row])
    xpix = reform(xarc[row,*])
    yflux = reform(arc[*,row] / max(arc[*,row]))

; measure the amplitude of every line in the final list

    dwave = mean(wset.coeff[1,*]*2/xrange)

    narclambda = n_elements(arclambda)
    amplitude = fltarr(narclambda)
    labellamp = strarr(narclambda)
    for ii = 0L, narclambda-1L do begin

       get_element, xwave, arclambda[ii]+3*dwave*[-1.0,1.0], xx
       get_element, xwave, arclambda[ii]+15*dwave*[-1.0,1.0], xx2
       amplitude[ii] = max(yflux[xx[0]:xx[1]])
;      plot, xwave[xx2[0]:xx2[1]], yflux[xx2[0]:xx2[1]], xsty=3, ysty=3, ps=10
;      djs_oplot, arclambda[ii]*[1,1], !y.crange, line=2, color='red', thick=2
;      cc = get_kbrd(1)
       
    endfor

    srt = reverse(sort(amplitude))
    labellamp = lamp[srt[0L:24L<(nlamp-1L)]]
;   labellamp = lamp
    
; OLD CODE - jm06jun28uofa
;   amplitude = interpol(yflux,xwave,arclambda)
;   srt = reverse(sort(amplitude))
;   labellamp = lamp[srt[0L:19L<(nlamp-1L)]]
    
; quality assurance plots

    ifitname = 'arc_'+strmid(arcfile,0,strpos(arcfile,'.fits'))
    
; plot 1
    
;   bigpeaks = trace_crude(yflux # (fltarr(10)+1),thresh=0.05,yset=yset)
;   peaks = reform(bigpeaks[0,*])
;
;   match, round(xpix), round(peaks), indx, suba
;
;   labellamp = lamp[indx]
    
    if keyword_set(doplot) or keyword_set(write) then begin

       if keyword_set(doplot) then iter = 1L else iter = 0L

       for i = 0L, iter do begin
          
          if (i eq 0L) then begin
             dfpsplot, datapath+'qaplot_'+ifitname+'_lineid.ps', /square, /color
             postthick = 4.0
          endif else begin
             im_window, 2, xratio=0.49, /square
             postthick = 2.0
          endelse

          srt = sort(xwave)
          lineid_plot, xwave[srt], yflux[srt], labellamp.lambda, labellamp.element, thick=postthick-2, $
            string(labellamp.lambda,format='(F6.1)'), xthick=postthick, ythick=postthick, $
            charsize=1.3, charthick=postthick, xsty=3, ysty=3, ps=10, yrange=[-0.05,1.0], /extend, $
            xtitle=textoidl('Wavelength [\AA]'), ytitle='Normalized Flux', $
            lcharsize=1.3, lcharthick=postthick-2, position=[0.15,0.12,0.95,0.75], $
            subtitle='['+ifitname+']'

          if (i eq 0L) then dfpsclose
          
       endfor
       cleanplot, /silent

    endif

; plot 2
    
; compute offset and stddev for each line center in pixels and in Angstroms

    meanpix = fltarr(nlamp)
    sigmapix = fltarr(nlamp)

    meanlam = fltarr(nlamp)
    sigmalam = fltarr(nlamp)
    
    for k = 0L, nlamp-1L do begin

; pixels
       
       djs_iterstat, xdiff[*,k], mean=mn, sigma=sg
       meanpix[k] = mn
       sigmapix[k] = sg

; wavelength

       djs_iterstat, yfit[k,*]-ytmp[k,*], mean=mn, sigma=sg
       meanlam[k] = mn
       sigmalam[k] = sg

    endfor

; determine the starting and ending wavelengths.  set the plot range

    traceset2xy, wset, transpose(replicate(wset.xmin,nrows)), minwave
    traceset2xy, wset, transpose(replicate(wset.xmax,nrows)), maxwave
   
    waverange = [min([minwave,maxwave]),max([minwave,maxwave])]
    splog, 'Wavelength range (Angstrom): ', strn(waverange[0]), ', ', strn(waverange[1])
    splog, 'Mean disperion (Angstrom/pixel): ', strn((waverange[1]-waverange[0])/ncols)

; plot the residuals in pixel position as a function of wavelength

    if keyword_set(doplot) or keyword_set(write) then begin

       if keyword_set(doplot) then iter = 1L else iter = 0L

       for i = 0L, iter do begin
      
          if (i eq 0L) then begin
             dfpsplot, datapath+'qaplot_'+ifitname+'.ps', /square, /color
             postthick = 4.0
          endif else begin
             im_window, 0, xratio=0.49, /square
             postthick = 2.0
          endelse

          yminmax = max(abs([max(meanlam+sigmalam),min(meanlam-sigmalam)]))
       
          ploterror, arclambda, meanlam, sigmalam, xrange=waverange, yrange=yminmax*[-1,1], ps=6, $
            xstyle=3, ystyle=3, xtitle=textoidl('Wavelength [\AA]'), ytitle=textoidl('Mean Residuals [\AA]'), $
            subtitle='['+ifitname+']', charsize=1.3, charthick=postthick, xthick=postthick, $
            ythick=postthick, thick=postthick, errthick=postthick, ymargin=[5,2]
          oplot, waverange, [0,0], thick=postthick
;         errplot, lamp.lambda, meanlam-sigmalam, meanlam+sigmalam, thick=postthick

          djs_iterstat, meanlam, sigrej=5.0, mean=rmean, sigma=rsig, median=rmedian

          stats = strtrim(string(rmean,format='(F12.3)'),2)+' \pm '+strtrim(string(rsig,format='(F12.3)'),2)+$
            ' ('+strtrim(string(rmedian,format='(F12.3)'),2)+')'
          legend, textoidl(stats), /right, /top, box=0, charsize=1.3, charthick=postthick
          
;         legend, [textoidl('\Delta')+' = '+string(mean(meanlam),format='(F8.5)'),$
;           textoidl('\sigma_{\Delta}')+' = '+string(stddev(meanlam),format='(F8.5)')], $
;           /right, /top, box=0, charsize=1.7, charthick=postthick
       
          if (i eq 0L) then dfpsclose
             
       endfor 

    endif

    if keyword_set(write) then begin

       if n_elements(wmapname) eq 0L then wmapname = 'wmap_'+ifitname+'.idlsave'
       if strcompress(wmapname,/remove) eq '' then wmapname = 'wmap_'+ifitname+'.idlsave'
         
       splog, 'Writing '+datapath+wmapname+'.'
       cmsave, filename=datapath+wmapname, wset, arclambda, xcen, archead, $
         names=['wset','arclambda','xarc','header']

       openw, lun, datapath+'qalog_'+ifitname+'.log', /get_lun
       printf, lun, '# Two dimensional wavelength map derived from ', arcfile
       printf, lun, '# '+strn(ncols), ' columns and ', strn(nrows), ' rows'
       printf, lun, '# Wavelength range (Angstrom): ', strn(waverange[0]), ', ', strn(waverange[1])
       printf, lun, '# Mean dispersion (Angstrom/pixel): ', strn(mean(wset.coeff[1,*]*2/xrange))
       printf, lun, '# '+strn(nlamp), ' arc lines used in the final solution'
       printf, lun, '# Minimum wavelength tolerance for a good arc line: ', strn(mintol), ' Angstrom'
       printf, lun, '# Mean residuals and RMS of the fit (Angstrom): ', strn(mean(meanlam)), ', ', strn(stddev(meanlam))
       printf, lun, '# Fitting function and order: ', strn(func_arc), ', ', strn(norder_arc)
       printf, lun, '#  '
       printf, lun, '# Nominal wavelength, mean measured wavelength, pixel position:'
       printf, lun, '#  '

       for i = 0L, nlamp-1L do printf, lun, arclambda[i], mean(yfit[i,*]), mean(xarc[*,i])
       
       free_lun, lun
       
    endif
    
    icleanup, cube
    
return
end
