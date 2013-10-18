;+
; NAME:
;       ISKYSUBTRACT2D
;
; PURPOSE:
;       Sky subtract a two-dimensional sky spectrum.
;
; INPUTS:
;       skyapfile - sky aperture file as written by ISKYSELECT
;
; OPTIONAL INPUTS:
;       datapath    - I/O path
;	skymethod   - a scalar integer indicating the type of sky
;                     subtraction desired:
;		        0 - no sky subtraction
;		        1 - mean sky value (with rejection)
;		        2 - median sky value (with rejection)
;		        3 - robust NORDER_SKY polynomial fit to each
;                           column independently
;                       4 - simultaneous two-dimensional b-spline fit  
;       sigclip_sky - iteratively reject SIGCLIP_SKY outliers 
;       skyfirst    - index of first object to sky subtract 
;       skylast     - index of last object to sky subtract
;       norder_sky  - order of the sky fit in the spatial dimension
;                     (only used if SKYMETHOD=3 or 4); overwrites the
;                     value of NORDER_SKY in SKYAPFILE
;       wmapname    - wavelength map name outputed by IARCFIT  
;       tracename   - spatial distortion map outputed by
;                     IFITDISTORTION  
;       extra       - keywords for BSPLINE_ITERFIT()
;
; KEYWORD PARAMETERS:
;       fullsky     - use all the spatial rows to model the 2D sky
;                     spectrum (supercedes SKYAPFILE, SKYAPERTURE,
;                     SKYLOWER, SKYUPPER & SKYPROMPT); appropriate for
;                     spectra of point sources such as standard stars
;       checkskysub - if SKYMETHOD=3 then plot the fit to each column
;                     and wait for a keystroke
;       gzip        - compress the output FITS file and quality
;                     assurance plot 
;       wfits       - write the sky-subtracted image to DATAPATH 
;
; OUTPUTS:
;       If WFITS=1 then each file in SKYAPFILE is written out with an
;       's' prepended and, optionally, compressed.
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SPLOG, IFORAGE(), IRESTORE_WMAP(), ICHOOSE_WMAP(), CMRESTORE,
;       READCOL, RD2DSPEC(), BSPLINE_ITERFIT(), SXADDPAR, SXADDHIST,
;       WRT2DSPEC, ICLEANUP, DJS_MEDIAN, DJS_OPLOT, ARM_PLOTCONFIG,
;       BSPLINE_VALU(), DJS_ITERSTAT, SXPAR(), ROBUST_POLY_FIT(),
;       IM_WINDOW 
;
; COMMENTS: 
;       Note:  The following comments pertain to backwards
;       compatibility with an old sky subtraction code and should not
;       impact any typical user.
;
;       There are two aperture schemes employed by ISKYSUBTRACT2D.  1)
;       2 apertures defined by SKYLOWER and SKYUPPER referenced to row
;       zero and row NROWS-1 respectively; and 2) N apertures which
;       are all referenced to row zero Scheme '1' is used if
;       'REFZERO=1' in SKYAPFILE is not written in the header.  In
;       this case the first aperture listed is treated as SKYLOWER,
;       the last aperture is treated as SKYUPPER and any other
;       apertures are ignored.  Scheme '2' is used 'REFZERO=1' is
;       written in the header of SKYAPFILE, as done automatically with
;       ISKYSELECT.  All ISKYSUBTRACT2D output (e.g., image headers)
;       adopt scheme '2'.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 December 8, U of A, with significant
;          contributions and testing by C. Tremonti
;       05-Feb-2004  Generalized for N sky apertures - A.R.Marble
;       17-Feb-2004  Interactive sky window selection - A.R.Marble
;       jm04jun23uofa - interactive selection of sky apertures was
;                       written out to ISKYSELECT; this routine was
;                       streamlined and re-documented
;       jm04sep15uofa - some major changes; added SKYMETHOD and
;                       CHECKSKYSUB inputs; the default now is to do
;                       polynomial fitting to each column
;                       independently; the b-spline fitting is now
;                       accessible by SKYMETHOD=4, and assumes no
;                       error in the fitting
;       jm05jan14uofa - NORDER_SKY optional input added
;       jm05jun21uofa - use IM_WINDOW to spawn monitor resolution 
;                       independent windows; write the wavelength map
;                       trace set as a fourth FITS extension; bug fix:
;                       the inverse variance input to
;                       BSPLINE_ITERFIT() was commented out, which
;                       prevents iterative rejection (thanks to
;                       Christy Tremonti) 
;       jm05jun27uofa - bug fix; floor the sky aperture to be at least
;                       one pixel (row) wide;
;       jm05jun29uofa - bug fix in the way the error in the sky
;                       subtraction was being computed; bug fix:
;                       FORAGE was not being properly indexed within
;                       the FOR loop on SKYLIST
;       jm06jan24uofa - added support: the input image(s) can have
;                       been wavelength-calibrated previously; added
;                       SPLINE_BKSPACE and SPLINE_NORD optional inputs 
;
; Copyright (C) 2003-2005, John Moustakas
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

pro iskysubtract2d, skyapfile, datapath=datapath, skymethod=skymethod, $
  sigclip_sky=sigclip_sky, skyfirst=skyfirst, skylast=skylast, $
  norder_sky=norder_sky1, spline_nord=spline_nord, spline_bkspace=spline_bkspace, $
  wmapname=wmapname, tracename=tracename, _extra=extra, fullsky=fullsky, $
  checkskysub=checkskysub, gzip=gzip, wfits=wfits

    nskyapfile = n_elements(skyapfile)

    if (nskyapfile eq 0L) then begin
       doc_library, 'iskysubtract2d'
       return
    endif 

    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(sigclip_sky) eq 0L then sigclip_sky = 3.0
    if n_elements(skymethod) eq 0L then skymethod = 4L ; b-spline!

    if (skymethod lt 0L) or (skymethod gt 4L) then skymethod = 3L
    case skymethod of
       0L: splog, 'SKYMETHOD = '+string(skymethod,format='(I0)')+' - No sky subtraction.'
       1L: splog, 'SKYMETHOD = '+string(skymethod,format='(I0)')+' - Mean sky subtraction.'
       2L: splog, 'SKYMETHOD = '+string(skymethod,format='(I0)')+' - Median sky subtraction.'
       3L: splog, 'SKYMETHOD = '+string(skymethod,format='(I0)')+' - Polynomial sky subtraction.'
       4L: splog, 'SKYMETHOD = '+string(skymethod,format='(I0)')+' - B-spline sky subtraction.'
    endcase

; SKYMETHOD = 4 preliminaries
    
    if (skymethod eq 4L) then begin
       
       tcalib = n_elements(tracename)
       wcalib = n_elements(wmapname)

       if (tcalib[0] eq 1L) then if (tracename[0] eq '') then tcalib = 0L
       if (wcalib[0] eq 1L) then if (wmapname[0] eq '') then wcalib = 0L

       if (tcalib eq 0L) and (wcalib eq 0L) then begin
          splog, 'Please specify WMAPNAME and (optionally) TRACENAME.'
          return
       endif

; restore the trace/distortion structure

       if (tcalib ne 0L) then begin

          if file_test(datapath+tracename,/regular) then begin

             splog, 'Restoring '+tracename+'.'
             cmrestore, datapath+tracename, traceinfo, /quiet

             npad = traceinfo.npad
             dmap = traceinfo.dmap ; distortion map padded by 2*NPAD pixels spatially 
             
          endif else begin

             splog, 'Unable to find trace structure '+tracename+'.'
             return

          endelse

       endif
       
; restore the wavelength map(s)

       awavemap = irestore_wmap(wmapname,datapath=datapath,$
         dwave=adwave,wheader=awheader,wset=awset)
       if (awavemap[0] eq -1L) then return

    endif
       
; read in SKYAPFILE

    if (file_test(datapath+skyapfile,/regular) eq 0L) then begin
       
       splog, 'Sky aperture file '+datapath+skyapfile+' not found.'
       return

    endif

    splog, 'Reading sky aperture file '+datapath+skyapfile+'.'

    nlines = numlines(datapath+skyapfile)

    skylist = strarr(nlines)
    skyfileobject = strarr(nlines)
    skypositions = strarr(nlines)
    skyapertures = strarr(nlines)
    skyorder = lonarr(nlines)
    lastposition = fltarr(nlines) - 1
    lastaperture = fltarr(nlines) - 1

    openr, lun, datapath+skyapfile, /get_lun
    line = ''

    refzeroeq1 = 0L
    for i = 0L, nlines-1L do begin

       readf, lun, line
       
       if strmid(line,0,1) eq '#' then begin

          if (strmatch(line,'*REFZERO=1*') eq 1L) then refzeroeq1 = 1L

       endif else begin

          words = strsplit(strcompress(line),' ',/extract)
          nwords = n_elements(words)
          skylist[i] = words[0]
          skyfileobject[i] = words[1]
          skyorder[i] = long(words[nwords-1])

          if not refzeroeq1 then begin

             skypositions[i] = strjoin(words[2], ' ')
             lastposition[i] = float(words[(nwords-1)/2])
             skyapertures[i] = strjoin(words[(nwords-1)/2+1], ' ')
             lastaperture[i] = float(words[nwords-2])

          endif else begin

             skypositions[i] = strjoin(words[2:(nwords-1)/2], ' ')
             skyapertures[i] = strjoin(words[(nwords-1)/2+1:nwords-2], ' ')

          endelse

       endelse

    endfor

    free_lun, lun

    good = where(skylist ne '',nimage)
    skylist = skylist[good]
    skyfileobject = skyfileobject[good]
    skypositions = skypositions[good]
    skyapertures = skyapertures[good]
    skyorder = skyorder[good]
    lastposition = lastposition[good]
    lastaperture = lastaperture[good]

; optionally crop the list of objects    

    if (n_elements(skyfirst) eq 0L) then skyfirst = 0L else skyfirst = skyfirst > 0L
    if (n_elements(skylast) eq 0L) then skylast = n_elements(skylist)-1L else $
      skylast = skylast < (n_elements(skylist)-1L)

    if (skyfirst gt n_elements(skylist)-1L) then begin
       splog, 'SKYFIRST exceeds the number of objects in SKYLIST!'
       return
    endif
    
; forage header information

    forage = iforage(datapath+skylist)
    if (size(forage,/type) ne 8L) then begin
       splog, 'Error in IFORAGE().'
       return
    endif
    
; reference last position (if needed) to row zero

    if not refzeroeq1 then begin
       
       wh = where(lastposition ne -1, count)
       if count gt 0L then begin
          skypositions[wh] = skypositions[wh]+' '+ $
            strcompress((forage.naxis2*forage.cd2_2- $
            lastposition-lastaperture)[wh],/remove_all)
          skyapertures[wh] = skyapertures[wh]+' '+strcompress(lastaperture,/remove_all)
          
       endif             

    endif 

    outname = 's'+repstr(skylist,'.gz','') ; root output file names

; loop on each image    

    for i = skyfirst, skylast do begin

       splog, 'Sky-subtracting '+skylist[i]+'.'

; the WWSET output from RD2DSPEC() is needed to check if the spectrum
; has been wavelength calibrated already (see, e.g., SINGS2D_STITCH) 
          
       cube = rd2dspec(skylist[i],datapath=datapath,wset=wwset,/silent)

       image = cube.image
       mask = cube.mask
       sigmamap = cube.sigmamap
       invvar = (1.0/sigmamap^2.0) ; * (mask eq 0)
       header = *cube.header

       ncols = cube.naxis1
       nrows = cube.naxis2
       rowaxis = findgen(nrows)
       colaxis = findgen(ncols)

       if (n_elements(norder_sky1) eq 0L) then $
         norder_sky = skyorder[i] else norder_sky = norder_sky1
       
; if FULLSKY=1 then use all the spatial rows, irrespective of the
; apertures that have been set
       
       if keyword_set(fullsky) then begin

          splog, 'Setting FULLSKY.'
          skyrows = lindgen(nrows) 
          skyposition = 0
          skyaperture = forage[i].naxis2

       endif else begin
          
          skyposition = strsplit(skypositions[i], ' ', /extract)
          skyaperture = strsplit(skyapertures[i], ' ', /extract)
          
; round to one decimal point    
          
          skyposition = round(10D0*skyposition)/10.0
          skyaperture = round(10D0*skyaperture)/10.0
          
; convert to pixels    
          
          skyposition = long(skyposition/forage[i].cd2_2)        ; [pixel]
          skyaperture = (long(skyaperture/forage[i].cd2_2)) > 1L ; jm05jun27uofa 

          skyrows = lindgen(skyaperture[0])+skyposition[0]
          if n_elements(skyposition) gt 1L then $
            for j=1,n_elements(skyposition)-1 do $
            skyrows = [skyrows, lindgen(skyaperture[j])+skyposition[j]]
          
       endelse 

       skyrows = skyrows[uniq(skyrows,sort(skyrows))] ; don't count rows twice
       nskyrows = n_elements(skyrows)

       skyimage = image[*,skyrows]
       skyinvvar = invvar[*,skyrows]

       skyfit = image*0.0
       skyfit_sigma = image*0.0

; --------------------------------------------------       
; MEAN or MEDIAN sky subtraction
; --------------------------------------------------       

       if (skymethod eq 1L) or (skymethod eq 2L) then begin

          wavemap = cube.wavemap
          wset = wwset

          for icol = 0L, ncols-1L do begin

; compute the sky value             
             
             skyvalues = reform(skyimage[icol,*])
             skysigmas = reform(1.0/sqrt(skyinvvar[icol,*]))
             
             good = where(finite(skyvalues),ngood)
             if (ngood eq 0L) then message, 'Problem computing the sky value at column '+$
               string(icol,format='(I0)')+'.'

             djs_iterstat, skyvalues, sigrej=sigclip_sky, mean=skymean, $
               median=skymedian, sigma=skysig, mask=skymask

             goodsky = where(skymask,ngoodsky)
             if (ngoodsky eq 0L) then message, 'Problem computing the sky value at column '+$
               string(icol,format='(I0)')+'.'

             if (skymethod eq 1L) then skyfit[icol,*] = skymean   ; mean
             if (skymethod eq 2L) then skyfit[icol,*] = skymedian ; median

; compute the error             
             
             skyresid = skyimage[icol,*]-skyfit[icol,skyrows]
             djs_iterstat, skyresid, sigrej=sigclip_sky, $
               mean=skymean, sigma=skysig, mask=skymask

             goodsky = where(skymask,ngoodsky)
             if (ngoodsky eq 0L) then message, 'Problem computing the '+$
               'sky error at column '+string(icol,format='(I0)')+'.'
             
             skyfit_sigma[icol,*] = skysig/sqrt(ngoodsky) ; error in the mean

             if keyword_set(checkskysub) then begin

                if !d.window ne 3L then im_window, 3, xratio=0.5, /square
                ploterror, skyrows[good], skyvalues[good], skysigmas[good], $
                  xsty=3, ysty=3, ps=3, charsize=2.0, charthick=2.0, xtitle='Row Number', $
                  ytitle='Counts', /nohat
                djs_oplot, rowaxis, skyfit[icol,*], color='green', thick=3.0, line=0
                djs_oplot, rowaxis, (skyfit[icol,*]+skyfit_sigma[icol,*]), color='green', thick=3.0, line=2
                djs_oplot, rowaxis, (skyfit[icol,*]-skyfit_sigma[icol,*]), color='green', thick=3.0, line=2
                legend, 'Column '+strn(icol,format='(I0)')+'/'+strn(ncols,format='(I0)'), $
                  /left, /top, box=0, charsize=2.0, charthick=2.0
                cc = get_kbrd(1)

             endif

          endfor 
          
       endif

; --------------------------------------------------       
; POLYNOMIAL sky subtraction
; --------------------------------------------------       

       if (skymethod eq 3L) then begin

          wavemap = cube.wavemap
          wset = wwset

          for icol = 0L, ncols-1L do begin

             skyvalues = reform(skyimage[icol,*])
             skysigmas = reform(1.0/sqrt(skyinvvar[icol,*]))
             
             good = where(finite(skyvalues),ngood)
             if (ngood eq 0L) then message, 'Problem computing the sky value at column '+$
               string(icol,format='(I0)')+'.'

; fit the sky - unweighted; set the error in the sky fit equal to the
; rms scatter of the sky
             
             skycoeff = poly_fit(skyrows[good],skyvalues[good],norder_sky,$
               yfit=skyfitcol,sigma=skycoeff_err,covar=covar,yband=skyband)
;              measure_errors=skysigmas[good])

             skyfit[icol,*] = poly(rowaxis,skycoeff)

; compute the error             
             
             skyresid = skyimage[icol,*]-skyfit[icol,skyrows]
             djs_iterstat, skyresid, sigrej=sigclip_sky, $
               mean=skymean, sigma=skysig, mask=skymask
             
             goodsky = where(skymask,ngoodsky)
             if (ngoodsky eq 0L) then message, 'Problem computing the '+$
               'sky error at column '+string(icol,format='(I0)')+'.'

             skyfit_sigma[icol,*] = skysig/sqrt(ngoodsky) ; error in the mean

             if keyword_set(checkskysub) then begin
                
                if !d.window ne 3L then im_window, 3, xratio=0.5, /square
                ploterror, skyrows[good], skyvalues[good], skysigmas[good], $
                  xsty=3, ysty=3, ps=3, charsize=2.0, charthick=2.0, xtitle='Row Number', $
                  ytitle='Counts', /nohat
                djs_oplot, rowaxis, skyfit[icol,*], color='green', thick=3.0, line=0
                djs_oplot, rowaxis, (skyfit[icol,*]+skyfit_sigma[icol,*]), color='green', thick=3.0, line=2
                djs_oplot, rowaxis, (skyfit[icol,*]-skyfit_sigma[icol,*]), color='green', thick=3.0, line=2
                legend, 'Column '+strn(icol,format='(I0)')+'/'+strn(ncols,format='(I0)'), $
                  /left, /top, box=0, charsize=2.0, charthick=2.0
                cc = get_kbrd(1)

             endif

          endfor
             
       endif

; --------------------------------------------------       
; 2D B-SPLINE sky subtraction
; --------------------------------------------------       

       if (skymethod eq 4L) then begin
       
; unpad the distortion structure, if it exists
          
          ymap_nodistort = (fltarr(ncols)+1.0) # findgen(nrows)
          if (tcalib eq 0L) then $
            ymap = (fltarr(ncols)+1.0) # findgen(nrows) else $
            ymap = dmap[*,npad:npad+nrows-1L]

; minimize the difference in RA, DEC to select the appropriate
; wavelength map, or use WWSET

          if (n_elements(wwset) eq 0L) then begin
             indx = ichoose_wmap(header,awheader,wmapname=wmapname)
             wmapname_out = wmapname[indx]
             wavemap = reform(awavemap[*,*,indx])
             dwave = (adwave[indx])[0]
             wset = (awset[indx])[0]
          endif else begin
             wavemap = cube.wavemap
             wset = wwset
             dwave = fix(100.0*(max(wavemap[*,nrows/2L])-min(wavemap[*,nrows/2L]))/(ncols-1))/100.0 ; [Angstrom/pixel]
          endelse

          if (n_elements(spline_nord) eq 0L) then spline_nord = 4L
          if (n_elements(spline_bkspace) eq 0L) then spline_bkspace = dwave

          file = strmid(skylist[i],0,strpos(skylist[i],'.fits'))
          
          skywavemap = wavemap[*,skyrows]
          skyymap = ymap[*,skyrows]
          skyymap_nodistort = ymap_nodistort[*,skyrows]

; use SKYYMAP *not* SKYYMAP_NODISTORT to model the distortions in the
; spatial dimension
          
          sset = bspline_iterfit(skywavemap,skyimage,nord=spline_nord,x2=skyymap,$
            npoly=norder_sky,bkspace=spline_bkspace,yfit=skyfit,upper=sigclip_sky,$
            lower=sigclip_sky,invvar=skyinvvar,_extra=extra,/silent)

;         bkpt = fltarr(ncols)
;         for icol = 0L, ncols-1L do bkpt[icol] = mean(skywavemap[icol,*])

;         bkpt = fltarr(3*ncols)
;         for icol = 0L, ncols-1L do bkpt[icol*3:(icol+1)*3-1L] = $
;           [min(skywavemap[icol,*]),mean(skywavemap[icol,*]),max(skywavemap[icol,*])]

;         sset = bspline_iterfit(skywavemap,skyimage,nord=spline_nord,x2=skyymap_nodistort,$
;           npoly=norder_sky,bkpt=bkpt,yfit=skyfit,upper=sigclip_sky,$
;           lower=sigclip_sky,invvar=skyinvvar,_extra=extra,/silent)

          skyfit = bspline_valu(wavemap,sset,x2=ymap)

;         djs_plot, skywavemap[*,10], skyimage[*,10], ps=4, xsty=3, ysty=3, xr=[5572,5578]; xr=[5560,5590]
;         djs_oplot, skywavemap[*,120], skyimage[*,120], ps=5, color=djs_icolor('green')
;         djs_oplot, skywavemap[*,10], skyfit[*,10], color='red'
;         djs_oplot, skywavemap[*,120], skyfit[*,120], color='cyan'
;         for ibk = 0L, n_elements(sset.fullbkpt)-1L do plots, sset.fullbkpt[ibk], !y.crange[0], ps=7

;         s = sort(skywavemap)
;         djs_plot, skywavemap[s], skyimage[s], ps=3, xsty=3, ysty=3, xr=[5400,6000] ; xr=[5572,5582];xr=[5560,5590]
;         djs_oplot, skywavemap[s], skyfit[s], color='red', ps=10
;         for ibk = 0L, n_elements(sset.fullbkpt)-1L do plots, sset.fullbkpt[ibk], !y.crange[0], ps=7

;         if keyword_set(checkskysub) then begin
;
;         endif
          
; compute the "simple" uncertainty in the sky subtraction, at each
; column 

;         skyfit_sigma = skyfit*0.0
;         skyfit_sigma = sqrt( (skyfit/gain)>0.0 + (rdnoise/gain)^2.0 )

          for icol = 0L, ncols-1L do begin

             skyresid = skyimage[icol,*]-skyfit[icol,skyrows]
             djs_iterstat, skyresid, sigrej=sigclip_sky, $
               mean=skymean, sigma=skysig, mask=skymask

             goodsky = where(skymask,ngoodsky)
             if (ngoodsky eq 0L) then message, 'Problem computing the '+$
               'sky error at column '+string(icol,format='(I0)')+'.'

             skyfit_sigma[icol,*] = skysig/sqrt(ngoodsky) ; error in the mean
;            skyfit_sigma[icol,*] = sqrt(total((skyimage[icol,goodsky]-skymean)^2.0)/(ngoodsky*(ngoodsky-1.0)))

          endfor

          if (n_elements(wmapname_out) ne 0L) then $
            sxaddpar, header, 'WMAPNAME', wmapname_out, ' wavelength map name', before='HISTORY'

       endif

       imnosky = image - skyfit
       signosky = sqrt(sigmamap^2.0 + skyfit_sigma^2.0)

; ---------------------------------------------------------------------------
; generate quality assurance plot
; ---------------------------------------------------------------------------
;
;      oldpmulti = !p.multi
;      arm_plotconfig, psfile=datapath+'qaplot_sky_s'+file+'.ps', $
;        nx=2, ny=8, coords=coords, xmargin=[1.25, 0.5], ymargin=[1.5,0.75], $
;        xspace=.75, yspace=[0.0,0.75,0.0,0.75,0.0,0.75,0.0], $
;        height=[1.0,0.625,1.0,0.625,1.0,0.625,1.0,0.625]
;      
;      xyouts, 0.5, 0.96, align=.5, /normal, charsize=1.5, charthick=2, $
;        '2D SKY SUBTRACTION QUALITY ASSURANCE PLOT'
;      xyouts, 0.5, 0.94, align=.5, charsize=1.0, charthick=2, /normal, $
;        's'+file+'.fits (norder_sky='+strn(norder_sky,f='(i)')+', shaded sky windows)'
;      xyouts, 0.5, 0.0, align=.5, /normal, charsize=1.0, charthick=2, $
;        'postscript file generated: '+systime()
;      xyouts, (coords[0,0]+coords[2,0])/2, 0.9, align=0.5, charsize=1.25, $
;        charthick=2, 'before', /normal
;      xyouts, (coords[0,1]+coords[2,1])/2, 0.9, align=0.5, charsize=1.25, $
;        charthick=2, 'after', /normal
;      
;      for j=0L, 3L do begin
;         
;         imgvals = total(image[j/4.*ncols:(j+1)/4.*ncols-1,*], 1) / (ncols/4.)
;         fitvals = total(skyfit[j/4.*ncols:(j+1)/4.*ncols-1,*], 1) / (ncols/4.)
;         subvals = total(imnosky[j/4.*ncols:(j+1)/4.*ncols-1,*], 1) / (ncols/4.)
;
;         djs_iterstat, imgvals[skyrows], sigrej=4, mask=msk
;         ok = where(msk)
;         
;         plot, rowaxis, /nodata, xstyle=7, yrange=[0,1], ystyle=5, $
;           position=[coords[0:1, 4*j+2], coords[2:3, 4*j]], /normal, /noerase
;         for k=0,n_elements(skyposition)-1 do begin
;            x = [skyposition[k], skyposition[k] + skyaperture[k]]
;            polyfill, [x[0], x[1], x[1], x[0], x[0]], [0, 0, 1, 1, 0], color=djs_icolor('grey')
;         endfor
;         
;         yrange1 = [min(imgvals)-0.1*(max(imgvals)-min(imgvals)), $
;                    max(imgvals)+0.1*(max(imgvals)-min(imgvals))]
;         yrange2 = [min(imgvals)-0.1*(max(imgvals[skyrows[ok]])-min(imgvals[skyrows[ok]])), $
;                    max(imgvals[skyrows[ok]])+0.1*(max(imgvals[skyrows[ok]])-min(imgvals[skyrows[ok]]))]
;
;         plot, rowaxis, imgvals, psym=10, thick=2, xthick=2, ythick=2, xstyle=3, $
;           yrange=yrange1, ystyle=1, $
;           position=coords[*,4*j], /normal, /noerase, charsize=1.5, charthick=2, $
;           xtickname=replicate(' ', 30), $
;           title='mean of columns: '+strn(j/4.*ncols, $
;                                          f='(i)')+' - '+strn((j+1)/4.*ncols, f='(i)')
;         djs_oplot, fitvals, thick=4, color='white'
;         oplot, fitvals, thick=2
;
;         plot, rowaxis, imgvals, psym=10, thick=2, xthick=2, ythick=2, xstyle=3, $
;           yrange=yrange2, ystyle=1, $
;           position=coords[*,4*j+2], /normal, /noerase, charsize=1.5, charthick=2, $
;           xtitle='row (pixels)'
;         djs_oplot, fitvals, thick=4, color='white'
;         oplot, fitvals, thick=2
;         xyouts, 0.075, (coords[3,4*j]+coords[1,4*j+2])/2., orientation=90, align=0.5, $
;           charthick=2, charsize=1., 'counts', /normal
;
;         plot, rowaxis, /nodata, xstyle=7, yrange=[0,1], ystyle=5, $
;           position=[coords[0:1,4*j+3], coords[2:3,4*j+1]], /normal, /noerase
;         for k=0,n_elements(skyposition)-1 do begin
;            x = [skyposition[k], skyposition[k] + skyaperture[k]]
;            polyfill, [x[0], x[1], x[1], x[0], x[0]], [0, 0, 1, 1, 0], color=djs_icolor('grey')
;         endfor
;         
;         yrange1 = [min(subvals)-0.1*(max(subvals)-min(subvals)), $
;                    max(subvals)+0.1*(max(subvals)-min(subvals))]
;         yrange2 = [min(subvals)-0.1*(max(subvals[skyrows[ok]])-min(subvals[skyrows[ok]])), $
;                    max(subvals[skyrows[ok]])+0.1*(max(subvals[skyrows[ok]])-min(subvals[skyrows[ok]]))]
;
;         plot, rowaxis, subvals, psym=10, thick=2, xthick=2, ythick=2, xstyle=3, $
;           yrange=yrange1, ystyle=1, $
;           position=coords[*,4*j+1], /normal, /noerase, charsize=1.5, charthick=2, $
;           xtickname=replicate(' ', 30), $
;           title='mean of columns: '+strn(j/4.*ncols, $
;                                          f='(i)')+' - '+strn((j+1)/4.*ncols, f='(i)')
;         djs_oplot, rowaxis, rowaxis*0, thick=4, color='white'
;         oplot, rowaxis, rowaxis*0, thick=2
;
;         plot, rowaxis, subvals, psym=10, thick=2, xthick=2, ythick=2, xstyle=3, $
;           yrange=yrange2, ystyle=1, $
;           position=coords[*,4*j+3], /normal, /noerase, charsize=1.5, charthick=2, $
;           xtitle='row (pixels)'
;         djs_oplot, rowaxis, rowaxis*0, thick=4, color='white'
;         oplot, rowaxis, rowaxis*0, thick=2
;         
;      endfor       
;      
;      device, /close
;      set_plot, 'x'
;      !p.multi = oldpmulti
;      if keyword_set(gzip) then spawn, ['gzip -f '+datapath+'qaplot_sky_s'+file+'.ps'], /sh
; ---------------------------------------------------------------------------

; update the header       

       sxaddpar, header, 'ISKYSB', skylist[i], ' input to ISKYSUBTRACT2D', $
         before='HISTORY'

       sxaddpar, header, 'SKYMETHD', skymethod, $
         ' sky method', before='HISTORY'
       sxaddpar, header, 'SKYNPOLY', fix(norder_sky), $
         ' sky spatial polynomial order', before='HISTORY'
       if (skymethod eq 4L) then $
         sxaddpar, header, 'SKYNORD', fix(spline_nord), $
         ' sky wavelength b-spline order', before='HISTORY'

       for j = 0L, n_elements(skyposition)-1L do begin
          
          sxaddpar, header, 'SKYROW'+string(j+1L,format='(I2.2)'), skyposition[j], $
            ' sky window starting row [pixel]', before='HISTORY'
          sxaddpar, header, 'SKYAP'+string(j+1L,format='(I2.2)'), skyaperture[j], $
            ' sky window aperture [pixel]', before='HISTORY'

       endfor

       sxaddhist, "'Sky spectrum subtracted "+im_today()+"'", header

; write out

       if keyword_set(wfits) then begin
          splog, 'Writing '+datapath+outname[i]+'.'
          wrt2dspec, outname[i], float(imnosky), float(signosky), cube.mask, $
            header, skyimage=float(skyfit), wset=wset, datapath=datapath, gzip=gzip
       endif

       icleanup, cube
       
    endfor

    free_lun, lun

return
end
