;+
; NAME:
;       IMONTECARLO_SKYSUBTRACT
;
; PURPOSE:
;       Investigate the error in the b-spline sky subtraction method
;       using a Monte Carlo technique.
;
; CALLING SEQUENCE:
;       imontecarlo_skysubtract, skyapfile, wmapfilelist, result, datapath=, $
;          sigclip_sky=, skyfirst=, skylast=, nmonte=, /write, _extra=extra
;
; INPUTS:
;       skyapfile     - sky aperture file as written by ISKYSELECT
;       wmapfilelist  - file providing the list of wavelength maps
;                       outputed by IARCFIT   
;
; OPTIONAL INPUTS:
;       datapath    - I/O path
;       sigclip_sky - iteratively reject SIGCLIP_SKY outliers 
;       skyfirst    - index of first object to sky subtract 
;       skylast     - index of last object to sky subtract
;       nmonte      - number of Monte Carlo iterations
;       extra       - keywords for BSPLINE_ITERFIT()
;
; KEYWORD PARAMETERS:
;       write - write RESULT to DATAPATH
;
; OUTPUTS:
;       result - data structure containing the results of the Monte
;                Carlo fitting
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       READCOL, SPLOG, IRESTORE_WMAP(), NUMLINES(), IFORAGE(),
;       RD2DSPEC(), BSPLINE_ITERFIT(), DJS_ITERSTAT, DJS_MEAN(),
;       MWRFITS 
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jun 29, U of A - written
;
; Copyright (C) 2005, John Moustakas
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

pro imontecarlo_skysubtract, skyapfile, wmapfilelist, result, datapath=datapath, $
  sigclip_sky=sigclip_sky, skyfirst=skyfirst, skylast=skylast, nmonte=nmonte, $
  write=write, _extra=extra
    
    if (n_elements(skyapfile) eq 0L) or (n_elements(wmapfilelist) eq 0L) then begin
       print, 'Syntax - imontecarlo_skysubtract, skyapfile, wmapfilelist, datapath=, $'
       print, '   sigclip_sky=, skyfirst=, skylast=, nmonte=, /write, _extra=extra'
       return
    endif
    
    if (n_elements(datapath) eq 0L) then datapath = cwd()

    if (file_test(datapath+skyapfile,/regular) eq 0B) then begin
       splog, 'Sky aperture file '+datapath+skyapfile+' not found.'
       return
    endif

    if file_test(datapath+wmapfilelist,/regular) then begin
       readcol, datapath+wmapfilelist, wmapname, format='A', /silent, comment='#'
    endif else begin
       splog, 'Wavelength map list '+datapath+wmapfilelist+' not found.'
       return
    endelse
    
    if n_elements(sigclip_sky) eq 0L then sigclip_sky = 3.0
    if (n_elements(nmonte) eq 0L) then nmonte = 500L

       montefile = 'monte_carlo_skysub.fits' ; output file

; restore the wavelength map(s)

    awavemap = irestore_wmap(wmapname,datapath=datapath,$
      dwave=adwave,wheader=awheader,wset=awset)
    if (awavemap[0] eq -1L) then return

; read in SKYAPFILE - this algorithm was taken from ISKYSUBTRACT2D
    
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
       if (count gt 0L) then begin
          skypositions[wh] = skypositions[wh]+' '+ $
            strcompress((forage.naxis2*forage.cd2_2- $
            lastposition-lastaperture)[wh],/remove_all)
          skyapertures[wh] = skyapertures[wh]+' '+strcompress(lastaperture,/remove_all)          
       endif             

    endif 

; loop on each object

    for i = skyfirst, skylast do begin

       splog, 'Monte Carlo sky-subtracting '+skylist[i]+'.'

       cube = rd2dspec(skylist[i],datapath=datapath,/silent)
       image = cube.image
       mask = cube.mask
       sigmamap = cube.sigmamap
       invvar = (1.0/sigmamap^2.0) ; * (mask eq 0)
       header = *cube.header
       icleanup, cube

       ncols = cube.naxis1
       nrows = cube.naxis2
       rowaxis = findgen(nrows)
       colaxis = findgen(ncols)

; define the sky apertures       
       
       skyposition = strsplit(skypositions[i], ' ', /extract)
       skyaperture = strsplit(skyapertures[i], ' ', /extract)
       
; round to one decimal point    
       
       skyposition = round(10D0*skyposition)/10.0
       skyaperture = round(10D0*skyaperture)/10.0
       
; convert to pixels    
       
       skyposition = long(skyposition/forage[i].cd2_2) ; [pixel]
       skyaperture = (long(skyaperture/forage[i].cd2_2)) > 1L ; jm05jun27uofa 

       skyrows = lindgen(skyaperture[0])+skyposition[0]
       if n_elements(skyposition) gt 1L then $
         for j=1,n_elements(skyposition)-1 do $
         skyrows = [skyrows, lindgen(skyaperture[j])+skyposition[j]]
       
       skyrows = skyrows[uniq(skyrows,sort(skyrows))] ; don't count rows twice
       nskyrows = n_elements(skyrows)

; minimize the difference in RA, DEC to select the appropriate
; wavelength map

       indx = ichoose_wmap(header,awheader,wmapname=wmapname)
       wmapname_out = wmapname[indx]
       wavemap = reform(awavemap[*,*,indx])
       dwave = (adwave[indx])[0]
       wset = (awset[indx])[0]

       file = strmid(skylist[i],0,strpos(skylist[i],'.fits'))

       skyimage = image[*,skyrows]
       skysigmap = sigmamap[*,skyrows]
       skyinvvar = invvar[*,skyrows]
       skywavemap = wavemap[*,skyrows]

       ymap = (fltarr(ncols)+1.0) # findgen(nrows)
       skyymap = ymap[*,skyrows]

; initialize the output data structure; do it within the loop because
; we need NCOLS; assume that the images in SKYLIST all have the same
; number of columns

       result1 = {$
         file:               skylist[i],    $
         object:             forage[i].object, $
         ra:                 forage[i].ra,     $
         dec:                forage[i].dec,    $
         skymean:            fltarr(ncols), $
         simple_sigma:       fltarr(ncols), $
         monte_bspline_mean: fltarr(ncols,nmonte), $
         monte_sigma:        fltarr(ncols)}
       
; Monte Carlo this bad boy

       t1 = systime(1)
       for imonte = 0L, nmonte-1L do begin
    
          if (((imonte+1L) mod 5L) eq 0L) then $
            print, format='("Monte Carlo iteration = ",I4,"/",I4,".",A4,$)', $
            imonte+1L, nmonte, string(13b)
          
          if (imonte eq 0L) then skyimage_monte = skyimage else $
            skyimage_monte = skyimage + randomn(seed,ncols,nskyrows)*skysigmap
          
          sset = bspline_iterfit(skywavemap,skyimage_monte,nord=4,x2=skyymap,$
            npoly=norder_sky,bkspace=dwave,yfit=skyfit,upper=sigclip_sky,$
            lower=sigclip_sky,invvar=skyinvvar,_extra=extra,/silent)
          skyfit = bspline_valu(wavemap,sset,x2=ymap)

          for icol = 0L, ncols-1L do begin
          
; compute the "simple" uncertainty in the sky subtraction, just once          

             if (imonte eq 0L) then begin

                djs_iterstat, skyimage[icol,*], sigrej=sigclip_sky, mean=skymean, $
                  sigma=skysig, mask=skymask

                goodsky = where(skymask,ngoodsky)
                if (ngoodsky eq 0L) then message, 'Problem computing the sky error at column '+$
                  string(icol,format='(I0)')+'.'

                skysig = sqrt(total((skyimage[icol,goodsky]-skymean)^2.0)/(ngoodsky*(ngoodsky-1.0)))

                result1.skymean[icol] = skymean
                result1.simple_sigma[icol] = skysig

             endif

; store the mean sky level according to the b-spline fit

             result1.monte_bspline_mean[icol,imonte] = djs_mean(skyfit[icol,*])
             
          endfor 
          
       endfor
       splog, 'Total time for Monte Carlo = '+strtrim(string((systime(1)-t1)/60.0,$
         format='(F12.1)'),2)+' minutes.'

; compute the error in the b-spline sky subtraction by collapsing the
; Monte Carlo dimension

       for icol = 0L, ncols-1L do result1.monte_sigma[icol] = $
         stddev(result1.monte_bspline_mean[icol,*])
;      plot, result1.skymean, result1.simple_sigma/result1.monte_sigma, ps=4, /xlog

; write out after every object in case this routine crashes

       if (n_elements(result) eq 0L) then result = result1 else result = [ [result], [result1] ]

; write out the results    

       if keyword_set(write) then begin

          splog, 'Writing '+datapath+montefile+'.'
          mwrfits, reform(result), datapath+montefile, /create
          spawn, ['gzip -f '+datapath+montefile], /sh

       endif

    endfor
    result = reform(result)

stop    
    
return
end
