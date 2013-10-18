;+
; NAME:
;	ICCDPROC
;
; PURPOSE:
;	All-purpose CCD spectroscopy processing routine.
;
; CALLING SEQUENCE:
;       iccdproc, proclist, datapath=, biasfile=, domefile=, $'
;          skyfile=, darkfile=, trim=, overscan=, norder_overscan=, $'
;          gain=, rdnoise=, saturation=, flatname=, flatlist=, $'
;          illumname=, badpixfile=, prefix=, /makeflat, /clean, $'
;          /wfits, /doplot, /gzip, objlist=, _extra=extra'
;
; INPUTS:
;	proclist   - file list (relative to DATAPATH) of images to
;                    process, including object frames and arc lamps
;
; OPTIONAL INPUTS:
;	datapath   - path to the data to be processed (default PWD) 
;       biasfile   - FITS file name of the bias frame
;       domefile   - FITS file name of the dome flat
;       skyfile    - FITS file name of the sky flat
;       darkfile   - FITS file name of the dark count image
;	trim       - [x1,x2,y1,y2] defining the good data region 
;                    (zero-indexed) 
;	overscan   - [x1,x2,y1,y2] defining the overscan region
;                    (zero-indexed) 
;	norder_overscan - order of the fit to the overscan section
;                         (default 2); if negative, then do not
;                         subtract the overscan
;       gain       - CCD gain [electron per ADU] (default 1.0)
;       rdnoise    - CCD read noise [electron] (default 0.0)
;       saturation - CCD saturation level [electron] (if not equal to
;                    zero then check for saturated pixels)
;	flatname   - file name of the master flat to read (if MAKEFLAT
;                    is non-zero) or to write out (if WRITE is
;                    non-zero); default 'masterflat.fits'
;       flatlist   - Boolean:  0B if DOMEFILE is a single FITS file or
;                    1B if DOMEFILE is a file name of multiple flat
;                    fields 
;       illumname  - file name of the illumination flat to read (if
;                    MAKEFLAT is non-zero) or to write out (if WRITE
;                    is non-zero); default 'illumflat.fits'
;	badpixfile - bad pixel file name to read (default badpix.dat;
;                    IRAF format)
;	prefix     - word or letter to prepend to the name of the
;                    output FITS file(s) (default 'r')
;	extra      - keywords for IOVERSCAN, IHEADER_UPDATE(),
;                    IMASTERFLAT(), and IILLUMFLAT()
;	
; KEYWORD PARAMETERS:
;	makeflat   - create a master flat field
;	clean      - linearly interpolate over bad pixels
;	wfits      - write out the processed data as multi-extension
;                    FITS files
;       doplot     - generate QA plots
;       gzip       - write GZIPPED FITS files
;
; OPTIONAL OUTPUTS:
;	objlist    - returns the list of object files processed (arc
;                    lamps in PROCLIST have been removed).
;
; COMMENTS:
;	There may be memory problems if a very large file list is
;	passed due to a very large array of structures that is
;	initialized at the beginning of the routine.  The solution is
;	to pass smaller file lists or get a bigger computer.
;
;       Pixels with negative variance after flat-fielding and bias-
;       and dark-subtraction are flagged as replaced with the read
;       noise.  This occurs because the noise model is not treated
;       totally correctly (I think).
;
;       Dark subtraction is not treated correctly yet.
;
;       This routine initializes the bad pixel masks for all the
;       processed images.  Dead pixels and rows that were indicated in
;       BADPIXFILE are assigned an integer bit mask according to
;       IMASK_BITS().  
;    
; EXAMPLE:
;
; INTERNAL ROUTINES:
;       HEAD_APPEND(), INITIALPROC(), CHOOSE_FLAT()
;
; PROCEDURES USED:
;	SXADDPAR, IOVERSCAN, VMAP_INIT(), IMASTERFLAT(),
;	DJS_MASKINTERP(), CWD(), SPLOG, READCOL, READFITS(), ITRIM,
;	WRT2DSPEC, ICLEANUP, IM_FITS_CUBE(), DJS_ITERSTAT,
;	IHEADER_UPDATE(), DJS_DIFF_ANGLE(), IM_HMS2DEC(), SXPAR(),
;	IILLUMFLAT(), IM_FITS_CUBE(), DFPSPLOT, DFPSCLOSE, PLOTHIST,
;	LEGEND, IMASK_BITS(), IM_WINDOW
;
; TODO:
;       [1] Check for saturated pixels.  Need to add a parameter in 
;       IBATCH, I think.
;
;       [2] Treat dark-subtraction properly. 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August, U of A
;	jm01oct15uofa - major overhaul
;       jm02novuofa   - restructuring and re-organization
;       jm03apr09uofa - checked out ISPEC v1.0.0
;       jm03oct15uofa, EXTEND must follow the last NAXISi header entry
;          (thanks to Kris Eriksen)
;       jm03dec07uofa - added GZIP support; additional error checking 
;       jm05jun17uofa - removed QFLAG input from IOVERSCAN(); use
;                       IM_WINDOW to spawn monitor resolution
;                       independent windows
;       arm05nov30uofa - PROCLIST files can have pathways relative to
;                        DATAPATH allowing raw files to be elsewhere
;
; Copyright (C) 2001-2003, 2005, John Moustakas
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

function head_append, cube, field, value, comment, _extra=extra
; jm03jan10uofa

    h = cube.header
    w = where(strcompress(h,/remove) ne '',nw)
    if nw ne 0L then h = h[w]

    sxaddpar, h, field, value, comment, _extra=extra
    nh = n_elements(h)
    
    cube.header[0:nh-1L] = h
    
return, cube
end

function initialproc, proclist, datapath=datapath, trim=trim, overscan=overscan, $
  norder_overscan=norder_overscan, xybad=xybad, noheaderupdate=noheaderupdate, _extra=extra

    fcount = n_elements(proclist)

    for i = 0L, fcount-1L do begin ; loop on each image

       if file_test(datapath+proclist[i]) eq 0L then begin
          splog, 'File '+datapath+proclist[i]+' not found.'
          return, -1L
       endif

       image = readfits(datapath+proclist[i],header,/silent)
       if (size(image,/n_dimension) ne 2L) then begin
          splog, 'Image is not two-dimensional!'
          return, -1L
       endif

       image = float(image)
       imsz = size(image,/dimension)
       
; initialize the bad pixel mask.  smart error checking if the bad
; pixel file is only one line.  assign a bit number of 2 to dead
; pixels & columns

       mask = fix(image*0.0) ; 0 is good

       if n_elements(xybad) ne 0L then begin ; input bad pixels

          bit = imask_bits(/badpix) ; bad pixel bit value
          
          ndim = size(xybad,/n_dimension)
          dims = size(xybad,/dimension)

          if ndim[0] eq 1L then jmax = 0L else jmax = dims[1]-1L
          
          for j = 0L, jmax do begin

             xy = xybad[*,j]
             mask[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)] = bit
;            mask[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)] = 2
;            mask[xy[0]>0L:xy[1]<(imsz[0]-1L),xy[2]>0L:xy[3]<(imsz[1]-1L)] = 0B

          endfor

       endif

; overscan-subtract, and trim

       itrim, mask, trim=trim

;      splog, 'Subtracting the overscan region.'
       ioverscan, image, header, overscan=overscan, $
         norder_overscan=norder_overscan, _extra=extra

       itrim, image, header=header, trim=trim

; update the header

       if not keyword_set(noheaderupdate) then header = iheader_update(header,_extra=extra)
       
; ----------------------------------------------------------------------
; test code:  fix the bad pixels (works well)

;      fixim = djs_maskinterp(image,mask,iaxis=0) ; linearly interpolate 
;      fixim = djs_maskinterp(image,(mask eq 0B),iaxis=0) ; linearly interpolate 
; ----------------------------------------------------------------------
       
; initialize the structure cube for the first image and fill it

       if i eq 0L then begin

          imsize = size(image,/dimension)
          ncols = imsize[0]
          nrows = imsize[1]
          
          cube = {$
            fname:  '',        $ ; file name
            object: '',        $ ; object type (galaxy,standard,etc)
            type:   '',        $ ; file type (zero,flat,comp,object)
            naxis1: 0L,        $ ; x size
            naxis2: 0L,        $ ; y size
            header: strarr(1000L), $
            image:  fltarr(ncols,nrows), $ ; image
            vmap:   fltarr(ncols,nrows), $ ; variance map
            mask:   intarr(ncols,nrows)} ; bad pixel mask
          cube = replicate(cube,fcount)
          
       endif

       cube[i].fname = proclist[i]
       cube[i].object = strn(sxpar(header,'object'))
       cube[i].type = strn(sxpar(header,'imagetyp'))
       cube[i].naxis1 = strn(sxpar(header,'naxis1'))
       cube[i].naxis2 = strn(sxpar(header,'naxis2'))
       cube[i].header = header[0L:n_elements(header)-1L]
       cube[i].image = image
       cube[i].mask = mask

       cube[i] = head_append(cube[i],'IRAW',proclist[i],' raw FITS file',before='HISTORY')
    
    endfor

return, cube
end

function choose_flat, header, flatheader
; if NDIM is 1 then return the zeroth index
    
    ndim = size(flatheader,/n_dimension)
    if (ndim eq 2L) then nflat = (size(flatheader,/dimension))[1] else return, 0L

    ra = 15.0*im_hms2dec(sxpar(header,'RA'))
    de = im_hms2dec(sxpar(header,'DEC'))

    diff = fltarr(nflat)
    
    for i = 0L, nflat-1L do diff[i] = $
      djs_diff_angle(ra,de,15.0*im_hms2dec(sxpar(flatheader[*,i],'RA')),$
      im_hms2dec(sxpar(flatheader[*,i],'DEC')))

    mindiff = min(diff,indx)

return, indx
end

pro iccdproc, proclist, datapath=datapath, biasfile=biasfile, domefile=domefile, $
  skyfile=skyfile, darkfile=darkfile, trim=trim, overscan=overscan, $
  norder_overscan=norder_overscan, gain=gain, rdnoise=rdnoise, saturation=saturation, $
  flatname=flatname, flatlist=flatlist, illumname=illumname, badpixfile=badpixfile, $
  prefix=prefix, makeflat=makeflat, clean=clean, wfits=wfits, doplot=doplot, gzip=gzip, $
  objlist=objlist, _extra=extra

    if not keyword_set(datapath) then datapath = cwd()

    if not keyword_set(flatname) then flatname = 'flat.fits'
    if not keyword_set(illumname) then illumname = 'illumflat.fits'
    if not keyword_set(prefix) then prefix = 'r'

    if n_elements(gain) eq 0L then gain = 1.0
    if n_elements(rdnoise) eq 0L then rdnoise = 0.0
    
    if (n_elements(proclist) eq 0L) and (not keyword_set(makeflat)) then begin
       print, 'Syntax - iccdproc, proclist, datapath=, biasfile=, domefile=, $'
       print, '   skyfile=, darkfile=, trim=, overscan=, norder_overscan=, $'
       print, '   gain=, rdnoise=, saturation=, flatname=, flatlist=, illumname=, badpixfile=, $'
       print, '   prefix=, /makeflat, /clean, /wfits, /doplot, /gzip, objlist=, _extra=extra'
       return
    endif
    
    fcount = n_elements(proclist)
    nbias  = n_elements(biasfile) & if (nbias gt 0L) then if (biasfile[0] eq '') then nbias = 0L
    ndome  = n_elements(domefile) & if (ndome gt 0L) then if (domefile[0] eq '') then ndome = 0L
    nsky   = n_elements(skyfile)  & if (nsky  gt 0L) then if (skyfile[0]  eq '') then nsky  = 0L
    ndark  = n_elements(darkfile) & if (ndark gt 0L) then if (darkfile[0] eq '') then ndark = 0L
    
; read the bad pixel file
    
    if n_elements(badpixfile) ne 0L then begin

       badpixname = file_test(datapath+badpixfile,/regular)

       if (badpixname eq 0L) or (strcompress(badpixfile,/remove) eq '') then $
         splog, 'Bad pixel file '+badpixfile+' not found.' else begin

          readcol, datapath+badpixfile, x1, x2, y1, y2, format='L,L,L,L', /silent, comment='#'
          splog, 'Identified '+strn(n_elements(x1))+' rows in the bad pixel file '+badpixfile+'.'

; IDL is zero-indexed but the bad pixel file is one-indexed (IRAF
; convention)          
          
          xybad = long(transpose([[x1],[x2],[y1],[y2]]))-1L 

       endelse
       
    endif

; ----------------------------------------------------------------------
; basic processing; store the data, the bias frame, dome flat, sky
; flat, and dark frame in different data cubes
; ----------------------------------------------------------------------

    splog, 'Initializing the output data structure.'

    if fcount ne 0L then begin
       cube = initialproc(proclist,datapath=datapath,trim=trim,$
         overscan=overscan,norder_overscan=norder_overscan,$
         xybad=xybad,_extra=extra)
       if size(cube,/type) ne 8L then begin
          splog, 'Error processing PROCLIST.'
          return
       endif
    endif

    if nbias ne 0L then begin
       biascube = initialproc(biasfile,datapath=datapath,trim=trim,$
         overscan=overscan,norder_overscan=norder_overscan,$
         xybad=xybad,_extra=extra,/noheaderupdate)
       if size(biascube,/type) ne 8L then begin
          splog, 'Error initializing BIASCUBE.'
          return
       endif
    endif
    
    if ndome ne 0L then begin
       domecube = initialproc(domefile,datapath=datapath,trim=trim,$
         overscan=overscan,norder_overscan=norder_overscan,$
         xybad=xybad,_extra=extra,/noheaderupdate)
       if size(domecube,/type) ne 8L then begin
          splog, 'Error initializing DOMECUBE.'
          return
       endif
    endif

    if nsky ne 0L then begin
       skycube = initialproc(skyfile,datapath=datapath,trim=trim,$
         overscan=overscan,norder_overscan=norder_overscan,$
         xybad=xybad,_extra=extra)
       if size(skycube,/type) ne 8L then begin
          splog, 'Error initializing SKYCUBE.'
          return
       endif
    endif

    if ndark ne 0L then begin
       darkcube = initialproc(darkfile,datapath=datapath,trim=trim,$
         overscan=overscan,norder_overscan=norder_overscan,$
         xybad=xybad,_extra=extra,/noheaderupdate)
       if size(darkcube,/type) ne 8L then begin
          splog, 'Error initializing DARKCUBE.'
          return
       endif
    endif

; ----------------------------------------------------------------------
; bias subtraction (should the dome flat be bias-subtracted? [probably
; doesn't matter])
; ----------------------------------------------------------------------
    
    if nbias eq 0L then begin

       splog, 'No bias frame passed...skipping bias subtraction.'

    endif else begin

       splog, 'Subtracting the bias frame '+biascube.fname+'.'
       
       bias = biascube.image
       biasname = biascube.fname
       
       if nsky ne 0L then begin
          skycube.image = skycube.image - bias
          skycube = head_append(skycube,'IBIAS',biasname,' bias frame',before='HISTORY')
       endif

       if ndark ne 0L then begin
          darkcube.image = darkcube.image - bias
          darkcube = head_append(darkcube,'IBIAS',biasname,' bias frame',before='HISTORY')
       endif

       if fcount ne 0L then begin
          for i = 0L, fcount-1L do begin
             cube[i].image = cube[i].image - bias
             cube[i] = head_append(cube[i],'IBIAS',biasname,' bias frame',before='HISTORY')
          endfor
        endif

     endelse 

; ----------------------------------------------------------------------
; dark current subtraction
; ----------------------------------------------------------------------
    
    if ndark eq 0L then begin

       splog, 'No dark frame passed...skipping dark subtraction.'

    endif else begin

       splog, 'Subtracting the dark frame '+darkcube.fname+'.'

       dark = darkcube.image
       darkname = darkcube.fname
       dexptime = sxpar(darkcube.header,'DARKTIME')
       
       if nsky ne 0L then begin
          exptime = sxpar(skycube.header,'EXPTIME')
          skycube.image = skycube.image - dark * exptime / dexptime
          skycube = head_append(skycube,'IDARK',darkname,' dark frame',before='HISTORY')
       endif

       if fcount ne 0L then begin
          for i = 0L, fcount-1L do begin
             exptime = sxpar(cube[i].header,'EXPTIME')
             cube[i].image = cube[i].image - dark * exptime / dexptime
             cube[i] = head_append(cube[i],'IDARK',darkname,' dark frame',before='HISTORY')
          endfor
        endif

     endelse 

; ----------------------------------------------------------------------
; constructing the master flat field
; ----------------------------------------------------------------------

    if keyword_set(makeflat) then begin ; interactively generate a flat field
    
       splog, 'Generating the flat field '+flatname+'.'

       if ndome eq 0L then begin
          splog, 'No dome flat has been provided!'
          icleanup, cube, /main & icleanup, biascube, /main
          icleanup, domecube, /main & icleanup, skycube, /main
          icleanup, darkcube, /main
       endif

; repair bad pixels in the dome flat.  DJS_MASKINTERP will interpolate
; where all mask values are not equal to zero

       for k = 0L, ndome-1L do begin
       
          dmask = domecube[k].mask
;         dmask = domecube[k].mask eq 0B

          domecube[k].image = djs_maskinterp(domecube[k].image,dmask,iaxis=0L)
          domecube[k].vmap = djs_maskinterp(domecube[k].vmap,dmask,iaxis=0L)

       endfor

; create and write to disk the master flat field.         

       flat = imasterflat(domecube.image,domecube.header,datapath=datapath,$
         flatname=flatname,flatheader=flatheader,/write,doplot=doplot,_extra=extra)

; illumination flat.  flat-field the sky flat with the dome flat that
; is nearest in RA to the sky flat itself
       
       if nsky ne 0L then begin

          indx = choose_flat(skycube.header,flatheader)
          goodflat = flat[*,*,indx]

; repair bad pixels in the sky flat

          smask = skycube.mask
;         smask = skycube.mask eq 0B
          
          skycube.image = djs_maskinterp(skycube.image,smask,iaxis=0L)
          skycube.vmap = djs_maskinterp(skycube.vmap,smask,iaxis=0L)
          
          illumflat = iillumflat(skycube.image,goodflat,skycube.fname,datapath=datapath, $
            illumname=illumname,illumheader=illumheader,/write,doplot=doplot,_extra=extra)

       endif else begin

          illumflat = make_array(domecube[0].naxis1,domecube[0].naxis2,value=1)
          illumname = ''

       endelse

       nflat = ndome
       
    endif else begin

; read the flat fields

       splog, 'Reading the flat field(s).'
       flatcube = im_fits_cube(flatname,path=datapath)

       if size(flatcube,/type) ne 8L then begin
          splog, 'Flat field(s) not found...returning.'

          icleanup, cube, /main & icleanup, biascube, /main
          icleanup, domecube, /main & icleanup, skycube, /main
          icleanup, darkcube, /main
       endif
       
       nflat = n_elements(flatcube)

       flat = flatcube.image
       for i = 0L, nflat-1L do if i eq 0L then flatheader = flatcube[i].header else $ ; this might choke jm02oct25uofa
         flatheader = [ [flatheader], [flatcube[i].header] ]

       if file_test(datapath+illumname,/regular) then begin
          splog, 'Reading '+datapath+illumname+'.'
          illumflat = readfits(datapath+illumname,illumheader,/silent)
       endif else begin
          splog, 'Illumination flat '+datapath+illumname+' not found...skipping illumination correction.'
          illumflat = flat[*,*,0]*0.0+1.0
          illumname = ''
       endelse

    endelse

; generate the normalized master flat field(s), which includes the 
; illumination correction
    
    masterflat = flat*0.0
    for j = 0L, nflat-1L do begin
       masterflat[*,*,j] = flat[*,*,j] * illumflat
       masterflat[*,*,j] = masterflat[*,*,j] / mean(masterflat[*,*,j])
    endfor
       
; ----------------------------------------------------------------------
; initialize the variance maps
; ----------------------------------------------------------------------

    if fcount ne 0L then for i = 0L, fcount-1L do cube[i].vmap = $
      vmap_init(cube[i].image,gain=gain,rdnoise=rdnoise)

; ----------------------------------------------------------------------
; flat-fielding
; ----------------------------------------------------------------------

; divide the objects and the variance maps by the flat field

    if fcount ne 0L then begin

       splog, 'Flat-fielding the objects.'
       for j = 0L, fcount-1L do begin

          indx = choose_flat(cube[j].header,flatheader)
          goodflat = masterflat[*,*,indx]
          
;         print, j+1L, strupcase(cube[objindx[j]].object), format='(I0,A20)'
          cube[j].image = cube[j].image / goodflat   ; object image
          cube[j].vmap = cube[j].vmap / (goodflat^2) ; variance map

          cube[j] = head_append(cube[j],'IFLAT',flatname[indx],$
            ' dome/quartz flat field',before='HISTORY')
          if (illumname ne '') then cube[j] = head_append(cube[j],'ISFLAT',$
            illumname,' illumination flat field',before='HISTORY')

       endfor 

    endif else splog, 'No objects in the file list.'

; flat-field the sky image
    
    if nsky ne 0L then begin

       indx = choose_flat(skycube.header,flatheader)
       goodflat = flat[*,*,indx]
          
       skycube.image = skycube.image / goodflat   ; object image
       skycube.vmap = skycube.vmap / (goodflat^2) ; variance map

       skycube = head_append(skycube,'IFLAT',flatname[indx],$
         ' flat field',before='HISTORY')
       if (illumname ne '') then skycube = head_append(skycube,'ISFLAT',$
         illumname,' sky flat',before='HISTORY')
       
    endif
       
; ----------------------------------------------------------------------
; write out
; ----------------------------------------------------------------------

; QA plot - histogram of the overscan-subtracted and trimmed bias
;           frame

    if (nbias ne 0L) and keyword_set(makeflat) then begin
    
       qaplotname = 'qaplot_bias_'+strmid(biasname,0,strpos(biasname,'.fits'))+'.ps'

       if keyword_set(doplot) then iter = 1L else iter = 0L

       djs_iterstat, bias, sigrej=3.0, median=biasmed, sigma=biassig
       xrange = biasmed+[-5*biassig,+10*biassig]
          
       for i = 0L, iter do begin

          if (i eq 0L) then begin
             dfpsplot, datapath+qaplotname, /square, /color
             postthick = 8.0
          endif else begin
             im_window, 3, xratio=0.3, /square
             postthick = 2.0
          endelse

          djs_iterstat, bias, sigma=bsigma, mean=bmean, median=bmedian, sigrej=3.0

          stats = strtrim(string(bmean,format='(F12.3)'),2)+' \pm '+$
            strtrim(string(bsigma,format='(F12.3)'),2)+$
            ' ('+strtrim(string(bmedian,format='(F12.3)'),2)+')'
    
          plothist, bias, xbin, ybin, bin=bsigma, /noplot
          plothist, bias, bin=bsigma, xstyle=3, ystyle=1, xthick=postthick, ythick=postthick, $
            charthick=postthick, charsize=1.3, thick=postthick, xtitle='Counts', $
            ytitle='Number', xrange=xrange, ymargin=[5,2], /fill, /fline, fcolor=djs_icolor('grey'), $
            title='Overscan-Subtracted Bias', subtitle='['+repstr(biasfile,'.fits','')+']', $
            yrange=minmax(ybin)*[1.0,1.15], forientation=45, fspacing=0.05, color=djs_icolor('')
          legend, textoidl(stats), /right, /top, charsize=1.3, charthick=postthick, box=0
          
          if (i eq 0L) then dfpsclose
          
       endfor

    endif 

    if keyword_set(wfits) then begin

; write objects (objects and comps) as multi-extension FITS files

       if fcount ne 0L then begin

          objlist = strarr(fcount) ; output files
          
          for j = 0L, fcount-1L do begin
    
; linearly interpolate over dead pixels and bad columns along the
; x-axis.  

; also identify pixels with negative variance and replace them with
; the read noise

             outimage = cube[j].image
             outmask = cube[j].mask
             outvmap = cube[j].vmap
             slash = strsplit(cube[j].fname, '/')
             lastslash = slash[n_elements(slash)-1]
             if lastslash eq 0 then outname = prefix+cube[j].fname $
             else outname = prefix+strmid(cube[j].fname,lastslash,strlen(cube[j].fname)-lastslash+1)
             header = cube[j].header
             
             neg = where(outvmap lt float(0),nneg)
             if nneg ne 0L then begin
                outvmap[neg] = (rdnoise/gain)^2 ; [ADU^2]
;               outmask[neg] = 0B
             endif
             
             if keyword_set(clean) then begin

                outimage = djs_maskinterp(outimage,outmask,iaxis=0L) 
                outvmap = djs_maskinterp(outvmap,outmask,iaxis=0L) ; <-- optional!
;               outimage = djs_maskinterp(outimage,(outmask eq 0B),iaxis=0L) 
;               outvmap = djs_maskinterp(outvmap,(outmask eq 0B),iaxis=0L) ; <-- optional!
             
             endif

             objlist[j] = outname
             splog, 'Writing '+datapath+outname+'.'

; image contains extensions; note that the FITS standard requires
; EXTEND to follow the last NAXISi entry (thanks to Kris Eriksen,
; 2003-Oct-15)              
             
             sxaddpar, header, 'EXTEND', 'T', after='NAXIS2' 
             
             wrt2dspec, outname, float(outimage), float(sqrt(outvmap)), outmask, $
               header, datapath=datapath, gzip=gzip

          endfor 

; crop the output list to include only object names (remove arcs)
          
          w = where((strupcase(cube.type) eq 'OBJECT') and (strmatch(cube.object,'*sky*') ne 1B),nobjects)
          if nobjects ne 0L then objlist = objlist[w] else objlist = ''
          
       endif

    endif
       
    icleanup, cube & icleanup, biascube & icleanup, domecube
    icleanup, skycube & icleanup, darkcube
    
return
end
