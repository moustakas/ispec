;+
; NAME:
;	ICRCOMBINE
;
; PURPOSE:
;	Combine cosmic-ray split images using CR_REJECT.
;
; CALLING SEQUENCE:
;       icrcombine, crimlist, combined_image, combined_noise, $
;          mastermask, header, datapath=, nsig=, dfactor=, $
;          user_pixshift=, crimname=, _extra=extra, /simple, $
;          /gzip, /verbose, /find_pixshift, /debug, /wfits
;
; INPUTS:
;	crimlist   - FITS file list of cosmic-ray split images 
;
; OPTIONAL INPUTS:
;	datapath   - path to the data
;	nsig       - see CR_REJECT documentation
;	dfactor    - see CR_REJECT documentation
;       extra      - keywords for VMAP_INIT() and CR_REJECT 
;       user_pixshift - if FIND_PIXSHIFT=1 then pixel shifts that have
;                       been determined manually can be passed via
;                       this keyword [NIMAGE] 
;	
; KEYWORD PARAMETERS:
;       simple        - combine simple FITS files, not iSPEC2d format 
;                       files; the proper gain and rdnoise should be 
;                       passed through _EXTRA to VMAP_INIT()
;       gzip          - write GZIPPED FITS files
;	verbose       - (lots) of verbose output, mostly from CR_REJECT
;       find_pixshift - solve for an integer pixel shift between every
;                       image (relative to the zeroth image) using
;                       cross-correlation  
;       debug         - inspect QA plots showing the derived pixel
;                       shifts (used only if FIND_PIXSHIFT=1); this
;                       keyword overwrites WFITS, unless SIMPLE=1
;	wfits         - write the cleaned image, sigma map, sky image,
;                       and updated bad pixel mask to disk 
;
; OUTPUTS:
;	combined_image - cosmic-ray rejected combined image, scaled by
;                        the total exposure time
;	combined_noise - combined noise image
;	mastermask     - master bad pixel mask, including input bad
;                        pixel mask and cosmic-ray mask
;	header         - output FITS header for COMBINED_IMAGE
;
; OPTIONAL OUTPUTS:
;	crimname - name of the output FITS file if WFITS=1
;
; COMMENTS:
;       Cosmic-ray pixels are assigned a bad pixel mask value
;       according to IMASK_BITS().  Recall that 0 is a good pixel. 
;
;	If only one image is given then the routine returns that image
;	and exits.
;
;       Here are some general issues to consider when combining 2D
;       spectra.  The images should be sequential in time.  Images
;       that are spaced far apart in time should not be combined due
;       to variations in atmospheric transparency and pointing
;       differences.  This routine will print to the screen some
;       useful information to decide whether the images should really
;       be combined (i.e., date, time, parallactic angle, position
;       angle, airmass, etc.).  For example, if the position angle of
;       one object is different then it should not be combined.
;
;       When combining the images we use the zeroth header as a
;       template for the output image header.  All the internal
;       iSPEC2d header entries are handled properly (kind of).  We
;       assume that OBJECT, IMAGETYP, DATE-OBS, EPOCH, GAIN, RDNOISE,
;       and all other instrument configuration parameters are the same
;       for all the images.  The output EXPTIME is the *total*
;       exposure time of all the images.
;
;       The following quantities are "mean" quantities of all the
;       observations:  RA, DEC, AIRMASS, PARANGLE, and JD.  "Time"
;       quantities such as UT, UTMIDDLE, and ST are not properly
;       treated and are taken from the zeroth image.
;
;       The output FITS name is based on the zeroth image and the
;       number of images combined.  For example, if we combined
;       'spec1.fits', 'spec2.fits', and 'spec3.fits' then the output
;       name will be 'spec1_3.fits', where 3 is the number of images
;       combined. 
;
;       We assume that the sky subtraction parameters (from
;       ISKYSUBTRACT2D) are the same for both images.
;
;       If FIND_PIXSHIFT=1 and the pixel shifts are non-zero, then we
;       copy rows without data with the corresponding rows in the
;       zeroth image.  Therefore, these rows will not be "correct" in
;       the output image, but it should not matter for small pixel
;       shifts.  This procedure is preferable to cropping the output
;       image. 
;
; EXAMPLE:
;
; TODO:
;       Handle mean "time" quantities in the output header (see
;       COMMENTS, above).
;
; INTERNAL SUPPORT ROUTINES:
;       FIND_PIXSHIFT()
;
; PROCEDURES USED:
;	RD2DSPEC(), SXPAR(), MAKE_CUBE(), CR_REJECT, CWD(), SPLOG,
;	DJS_MEAN(), SXADDPAR, SXADDHIST, SXDELPAR, WRT2DSPEC,
;	ICLEANUP, IFORAGE(), IM_HMS2DEC(), IM_DEC2HMS(), IMASK_BITS(),
;	VMAP_INIT(), IM_WINDOW, DJS_PLOT, LEGEND, IM_ZCOMPUTE(),
;	DJS_MASKINTERP() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 July 8, U of A
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03dec07uofa - added GZIP keyword; properly handle the
;                       combined sky images 
;       jm03dec22uofa - added SIMPLE and _EXTRA keywords
;       jm05jun17uofa - read and average MJD-OBS keyword; changed the
;                       default NSIG parameter to increase the number
;                       of iterations
;       jm05jun28uofa - wrote FIND_PIXSHIFT(); added DEBUG keyword and
;                       USER_PIXSHIFT optional input
;       jm05jul20uofa - bug fix: if CR_REJECT rejects all pixels from
;                       an image stack, then interpolate smoothly over
;                       those pixels (previously, CR_REJECT was
;                       setting the error to 0.0, which gave
;                       infinities when sky subtracting, etc.) 
;      
; Copyright (C) 2002-2003, 2005, John Moustakas
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

function make_cube, incube, sigmamap=sigmamap, sky=sky, mask=mask
; jm0feb14uofa
; jm03dec8uofa - added SKY keyword

    nim = (size(incube))[3]
    imsize = size(incube[0].image,/dimension)
    xsize = imsize[0]
    ysize = imsize[1]

    if keyword_set(mask) then $
      outcube = make_array(xsize,ysize,nim,/integer) else $
      outcube = make_array(xsize,ysize,nim,/float)

    if keyword_set(sigmamap) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].sigmamap
    if keyword_set(sky) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].sky 
    if keyword_set(mask) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].mask 

    if (not keyword_set(sigmamap)) and (not keyword_set(sky)) and (not keyword_set(mask)) then $
      for k = 0L, nim-1L do outcube[*,*,k] = incube[k].image

return, outcube
end

function find_pixshift, input_cube, noise_cube, nlags=nlags, $
  user_pixshift=user_pixshift, crimlist=crimlist, debug=debug

    imsize = size(input_cube,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    nimage = imsize[2]
    midcol = ncols / 2L
    pixshift = fltarr(nimage)
    rowaxis = lindgen(nrows)

    nmed = 10L
    nfrac = fix(0.01*ncols)
    
    if (n_elements(user_pixshift) eq 0L) then user = 0L else begin
       if (n_elements(user_pixshift) ne nimage) then begin
          splog, 'Dimensions of USER_PIXSHIFT and CRIMLIST do not agree.'
          user = 0L
       endif else user = 1L
    endelse
    
    if (n_elements(nlags) eq 0L) then nlags = 50L
    lags = [-reverse(lindgen(nlags)),lindgen(nlags)]

    bigindx = lindgen(nrows+2*nlags)
    bigrefspec = fltarr(nrows+2*nlags) ; pad by zeros
    bigspec = bigrefspec*0.0
    bigivar = bigrefspec*0.0

    refspec = djs_median(input_cube[(midcol-nmed)>0L:(midcol+nmed)<(ncols-1L),*,0],1)
;   refspec = total(input_cube[*,*,0],1)/float(nrows)
;   refspec = im_normalize(total(input_cube[*,*,0],1),/max)
    bigrefspec[nlags:nrows+nlags-1L] = refspec

    objlist = repstr(crimlist,'.fits','')
    
    for ii = 1L, nimage-1L do begin
    
       spec = djs_median(input_cube[(midcol-nmed)>0L:(midcol+nmed)<(ncols-1L),*,ii],1)
       ivar = djs_median(noise_cube[(midcol-nmed)>0L:(midcol+nmed)<(ncols-1L),*,ii],1)
;      spec = total(input_cube[*,*,ii],1)/float(nrows)
;      ivar = sqrt(total(noise_cube[*,*,ii]^2,1)/float(nrows))
;      spec = im_normalize(total(input_cube[*,*,ii],1),/max,norm=norm)

       bigspec[nlags:nrows+nlags-1L] = spec
       bigivar[nlags:nrows+nlags-1L] = ivar

       if user then pixshift[ii] = user_pixshift[ii] else begin
          if keyword_set(debug) then im_window, 0, xratio=0.4, /square
          pans = im_zcompute(bigspec,bigivar,bigrefspec,nfind=1L,$
            pspace=1L,pmin=min(lags),pmax=max(lags),bestlag=bestlag,$
            doplot=debug)
          if (pans.z_err gt 0.0) then pixshift[ii] = bestlag[0]
       endelse

; compute the mean upper, lower, and central profiles; plot the mean
; ratio at each of these three positions

       profile_cen = total(input_cube[midcol-nfrac/2L:midcol+nfrac/2L,*,ii],1)/float(nfrac)
       profile_lo = total(input_cube[0L:nfrac-1L,*,ii],1)/float(nfrac)
       profile_up = total(input_cube[ncols-nfrac:ncols-1L,*,ii],1)/float(nfrac)

       refprofile_cen = total(input_cube[midcol-nfrac/2L:midcol+nfrac/2L,*,0],1)/float(nfrac)
       refprofile_lo = total(input_cube[0L:nfrac-1L,*,0],1)/float(nfrac)
       refprofile_up = total(input_cube[ncols-nfrac:ncols-1L,*,0],1)/float(nfrac)

       good_cen = where((profile_cen gt 0.0) and (refprofile_cen gt 0.0))
       ratio_cen = refprofile_cen[good_cen]/shift(profile_cen[good_cen],pixshift[ii])
       stats_cen = strtrim(string(median(ratio_cen),format='(F12.2)'),2)+' ('+$
         strtrim(string(mean(ratio_cen),format='(F12.2)'),2)+'\pm'+$
         strtrim(string(stddev(ratio_cen),format='(F12.2)'),2)+')'

       good_lo = where((profile_lo gt 0.0) and (refprofile_lo gt 0.0))
       ratio_lo = refprofile_lo[good_lo]/shift(profile_lo[good_lo],pixshift[ii])
       stats_lo = strtrim(string(median(ratio_lo),format='(F12.2)'),2)+' ('+$
         strtrim(string(mean(ratio_lo),format='(F12.2)'),2)+'\pm'+$
         strtrim(string(stddev(ratio_lo),format='(F12.2)'),2)+')'

       good_up = where((profile_up gt 0.0) and (refprofile_up gt 0.0))
       ratio_up = refprofile_up[good_up]/shift(profile_up[good_up],pixshift[ii])
       stats_up = strtrim(string(median(ratio_up),format='(F12.2)'),2)+' ('+$
         strtrim(string(mean(ratio_up),format='(F12.2)'),2)+'\pm'+$
         strtrim(string(stddev(ratio_up),format='(F12.2)'),2)+')'

; debug plot       
       
       if keyword_set(debug) then begin

          title = objlist[ii]+' shifted to '+objlist[0]

          im_window, 2, xratio=0.7, yratio=0.5
;         im_window, 2, xratio=0.4, /square

; center          

          yrange = [min(refprofile_cen[good_cen])<min(profile_cen[good_cen]),$
            max(refprofile_cen[good_cen])>max(profile_cen[good_cen])]

          plot, rowaxis, refprofile_cen, ps=10, xsty=1, ysty=11, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='', xtickname=replicate(' ',10), ytitle='Counts', $
            title=title, thick=1.0, position=[0.11,0.45,0.89,0.9]
          axis, /yaxis, yrange=yrange, charsize=1.5, charthick=2.0, ytitle='Counts'
          djs_oplot, rowaxis, shift(profile_cen,pixshift[ii]), ps=10, color='red', thick=1.0
          djs_oplot, rowaxis, profile_cen, ps=3, color='yellow', thick=0.1
          legend, 'Best Shift = '+string(pixshift[ii],format='(I0)'), /right, $
            /top, box=0, charsize=1.5, charthick=2.0
          legend, 'Center Column', /left, /top, box=0, charsize=1.5, charthick=2.0

; lower          

          yrange = [min(refprofile_lo[good_lo])<min(profile_lo[good_lo]),$
            max(refprofile_lo[good_lo])>max(profile_lo[good_lo])]

          plot, rowaxis, refprofile_lo, ps=10, xsty=1, ysty=3, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='Row Number', ytitle='Counts', $
            thick=1.0, position=[0.11,0.1,0.49,0.45], /noerase
          legend, 'Lower Column', /left, /top, box=0, charsize=1.5, charthick=2.0
          djs_oplot, rowaxis, shift(profile_lo,pixshift[ii]), ps=10, color='red', thick=1.0
          djs_oplot, rowaxis, profile_lo, ps=3, color='yellow', thick=0.1

; upper          

          yrange = [min(refprofile_up[good_up])<min(profile_up[good_up]),$
            max(refprofile_up[good_up])>max(profile_up[good_up])]

          plot, rowaxis, refprofile_up, ps=10, xsty=1, ysty=11, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='Row Number', ytitle='', $
            thick=1.0, position=[0.49,0.1,0.89,0.45], /noerase, ytickname=replicate(' ',10)
          axis, /yaxis, yrange=yrange, charsize=1.5, charthick=2.0, ytitle='Counts'
          legend, 'Upper Column', /left, /top, box=0, charsize=1.5, charthick=2.0
          djs_oplot, rowaxis, shift(profile_up,pixshift[ii]), ps=10, color='red', thick=1.0
          djs_oplot, rowaxis, profile_up, ps=3, color='yellow', thick=0.1

; ------------          
; SCALE FACTOR
; ------------          
          
          im_window, 1, xratio=0.7, yratio=0.5
          
; center          

          yrange = minmax(ratio_cen)*[1.0,1.1]

          plot, rowaxis, ratio_cen, ps=10, xsty=1, ysty=11, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='', xtickname=replicate(' ',10), ytitle='Ratio', $
            title=title, thick=1.0, position=[0.11,0.45,0.89,0.9]
          axis, /yaxis, yrange=yrange, charsize=1.5, charthick=2.0, ytitle='Ratio'
          legend, 'Center Column', /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, textoidl(stats_cen), /right, /top, box=0, charsize=1.5, charthick=2.0

; lower          

          yrange = minmax(ratio_lo)*[1.0,1.1]

          plot, rowaxis, ratio_lo, ps=10, xsty=1, ysty=3, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='Row Number', ytitle='Ratio', $
            thick=1.0, position=[0.11,0.1,0.49,0.45], /noerase
          legend, 'Lower Column', /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, textoidl(stats_lo), /right, /top, box=0, charsize=1.5, charthick=2.0

; upper          

          yrange = minmax(ratio_up)*[1.0,1.1]

          plot, rowaxis, ratio_up, ps=10, xsty=1, ysty=11, yrange=yrange, $
            xthick=2, ythick=2, charsize=1.5, xrange=[0,nrows]+pixshift[ii], $
            charthick=2.0, xtitle='Row Number', ytitle='', $
            thick=1.0, position=[0.49,0.1,0.89,0.45], /noerase, ytickname=replicate(' ',10)
          axis, /yaxis, yrange=yrange, charsize=1.5, charthick=2.0, ytitle='Ratio'
          legend, 'Upper Column', /left, /top, box=0, charsize=1.5, charthick=2.0
          legend, textoidl(stats_up), /right, /top, box=0, charsize=1.5, charthick=2.0
          
;         yrange = [min(bigrefspec[good])<min(bigspec[good]),max(bigrefspec[good])>max(bigspec[good])]
;         plot, bigindx, bigrefspec, ps=10, xsty=1, ysty=3, yrange=yrange, $ ;yrange=[0,1.1], $
;           xthick=2, ythick=2, charsize=1.5, xrange=[nlags,nlags+nrows]+pixshift[ii], $
;           charthick=2.0, xtitle='Row Number + '+string(nlags,format='(I0)'), $
;           ytitle='Counts', title=title, thick=1.0
;         djs_oplot, bigindx, shift(bigspec,pixshift[ii]), ps=10, color='red', thick=1.0
;         djs_oplot, bigindx, bigspec, ps=3, color='yellow', thick=0.1
                    
          splog, 'Press any key to continue.'
          cc = get_kbrd(1)

       endif
       
    endfor

return, pixshift
end

pro icrcombine, crimlist, combined_image, combined_noise, mastermask, header, $
  datapath=datapath, nsig=nsig, dfactor=dfactor, user_pixshift=user_pixshift, $
  crimname=crimname, _extra=extra, simple=simple, gzip=gzip, verbose=verbose, $
  find_pixshift=find_pixshift, debug=debug, wfits=wfits
    
    if (n_elements(crimlist) eq 0L) then begin
       print, 'Syntax - icrcombine, crimlist, combined_image, combined_noise, $'
       print, '   mastermask, header, datapath=, nsig=, dfactor=, user_pixshift=, $'
       print, '   crimname=, _extra=extra, /simple, /gzip, /verbose, /find_pixshift, $'
       print, '   /debug, /wfits'
       return
    endif
    
    nimage = n_elements(crimlist)
    case nimage of
       0L: begin
          splog, 'CRIMLIST not defined.'
          return
       end
       1L: begin
          splog, 'Only one image passed...returning.'
          return
       end
       else:
    endcase
    
    if (size(crimlist[0],/type) ne 7L) or (strmatch(crimlist[0],'*.fits*') eq 0B) then begin
       splog, 'Files must be type string FITS files.'
       return
    endif

; set defaults
    
    if not keyword_set(datapath) then datapath = cwd()

    if n_elements(nsig) eq 0L then nsig = [8,6,4,3,3,3] ; [8,5,4,3] ; [4,3,3,3]
    if n_elements(dfactor) eq 0L then dfactor = 0.3
    
; use the zeroth image to form the output name

    if n_elements(crimname) eq 0L then $
      crimname = strmid(crimlist[0],0,strpos(crimlist[0],'.fits'))+'_'+string(nimage,format='(I0)')+'.fits'

    if keyword_set(simple) then begin

       info = iforage(datapath+crimlist)
;      header = headfits(crimlist[0])
       exptime = info.exptime

       input_cube = make_array(info[0].naxis1,info[0].naxis2,nimage,/float)
       noise_cube = make_array(info[0].naxis1,info[0].naxis2,nimage,/float)

       for i = 0L, nimage-1L do begin
          input_cube[*,*,i] = readfits(datapath+crimlist[i],header,/silent)
          noise_cube[*,*,i] = vmap_init(input_cube[*,*,i],_extra=extra)
       endfor

       noise_cube = sqrt(noise_cube)

       cr_reject, input_cube, -1, 0, 0, 1.0, combined_image, combined_noise, $
         combined_npix, mask_cube=mask_cube, exptime=exptime, verbose=verbose, $
         nsig=nsig, dilation=dilation, dfactor=dfactor, /noclearmask, /restore_sky, $
         /noskyadjust, weighting=0, noise_cube=noise_cube

       mkhdr, header, combined_image
       sxaddpar, header, 'INCOMB', nimage, ' number of images stacked in ICRCOMBINE'

       if keyword_set(wfits) then begin
          
          splog, 'Writing '+datapath+crimname+'.'
          writefits, datapath+crimname, float(combined_image), header
          if keyword_set(gzip) then spawn, ['gzip -f '+datapath+crimname], /sh
          
       endif
       return
       
    endif

    if keyword_set(debug) then wfits = 0L ; NOTE!

; read the data cube
    
    cube = rd2dspec(crimlist,datapath=datapath)
    forage = iforage(datapath+crimlist)
    ncols = cube[0].naxis1
    nrows = cube[0].naxis2

    header = *cube[0].header ; output header

; image exposure time

    exptime = dblarr(nimage)+1.0
    for i = 0L, nimage-1L do if sxpar(*cube[i].header,'EXPTIME') ne 0L then $
      exptime[i] = sxpar(*cube[i].header,'EXPTIME')
    time = total(exptime)

    oexptime = sxpar(*cube[0].header,'OEXPTIME',count=oexptime_count) ; original exposure time
    if (oexptime_count ne 0L) then begin
       for i = 1L, nimage-1L do oexptime = [oexptime,sxpar(*cube[i].header,'OEXPTIME')]
       oexptime_total = total(oexptime)
    endif
    
; input data and noise cubes

    input_cube = make_cube(cube)
    noise_cube = make_cube(cube,/sigmamap)
    sky_cube = make_cube(cube,/sky)
    mask_cube = make_cube(cube,/mask)

    inputmask = mask_cube*0B+1B ; give CR_REJECT a mask where all pixels are good

;    if keyword_set(nonegative) then begin
;
;       for j = 0L, nimage-1L do begin
;          msk = inputmask[*,*,j]
;          neg = where(input_cube[*,*,j] lt float(0),nneg)
;          if nneg ne 0L then msk[neg] = 0B
;          outimage = djs_maskinterp(outimage,(outmask eq 0B),iaxis=0L) 
;
;          inputmask[*,*,j] = msk
;       endfor
;
;    endif

; solve for an integer pixel shift between images, if requested    
    
    if keyword_set(find_pixshift) then begin

       pixshift = find_pixshift(input_cube,noise_cube,nlags=nlags,$
         user_pixshift=user_pixshift,crimlist=crimlist,debug=debug)
       splog, 'Best pixel shifts: '+strjoin(string(pixshift,format='(I0)'),', ')+' pixels.'

       for ii = 1L, nimage-1L do begin
          
; replace rows without data with the pixel values in the zeroth image;
; if the pixel shift is zero then no replacement is necessary
          
          shift_image = shift(input_cube[*,*,ii],0,pixshift[ii])
          shift_noise = shift(noise_cube[*,*,ii],0,pixshift[ii])
          shift_sky   = shift(sky_cube[*,*,ii],0,pixshift[ii])
          shift_mask  = shift(mask_cube[*,*,ii],0,pixshift[ii])

          if (pixshift[ii] lt 0.0) then begin

             shift_image[*,pixshift[ii]+nrows-1L:nrows-1L] = input_cube[*,pixshift[ii]+nrows-1L:nrows-1L,0]
             shift_noise[*,pixshift[ii]+nrows-1L:nrows-1L] = noise_cube[*,pixshift[ii]+nrows-1L:nrows-1L,0]
             shift_sky[*,pixshift[ii]+nrows-1L:nrows-1L]   = sky_cube[*,pixshift[ii]+nrows-1L:nrows-1L,0]
             shift_mask[*,pixshift[ii]+nrows-1L:nrows-1L]  = mask_cube[*,pixshift[ii]+nrows-1L:nrows-1L,0]

          endif
             
          if (pixshift[ii] gt 0.0) then begin

             shift_image[*,0L:pixshift[ii]] = input_cube[*,0L:pixshift[ii],0]
             shift_noise[*,0L:pixshift[ii]] = noise_cube[*,0L:pixshift[ii],0]
             shift_sky[*,0L:pixshift[ii]]   = sky_cube[*,0L:pixshift[ii],0]
             shift_mask[*,0L:pixshift[ii]]  = mask_cube[*,0L:pixshift[ii],0]

          endif

; restore the data cubes          
          
          input_cube[*,*,ii] = shift_image
          noise_cube[*,*,ii] = shift_noise
          sky_cube[*,*,ii] = shift_sky
          mask_cube[*,*,ii] = shift_mask

       endfor

; the code below crops the rows without data in all NIMAGE images; we
; do not want to crop, however, because we would like all the images
; in each run to be the same size
       
;      neg = where(pixshift lt 0.0,nneg)
;      if (nneg ne 0L) then negshift = long(min(pixshift)) else negshift = 0L
;      pos = where(pixshift gt 0.0,npos)
;      if (npos ne 0L) then posshift = long(max(pixshift)) else posshift = 0L
       
;      input_cube = input_cube[*,posshift:nrows+negshift-1L,*]
;      noise_cube = noise_cube[*,posshift:nrows+negshift-1L,*]
;      sky_cube = sky_cube[*,posshift:nrows+negshift-1L,*]
;      mask_cube = mask_cube[*,posshift:nrows+negshift-1L,*]
;      nrows = long(nrows + negshift - posshift)

; test that the shifts worked
       
;      print, find_pixshift(input_cube,nlags=nlags,crimlist=crimlist,/debug)

       if keyword_set(debug) then begin
          icleanup, cube
          return ; NOTE!                
       endif
       
    endif else pixshift = fltarr(nimage)

; scale the combined sky spectrum simply with the exposure time

    for i = 0L, nimage-1L do sky_cube[*,*,i] = sky_cube[*,*,i]/exptime[i]
    combined_sky = time * total(sky_cube,3) / nimage

; cosmic-ray rejection: use the variance map with 4 sigma-clipping
; iterations.  do not allow the cosmic ray mask to be reset after each
; iteration (NOCLEARMASK).  scale the combined image to the total
; exposure time (WEIGHTING).  use dilation and do not sky-subtract.

    cr_reject, input_cube, -1, 0, 0, 1.0, combined_image, combined_noise, $
      combined_npix, mask_cube=inputmask, exptime=exptime, verbose=verbose, $
      nsig=nsig, dilation=dilation, dfactor=dfactor, /noclearmask, /restore_sky, $
      /noskyadjust, weighting=0, noise_cube=noise_cube, null_value=-999.0;, /xmedsky

; interpolate pixels that were rejected from all images

    null = where(combined_image eq -999.0,nnull)

    nullmask = combined_image eq -999.0
    combined_image = djs_maskinterp(combined_image,nullmask,iaxis=1L,/const)
    combined_noise = sqrt(djs_maskinterp(combined_noise^2,nullmask,iaxis=1L,/const))

; test code below:  we could constrain bad pixels to be a minimum size
; (e.g., two or more connected pixels)    
    
;   im = combined_npix ne nimage
;   imdisp, imgscl(im,/log), margin=0.01, /erase
;   out = dkboxstats(im,boxstat='total')

; COMBINED_NPIX = NIMAGE for good pixels and otherwise for flagged
; pixels

    crmask = combined_npix ne nimage ; 0 is good, 1 is bad
    splog, 'Combined image: '+string(long(total(crmask)),format='(I0)')+' bad pixels.'
    
; collapse the input mask cube into two dimensions.  note that we are
; assuming that the each bad pixel mask has either a value of 0 (zero)
; (good pixel), or a value of 2 (due to a dead pixel or column
; assigned in ICCDPROC).  if any other bad pixels (ILA_COSMIC) or
; (ICALIBRATE) have been assigned then the code below is not properly
; generalized! (jm03apr10uofa)

    newmask = fix(total(mask_cube,3)/nimage)
    
; assign the appropriate bit number to the bad pixel mask

    tempmask = crmask*fix(0)
    tempmask[where(crmask eq 1B)] = imask_bits(/crsplits)
;   tempmask[where(crmask eq 1B)] = fix(4)

    mastermask = newmask + tempmask ; 0 is good
;   mastermask = (byte(total(mask_cube eq 0B,3)) + (combined_npix ne nimage)) eq 0B

; generate the combined header using the zeroth image as a template 

    info = struct_trimtags(forage,select=['FILE','OBJECT','DATE',$
      'UT','RA','DEC','EXPTIME','PA','AIRMASS','PARANGLE']);,'SCANLEN'])
    struct_print, info
    
    mean_ra = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.ra))),2),' ',':')
    mean_dec = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(forage.dec))),2),' ',':')
    mean_airmass = djs_mean(forage.airmass)   ; mean airmass
    mean_zd = djs_mean(forage.zd)             ; zenith angle
    mean_parangle = djs_mean(forage.parangle) ; mean parallactic angle
    mean_jd = djs_mean(forage.jd)             ; mean Julian date
    mean_mjd = djs_mean(forage.mjd)           ; mean modified Julian date
    
; update the output header

    sxaddpar, header, 'EXPTIME', time, ' exposure time [seconds]', $
      format='(F12.3)'
    if (oexptime_count ne 0L) then sxaddpar, header, 'OEXPTIME', oexptime_total, $
      ' original exposure time [seconds]', format='(F12.3)'
    sxaddpar, header, 'RA', mean_ra, ' right ascension [HMS]'
    sxaddpar, header, 'DEC', mean_dec, ' declination [DMS]', after='RA'
    sxaddpar, header, 'AIRMASS', float(mean_airmass), ' airmass at mid-exposure', format='(F12.4)'
    sxaddpar, header, 'ZD', float(mean_zd), after='AIRMASS', $
      ' zenith distance at mid-exposure [degrees]', format='(F12.2)'
    sxaddpar, header, 'PARANGLE', float(mean_parangle), after='AIRMASS', $
      ' parallactic angle at mid-exposure [degrees]', format='(F12.2)'
    
; handle ISPEC header keywords
    
    for i = 1L, nimage do sxaddpar, header, 'IRAW'+string(i,format='(I2.2)'), $
      sxpar(*cube[i-1L].header,'IRAW'), ' raw FITS file', before='HISTORY'
    junk = sxpar(header,'ISKYSB',count=subcheck)
    if subcheck then for i = 1L, nimage do sxaddpar, header, 'ISKYSB'+string(i,format='(I2.2)'), $
      sxpar(*cube[i-1L].header,'ISKYSB'), ' input to ISKYSUBTRACT2D', before='HISTORY'
    for i = 1L, nimage do sxaddpar, header, 'ICOMB'+string(i,format='(I2.2)'), $
      crimlist[i-1L], ' input to ICRCOMBINE', before='HISTORY'
    for i = 1L, nimage do sxaddpar, header, 'IPIX'+string(i,format='(I2.2)'), $
      long(pixshift[i-1L]), ' pixel shift applied in ICRCOMBINE', before='HISTORY' ; NOTE LONG()!
       
    sxaddpar, header, 'INCOMB', nimage, ' number of images stacked in ICRCOMBINE', $
      before='HISTORY'

    sxaddhist, "'Cosmic ray split images combined "+im_today()+"'", header

    sxdelpar, header, 'IRAW'
    
    if keyword_set(wfits) then begin

       splog, 'Writing '+datapath+crimname+'.'
       wrt2dspec, crimname, float(combined_image), float(combined_noise), $
         mastermask, header, skyimage=float(combined_sky), datapath=datapath, $
         gzip=gzip
    
    endif
    
; clean up memory
    
    icleanup, cube

return
end
