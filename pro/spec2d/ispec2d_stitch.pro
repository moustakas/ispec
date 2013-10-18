;+
; NAME:
;       ISPEC2D_STITCH
;
; PURPOSE:
;       Stitch together 
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;       Only works on two images.  Shift the second image relative to
;       the first. 
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jun 27, U of A - written, based on
;                                           SINGS_STITCH 
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

pro ispec2d_stitch, spec1list, spec2list, stitchlist, objectlist, $
  nlags=nlags, user_pixshift=user_pixshift, datapath=datapath, $
  debug=debug, wfits=wfits, gzip=gzip

    ngalaxy = n_elements(spec1list)
    if (n_elements(spec2list) ne ngalaxy) or (n_elements(stitchlist) ne ngalaxy) then begin
       print, 'Syntax - ispec2d_stitch, spec1list, spec2list, stitchlist, /wfits'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(nlags) eq 0L) then nlags = 51L
    if (n_elements(user_pixshift) eq 0L) then user = 0L else user = 1L

; call this routine recursively    
    
    if (ngalaxy gt 1L) then begin
       for igalaxy = 0L, ngalaxy-1L do $
         ispec2d_stitch, spec1list[igalaxy], spec2list[igalaxy], $
           stitchlist[igalaxy], objectlist[igalaxy], nlags=nlags, $
           user_pixshift=user_pixshift, datapath=datapath, debug=debug, $
           wfits=wfits, gzip=gzip
       return
    endif

; forage header information and read the 2D spectra
    
    info = iforage(datapath+[spec1list,spec2list])
    spec1 = rd2dspec(spec1list,datapath=datapath,wset=wset)
    spec2 = rd2dspec(spec2list,datapath=datapath)

    ncols = spec1.naxis1
    nrows = spec1.naxis2
    midcol = ncols / 2L
    midrow = nrows / 2L
    nmed = 10L

    lags = lindgen(2*nrows)-nrows
    bigprofile1 = lags*0.0D
    bigprofile2 = lags*0.0
    
    profile1 = total(spec1.image[midcol-nmed/2L:midcol+nmed/2L,*],1)/float(nmed)
    profile1 = profile1/max(profile1)
    profile2 = total(spec2.image[midcol-nmed/2L:midcol+nmed/2L,*],1)/float(nmed)
    profile2 = profile2/max(profile2)

    bigprofile1[0L:nrows-1L] = profile1
    bigprofile2[0L:nrows-1L] = profile2
;   bigprofile2[nrows:2*nrows-1L] = profile2
    
    result = djs_correlate(bigprofile1,bigprofile2,lags)
    
    
    refspec = im_normalize(total(spec2.image,1),rowaxis,/max,normwave=nrows-15L,binsize=14)
    bigrefspec[0:nrows-1L] = refspec

    spec = im_normalize(total(spec1.image,1),rowaxis,/max,normwave=15L,binsize=14,const=normconst)
    bigspec[nrows:bignrows-1L] = spec
    bigivar[nrows:bignrows-1L] = 1.0/total((spec1.sigmamap/normconst)^2,1)

    lags = -reverse(lindgen(nlags))
    result = djs_correlate(bigspec,bigrefspec,lags)
    maxcorr = max(result,bestlag)
    pixshift = lags[bestlag]
       
    plot, bigrefspec, ps=10, xsty=3, ysty=3
    djs_oplot, bigspec, ps=10, color='red'
    djs_oplot, shift(bigspec,pixshift), ps=10, color='green'

    plot, bigrefspec, ps=10, xsty=3, ysty=3
    djs_oplot, bigspec, ps=10, color='red'
    djs_oplot, shift(bigspec,pixshift), ps=10, color='green'

stop    
    
; compute the mean upper, lower, and central profiles; plot the mean
; ratio at each of these three positions

    profile_cen = total(input_cube[midcol-nmed/2L:midcol+nmed/2L,*,ii],1)/float(nmed)
    profile_lo = total(input_cube[0L:nmed-1L,*,ii],1)/float(nmed)
    profile_up = total(input_cube[ncols-nmed:ncols-1L,*,ii],1)/float(nmed)

    refprofile_cen = total(input_cube[midcol-nmed/2L:midcol+nmed/2L,*,0],1)/float(nmed)
    refprofile_lo = total(input_cube[0L:nmed-1L,*,0],1)/float(nmed)
    refprofile_up = total(input_cube[ncols-nmed:ncols-1L,*,0],1)/float(nmed)

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
       
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)

    endif

stop
    
; HERE!       
    
; cross-correlate the two spatial profiles to determine the pixel
; shift; shift SPEC2 to match SPEC1

    ncols = spec1.naxis1
    nrows = spec1.naxis2
    rowaxis = lindgen(nrows)
    bignrows = 2*nrows
    bigindx = lindgen(bignrows)
    
    bigrefspec = bigindx*0.0
    bigspec = bigindx*0.0

    refspec = im_normalize(total(spec2.image,1),rowaxis,/max,normwave=nrows-15L,binsize=14)
    bigrefspec[0:nrows-1L] = refspec

    spec = im_normalize(total(spec1.image,1),rowaxis,/max,normwave=15L,binsize=14)
    bigspec[nrows:bignrows-1L] = spec

    lags = -reverse(lindgen(nlags))
;   nlags = 101L
;   lags = lindgen(nlags)-nlags/2L
    result = djs_correlate(bigspec,bigrefspec,lags)
    maxcorr = max(result,bestlag)
    pixshift = lags[bestlag]
       
; form the output arrays and co-add the spectra; construct the average
; image and sky spectrum; only divide by two in the overlapping rows

    outnrows = bignrows + pixshift
    
    outimage1 = dblarr(ncols,outnrows) & outimage2 = outimage1*0.0
    outimage1[*,0L:nrows-1L] = spec2.image
    outimage2[*,nrows+pixshift:outnrows-1L] = spec1.image

    outskyimage1 = dblarr(ncols,outnrows) & outskyimage2 = outskyimage1*0.0
    outskyimage1[*,0L:nrows-1L] = spec2.sky
    outskyimage2[*,nrows+pixshift:outnrows-1L] = spec1.sky

    denom = outimage1*0.0+1.0
    denom[where((abs(outimage1) gt 0.0) and (abs(outimage2) gt 0.0))] = 2.0
    
    outimage = float(total([ [ [outimage1] ], [ [outimage2] ] ],3)/denom)
    outskyimage = float(total([ [ [outskyimage1] ], [ [outskyimage2] ] ],3)/denom)

; error map
    
    outsigmap1 = dblarr(ncols,outnrows) & outsigmap2 = outsigmap1*0.0
    outsigmap1[*,0L:nrows-1L] = spec2.sigmamap
    outsigmap2[*,nrows+pixshift:outnrows-1L] = spec1.sigmamap
    outsigmap = float(sqrt(total([ [ [outsigmap1^2.0] ], [ [outsigmap2^2.0] ] ],3)/denom))

; bad pixel mask
    
    outmask1 = fltarr(ncols,outnrows) & outmask2 = outmask1*0.0
    outmask1[*,0L:nrows-1L] = spec2.mask
    outmask2[*,nrows+pixshift:outnrows-1L] = spec1.mask
    outmask = fix(total([ [ [outmask1] ], [ [outmask2] ] ],3))

; telluric spectrum    

    telluric_spec = spec1.telluric_spec
    telluric_head = spec1.telluric_head

; update the header; assume that EXPTIME, DATE-OBS, SCANLEN, PA, CCD
; noise parameters, wavelength parameters, RMAXERR, and RMEANERR, are
; the same; retain the time parameters of the zeroth exposure (UT,
; UTMIDDLE, JD, and ST) and all the ISPEC2D keywords (IRAW, etc.);
; compute the mean AIRMASS, ZD, PARANGLE, RA, and DEC; compute the new
; reference pixel CRPIX2; add a history note
    
    outheader = *spec1.header

    meanra = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(info.ra))),2),' ',':')
    meande = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(info.dec))),2),' ',':')
    
    sxaddpar, outheader, 'OBJECT', objectlist
    sxaddpar, outheader, 'RA', meanra, ' right ascension'
    sxaddpar, outheader, 'DEC', meande, ' declination'
    sxaddpar, outheader, 'ZD', float(djs_mean(info.zd)), ' mean zenith distance [degrees]'
    sxaddpar, outheader, 'AIRMASS', float(djs_mean(info.airmass)), ' effective airmass'
    sxaddpar, outheader, 'PARANGLE', float(djs_mean(info.parangle)), ' parallactic angle ([0-180] degrees)'
    sxaddpar, outheader, 'CRPIX2', outnrows/2L+1L, ' reference pixel number'
    sxaddhist, "'Images "+strjoin(repstr(info.file,'.fits.gz',''),' and ')+" stitched together "+im_today()+"'", outheader
    
; write out

    if keyword_set(wfits) then begin

       splog, 'Writing '+datapath+stitchlist+'.'

       wrt2dspec, stitchlist, outimage, outsigmap, outmask, outheader, $
         skyimage=outskyimage, wset=wset, telluric_spec=telluric_spec, $
         telluric_head=telluric_head, datapath=datapath, gzip=gzip

    endif

    icleanup, spec1
    icleanup, spec2

; debugging plots    
    
    if not keyword_set(wfits) then begin

       im_window, 0, /square

       plot, lags, result, xsty=3, ysty=3, xthick=2, ythick=2, charsize=2.0, $
         charthick=2.0, xtitle='Pixel Shift', ytitle='Cross-Correlation'
       oplot, pixshift*[1,1], !y.crange, line=2, thick=2

       im_window, 2, /square
       
;      plot, bigindx, bigrefspec, ps=10, xsty=3, ysty=3, yrange=[0,1.1], $
;        xrange=[50,200], xthick=2, ythick=2, charsize=2.0, $
;        charthick=2.0, xtitle='Pixel Shift', ytitle='Relative Intensity'
;      djs_oplot, bigindx, bigspec, ps=10, color='red'
;      djs_oplot, bigindx, shift(bigspec,pixshift), ps=10, color='yellow'

       scale = 1E15
       djs_plot, scale*outimage[nrows/2L,*], ps=10, xsty=3, ysty=3, thick=2.0, $
         xthick=2, ythick=2, charsize=2.0, charthick=2.0, xtitle='Row', $
         ytitle='Flux at Center Row [10^{-15} '+flam_units()+']'
       djs_oplot, scale*outimage1[nrows/2L,*], ps=10, color='red'
       djs_oplot, scale*outimage2[nrows/2L,*], ps=10, color='yellow'

       legend, ['Combined','Spectrum 1','Spectrum 2'], /left, /top, box=0, $
         charthick=2.0, charsize=1.8, textcolor=djs_icolor(['white','red','yellow'])

       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
       
    endif

return
end    
