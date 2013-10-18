;+
; NAME:
;	RKCRSPLIT
;
; PURPOSE:
;	Combine cosmic-ray split images.
;
; CALLING SEQUENCE:
;
; INPUTS:
;	crpairs - FITS file names of the CR-split images
;
; OPTIONAL INPUTS:
;	sigma      - cosmic ray threshold in sigma (default 3.0)
;	fudge      - additional threshold factor proportional to image
;                    counts to account for pointing error (default
;                    0.9, conservative)
;	datapath   - path to the data
;	gain       - detector gain [electrons/ADU; default 1.0]
;	rdnoise    - detector read noise [electron; default 0.0]
;	extra      - keywords for RKSPECREGISTER
;	
; KEYWORD PARAMETERS:
;	wfits      - write the cleaned image, sigma map, and updated
;                    bad pixel mask to disk
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;	sigmamap - new sigma map with masked pixels set to a very
;                  large number
;	header   - new header for the combined image
;	crimname - name of the output FITS file if WFITS=1
;
; COMMENTS:
;	The procedure first generates "minimum" images by finding the
;	minimum of a 3x3 box centered on each pixel in every image.
;	Next a "comparison" image is created by determining the
;	minimum value of each pixel in the "minimum" image stack.
;	Finally, for every image the "maximum deviation" is computed;
;	pixels greater than the comparison pixel value plus the
;	maximum deviation are flagged as bad.  Bad pixels are replaced
;	by good pixel counts and scaled to the total exposure time.
;	See Dolphin 2000, PASP, 12, 1383 for more details.
;
;	The object headers must have an exposure time entry, otherwise the
;       exposure time is set to one second, and an airmass entry.  The
;       gain and read ; noise are assumed to be the same for all
;       images.  All the ; images are assumed to be photometric.
;
; EXAMPLE:
;
;
; PROCEDURES USED:
;	RD2DSPEC(), SXPAR(), MAKE_CUBE(), DKBOXSTATS(), CWD(),
;	RKSPECREGISTER 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 27, U of A
;	jm01dec18uofa, generalized to an image stack and cleaned up 
;-

pro dev_rkcrsplit, crimlist, cleanimage, cleansigmap, mastermask, header, sigma=sigma, $
               fudge=fudge, datapath=datapath, gain=gain, rdnoise=rdnoise, $
               crimname=crimname, wfits=wfits

    if (n_elements(crimlist) eq 0L) then begin
       print, 'Syntax - rkcrsplit, crimlist, cleanimage, cleansigmap, mastermask, $'
       print, '           header, [sigma=, fudge=, datapath=, gain=, rdnoise=, crimname=], $'
       print, '           wfits=wfits'
       return
    endif
    
    if (size(crimlist[0],/type) ne 7L) or (strmatch(crimlist[0],'*.fits') eq 0B) then begin
       splog, 'Files must be type string FITS files.'
       return
    endif

    nimage = n_elements(crimlist)
    if nimage lt 2L then begin
       splog, 'Need at least two images to combine.'
       return
    endif
    
; set defaults
    
    if not keyword_set(datapath) then datapath = cwd()
    if not keyword_set(sigma) then sigma = 3.0
    if not keyword_set(fudge) then fudge = 0.9
    if not keyword_set(gain) then gain = 1.0
    if not keyword_set(rdnoise) then rdnoise = 0.0

    cube = rd2dspec(crimlist,datapath=datapath)

    imsize = size(cube[0].image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]

    rowaxis = findgen(nrows)
    
; image exposure time

    exptime = dblarr(nimage)+1.0
    for i = 0L, nimage-1L do if sxpar(*cube[i].header,'EXPTIME') ne 0L then $
      exptime[i] = sxpar(*cube[i].header,'EXPTIME')

    time = total(exptime,/double) ; total integration time [s]

; generate an image and sigma cube in counts per second

    imcube = make_cube(cube)
    sigmacube = make_cube(cube,/sigmamap)

    for i = 0L, nimage-1L do imcube[*,*,i] = imcube[*,*,i]/exptime[i]
    for i = 0L, nimage-1L do sigmacube[*,*,i] = sigmacube[*,*,i]/exptime[i]

; register each images using their spatial profiles

;   rkspecregister, imcube, errcube=sigmacube, maskcube=cube.mask, $
;     pixshift=pixshift, _extra=extra
    pixshift = replicate(0.0,nimage)
    
; combined image and sigma map, including bad pixels (replaced below),
; in mean counts
    
    combined_image = time * total(imcube,3) / nimage
    combined_sigmamap = time * sqrt(total(sigmacube^2.0,3)) / nimage

    splog, 'Generating the minimum and comparison images.'
    minimcube = imcube*0.0
    for j = 0L, nimage-1L do minimcube[*,*,j] = dkboxstats(reform(imcube[*,*,j]),xwidth=3,boxstat='min')

; create the comparison image of minimum pixel values    

    compim = fltarr(ncols,nrows)
    for j = 0L, nrows-1L do for k = 0L, ncols-1L do compim[k,j] = min(minimcube[k,j,*])
    
    crmask = bytarr(ncols,nrows,nimage) ; bad pixel mask 

; compute the maximum and minimum deviation for each image and compare 
    
    splog, 'Flagging and rejecting cosmic rays.'
    for i = 0L, nimage-1L do begin

       maxdev = sqrt( sigma^2.0 * ( (rdnoise/gain)^2.0 + imcube[*,*,i] * exptime[i] / gain ) + $ 
                      (fudge * imcube[*,*,i] * exptime[i])^2.0 ) / exptime[i]

       crmask[*,*,i] = (imcube[*,*,i] lt (compim + maxdev)) and (imcube[*,*,i] gt -(compim + maxdev)) ; 1B is good

    endfor

; if a pixel is bad in the whole image stack, put it in the final bad
; pixel mask, otherwise scale the good pixel value by the total
; exposure time and replace that pixel in the combined image
        
    splog, 'Scaling the counts in bad pixels.'

    image_mask = bytarr(ncols,nrows,nimage) ; new image mask (bad pixels and CRs)
    for k = 0L, nimage-1L do image_mask[*,*,k] = (cube[k].mask eq 0B) + (crmask[*,*,k] eq 0B)

    image_mask = image_mask ge 1B          ; 0B is good
    mastermask = byte(total(image_mask,3)) ; combined image mask

    for k = 0L, nimage-1L do splog, 'Image '+string(k+1,format='(I0)')+': '+$
      string(long(total(image_mask[*,*,k])),format='(I0)')+' bad pixels.'
    
    vbad = where(mastermask ge nimage,nvbad)                        ; pixels bad in all images
    bad = where((mastermask lt nimage) and (mastermask gt 0B),nbad) ; repairable pixels

; scale "bad" pixels    
    
    if nbad ne 0L then begin

       goodcounts = fltarr(nbad,nimage)
       goodsigcounts = fltarr(nbad,nimage)

       for k = 0L, nimage-1L do begin

          masksub = (reform(image_mask[*,*,k]))[bad] ; 0B is good
          good = where(masksub eq 0B,ngood)

          imsub = (reform(imcube[*,*,k]))[bad]
          sigimsub = (reform(sigmacube[*,*,k]))[bad]

          if ngood ne 0L then begin
             goodcounts[good,k] = imsub[good]
             goodsigcounts[good,k] = sigimsub[good]
          endif

       endfor

       combined_image[bad] = time * total(goodcounts,2) / (nimage - mastermask[bad])                 ; mean counts
       combined_sigmamap[bad] = time * sqrt(total(goodsigcounts^2.0,2)) / (nimage - mastermask[bad])

       mastermask[bad] = 0B ; repaired! (i.e. now good)
       
    endif

; in the combined image interpolate over "very bad" pixels (across the x-axis)

    cleanimage = djs_maskinterp(combined_image,mastermask,iaxis=0L)
    cleansigmamap = djs_maskinterp(combined_sigmamap,mastermask,iaxis=0L)

    splog, 'Combined image: '+string(long(total(mastermask)),format='(I0)')+' bad pixels.'
    
    mastermask = mastermask eq 0B ; 1B is good

; additional rejection of hot pixels:  (1) the sky-subtracted pixel
; value is more than 10 times the average of adjacent sky-subtracted
; pixel values; (2) the sky-subtracted pixel value is more than 7
; standard deviations above the average value of adjacent pixels

;   sky, crrejim, smode, ssig

; update the header using the zeroth image as a template

    header = *cube[0].header
    
;   ut0 = im_hms2dec([sxpar(*cube[0].header,'UTMIDDLE'),sxpar(*cube[0].header,'UT')])
;   ut1 = im_hms2dec([sxpar(*cube[1].header,'UTMIDDLE'),sxpar(*cube[1].header,'UT')])
;   utmiddle = im_dec2hms((2.0*(ut0[0]-ut0[1])+2.0*(ut1[0]-ut1[1]))/2.0)

; possibly want to check for the airmass header entry
    
    airmass_vector = fltarr(nimage)
    for i = 0L, nimage-1L do airmass_vector[i] = sxpar(*cube[i].header,'AIRMASS')
    mean_airmass = djs_mean(airmass_vector)

; compute the mean parallactic angle of the observation

    parangle_vector = fltarr(nimage)
    for i = 0L, nimage-1L do parangle_vector[i] = sxpar(*cube[i].header,'PARANGLE')
    mean_parangle = djs_mean(parangle_vector)
    
    sxaddpar, header, 'EXPTIME', time, format='(F12.3)'
    sxaddpar, header, 'AIRMASS', mean_airmass, format='(F6.4)'
    sxaddpar, header, 'PARANGLE', mean_parangle, format='(F12.3)'
    sxaddhist, "'CR split images ...", header
    for j = 0L, nimage-1L do sxaddhist, '    '+crimlist[j], header
    sxaddhist, "...combined "+im_today()+"'", header
    sxaddhist, "'Pixel shifts = "+strjoin(strcompress(pixshift,/remove),', ')+"'", header

    if keyword_set(wfits) then begin

; use the zeroth image to form the output name
       
       crimname = strmid(crimlist[0],0,strpos(crimlist[0],'.fits'))+'_'+$
         string(nimage,format='(I0)')+'.fits'
;      crimname = 'test.fits'
       
       splog, 'Writing '+datapath+crimname+'.'
       wrt2dspec, crimname, float(cleanimage), float(cleansigmamap), byte(mastermask), $
         header, datapath=datapath
    
    endif
    
; clean up memory
    
    rkcleanup, cube

return
end
