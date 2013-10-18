;+
; NAME:
;	RKSPECREGISTER
;
; PURPOSE:
;	Register (by partial-pixels) two or more two-dimensional
;	spectral images using their spatial profiles.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;	imcube   - input image cube [NX,NY,NIMAGE]
;
; OPTIONAL INPUTS:
;	errcube  - error cube for IMCUBE
;	maskcube - bad pixel mask for IMCUBE (1B is good)
;	nsamp    - resampling factor [default 10]
;	maxshift - maximum allowable pixel shift [default 10]
;	method   - long integer specifying the technique to use for
;                  registering the images by partial pixels:  (NOT
;                  IMPLEMENTED)
;	  0:  linear interpolation
;	  1:  spline interpolation (default)
;	  2:  damped sinc function
;	
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;	pixshift - vector of pixel shifts for each image relative to
;                  the zeroth image in IMCUBE
;
; COMMENTS:
;	The images are only registered in the spatial dimension,
;	assuming the second dimension of IMCUBE is the spatial one.
;	If MASKCUBE is passed then any pixels effected by the
;	registering process will be flagged as bad.  All the data
;	cubes are returned modified.
;
; EXAMPLE:
;
;
; PROCEDURES USED:
;	DJS_MEDIAN()
;
; INTERNAL ROUTINES:
;	
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Feb 28, U of A, excised code from RKCRSPLIT 
;	jm02jun7uofa - added padding with zeros rather than
;                      extrapolation
;
;-

function imageshift, image, nsamp, pixshift, maxshift
; shift an image by partial pixels

    imsize = size(image,/dimension)
    nx = imsize[0]
    ny = imsize[1]

    saxis = findgen(ny)           ; spatial axis
    finesaxis = findgen(ny*nsamp) ; sampled spatial axis

stop
    
; resample the image in the spatial dimension and pad with zeros

;   fineimage = frebin(image,nx,ny*nsamp,/total)
    fineimage = congrid(image,nx,ny*nsamp,cubic=-0.5,/minus_one)

    shiftimage = fltarr(nx,(ny+2*pixshift)*nsamp)
    shiftimage[*,pixshift*nsamp:(ny+pixshift)*nsamp-1] = fineimage

; shift and crop
    
;   outimage = shift(fineimage,0,pixshift*nsamp) 
    outimage = shift(shiftimage,0,pixshift*nsamp)
;   outimage = sshift2d(fineimage,[0,pixshift*nsamp]) ; shift using a damped sinc function

; since SHIFT wraps around the edges, replace the edges with spline
; interpolated points

    arows = long(maxshift*nsamp)
    if pixshift lt 0.0 then srows = lindgen(arows)+long(ny*nsamp)-arows else srows = lindgen(arows)
    
;   arows = long(pixshift*nsamp) ; number of affected rows
;   if arows lt 0L then srows = [lindgen(-arows)+long(ny*nsamp)+arows,erows] else srows = lindgen(arows)  

    for k = 0L, nx-1L do outimage[k,srows] = interpol(fineimage[k,srows],srows,srows+abs(pixshift*nsamp),/spline) ; extend
;   for k = 0L, nx-1L do outimage[k,*] = interpol(fineimage[k,*],finesaxis,finesaxis+pixshift*nsamp,spline=spline)

;   djs_plot, fineimage[550,*], xsty=3, ysty=3, ps=10, xr=[1240,1265]
;   djs_oplot, outimage[550,*], ps=10, thick=2.0

; now rebin back to the original size

    outimage = congrid(outimage,nx,ny,cubic=-0.5)

return, outimage
end
    
pro dev_rkspecregister, imcube, errcube=errcube, maskcube=maskcube, nsamp=nsamp, $
  maxshift=maxshift, method=method, pixshift=pixshift, doplot=doplot

    if n_elements(imcube) eq 0L then begin
       print, 'Syntax - rkspecregister, '
       return
    endif

    if not keyword_set(nsamp) then nsamp = 10.0
    if not keyword_set(maxshift) then maxshift = 10.0
    if n_elements(method) eq 0L then method = 1L else method = long(method[0])
   
; error checking
    
    ndim = size(imcube,/n_dimension)
    if ndim ne 3L then message, 'IMCUBE must be three-dimensional!'

    imsize = size(imcube,/dimension)
    nx = imsize[0]
    ny = imsize[1]
    nimage = imsize[2]

    nerr = n_elements(errcube)
    if nerr ne 0L then begin

       ndimerr = size(errcube,/n_dimension)
       if ndimerr ne 3L then message, 'ERRCUBE must be three-dimensional!'
       errsize = size(errcube,/dimension)
       if (errsize[0] ne nx) and (errsize[1] ne ny) and (errsize[2] ne nimage) then $
         message, 'Dimensions of IMCUBE and ERRCUBE do not agree!'

    endif
       
    nmask = n_elements(maskcube)
    if nmask ne 0L then begin

       ndimmask = size(maskcube,/n_dimension)
       if ndimmask ne 3L then message, 'MASKCUBE must be three-dimensional!'
       masksize = size(maskcube,/dimension)
       if (masksize[0] ne nx) and (masksize[1] ne ny) and (masksize[2] ne nimage) then $
         message, 'Dimensions of IMCUBE and MASKCUBE do not agree!'
       mask = maskcube

    endif

; generate the spatial profiles by collapsing along the NX dimension
; then resample by a factor of NSAMP.  do not conserve flux

    sprofile = total(imcube,1)
    sprofine = congrid(sprofile,ny*nsamp,nimage,cubic=-0.5)

; construct a chi2 statistic for the square deviations of each
; spatial profile relative to the zeroth spatial profile.  to speed
; the computation, constrain the maximum pixel shift to be MAXSHIFT
    
    npix = maxshift*nsamp
    lags = lindgen(2*npix)-npix
    nlags = n_elements(lags)
    chi2array = fltarr(nlags,nimage-1L)

    padvals = replicate(djs_median(sprofine[*,0]),npix)             ; pad with the median
    refprofile = [padvals,sprofine[*,0],padvals]/max(sprofine[*,0]) ; reference profile (zeroth image)

    pixshift = fltarr(nimage)   ; solve for the output pixel shifts

    testcube = imcube
    
    for i = 1L, nimage-1L do begin

       padvals = replicate(djs_median(sprofine[*,i]),npix)
       profile = [padvals,sprofine[*,i],padvals]/max(sprofine[*,i])

       for j = 0L, nlags-1L do chi2array[j,i-1] = total((refprofile-shift(profile,lags[j]))^2.0)
       minchi2 = min(chi2array[*,i-1],lagbest) ; minimum chi2

       if keyword_set(doplot) then begin
       
          djs_plot, lags, chi2array, xsty=3, ysty=3       
          djs_plot, refprofile, ysty=3, ps=10, xsty=3
          djs_oplot, profile, ps=10, color='red'
          djs_oplot, shift(profile,lags[lagbest]), ps=10, color='green'

       endif
          
       pixshift[i] = lags[lagbest]/nsamp
       print, 'Image '+strn(i)+': PIXSHIFT [pixels] = '+strn(pixshift[i])

       if pixshift[i] ge maxshift then message, 'WARNING: PIXSHIFT exceeds MAXSHIFT.', /info

       intshift = fix(pixshift[i])         ; nearest integer shift
       w1 = abs(pixshift[i])-abs(intshift) ; weight 1
       w2 = 1.0-w1                         ; weight 2

;      for j = 0L, nx-1L do for k = intshift+1L, ny-intshift-2L do $
;        test[j,k,i] = w1*test[j,k-1L,i]+w2*test[j,k,i]

;      plot, imcube[610,*,1]-test[610,*,1], xsty=3, ysty=3, ps=10
;      djs_oplot, test[610,*,1], color='red'

       if abs(pixshift[i]) gt 0.0 then begin ; non-zero pixel shift
          
          for j = 0L, nx-1L do for k = 1L, ny-1L do $
            imcube[j,k,i] = w1*imcube[j,k-1L,i]+w2*imcube[j,k,i]

stop
          
;         for j = 0L, nx-1L do for k = intshift+1L, ny-1L do $
;           imcube[j,k,i] = w1*imcube[j,k-1L,i]+w2*imcube[j,k,i]

;          if pixshift[i] lt 0.0 then begin
;             for j = 0L, nx-1L do for k = intshift+1L, ny-intshift-1L do $
;               imcube[j,k-intshift-1L,i] = w1*imcube[j,k-1L,i]+w2*imcube[j,k,i]
;          endif else begin
;             for j = 0L, nx-1L do for k = intshift+1L, ny-intshift-1L do $
;               imcube[j,k-intshift-1L,i] = w1*imcube[j,k-1L,i]+w2*imcube[j,k,i]
;          endelse

       endif
       
;      if pixshift[i] lt 0.0 then subshift = -(abs(pixshift[i])-abs(intshift)) else subshift = pixshift[i]-intshift
       
;       if pixshift[i] ne 0.0 then begin ; register the images
;          
;          imcube[*,*,i] = imageshift(imcube[*,*,i],nsamp,pixshift[i],maxshift)
;          if nerr ne 0L then errcube[*,*,i] = imageshift(errcube[*,*,i],nsamp,pixshift[i],maxshift)
;          if nmask ne 0L then maskcube[*,*,i] = imageshift(mask[*,*,i] eq 0B,nsamp,pixshift[i],maxshift) eq 0B
;          
;       endif
       
    endfor    

;   imcube = imcube[*,0L:ny-1L-(intshift+1L),*]
;   imcube = imcube[*,intshift+1L:ny-1L,*]
;   imcube = imcube[*,1L:ny-1L,*]
    imcube = imcube[*,0L:ny-2L,*]

    plot, testcube[510,*,0], ps=10, xsty=3, ysty=3
    djs_oplot, testcube[510,*,1], ps=10, color='green'

    plot, imcube[510,*,0], ps=10, xsty=3, ysty=3
    djs_oplot, imcube[510,*,1], ps=10, color='green'

stop
    
; check that it worked

    rkspecregister, imcube, /doplot
    
stop

    
return
end    
