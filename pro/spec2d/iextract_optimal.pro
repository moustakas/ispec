;+
; NAME:  IEXTRACT_OPTIMAL_SPROF
;       
; PURPOSE:  Construct a 2D spatial profile of a CCD image.
;
; CALLING SEQUENCE:
;       result = iextract_opt_sprof(image, spmask, opt_ncol,
;                                   [/opt_check, /opt_gauss]) 
;
; INPUTS:
;	image       - two-dimensional wavelength calibrated,
;                     sky-subtracted image frame
;	spmask      - two-dimensional bad pixel mask (use=0, ignore=1) 
;       opt_ncol    - number of columns to average when fitting in
;                     spatial direction, see discussion in comments
;                     section 
;
; KEYWORDS:
;       opt_check   - display spatial profile fitting steps
;       opt_gauss   - replace each column profile with fitted gaussian
;                     (not advisable if the spatial profile is not
;                     gaussian!) 
;
; OUTPUTS: returns the 2d spatial profile image
;
; COMMENTS: A 2D image is created which models the spatial
;           illumination pattern from a point source.  This is done by
;           calculating the mean column profile at specified intervals
;           utilizing sigma-clipping and ignoring masked pixels (see
;           OPT_IGNORE).  In the unlikely event that none of the
;           columns in a given interval have usable entries for a
;           certain row (eg, a bad pixel row on the detector), an
;           interpolated value from a high-order spline fit which
;           preserves all other values is used. 
;
;           The number of columns averaged is determined by OPT_NCOL.
;           OPT_NCOL should be small enough that any changes in the
;           illumination pattern (eg, travelling peak center or
;           widening of the width) vary approximately linearly across
;           the columns of interest.  However, OPT_NCOL must be
;           sufficiently large that the resulting mean profile is free
;           of noise.  [Note: The reliability of the mean profile can
;           be checked visually by setting the OPT_CHECK keyword.
;           Each profile is then plotted with error bars corresponding
;           to the 1-sigma uncertainties in the determination of each
;           mean value.  If these error bars are uncomfortably large,
;           OPT_NCOL should be increased.)
;
;           The mean profiles are normalized to the amplitude of a
;           gaussian fit in order to remove spectral shape in the 
;           illumination pattern.  They are then fit along each row
;           with a high-order spline that preserves all values without
;           smoothing.  This is done in order to evaluate the
;           illumination pattern at each column (which is then
;           replaced with a gaussian fit if the OPT_GAUSS keyword if
;           set).  Any negative values in the resulting 2D image are
;           redefined to be zero, and each column is normalized to sum
;           to one. 
;          
;           The resulting 2D spatial profile will likely look "choppy"
;           (ie, undulating peak) due to the normally sparse sampling
;           in the spatial dimension.
;
; PROCEDURES USED: MPFITPEAK(), DJS_PLOT, DJS_OPLOT,
;                  LEGEND, BSPLINE_ITERFIT(), BSPLINE_VALU()
;
; MODIFICATION HISTORY:
;       written by A.R.Marble, Steward Observatory, June 2003
;       streamlined by A.R.Marble, 2004 March 4
;       OPT_GAUSS keyword added by A.R.Marble, 2004 May 21
;-

function iextract_optimal_sprof, image, spmask, opt_ncol, opt_skip=opt_skip, $
                                 opt_check=opt_check, opt_gauss=opt_gauss

    if n_elements(opt_skip) eq 0L then opt_skip = 10 else opt_skip = fix(opt_skip)

    ncols = n_elements(image[*,0]) ; number of columns
    nrows = n_elements(image[0,*]) ; number of rows
    nbins = ncols / opt_ncol       ; number of column bins

    if n_elements(spmask[*,0]) ne ncols or n_elements(spmask[0,*]) ne nrows then $
      message, 'IMAGE and MASK have incompatible dimensions' $
    else spmask = byte(spmask)

    splog, 'constructing spatial profile...'    

    sprof   = image * 0.0
    nimage  = image * 0.0
    binsfit = fltarr(nbins, nrows)

    binsaxis = lindgen(nbins) * opt_ncol + opt_ncol / 2L
    rowaxis  = lindgen(nrows)
    colaxis  = lindgen(ncols)

    if keyword_set(opt_check) then window, 0, xs=500, ys=500

    binmean = fltarr(nrows)
    binmeanerr = fltarr(nrows)

    for i = 0L, ncols-1L do nimage[i,*] = image[i,*] / total(image[i,*])

    for i = 0L, nbins-1L do begin

       bincols = lindgen(opt_ncol) + i * opt_ncol ; column indices

       okay = lonarr(nrows)

; form mean profile of binned columns, disregarding masked pixels

       for j = 0L, nrows-1L do begin
          
          wh = where(spmask[bincols, j] eq 0L and finite(nimage[bincols,j]) eq 1L, count)
          if count gt 1L then begin
             djs_iterstat, (nimage[bincols,j])[wh], mask=djsmask, mean=djsmean, sigrej=3
             if total(djsmask) gt 1L then begin
                use = (bincols[wh])[where(djsmask eq 1)]
                binmean[j] = djsmean
                stdev = sqrt((1d0/(total(djsmask)-1))*total((nimage[use,j]-djsmean)^2))
                binmeanerr[j] = stdev / sqrt(total(djsmask)) ; error in the mean
                okay[j] = 1L
             endif
          endif 

       endfor

       if total(okay) eq 0 then message, 'no valid pixels!'

; normalize profile to gaussian amplitude

       gaussfit = mpfitpeak(rowaxis[where(okay)], binmean[where(okay)], $
                          a, nterms=3, /positive)
       binmean = binmean / a[0]
       binmeanerr = binmeanerr / a[0]

       for k = 0L, opt_ncol-1L do nimage[bincols[k],*] = nimage[bincols[k],*] / a[0]

; evaluate missing values with high-order spline

       if total(okay) eq nrows then binsfit[i,*] = binmean else begin
          
          sset = bspline_iterfit(rowaxis[where(okay)], binmean[where(okay)], $
                                nbkpts=nrows, nord=3, /silent)
          binsfit[i,*] = bspline_valu(rowaxis, sset)
          
       endelse

       binsfit = binsfit > 0    ; force positive values

       if keyword_set(opt_check) then begin

          djs_plot, rowaxis, binsfit[i,*], xsty=3, ysty=3, ps=4, charthick=2.0, $
            charsize=1.3, xthick=2.0, ythick=2.0, xtitle='Row Number', $
            ytitle='Normalized Counts', thick=2.0, color='red'
          if total(okay) gt 0L then begin
             oplot, rowaxis[where(okay)], binsfit[i,where(okay)], psym=4, thick=2.0
             errplot, rowaxis[where(okay)], (binsfit[i,*]-binmeanerr)[where(okay)], $
               (binsfit[i,*]+binmeanerr)[where(okay)]
          endif
          text = 'Mean column '+strn(i,format='(I0)')+'/'+strn(nbins,format='(I0)')
          if total(okay) lt nrows then text = [text, 'interpolated values in red']
          legend, text, /left, /top, box=0, charsize=1.3, charthick=2.0
          cc = get_kbrd(1)

       endif
       
    endfor

; fit rows with high order spline in order to "fill-in" columns

    for i = 0L, nrows-1L do begin

       sset = bspline_iterfit(binsaxis, reform(binsfit[*,i]), nbkpts=nbins-1, $
                              nord=3, /silent)
       sprof[*,i] = bspline_valu(colaxis, sset) > 0 ; force positive values

       if keyword_set(opt_check) then begin

          djs_plot, binsaxis, binsfit[*,i], xsty=3, ysty=3, charthick=2.0, $
            charsize=1.3, xthick=2.0, ythick=2.0, xtitle='Column Number', $
            ytitle='Normalized Counts', thick=2.0, ps=4, yrange=minmax(binsfit)
          djs_oplot, colaxis, sprof[*,i], color='green', thick=2.0
          legend, 'Row '+strn(i,format='(I0)')+'/'+strn(nrows,format='(I0)'), $
            /left, /top, box=0, charsize=1.3, charthick=2.0
          cc = get_kbrd(1)

       endif
       
    endfor

; fit gaussian to each column (if desired)

    if keyword_set(opt_gauss) then begin

;;        len = strlen(string(ncols, f='(i0)'))
;;        fmt = 'i'+string(len,f='(i0)')+'.'+string(len,f='(i0)')

       model = mpfitpeak(rowaxis, djs_median(sprof[ncols/2.-50:ncols/2.+50,*],1), $
                         modelpars, nterms=3, /positive)

       parinfo = replicate({value: 0.0D, fixed: 0B},3)
       parinfo.value = modelpars
       parinfo.fixed = [0,0,1]

       for i=0,ncols-1 do begin
          
;;           print, format='("--> fitting gaussian curve to column ",'+fmt+'," / ",'+ $
;;             fmt+',a1,$)', i, ncols, string(13b)

          gaussfit = mpfitpeak(rowaxis, sprof[i,*], a, nterms=3, /positive, parinfo=parinfo)
          curve = gaussian(findgen(nrows*10.)/10., a )

          if keyword_set(opt_check) and fix(i)/opt_skip eq float(i)/opt_skip then begin
          
             djs_plot, rowaxis, sprof[i,*], xsty=3, ysty=3, charthick=2.0, $
               charsize=1.3, xthick=2.0, ythick=2.0, xtitle='Row Number', $
               ytitle='Normalized Counts', thick=2.0, ps=4, yrange=[0,1]
             djs_oplot, rowaxis, nimage[i,*], thick=2.0, ps=4, color='red'
             djs_oplot, findgen(nrows*10.)/10., curve, color='green', thick=2.0
             legend, ['Column '+strn(i,format='(I0)')+'/'+strn(ncols,format='(I0)'), $
                      'red --> raw image'], /left, /top, box=0, charsize=1.3, charthick=2.0
             cc = get_kbrd(1)

          endif
       
          sprof[i,*] = gaussfit
          
       endfor
       
      
    endif

; normalize column totals to one

    for i = 0L, ncols-1L do sprof[i,*] = sprof[i,*] / total(sprof[i,*])

    if keyword_set(opt_check) then begin

       !p.multi=[2,1,2]
       erase
       shade_surf, image, zr=minmax(image), zst=3, title='image', charsize=1.3
       shade_surf, sprof, zr=minmax(sprof), zst=3, title='spatial illumination pattern', $
         charsize=1.3
       cc = get_kbrd(1)
       !p.multi=[1,1,1]
       erase

    endif
    
    return, sprof
    
 end
 
;+
; NAME:  IEXTRACT_OPTIMAL_1D
;       
; PURPOSE:  Optimally extract a one-dimensional spectrum 
;           (weight by variance).
;
; CALLING SEQUENCE:
;       result = iextract_optimal_1D(image, sprof, mask, var, $ 
;                                    opt_nsig, [sig1d=, mask1d=] 
;
; INPUTS:
;	image       - two-dimensional wavelength calibrated,
;                     sky-subtracted image frame
;       sprof       - two-dimensional spatial profile image; each
;                     column should sum to one
;	mask        - two-dimensional bad pixel mask
;       var         - two-dimensional variance map
;       sky         - two-dimensional sky image
;       opt_nsig    - multiplicative factor for cosmic-ray sigma 
;                     clipping (default = 4)
;
; OUTPUTS: returns one-dimensional optimally extracted spectrum
;
; OPTIONAL OUTPUTS:
;       sig1d       - corresponding one-dimensional error array
;       mask1d      - composite mask of combined pixels
;       sky1d       - one-dimensional sky array 
;
; PROCEDURES USED: 
;
; MODIFICATION HISTORY:
;       written by A.R.Marble, Steward Observatory, June 2003
;       streamlined by A.R.Marble, 2004 March 4
;-

function iextract_optimal_1d, image, sprof, mask, var, sky, $
                 opt_nsig, sig1d=sig1d, mask1d=mask1d, sky1d=sky1d

    splog,'optimally extracting 1D spectrum...'

    ncols = n_elements(image[*,0])
    nrows = n_elements(image[0,*])

    if n_elements(mask[*,0]) ne ncols or n_elements(mask[0,*]) ne nrows then $
      message, 'IMAGE and MASK have incompatible dimensions' $
    else mask = byte(mask)

    spec1d = fltarr(ncols)
    sig1d  = fltarr(ncols)
    sky1d  = fltarr(ncols)
    mask1d = fltarr(ncols)

    ttlflx = total(image, 2) ; includes any residual cosmic rays

; - form cosmic-ray mask, include bad pixels

    for i = 0L, ncols-1L do begin          
       
       old = -1L       ; number of bad pixels from previous iteration
       
; - identify most deviant pixel per iteration until no new ones found

       while old lt total(mask[i,*]) do begin

          old = total(mask[i,*])

; - calculate deviations of pixels from spatial profile expectation

          deviations = (1-mask[i,*]) * (image[i,*] - ttlflx[i] * sprof[i,*])^2. / var[i,*] 

; - exclude only most deviant pixel per iteration

          if max(deviations) gt opt_nsig^2 then $
            mask[i, (where(deviations eq max(deviations)))[0]] = 1L


; - optimally extract flux

          good = where(mask[i,*] eq 0, count)
          if count gt 0L then spec1d[i] = $
            total(sprof[i,good] * image[i,good] / var[i,good]) / $
            total(sprof[i,good]^2. / var[i,good]) else spec1d[i] = 0
          
          ttlflx[i] = spec1d[i]
            
       endwhile   

; - "optimally" extract sky

       if count gt 0L then sky1d[i] = $
         total(sprof[i,good] * sky[i,good] / var[i,good]) / $
         total(sprof[i,good]^2. / var[i,good]) else sky1d[i] = 0
       
; - variance of optimal extraction

       if count gt 0L then sig1d[i] = $
         sqrt(total(sprof[i,good]) / $
              total(sprof[i,good]^2. / var[i,good]))

; - collapse mask image (require at least 10% of uncorrupted flux)

       if total(sprof[i,*] * (1-mask[i,*])) lt 0.1 then mask1d[i] = fix(1)
       
    endfor

    return, spec1d

 end

;+
; NAME:  IEXTRACT_OPTIMAL
;       
; PURPOSE:  Optimally extract a one-dimensional spectrum 
;           (weight by variance).
;
; CALLING SEQUENCE:
;       iextract_optimal, imnosky, sigmamap, mask, apinfo, $
;       [skyimage=skyimage, opt_ncol=, opt_ignore=, opt_nsig=], $
;       [spec1d=, sig1d=, mask1d=, sprof=, /opt_check]
;
; INPUTS:
;	imnosky     - two-dimensional wavelength calibrated,
;                     sky-subtracted image frame
;	sigmamap    - corresponding two-dimensional error map
;	mask        - corresponding two-dimensional bad pixel mask
;	apinfo      - aperture extraction parameters
;
; OPTIONAL INPUTS:
;       skyimage    - two-dimensional sky image
;       
; OPTIONAL INPUTS:
;       opt_ncol    - number of columns to average when fitting in
;                     spatial direction (default=41); see discussion
;                     in comments section of IEXTRACT_OPTIMAL_SPROF
;       opt_ignore  - array of mask values that should not be used
;                     when fitting the spatial profile or summing the
;                     1D spectrum (default=[2,8,16,32]) 
;       opt_nsig    - multiplicative factor for cosmic-ray sigma 
;                     clipping (default = 4)
;
; KEYWORDS:
;       opt_check   - display spatial profile fitting steps
;
; OPTIONAL OUTPUTS:
;       spec1d      - one-dimensional optimally extracted spectrum
;       sig1d       - corresponding one-dimensional error array
;       mask1d      - composite mask of combined pixels
;       sky1d       - one-dimensional "optimally extracted" sky spectrum
;       sprof       - two-dimensional spatial profile
;
; COMMENTS:  Modeled after the optimal extraction algorithm in Table 1
;            of 'AN OPTIMAL EXTRACTION ALGORITHM FOR CCD SPECTROSCOPY'
;            by Keith Horne (PASP, 1986, 98, 609).
;
;            Optimal extraction is only ever intended for point
;            sources where all the light originates from a common
;            origin.  In the case of extended objects, weighting by
;            the variance will discount light from regions which are
;            less luminous.
;
; PROCEDURES USED: IEXTRACT_OPTIMAL_SPROF(), IEXTRACT_OPTIMAL_1D()
;
; TODO: generate quality assurance plots
;
; MODIFICATION HISTORY:
;       written by A.R.Marble, Steward Observatory, June 2003
;       streamlined by A.R.Marble, 2004 Feb 9
;       SKYIMAGE converted to optional input by A.R.Marble, 2004 May 12
;-

pro iextract_optimal, imnosky, sigmamap, mask, apinfo, skyimage=skyimage, $
            opt_ncol=opt_ncol, $
            opt_ignore=opt_ignore, opt_nsig=opt_nsig, opt_check=opt_check, $
            spec1d=spec1d, sig1d=sig1d, mask1d=mask1d, sky1d=sky1d, sprof=sprof, $
            opt_gauss=opt_gauss

    up = floor(apinfo.upper) ; upper boundary index
    lo = ceil(apinfo.lower)  ; lower boundary index

    image = imnosky[*,min(lo):max(up)] ; extract 2D box around trace
    var   = (sigmamap[*,min(lo):max(up)])^2.
    msk   = mask[*,min(lo):max(up)] 
    if n_elements(skyimage) ne 0L then sky = skyimage[*,min(lo):max(up)] $
    else sky = image * 0.       ; extract 2D box around trace

    ncols = n_elements(image[*,0]) ; number of columns
    nrows = n_elements(image[0,*]) ; number of rows

; defaults

    if n_elements(opt_nsig)    eq 0L then opt_nsig    = 4.0
    if n_elements(opt_ignore)  eq 0L then opt_ignore  = [2,8,16,32] $
                                     else opt_ignore  = fix(opt_ignore)
    if n_elements(opt_ncol)    eq 0L then opt_ncol    = 41L

    nbins = ncols / opt_ncol    ; number of column bins

    colaxis = findgen(ncols)    ; index array for columns
    rowaxis = findgen(nrows)    ; index array for rows

; fit spatial profile
    
    badmask = bytarr(ncols, nrows) * 0L
    for j = 0L, n_elements(opt_ignore)-1L do begin
       wh = where(msk/opt_ignore[j] mod 2, count)
       if count gt 0 then badmask[wh] = 1L 
    endfor
    sprof = iextract_optimal_sprof(image, badmask, opt_ncol, opt_skip=opt_skip, $
                                   opt_check=0L, opt_gauss=opt_gauss)

; optimally extract 1d spectrum
    
    spec1d = iextract_optimal_1d(image, sprof, badmask, var, sky, opt_nsig, $
                                 sig1d=sig1d, mask1d=mask1d, sky1d=sky1d)

; refit spatial profile with improved cosmic-ray mask

    sprof = iextract_optimal_sprof(image, badmask, opt_ncol, opt_skip=opt_skip, $
                                   opt_check=opt_check, opt_gauss=opt_gauss)

; optimally re-extract 1d spectrum using improved spatial profile image

    spec1d = iextract_optimal_1d(image, sprof, badmask, var, sky, $
      opt_nsig, sig1d=sig1d, mask1d=mask1d, sky1d=sky1d)
 end
 









