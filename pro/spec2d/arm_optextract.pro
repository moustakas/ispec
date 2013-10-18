;+
; NAME:  
;
;   ARM_OPTEXTRACT
;       
; PURPOSE:  
;
;   Optimally extract the 1D spectrum of a point source from a 2D
;   image.
;
; CALLING SEQUENCE:
;
;   arm_optextract, spec2d, err2d, [mask2d=, psfile=, sky2d=, $
;     lower=, upper=, maskrej=, nsigma=, column_frac=, trace_ord=, $
;     spec1d=, err1d=, mask1d=, sky1d=, model=, spline_int=, $
;     spline_ord=, spline_nsig=, spline_sm=, nlevel=, outmask=,
;     outtrace=, trace=, column_nmin=, _extra=, floor=, /column_med,
;     /spline_eq, /primary, /silent]
;
; INPUTS:
;
;   spec2d - two-dimensional, sky-subtracted image frame with the
;            dispersion dimension coinciding with the column axis 
;            (wavelength changing with columns) 
;   err2d  - corresponding two-dimensional 1-sigma uncertainty map
;   
; OPTIONAL GENERAL INPUTS:
;
;   lower   - lower boundary of image rows to use (default is
;             smallest row in image); the rows bounded by LOWER and
;             UPPER should only include light from the single object for
;             which a spectrum is being extracted  
;   upper   - upper boundary of image rows to use (default is
;             largest row in image) 
;   mask2d  - two-dimensional bad pixel mask (default is all zeroes,
;             see MASKREJ)
;   maskrej - array of mask values corresponding to MASK2D values
;             which should not be used in extraction (default > 0)
;   sky2d   - two-dimensional sky image corresponding to SPEC2D
;   nsigma  - pixels deviating from the illumination model by NSIGMA
;             (or more) sigma will not be used when optimally
;             extracting the 1D arrays
;   floor   - the spatial profile is redefined to be zero where it
;             falls below FLOOR times the profile maximum (default = 0.01)
;   
; OPTIONAL TRACE INPUTS:
;
;   trace     - array of row positions for spatial profile peak in each
;               column (with respect to spec2d), overrides routine's
;               trace fit
;   trace_ord - order of polynomial to fit trace (position of
;               illumination peak as a function of column) with
;               (default = 3, linear is 1) 
;
; OPTIONAL SPATIAL PROFILE FIT INPUTS:
;
;   column_frac - two-element array indicating which fraction of columns
;                 to ignore when forming 1D spatial profile: fraction of
;                 columns with lowest flux followed by one minus the
;                 fraction of columns with highest flux 
;                 (default = [0.05, 0.95]) 
;   column_nmin - minimum number of medianed/averaged values required
;                 for the result to be used in the spatial profiel fit
;                 (default = 1)
;   spline_int  - breakpoint interval used by BSPLINE_ITERFIT() when
;                 fitting the 1D spatial profile; note that a value of
;                 100 corresponds to 1 row in the original image due to
;                 oversampling supplied by curvature in the trace
;                 (e.g., a value of 2 means that a breakpoint is placed
;                 at every other pixel in the 1D spatial profile)
;   spline_ord  - order for spline fit (default = 3)
;   spline_nsig - deviations of SPLINE_NSIG (or more) sigma will be
;                 rejected during spline fitting (default=4); to turn
;                 off rejection, set SPLINE_NSIG=0
;   spline_sm   - pixel width (corresponding to oversampled
;                 resolution) of smoothing kernel used by POLY_SMOOTH()
;                 to smooth the spline fit (default=0)     
;   _extra      - any inputs to BSPLINE_ITERFIT() will be passed to
;                 that routine, overiding those specified by
;                 ARM_OPTEXTRACT 
;
; OPTIONAL PLOT INPUTS:
;
;   psfile - filename (including path) of postscript file to generate
;            for diagnostic plots (default = none)  
;   nlevel - number of evenly spaced intervals to use for contour plot 
;
; KEYWORDS:
;
;   column_med - use median values when forming 1D spatial profile,
;                rather than error-weighted means 
;   spline_eq  - weight 1D spatial profile values equally when
;                performing spline fit (instead of weighting them by
;                their variance)
;   primary    - redefine the 1D spatial profile to be zero beyond the
;                point where the profile first falls below FLOOR times
;                the profile maximum (on either side of the peak),
;                this can be useful for elliminating contamination
;                from spurious emission beyond the illumination
;                profile 
;   silent     - suppress messages
;
; OPTIONAL OUTPUTS:
;
;   spec1d   - one-dimensional optimally extracted spectrum
;   err1d    - corresponding one-dimensional error array
;   mask1d   - corresponding one-dimensional mask array (0:good, 1:bad)  
;   sky1d    - one-dimensional optimally extracted sky spectrum
;   model    - two-dimensional illumination model
;   outmask  - output mask with dimensions equal to SPEC2D and the
;              following values: 0 - good
;                                1 - N/A (row < LOWER or row > UPPER)
;                                2 - already rejected by MASK2D
;                                4 - rejected by ARM_SPECEXTRACT
;   outtrace - trace fit by ARM_OPTEXTRACT, coordinates refer to SPEC2D
;              (not region defined by LOWER/UPPER)
;
; PROCEDURES USED: 
;
;   ARM_PLOTCONFIG, ARM_ORDEROFMAG(), ARM_OPLOTLINE,
;   BSPLINE_ITERFIT(), BSPLINE_VALU(), DJS_PLOT, DJS_OPLOT, POLY(),
;   POLY_SMOOTH(), TEXTOIDL(), TRACE_FWEIGHT() 
;
; OVERVIEW:
;
;   Unlike traditional spectral extraction methods, this routine
;   models the 2D illumination pattern of a ccd image in order to most
;   reliably reconstruct 1D spectral information.  This has two
;   primary benefits.  First, such a model allows for robust
;   identification of outlier pixels (eg, cosmic rays) including
;   residual effects not completely elliminated by other means (eg,
;   combining multiple images).  For this purpose, the higher the data
;   quality, the better the results.  Secondly, properly modeling the
;   2D illumination pattern allows for unequal weighting of pixels
;   without degradation of spectral information (only appropriate for
;   point sources observations!).  Variance-weighting (or so called
;   "optimal extraction") results in an improvement in signal-to-noise
;   of the extracted spectrum.  For this purpose, the lower the data
;   quality (within reason), the greater the improvement.
;
;   Reliable results are contingent upon construction of a successful
;   model (which in turn requires accurate error information).  The
;   two components of this model are the 1D spatial profile
;   (illumination pattern along the pixels of a single row) and the 
;   trace (changing position of the 1D spatial profile as a function
;   of column number).  Because the nature of these components can
;   vary widely, this routine necessarily incorporates many adjustable
;   parameters.  The user can best gauge the reliability of the
;   trace/profile fitting using the diagnostic plots printed to
;   postscript (see PSFILE).  In the event that parameter defaults are
;   not sufficient, a full description of the algorithm is available
;   for reference below.  
;
; FEATURES:
;
;   ARM_OPTEXTRACT has a number of options which may be helpful in
;   various circumstances.  In addition to being documented above
;   as input and output parameters, some are highlighted here.
;
;   * optimal extraction
;   * identification of outlier pixels (see OUTMASK)
;   * flexible interpretation of input masks (see MASKREJ)
;   * user-supplied trace for instances where it is difficult to fit
;     (eg, too noisy, heavily absorbed spectrum) but known from other
;     observations (see TRACE and OUTTRACE)
;   * improved sampling of the 1D spatial profile due to trace removal 
;   * ability to affect which columns are used to construct 1D spatial
;     profile (see COLUMN_FRAC, COLUMN_NMIN and COLUMN_MED)
;   * ability to customize spline fit to 1D spatial profile as needed
;     (see SPLINE_ parameters; also parameters for BSPLINE_ITERFIT()
;     routine can be passed through _EXTRA)
;
; DETAILS:
;
;   Standard methods for extracting a 1D spectrum from a 2D ccd image
;   weight all illuminated pixels equally, despite the fact that they
;   are not illuminated equally and thus have unequal signal to noise 
;   ratios.  This is done in order to include all light at a given
;   wavelength, preserving spectrophotometric accuracy.  When summing
;   the light contributions of pixels illuminated in the spatial
;   direction, one would prefer to weight them by their uncertainty,
;   reducing the combined noise.  However, since they are not
;   illuminated equally, this would result in a loss of spectroscopic
;   information. Optimal extraction is a methodology which allows you
;   to weight pixels by their variance while avoiding any spectral
;   degradation.  This is made possible by modeling the pattern of
;   changing illumination across the chip.  
;
;   It is important to note that optimal extraction is only intended
;   for point sources where all the light originates from a common
;   origin.  In the case of extended objects, weighting pixels by
;   variance may discount light from less luminous regions.  Also,
;   proper modelling of the illumination pattern requires that the
;   spectral signature remain constant for all illuminated rows.

;   The majority of this routine is dedicated to modeling the
;   illumination of the 2D ccd image.  Often, the spatial profile is a
;   simple gaussian curve.  However, factors such as variable seeing,
;   combining somewhat incoincident images or observing extended
;   objects can result in a spatial profile which is not characterized 
;   well by any function.  Similarly, the position of the spatial
;   profile as a function of column number (the trace) can migrate
;   significantly across the image, resulting in distinctly non-linear
;   profiles in the dispersion direction.  This migration is expected
;   to vary smoothly from pixel to pixel; however, the net effect can
;   range from flat lines to double peaked profiles depending on where
;   the illumination pattern intersects a given row.  Thus, reliable
;   modeling of the 2D illumination pattern requires a high degree of
;   generality, which in turns necessitates a high degree of
;   flexibility in customizing parameters as needed.  
;
;   The illumination model is essentially constructed by first fitting
;   and removing the trace, fitting the 1D spatial profile resulting
;   from collapsing the column axis, and resampling that fit for each
;   column with the trace restored.  The actual algorithm is described
;   in detail (and in order of operation for the purpose of
;   documenting the code) below. 
;
;   The 2D arrays are cropped to retain only the relevant rows
;   described by LOWER and UPPER.  This region of interest should only
;   include light from the single object for which a spectrum is being
;   extracted.  Additional objects with different spectral signatures
;   can adversely affect the construction of the illumination model.
;   A modified version of MASK2D is created, comprised of only ones
;   (ignore) and zeroes (use).  MASKREJ allows the user to specify how
;   different MASK2D values will be interpreted, accomodating mask
;   values which may indicate properties that do not reflect unusable
;   data.  The error input array ERR2D is checked for unmasked pixels
;   with [unphysical] zero entries, and the mask array is modified to
;   exclude any found. 
;
;   The trace (treated here as the running flux-weighted center of the
;   spatial profile) is fit by a TRACE_ORD order polynomial.  The
;   reliability of this fit can be checked visually via inspection of
;   the output postscript page defined by PSFILE (NLEVEL contours from
;   a normalized version of SPEC2D are plotted in the first window
;   with the trace represented by a dashed line in both the first and
;   second windows).  Note that the fitted trace can be returned using
;   the OUTTRACE parameter.  Also, for cases where the user is
;   convinced that the trace is already known, it can be supplied
;   using TRACE rather than being fit by ARM_OPTEXTRACT.  This may be
;   especially useful if the trace is difficult to fit reliably due to
;   noise or heavily absorbed spectra.
;
;   Spectral variations are removed by normalizing the area under each
;   column to unity.  Then a subset of columns suitable for forming
;   the 1D spatial profile are selected.  The first criterion for this
;   selection is defined by COLUMN_FRAC, which allows the user to
;   deselect columns with the lowest and highest cumulative values
;   which might be affected by increased noise or spurious pixels (eg,
;   cosmic ray hits).  Any unmasked columns with zero flux will be
;   deselected regardless of how COLUMN_FRAC is defined.  The standard
;   deviation is calculated for each of the remaining rows.  Any
;   column deviating from the resulting distribution by 3-sigma (or
;   greater) is deselected as a further precaution against bias due to
;   cosmic rays or ccd defects.
;
;   The selected subset of columns are individually shifted to remove
;   the fitted trace to one hundredth of a pixel.  This aligns the
;   spatial profiles and allows them to be collapsed, forming an
;   oversampled (if the trace is not a constant value) 1D spatial
;   profile.  Note that this assumes that the shape of the
;   illumination pattern does not change from column to column, an
;   assumption which is expected to always be valid.  Columns shifted
;   by the same value are represented in the collapsed profile by
;   their error-weighted mean (unless the COLUMN_MED keyword is set,
;   in which case the median is used).  The result is shown in the
;   third window of the postscript page as black dots with grey bars
;   reflecting the 1-sigma uncertainties.  Depending on the shape of
;   the trace, different average/median values in the 1D spatial
;   profile will reflect different numbers of columns.  This number,
;   N, is plotted below the 1D spatial profile in the fourth window.
;   To reduce the affects of noisy or biased values, the user can use
;   COLUMN_NMIN to specify a minimum value of N required for the
;   average/median to be included in the 1D spatial  profile.  Values
;   corresponding to N < COLUMN_NMIN are plotted in orange as opposed
;   to black. 
;
;   A spline curve is fit to the 1D spatial profile using the
;   BSPLINE_ITERFIT() routine.  Breakpoints between which SPLINE_ORD
;   order polynomials are fit are placed at intervals of SPLINE_INT
;   oversampled pixels (recall that 100 oversampled pixels is the
;   equivalent of 1 pixel at the data's original resolution).  Profile
;   values are weighted by their propogated variance (unless SPLINE_EQ
;   is used to specify equal weighting).  Outliers deviating from
;   the fit by SPLINE_NSIG sigma (or greater) are iteratively
;   rejected and indicated on the plot by red X's.  Additional
;   customization of the spline fit can be  achieved by passing
;   parameters used by BSPLINE_ITERFIT().  These values will override
;   any defined by ARM_OPTEXTRACT.  Finally, if SPLINE_SM is non-zero,
;   the spline fit will be smoothed with a kernel width of SPLINE_SM
;   oversampled pixels.  The fitted profile is forced to be greater
;   than or equal to zero.  Values falling below FLOOR times the
;   profile maximum are redefined as zero.  If the PRIMARY keyword is
;   set, points beyond where the profile first reaches such a value
;   (FLOOR times the profile maximum) on either side of the profile
;   peak are all redefined as zero.  This can prevent the spatial
;   profile from being adversely affected by spurious light which was
;   not removed by other means.  The resulting profile is overplotted
;   with a green line.   
;
;   The 2D illumination model is constructed by normalizing the area
;   under the supersampled spatial profile to unity and resampling it
;   (at the original resolution) for each image column according to
;   the corresponding trace values.
;
;   A 1D spectrum is iteratively extracted (weighted by variance) from
;   the ccd image.  With each iteration, the most deviant (a greater
;   than NSIGMA sigma outlier from the expected shape of the
;   illumination model) un-masked pixel per column is identified and
;   masked until none remain, and the spectrum (SPEC1D) is
;   recalculated.  Once this is completed, the corresponding error
;   (ERR1D) and sky (SKY1D, if SKY2D was passed) arrays are similarly
;   extracted.  Columns for which no unmasked pixels are available are
;   give a mask (MASK1D) value of 1, otherwise 0.
;
;   For more information on optimal extraction, refer to 'AN OPTIMAL
;   EXTRACTION ALGORITHM FOR CCD SPECTROSCOPY' by Keith Horne (PASP,
;   1986, 98, 609), which describes an alternate algorithm.
;
; MODIFICATION HISTORY:
;
;   written by A.R.Marble, Steward Observatory, 2005 May 4
;-

pro arm_optextract, spec2d, err2d, mask2d=mask2d, sky2d=sky2d, lower=lower, $
       upper=upper, column_frac=column_frac, nsigma=nsigma, psfile=psfile, $
       trace=trace, trace_ord=trace_ord, column_nmin=column_nmin, floor=floor, $
       spec1d=spec1d, err1d=err1d, mask1d=mask1d, sky1d=sky1d, model=model, $
       maskrej=maskrej, spline_int=spline_int, spline_ord=spline_ord, $
       spline_sm=spline_sm, _extra=extra, nlevel=nlevel, outmask=outmask, $
       spline_nsig=spline_nsig, COLUMN_MED=column_med, $
       SPLINE_EQ=spline_eq, PRIMARY=primary, outtrace=outtrace, SILENT=silent

; defaults, initialization & error-checking

    if N_ELEMENTS(err2d) eq 0L then begin
       PRINT, 'Syntax - ARM_OPTEXTRACT, spec2d, err2d, [mask2d=, sky2d=, trace=, outtrace=, $'
       PRINT, '         psfile=, lower=, upper=, column_frac=, nsigma=, psfile=, $' 
       PRINT, '         trace_ord=, spec1d=, err1d=, mask1d=, sky1d=, model=, outmask=, $'
       PRINT, '         maskrej=, spline_int=, spline_ord=, spline_sm=, nlevel=, $'
       PRINT, '         column_nmin=, spline_nsig=, _extra=, /column_med, /spline_eq, /silent]'
       return
    endif
    
    if N_ELEMENTS(spline_int)  eq 0L then spline_int  = 100
    if N_ELEMENTS(spline_ord)  eq 0L then spline_ord  = 3
    if N_ELEMENTS(spline_sm)   eq 0L then spline_sm   = 0
    if N_ELEMENTS(nsigma)      eq 0L then nsigma      = 4.0
    if N_ELEMENTS(trace_ord)   eq 0L then trace_ord   = 3
    if N_ELEMENTS(nlevel)      eq 0L then nlevel      = 10 else nlevel = ROUND(nlevel) > 1
    if N_ELEMENTS(column_nmin) eq 0L then column_nmin = 1  else column_nmin = ROUND(column_nmin) > 1
    if N_ELEMENTS(spline_nsig) eq 0L then spline_nsig = 4 
    if N_ELEMENTS(floor)       eq 0L then floor       = 0.01

    if N_ELEMENTS(mask2d) eq 0L then begin
       mask2d = spec2d * 0
       maskrej = 1L
    endif

    if N_ELEMENTS(spec2d[*,0]) ne N_ELEMENTS(err2d[*,0]) or $
      N_ELEMENTS(spec2d[*,0]) ne N_ELEMENTS(mask2d[*,0]) then begin
       MESSAGE, 'SPEC2D, ERR2D and/or MASK2D dimensions do not agree.', /continue
       return
    endif

    def_frac = [0.05, 0.95]
    if N_ELEMENTS(column_frac) eq 0L then frac = def_frac else frac = column_frac
    if frac[1] lt frac[0] then begin
       MESSAGE, 'Invalid COLUMN_FRAC values, using default.', /continue
       frac = def_frac
    endif

    if KEYWORD_SET(lower) then begin
       lo = FLOOR(lower[0])
       if lo ne lower[0] then MESSAGE, 'LOWER row boundary changed from '+ $
         STRTRIM(lower[0],2)+' to '+STRTRIM(lo,2)+'.', /continue
    endif else lo = 0L

    if KEYWORD_SET(upper) then begin
       up = CEIL(upper[0]) > 0 
       if up ne upper[0] then MESSAGE, 'UPPER row boundary changed from '+ $
         STRTRIM(upper[0],2)+' to '+STRTRIM(up,2)+'.', /continue
    endif else up = N_ELEMENTS(spec2d[0,*])-1 

; crop desired region of image

    im    = spec2d[*,lo:up]     
    var   = (err2d[*,lo:up])^2d
    msk   = mask2d[*,lo:up] 
    if N_ELEMENTS(sky2d) ne 0L then sky = sky2d[*,lo:up] else sky = im * 0.          

    ncols = N_ELEMENTS(im[*,0]) ; number of columns
    nrows = N_ELEMENTS(im[0,*]) ; number of rows
    
; interpret mask

    mask = BYTARR(ncols, nrows) * 0L
    if N_ELEMENTS(maskrej) gt 0L then begin
       for j = 0L, N_ELEMENTS(maskrej)-1L do begin
          wh = WHERE(msk/maskrej[j] mod 2, count)
          if count gt 0 then mask[wh] = 1L 
       endfor
    endif else begin
       wh = WHERE(msk ne 0, count)
       if count gt 0 then mask[wh] = 1L
    endelse
    origmask = mask

; check error image for negative and zero values

    zero = WHERE(err2d[*,lo:up] le 0 and mask eq 0, nzero)
    if nzero gt 0 then begin
       mask[zero] = 1L
       MESSAGE, 'WARNING: '+STRTRIM(STRING(nzero,f='(i)'),2)+' unmasked '+ $
         'negative/null values in ERR2D have been found and masked.', /continue
    endif

; trace image

    if not KEYWORD_SET(silent) then MESSAGE, 'constructing 2D illumination model...', /cont

    if N_ELEMENTS(trace) ne 0L then begin
       if N_ELEMENTS(trace) ne nrows then MESSAGE, $
         'TRACE has incompatible number of elements, using fitted trace.', /continue $
       else if MIN(trace-lo) lt 0 then MESSAGE, $
         'TRACE does not lie between LOWER and UPPER, using fitted trace.', /continue $
       else tr = trace
    endif

    if N_ELEMENTS(tr) ne nrows then begin
       colaxis = findgen(ncols)
       cen = FLTARR(ncols)
       for j = 0L, ncols-1L do begin
          dummy = MAX(im[j,*]*(1d0-mask[j,*]), cenmax)
          cen[j] = cenmax
       endfor
       for i = 0L, 2L do cen = TRACE_FWEIGHT(TRANSPOSE(im), cen, colaxis)
       DJS_ITERSTAT, cen, sigrej=2.5, sigma=censig, median=cenmed, mask=cenmask
       cengood = WHERE(cenmask eq 1B, ngood)
       if ngood ne 0L then yc = $
         POLY(colaxis, POLY_FIT(colaxis[cengood], cen[cengood], trace_ord)) $
       else begin
          MESSAGE, 'ABORTED: No good points to fit to trace.', /continue
          return
       endelse       
    endif else yc = trace - lo

; form illumination model

    model = im * 0.0
    nim = DBLARR(ncols, nrows)
    nerr = DBLARR(ncols, nrows)
    colttl = TOTAL(im*(1d0-mask),2)
    goodcol = WHERE(colttl gt 0, ngoodcol)
    if ngoodcol gt 0L then begin
       nim[goodcol,*] = im[goodcol,*] / (colttl[goodcol] # (DBLARR(nrows)+1))
       nerr[goodcol,*] = sqrt(var[goodcol,*]) / (colttl[goodcol] # (DBLARR(nrows)+1))
    endif else begin
       MESSAGE, 'ABORTED: No unmasked pixels with non-zero flux!', /continue
       return
    endelse

;   determine which columns to use for creating 1D spatial profile

    colttl = TOTAL(im, 2)
    use1 = (SORT(colttl))[((frac[0]*ncols)>0)<(ncols-1):(frac[1]*ncols)<(ncols-1)]
    nuse1 = N_ELEMENTS(use1)

    zerottl = WHERE(colttl[use1] le 0 and TOTAL(mask[use1,*],2) lt nrows, nzerottl)
    if nzerottl gt 0L then begin
       if nzerottl eq nuse1 then begin
          MESSAGE, 'ABORTED: No unmasked pixels with non-zero flux!', /continue
          return
       endif
       REMOVE, zerottl, use1
       nuse1 = N_ELEMENTS(use1)
    endif

    std = FLTARR(nuse1)
    for i = 0L, nuse1-1 do std[i] = SQRT(VARIANCE(nim[use1[i],*]))
    DJS_ITERSTAT, std, mean=stdavg, sigma=stdstd
    outlier = WHERE(ABS(std-stdavg)/stdstd gt 3, noutlier, complement=ok)
    if noutlier lt nuse1 then begin
       if noutlier gt 0L then begin
          use2 = use1[ok]
          nuse2 = N_ELEMENTS(use2)
       endif else begin
          use2 = use1
          nuse2 = nuse1
       endelse
    endif

;   create 1D spatial profile
    
    nint = 100.
    ycr = ROUND(yc*nint)
    refycr = ROUND(MEAN(yc*nint))
    shift = refycr - ycr

    sprof1d = DBLARR((nrows-1)*nint+1)
    errsprof1d = DBLARR((nrows-1)*nint+1)
    nsprof1d = LONARR((nrows-1)*nint+1)

    ushift = shift[UNIQ(shift, SORT(shift))]

    need = BYTARR((nrows-1)*nint+1)
    for i = 0L, N_ELEMENTS(ushift)-1 do begin
       match = WHERE(shift[use2] eq ushift[i], nmatch)
       pxl = FINDGEN(nrows)*nint+ushift[i]
       need[pxl[WHERE(pxl ge 0 and pxl le (nrows-1)*nint)]] = 1
       if nmatch gt 0L then begin
          for j = 0, nrows-1 do begin
             if j+ushift[i]/nint ge 0 and j+ushift[i]/nint le (nrows-1) then begin
                valid = WHERE(mask[use2[match],j] eq 0, nvalid)
                if nvalid gt 0L then begin
                   nsprof1d[j*nint+ushift[i]] = nvalid
                   vals = nim[use2[match[valid]],j]
                   errs = nerr[use2[match[valid]],j]
                   if KEYWORD_SET(column_med) then sprof1d[j*nint+ushift[i]] = median(vals) $
                   else sprof1d[j*nint+ushift[i]] = (TOTAL(vals/errs^2) / TOTAL(1d/errs^2))  
                   errsprof1d[j*nint+ushift[i]] = SQRT(1d/TOTAL(1d/errs^2))
                endif
             endif
          endfor
       endif
    endfor

    valid = WHERE(nsprof1d ge column_nmin, nvalid)
    if nvalid eq 0L then begin
       column_nmin = 1
       MESSAGE, 'COLUMN_NMIN is too large... changed to 1.', /continue
       valid = WHERE(nsprof1d ge column_nmin, nvalid)
    endif

;   fit spline to 1D spatial profile

    if KEYWORD_SET(spline_eq) then invvar = valid * 0 + 1. $
    else invvar = 1d0 / (errsprof1d[valid]^2d)

    bkpt=LINDGEN(((nrows-1)*nint+1)/spline_int)*spline_int
    if bkpt[n_elements(bkpt)-1] ne (nrows-1)*nint then bkpt = [bkpt, (nrows-1)*nint]

    sset = BSPLINE_ITERFIT(valid, sprof1d[valid], nord=spline_ord, $
                           bkpt=bkpt, $
                           invvar=invvar, upper=spline_nsig, $
                           lower=spline_nsig, _extra=extra, outmask=bsplinemask)
    sprof1dfit = BSPLINE_VALU(DINDGEN(N_ELEMENTS(sprof1d)), sset) > 0

    if spline_sm ne 0 then sprof1dfitsm = POLY_SMOOTH(sprof1dfit, spline_sm*nint) > 0 $
      else sprof1dfitsm = sprof1dfit

;   floor spatial profile to zero where it is < FLOOR times profile max

    floored = WHERE(sprof1dfitsm le floor*MAX(sprof1dfitsm), nfloored)
    if nfloored gt 0L then begin
       floorfunc = LONARR((nrows-1)*nint+1) + 1
       floorfunc[floored] = 0
       if KEYWORD_SET(primary) then begin
          left = WHERE(floored lt refycr, nleft)
          if nleft gt 0L then floorfunc[0:MAX(floored[left])] = 0
          right = WHERE(floored gt refycr, nright)
          if nright gt 0L then floorfunc[MIN(floored[right]):(nrows-1)*nint] = 0
       endif
       if MAX(floorfunc) eq 1 then sprof1dfitsm = sprof1dfitsm * floorfunc $
       else MESSAGE, 'spatial profile error, disregarding PRIMARY keyword.', /continue
    endif

;   normalize area under spatial profile curve to unity

    area = TOTAL(sprof1dfitsm) / nint
    if area eq 0 then begin
       MESSAGE, 'area under the spatial profile is zero! ABORTING.', /continue
       return
    endif
    sprof1dfitsm = sprof1dfitsm / area
    sprof1d      = sprof1d      / area
    errsprof1d   = errsprof1d   / area

;   resample fitted profile to form 2D illumination model

    for i = 0L, ncols-1 do begin
       xshift = LINDGEN(nrows)*nint + shift[i]
       valid = WHERE(xshift ge 0 and xshift le (nrows-1)*nint)
       model[i,valid] = sprof1dfitsm[xshift[valid]]
    endfor

; output trace (if desired)

    outtrace = yc + lo

; extract spectra
       
    if not KEYWORD_SET(silent) then MESSAGE, 'optimally extracting spectra...', /continue

    err1d  = FLTARR(ncols)
    sky1d  = FLTARR(ncols)
    mask1d = FLTARR(ncols)

    spec1d = TOTAL(im, 2) ; includes any residual cosmic rays

    for i = 0L, ncols-1L do begin          

       old = -1L       ; number of bad pixels from previous iteration
       
       while old lt TOTAL(mask[i,*]) do begin

          old = TOTAL(mask[i,*])

;   exclude most deviant pixel until none remain

          deviations = (1-mask[i,*]) * (im[i,*] - spec1d[i] * model[i,*])^2. / var[i,*] 
          if MAX(deviations) gt nsigma^2 then $
            mask[i, (WHERE(deviations eq MAX(deviations)))[0]] = 1L

;   optimally extract flux

          good = WHERE(mask[i,*] eq 0, count)
          if count gt 0L then spec1d[i] = $
            TOTAL(model[i,good] * im[i,good] / var[i,good]) / $
            TOTAL(model[i,good]^2 / var[i,good]) else spec1d[i] = 0

       endwhile   

;   optimally extract sky

       if count gt 0L then sky1d[i] = $
         TOTAL(model[i,good] * sky[i,good] / var[i,good]) / $
         TOTAL(model[i,good]^2. / var[i,good]) else sky1d[i] = 0

;   variance of optimal extraction

       if count gt 0L then err1d[i] = $
         SQRT(TOTAL(model[i,good]) / $
              TOTAL(model[i,good]^2. / var[i,good]))

;   collapse mask image 

       illum = WHERE(model[i,*] gt MAX(model)*floor, nillum)
if nillum eq 0 then stop
       if MIN(mask[i,illum]) gt 0 then mask1d[i] = FIX(1)

       outmask = BYTE(mask2d * 0) + 1
       outmask[*, lo:up] = origmask * 2
       newmask = WHERE(mask gt 0 and origmask eq 0, nnewmask)
       if nnewmask gt 0 then outmask[newmask] = 4

    endfor

; prepare plot information

    if N_ELEMENTS(psfile) gt 0L then begin
       good = WHERE(sprof1d gt MAX(sprof1d)/1d2, ngood)
       xr = [(MIN(good)/nint)>0,(MAX(good)/nint)<(N_ELEMENTS(nim[0,*])-1)]
       yr = [0-.1*MAX(sprof1d), 1.2*MAX(sprof1d)]
    endif

;; ; plot to screen (if desired)
;; 
;;     if KEYWORD_SET(doplot) then begin
;; 
;;        WINDOW, 0, xsize=600, ysize=400
;;        DJS_PLOT, [0,0], /nodate, /noerase, xthick=2, ythick=2, charth=2, xr=xr, chars=1.5, $
;;          yr=yr, xst=3, yst=1, ytit='normalized 1D spatial profile', $
;;          xtickname=REPLICATE(' ',60)
;;        valid = WHERE(nsprof1d ge column_nmin, nvalid)
;;        discarded = WHERE(bsplinemask eq 0, ndiscarded)
;;        if ndiscarded gt 0L then DJS_OPLOT, valid[discarded]/nint, $
;;          sprof1d[valid[discarded]], psym=7, color='red', thick=2, symsize=.5
;; 
;;        if nvalid gt 0L then begin
;;           ymn = sprof1d[valid] - 1*errsprof1d[valid]
;;           ypl = sprof1d[valid] + 1*errsprof1d[valid]
;;           for i = 0L, nvalid-1 do DJS_OPLOT, valid[i]/nint*[1,1], [ymn[i],ypl[i]], $
;;             thick=2, color='grey'
;;           DJS_OPLOT, valid/nint, sprof1d[valid], ps=2, symsize=.2
;;        endif
;;        reject = WHERE(nsprof1d gt 0 and nsprof1d lt column_nmin, nreject)
;;        if nreject gt 0L then begin
;;           ymn = sprof1d[reject] - 1*errsprof1d[reject]
;;           ypl = sprof1d[reject] + 1*errsprof1d[reject]
;;           for i = 0L, nreject-1 do DJS_OPLOT, reject[i]/nint*[1,1], [ymn[i],ypl[i]], $
;;             thick=2, color='orange'
;;           DJS_OPLOT, reject/nint, sprof1d[reject], ps=2, symsize=.2, color='orange'
;;        endif
;;        needed = WHERE(need eq 1, nneeded)
;;        cont = needed[1:nneeded-1] - needed[0:nneeded-2]
;;        x1 = needed[0]
;;        x2 = needed[nneeded-1]
;;        jump = WHERE(cont ne 1, njump)
;;        if njump gt 0L then begin
;;           x1 = [x1, needed[jump+1]]
;;           x2 = [needed[jump], x2]
;;        endif
;;        for i = 0L, n_elements(x1)-1 do begin
;;           xx = LINDGEN(x2[i]-x1[i]+1)+x1[i]
;;           DJS_OPLOT, xx/nint, sprof1dfitsm[xx], color='green', thick=4
;;        endfor
;;         
;;     endif

; plot to postscript file (if desired)

    if N_ELEMENTS(psfile) gt 0L then begin

       if not KEYWORD_SET(silent) then MESSAGE, 'plotting to postcript...', /cont
    
       oldpmulti = !p.multi
       ARM_PLOTCONFIG, psfile=psfile, /writeover, nx=1, ny=5, ymargin=[1.0,1.0], $
         xmargin=[1.0,0.5], yspace=[0.05,0.5,0.05,0.5], height=[1.2,0.8,2.55,0.8,2.5], coord=crd
       title = 'ARM_OPTEXTRACT diagnostic plots ('+ $
         (REVERSE(STRSPLIT(psfile,'/',/extract)))[0]+')!Csee documentation in routine header'
       XYOUTS, 0.5, 0.95, charth=2, chars=1, title, align=.5, /normal

;   trace plots
       
       dummy = MIN(yc[[0,ncols/2,ncols-1]], mn)
       CONTOUR, nim, pos=crd[*,0], xst=1, yst=1, xthick=2, ythick=2, $
         charthick=2, levels = (FINDGEN(nlevel)/(nlevel)+(0.5/nlevel)) * $
         MAX((nim[use2,*])[WHERE(FINITE(nim[use2,*]))]), $
         chars=1.5, xtickname=REPLICATE(' ',60), yminor=1, ytit='row'
       DJS_OPLOT, yc, color='white', thick=3
       DJS_OPLOT, yc, color='black', thick=3, line=2

       PLOT, yc, thick=3, line=2, xst=1, yst=3, xthick=2, ythick=2, /noerase, $
         charthick=2, pos=crd[*,1], chars=1.5, yr=MINMAX(yc), xtitle='column'
       PLOT, [0,1], [0,1], /nodata, /noerase, xst=5, yst=5, pos=crd[*,1]
       XYOUTS, ([0.03,0.5,.97])[mn], 0.8, chars=.8, charth=2, align=([0,.5,1])[mn], $
         'TRACE_ORD='+STRTRIM(STRING(trace_ord,f='(i)'),2)+', flux-weighted trace'
       XYOUTS, ([0.03,0.5,.97])[2-mn], 0.1, chars=.8, charth=2, align=([0,.5,1])[2-mn], $
         '(enlarged view without NLEVEL='+STRTRIM(STRING(nlevel,f='(i)'),2)+' contours)'

;   spatial profile plot

       DJS_PLOT, [0,0], /nodate, /noerase, xthick=2, ythick=2, charth=2, xr=xr, chars=1.5, $
         yr=yr, xst=3, yst=1, pos=crd[*,2], ytit='normalized 1D spatial profile', $
         xtickname=REPLICATE(' ',60)
       valid = WHERE(nsprof1d ge column_nmin, nvalid)
       discarded = WHERE(bsplinemask eq 0, ndiscarded)
       if ndiscarded gt 0L then DJS_OPLOT, valid[discarded]/nint, $
         sprof1d[valid[discarded]], psym=7, color='red', thick=2, symsize=.5
       reject = WHERE(nsprof1d gt 0 and nsprof1d lt column_nmin, nreject)
       if nreject gt 0L then begin
          ymn = sprof1d[reject] - 1*errsprof1d[reject]
          ypl = sprof1d[reject] + 1*errsprof1d[reject]
          for i = 0L, nreject-1 do DJS_OPLOT, reject[i]/nint*[1,1], [ymn[i],ypl[i]], $
            thick=2, color='orange'
          DJS_OPLOT, reject/nint, sprof1d[reject], ps=2, symsize=.2, color='orange'
       endif
       if nvalid gt 0L then begin
          ymn = sprof1d[valid] - 1*errsprof1d[valid]
          ypl = sprof1d[valid] + 1*errsprof1d[valid]
          for i = 0L, nvalid-1 do DJS_OPLOT, valid[i]/nint*[1,1], [ymn[i],ypl[i]], $
            thick=2, color='grey'
          DJS_OPLOT, valid/nint, sprof1d[valid], ps=2, symsize=.2
       endif
       needed = WHERE(need eq 1, nneeded)
       cont = needed[1:nneeded-1] - needed[0:nneeded-2]
       x1 = needed[0]
       x2 = needed[nneeded-1]
       jump = WHERE(cont ne 1, njump)
       if njump gt 0L then begin
          x1 = [x1, needed[jump+1]]
          x2 = [needed[jump], x2]
       endif
       for i = 0L, n_elements(x1)-1 do begin
          xx = LINDGEN(x2[i]-x1[i]+1)+x1[i]
          DJS_OPLOT, xx/nint, sprof1dfitsm[xx], color='green', thick=4
       endfor
       PLOT, [0,0], /nodata, /noerase, pos=crd[*,2], xst=5, yst=5, xr=[0,1], yr=[0,1]
       if spline_int eq 100 then rowtxt = ' row' else rowtxt = ' rows'
       XYOUTS, 0.03, 0.925, align=0, charth=2, chars=0.8, $
         'SPLINE_INT = '+STRTRIM(STRING(spline_int,f='(i)'),2) + $
         ' ('+STRTRIM(STRING(spline_int/1d2,f='(f4.2)'),2)+rowtxt+')'
       XYOUTS, 0.03, 0.875, align=0, charth=2, chars=0.8, $
         'SPLINE_ORD = '+STRTRIM(STRING(spline_ord,f='(i)'),2)
       XYOUTS, 0.03, 0.825, align=0, charth=2, chars=0.8, $
         'SPLINE_NSIG = '+STRTRIM(STRING(spline_nsig,f='(i)'),2)
       XYOUTS, 0.03, 0.775, align=0, charth=2, chars=0.8, $
         'SPLINE_SM   = '+STRTRIM(STRING(spline_sm,f='(i)'),2)
       XYOUTS, 0.72, 0.925, align=0, charth=2, chars=0.8, $
         'COLUMN_FRAC = ['+STRTRIM(STRING(frac[0],f='(f4.2)'),2)+', '+$
         STRTRIM(STRING(frac[1],f='(f4.2)'),2)+']'
       XYOUTS, 0.72, 0.87, align=0, charth=2, chars=0.8, $
         'COLUMN_NMIN = '+STRTRIM(STRING(column_nmin,f='(i)'),2)
       XYOUTS, 0.72, 0.815, align=0, charth=2, chars=0.8, $
         'FLOOR = '+STRTRIM(STRING(floor,f='(f5.3)'),2)
       PLOT, [0,0], /nodata, /noerase, xthick=2, ythick=2, pos=crd[*,3], xr=xr, xst=3, $
         yr=[0.11,2*MAX(nsprof1d[valid])], yst=1, charsize=1.5, charthick=2, xtit='row', $
         ytit=TEXTOIDL('N'), /ylog
       if nvalid gt 0 then OPLOT, valid/nint, nsprof1d[valid], ps=2, symsize=.2
       if nreject gt 0 then DJS_OPLOT, reject/nint, nsprof1d[reject], ps=2, symsize=.2, $
         color='orange'
        
       good = WHERE(mask1d eq 0, ngood, complement=bad)
       nonopt = TOTAL(im,2)

       nonoptgood = WHERE(TOTAL(origmask, 2) eq 0, nnonoptgood, complement=nonoptbad)
       nonoptmask = INTARR(ncols, nrows)
       illum = WHERE(model gt MAX(model)/1d2, nillum)
       if nillum gt 0L then nonoptmask[illum] = 1
       nonopt = TOTAL(im*nonoptmask, 2)

       power = ARM_ORDEROFMAG(MEDIAN(spec1d))

       if ngood gt 0L then yr = [0, 1.33*MAX(spec1d[good])/1d1^power] else $
         if nnonoptgood gt 0L then yr = [0,1.33*MAX(nonopt)] else yr = [0,1]

       PLOT, [0,0], /nodata, /noerase, xr=[0,ncols-1], yr=yr, xst=7, yst=5, pos=crd[*,4]
       if nnonoptgood lt ncols then ARM_OPLOTLINE, nonoptbad, color='green', line=2
       if ngood lt ncols then ARM_OPLOTLINE, bad, color='black', line=1
       if nnonoptgood gt 0L then DJS_OPLOT, nonoptgood, nonopt[nonoptgood]/1d1^power, $
         psym=10, color='green'
       if ngood gt 0L then OPLOT, good, spec1d[good]/1d1^power, psym=10
       PLOT, [0,0], /nodata, /noerase, pos=crd[*,4], xst=3, yst=1, xr=[0,ncols-1], $
         yr=yr, xtit='column', xthick=2, ythick=2, charth=2, charsize=1.5, $
         ytit=TEXTOIDL('[extracted flux] / 10^{'+ STRTRIM(STRING(power,f='(i)'),2)+'}')

       PLOT, [0,1], [0,1], /nodata, /noerase, xst=5, yst=5, pos=crd[*,4]
       XYOUTS, 0.03, 0.90, align=0, charth=2, charsize=0.8, $
         'GREEN: without optimal extraction (dashed lines indicate masked pixels)'
       XYOUTS, 0.03, 0.85, align=0, charth=2, charsize=0.8, $
         'BLACK: with optimal extraction (dotted lines indicate masked pixels)'
       XYOUTS, 0.97, 0.9, align=1, charth=2, chars=0.8, $
         'NSIGMA = '+STRTRIM(STRING(nsigma,f='(i)'),2)

       DEVICE, /close
       SET_PLOT, 'x'
       !p.multi = oldpmulti
    endif

 end
 









