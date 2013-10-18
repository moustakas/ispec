;+
; NAME:
;	ISKYSUBTRACT
;
; PURPOSE:
;       Sky subtract a two-dimensional spectrum.
;
; CALLING SEQUENCE:
;       iskysubtract, image, sigmamap, header, apinfo, skyinfo, $
;          skymethod=, norder_skysub=, skyinput=, /noskysub, $
;          /checkskysub
;
; INPUTS:
;	image       - two-dimensional image
;	sigmamap    - corresponding two-dimensional error map
;                     (modified) 
;	header      - header associated with IMAGE
;       apinfo      - aperture extraction info structure defined in
;                     IEXTRACT
;
; OPTIONAL INPUTS:
;	skymethod     - a scalar integer indicating the type of sky
;                       subtraction desired:
;		          0 - no sky subtraction
;		          1 - mean sky value (with rejection)
;		          2 - median sky value (with rejection)
;		          3 - polynomial fit (see NORDER_SKYSUB)
;       norder_skysub - order of the sky-subtraction polynomial for
;                       SKYMETHOD = 3 (default NORDER_SKYSUB = 1) 
;	skyinput      - input two dimensional sky model (not
;                       supported yet) 
;
; KEYWORD PARAMETERS:
;	noskysub    - do not sky subtract
;       checkskysub - generate a plot of the sky model fit at every
;                     column for verification purposes
;
; OUTPUTS:
;	skyinfo     - output structure with the following fields:
;	  skymethod      - SKYMETHOD (above)
;	  skyfit         - model two dimensional sky spectrum
;         imnosky        - sky-subtracted image (IMAGE - SKYFIT)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       SKYMETHOD = 4 and SKYINPUT are not fully supported and 
;       generalized yet.
; 
;       SIGMAMAP is modified on output to include the error in the sky
;       subtraction unless NOSKYSUB=1.
;
; EXAMPLE:
;
; PROCEDURES USED:
;	DJS_ITERSTAT, FUNC_FIT(), MAKE_WAVE(), DJS_OPLOT, LEGEND 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 December 27, U of A, excised from IEXTRACT
;       jm02oct16uofa - cleaned up documentation and code
;       jm03apr11uofa - checked out ISPEC v1.0.0
;       jm03apr28uofa - added some error checking.  if SKYVALUES had
;                       values of zero then POLY_FIT() was crashing
;       jm03dec07uofa - major bug fix: propagated the error in the sky
;                       subtraction into SIGMAMAP
;
; Copyright (C) 2001-2003, John Moustakas
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

pro iskysubtract, image, sigmamap, header, apinfo, skyinfo, skymethod=skymethod, $
  norder_skysub=norder_skysub, skyinput=skyinput, noskysub=noskysub, $
  checkskysub=checkskysub

    if n_elements(image) eq 0L then begin
       print, 'Syntax - iskysubtract, image, sigmamap, header, apinfo, skyinfo, $ '
       print, '   skymethod=, norder_skysub=, skyinput=, /noskysub, /checkskysub'
       return
    endif
    
    imsize = size(image,/dimension)
    ncols = imsize[0]
    nrows = imsize[1]
    rowaxis = findgen(nrows)
    colaxis = findgen(ncols)

    skyinfo = { $
      skymethod:   -1L, $
      skyfit:       make_array(ncols,nrows,/float), $
      skyfit_sigma: make_array(ncols,nrows,/float), $
      imnosky:      make_array(ncols,nrows,/float) $
      }

    if (n_elements(norder_skysub) eq 0L) then norder_skysub = 1 else $ ; order of the sky fit
      norder_skysub = long(norder_skysub)

    if (n_elements(skyinput)) eq 0L then begin ; no input sky spectrum

       if keyword_set(noskysub) then begin
          
          imnosky = image
          skyinfo.imnosky = imnosky
          splog, 'No sky subtraction.'
          return

       endif else begin
    
          if n_elements(skymethod) eq 0L then skymethod = 3L else skymethod = long(skymethod)
          skyinfo.skymethod = skymethod

          skyfit = image*0.0       ; sky background model
          skyfit_sigma = image*0.0 ; sky background error map

          for i = 0L, ncols-1L do begin

             skyrowslow = lindgen(apinfo.skyap[0])+apinfo.skylow[i]
             skyrowshigh = reverse(apinfo.skyup[i]-lindgen(apinfo.skyap[1]))
             skyrows = [skyrowslow,skyrowshigh]

;            plot, skyrows, image[ncols-1,skyrows], xsty=3, ysty=3, ps=4

             skyvalues = reform(image[i,skyrows])
             skyerrors = reform(sigmamap[i,skyrows])
             zero = where(sigmamap[i,skyrows] ne 0,nzero)
             if nzero eq 0L then skyinvvar = sigmamap[i,skyrows]*0.0 else begin
                skyinvvar = sigmamap[i,skyrows]*0.0
                skyinvvar[zero] = 1.0/reform(sigmamap[i,skyrows[zero]]^2.0D)
             endelse

; reject outliers (4-sigma)
             
             djs_iterstat, skyvalues, sigrej=4.0, mean=skymean, median=skymedian, $
               sigma=skysigma, mask=skymask
             good = where((skymask eq 1B) and (skyvalues ne float(0)),ngood)

; if there are no good sky values (e.g., if they are all zero) then
; set the sky model to zero, otherwise fit it
             
             if ngood eq 0L then begin

                message, 'No good sky values in column '+string(i,format='(I5)')+'.', /info
                skyfit[i,*] = 0.0

             endif else begin
                
                case skymethod of
                   0:           ; no sky subtraction
                   1: begin
                      skyfit[i,*] = skymean
                      skyfit_sigma[i,*] = skysigma
                   end
                   2: begin
                      skyfit[i,*] = skymedian
                      skyfit_sigma[i,*] = skysigma
                   end
                   3: begin
                   
;                    skycoef = func_fit(skyrows[good],skyvalues[good],norder_skysub,$
;                      invvar=skyinvvar[good],function_name='poly')
;                     skyfit[i,*] = poly(rowaxis,skycoef)
;                  
;                    skycoef = func_fit(skyrows[good],skyvalues[good],norder_skysub,$
;                       invvar=skyinvvar[good],function_name='poly')
;                     skyfit[i,*] = polyleg(rowaxis,skycoef)

                      skycoef = poly_fit(skyrows[good],skyvalues[good],norder_skysub,$
                        measure_errors=skyerrors[good],sigma=skycoef_err,/double,$
                        chisq=chisq,covar=covar)
                      skyfit[i,*] = poly(rowaxis,skycoef)
                      skyfit_sigma[i,*] = sqrt(poly(rowaxis^2,skycoef_err^2))

                      if keyword_set(checkskysub) then begin

                         norm = max(skyvalues[good])
                         
                         if !d.window ne 3L then window, 3, xs=400, ys=400
                         ploterror, skyrows[good], skyvalues[good]/norm, skyerrors[good]/norm, $
                           xsty=3, ysty=3, ps=4, charsize=2.0, $
                           charthick=2.0, xtitle='Row Number', ytitle='Normalized Flux', yrange=yrange
                         djs_oplot, rowaxis, skyfit[i,*]/norm, color='green', thick=3.0, line=0
                         djs_oplot, rowaxis, (skyfit[i,*]+skyfit_sigma[i,*])/norm, color='green', thick=3.0, line=2
                         djs_oplot, rowaxis, (skyfit[i,*]-skyfit_sigma[i,*])/norm, color='green', thick=3.0, line=2
                         legend, 'Column '+strn(i,format='(I0)')+'/'+strn(ncols,format='(I0)'), $
                           /left, /top, box=0, charsize=2.0, charthick=2.0
                         cc = get_kbrd(1)

                      endif

                   end

                   else: message, 'Invalid sky method selection!'

                endcase

             endelse

          endfor

          skyinfo.skyfit = skyfit
          skyinfo.skyfit_sigma = skyfit_sigma

       endelse ; close sky-subtraction condition

    endif else begin ; SKYINPUT defined (not supported yet)

       splog, 'SKYINPUT not supported yet.'
       
;       wave = make_wave(header)
;       exptime = sxpar(header,'EXPTIME')
;          
;       splog, 'Two-dimensional sky spectrum input.'
;
;       skyinsize = size(skyinput.skyflux,/dimension)
;       if (skyinsize[0] ne ncols) and (skyinsize[1] ne nrows) then $
;         message, 'Input sky and image dimensions do not match!'
;
;; cross-correlate the middle row sky vectors to match the wavelengths
;
;; resample the wavelength grid
;
;;      skyfit = image*0.0
;;      invvar = skyfit+1.0
;;      for j = 0L, nrows-1L do begin
;;          
;;         combine1fiber, alog10(skyinput.wave), skyinput.skyflux[*,j], invvar[*,j], $
;;           newloglam=alog10(wave), newflux=newflux
;;         skyfit[*,j] = newflux
;;               
;;      endfor
;             
;; normalization for the image
;             
;       get_element, wave, [5475.0,5525.0], xx ; <-- not general!
;       imnorm = djs_mean(image[xx[0]:xx[1],nrows/2L])
;
;; normalization for the input sky spectrum
;
;       get_element, skyinput.wave, [5475.0,5525.0], xx
;       innorm = djs_mean(skyinput.skyflux[xx[0]:xx[1],nrows/2L])
;
;       skyfit = skyfit*imnorm/innorm ; scale
             
    endelse 

    imnosky = image - skyfit                         ; subtract the sky model
    sigmamap = sqrt(sigmamap^2.0 + skyfit_sigma^2.0) ; propagate the error in the sky subtraction
    skyinfo.imnosky = imnosky

return
end
