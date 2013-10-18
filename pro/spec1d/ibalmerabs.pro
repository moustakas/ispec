;+
; NAME:
;       IBALMERABS()
;
; PURPOSE:
;       Measure the Balmer absorption line equivalent widths. 
;
; INPUTS:
;       wave     - *observed* wavelength vector [NPIX]
;       flux     - *observed* spectrum in erg/s/cm2/Angstrom [NPIX] 
;
; OPTIONAL INPUTS:
;       ferr        - *oberved* error spectrum corresponding to FLUX
;                     [NPIX]  
;       balmerwaves - *rest* wavelengths of the Balmer absorption
;                     lines to measure (default is all of them)
;       balmersigma - *rest* line-width of the Balmer absorption
;                     lines; measure the equivalent width in the
;                     interval BALMERWAVEs +/- GSIGMA*BALMERSIGMA
;                     [Angstrom] (default 5.0)
;       gsigma      - see BALMERSIGMA (default 5.0)
;       z           - object redshift; if non-zero, then everything is
;                     shifted into the restframe; all quantities are
;                     returned in the rest frame
;       extra       - extra keywords for IABSLINEEW()
;
; KEYWORD PARAMETERS:
;       return_balmerabs - simply returns the initialized BALMERABS
;                          data structure without attempting to make
;                          any measurements; used to initialize this
;                          structure in IFITSPEC
;       debug            - generate a plot useful for debugging
;       postscript       - keyword for IABSLINEEW()
;
; OUTPUTS:
;       balmerabs - data structure with the results (all measurements
;                   in the rest frame) 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       If Z is non-zero then be sure to specify BALMERSIGMA.   
;
; PROCEDURES USED:
;       BALMER_DATA(), IABSLINEEW()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 November 26, U of A
;       jm03dec01uofa - added POSTSCRIPT keyword
;       jm04apr06uofa - added Z optional input
;       jm04sep13uofa - added RETURN_BALMERABS keyword
;       jm05jul25uofa - output changed to floating-point precision 
;       jm05jul28uofa - call IABSLINEEW() with ALLOW_PARTIAL=1 
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

function ibalmerabs, wave, flux, ferr=ferr, balmerwaves=balmerwaves, $
  balmersigma=balmersigma, gsigma=gsigma, z=z, return_balmerabs=return_balmerabs, $
  debug=debug, postscript=postscript, _extra=extra

    npix = n_elements(wave)
    nflux = n_elements(flux)

    if (npix eq 0L) or (nflux eq 0L) then begin
       doc_library, 'ibalmerabs'
       return, -1L
    endif

    if n_elements(ferr) eq 0L then ferr = flux*0.0+1.0
    nferr = n_elements(ferr)
    
    if (npix ne nflux) or (npix ne nferr) then begin
       print, 'Dimensions of WAVE, FLUX, and FERR do not agree.'
       return, -1L
    endif

    if n_elements(gsigma) eq 0L then gsigma = 3.0 ; 3-sigma
    if n_elements(z) eq 0L then z = 0.0D

; initialize the Balmer line data and then crop to include a possible
; subset of all the supported Balmer lines
    
    bdata = balmer_data(nbdata=nbdata)

    if n_elements(balmerwaves) eq 0L then balmerwaves = bdata.wave
    nbalmer = n_elements(balmerwaves)

    if n_elements(balmersigma) eq 0L then balmersigma = replicate(5.0,nbalmer)

    if (z gt 0.0) then begin

       restwave = wave/(1.0+z)
       restflux = flux*(1.0+z)
       restferr = ferr*(1.0+z)
       
    endif else begin

       restwave = wave
       restflux = flux
       restferr = ferr
       
    endelse
       
; initialize the output data structure
    
    balmerabs = {$
      babs_line:            '', $
      babs_wave:           0.0, $
      babs:                0.0, $
      babs_err:           -2.0, $
      babs_ew:             0.0, $
      babs_ew_err:        -2.0, $
      babs_continuum:      0.0, $ ; continuum at line center
      babs_continuum_err: -2.0, $ ; continuum error at line center
      babs_lline:          0.0, $ ; lower bound of the line
      babs_uline:          0.0}   ; upper bound of the line
    balmerabs = replicate(balmerabs,nbdata)
    balmerabs.babs_line = bdata.line
    balmerabs.babs_wave = bdata.wave

    if keyword_set(return_balmerabs) then return, balmerabs

    for i = 0L, nbalmer-1L do begin

       match = where(fix(balmerwaves[i]) eq fix(bdata.wave),nmatch)
       if (nmatch ne 1L) then message, 'This should never happen.'

; define the width of each line and then compute the equivalent widths

       lline = bdata[match].wave - gsigma*balmersigma[i]
       uline = bdata[match].wave + gsigma*balmersigma[i]

       absline = iabslineew(restwave,restflux,bdata[match].wave,ferr=restferr,$
         llimit=bdata[match].llimit,lwidth=bdata[match].lwidth,$
         ulimit=bdata[match].ulimit,uwidth=bdata[match].uwidth,$
         lline=lline,uline=uline,ncoeff=2,wplot=wplot,label=bdata[match].label,$
         _extra=extra,absplot=absplot,/allow_partial,debug=debug,$
         postscript=postscript,keystroke=debug)

; fill the output data structure and return    

       balmerabs[match].babs = absline.lineflux
       balmerabs[match].babs_err = absline.lineflux_err
       balmerabs[match].babs_ew = absline.lineew
       balmerabs[match].babs_ew_err = absline.lineew_err
       balmerabs[match].babs_continuum = absline.linec
       balmerabs[match].babs_continuum_err = absline.linec_err
       balmerabs[match].babs_lline = absline.lline
       balmerabs[match].babs_uline = absline.uline

    endfor

return, balmerabs
end    
