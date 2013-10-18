;+
; NAME:
;	DARK_CHECK
;
; PURPOSE:
;	Compute the number of dark current electrons.
;
; INPUTS:
;	darklist  - file list of bias-subtracted, overscan-subtracted 
;                   dark frames
;
; OPTIONAL INPUTS:
;	datapath  - path to the data
;	timescale - scale the dark frames to this time scale [s] 
;	
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;	darkmean   - mean dark counts
;	darkmedian - median dark counts
;	darksigma  - standard deviation of the dark counts
;
; COMMENTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;	RD2DSPEC(), ICLEANUP, SXPAR, QZAP, DJS_ITERSTAT, DJS_AVSIGCLIP
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 17, U of A
;
; Copyright (C) 2001, John Moustakas
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

pro dark_check, darklist, datapath=datapath, timescale=timescale, darkmean=darkmean, $
  darkmedian=darkmedian, darksigma=darksigma

    darkcube = rd2dspec(darklist,datapath=datapath)
    darkindx = where(strlowcase(darkcube.type) eq 'dark',ndark)

    if ndark eq 0L then begin
       print, 'Please pass at least one dark frame.'
       icleaup, darkcube
    endif else darkcube = darkcube[darkindx]

    if not keyword_set(timescale) then timescale = 900.0 ; [s]

    exptime = fltarr(ndark)
    for i = 0L, ndark-1L do exptime[i] = float(sxpar(*darkcube[i].header,'darktime'))

; zap cosmic rays

    imsize = size(darkcube[0].image,/dimension)
    naxis1 = imsize[0]
    naxis2 = imsize[1]

    imcube = make_array(naxis1,naxis2,ndark,/float)
    crmask = make_array(naxis1,naxis2,ndark,/byte)
    
    for i = 0L, ndark-1L do begin
       qzap, darkcube[i].image, newim, newmask, nsigma=3.0, $
         boxsize=5, skyfiltsize=10.0, nrings=2.0, fluxratio=0.1
       imcube[*,*,i] = newim
       crmask[*,*,i] = newmask
    endfor

    crmask = crmask eq 0B ; 1B is good
    
; convert the images to electrons/second and scale to TIMESCALE

    gain = sxpar(*darkcube[0].header,'GAIN')

    scalecube = imcube*0.0
    for i = 0L, ndark-1L do scalecube[*,*,i] = imcube[*,*,i]*gain*timescale/exptime[i]

    darkmaster = djs_avsigclip(scalecube,3,sigrej=2.0,maxiter=10.0,inmask=crmask,outmask=outmask)
    darkmask = byte(total(outmask,3))
    
; compute statistics

    djs_iterstat, darkmaster, sigrej=2.0, mean=darkmean, median=darkmedian, sigma=darksigma, mask=darkmask

    icleanup, darkcube

return
end
