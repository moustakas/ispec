;+
; NAME:
;	GAIN_RDNOISE
;
; PURPOSE:
;	Measure the CCD gain and read noise from two flat field images 
;	and two bias frames.
;
; CALLING SEQUENCE:
;	gain_rdnoise, flatfile1, flatfile2, biasfile1, biasfile2,
;	   xskip=, yave=, gain=, rdnoise=
;
; INPUTS:
;	flatfile1  - File name for flat #1
;	flatfile2  - File name for flat #2
;	biasfile1  - File name for bias #1
;	biasfile2  - File name for bias #2
;
; OPTIONAL INPUTS:
;	xskip - number of columns to ignore at the beginning and end
;               of each image (default to 50) 
;	yave  - number of rows to analyze for each calculation of gain 
;               and read noise (default to 15) 
;
; OUTPUTS:
;	gain   - gain in electrons/ADU 
;	rnoise - read noise in electrons 
;
; COMMENTS:
;	The appropriate formulae can be found in Howell's Handbook of
;	CCD Astronomy.  This routine is taken almost entirely from
;	David Schlegel's SPGAIN routine. 
;
; PROCEDURES USED:
;	READFITS(), DJS_ITERSTAT
;
; MODIFICATION HISTORY:
;	2001 July 23, J. Moustakas, U of A, written
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

pro gain_rdnoise, flatfile1, flatfile2, biasfile1, biasfile2, xskip=xskip, $
        yave=yave, gain=gain, rdnoise=rdnoise

    if (n_params() ne 4L) then begin
       print, 'Syntax - gain_rdnoise, flatfile1, flatfile2, biasfile1, $'
       print, '   biasfile2, xskip=, yave=, gain=, rdnoise='
       return
    endif

    if (n_elements(xskip) eq 0L) then xskip = 50L
    if (n_elements(yave) eq 0L) then yave = 15L

; read in the flat fields and the bias frames

    flatimg1 = readfits(flatfile1,flathead1,/silent)
    flatimg2 = readfits(flatfile2,flathead2,/silent)

    biasimg1 = readfits(biasfile1,biashead1,/silent)
    biasimg2 = readfits(biasfile2,biashead2,/silent)

    dims = size(flatimg1,/dimen)
    nx = dims[0]
    ny = dims[1]

    nyblock = fix(ny/yave)          ; row subsection
    gainarr = fltarr(nx,nyblock)
    rdnoisearr = fltarr(nx,nyblock)
    corrfac = sqrt(yave/(yave-1.0))   ; correction factor for measured sigmas

    for ix = 0, nx-1 do begin

       for iy = 0, nyblock-1 do begin

          y1 = iy * yave
          y2 = y1 + yave - 1
          
          flatsub1 = float(flatimg1[ix,y1:y2])
          flatsub2 = float(flatimg2[ix,y1:y2])
          
; compute statistics for flats

          djs_iterstat, flatsub1, sigrej=sigrej, maxiter=maxiter, mean=flatmean1
          djs_iterstat, flatsub2, sigrej=sigrej, maxiter=maxiter, mean=flatmean2
          djs_iterstat, flatsub2 - flatsub1, sigrej=sigrej, maxiter=maxiter, sigma=flatdifsig

          flatdifsig = corrfac * flatdifsig

; compute statistics for biases

          biassub1 = biasimg1[ix,y1:y2]
          biassub2 = biasimg2[ix,y1:y2]

          djs_iterstat, biassub1, sigrej=sigrej, maxiter=maxiter, mean=biasmean1
          djs_iterstat, biassub2, sigrej=sigrej, maxiter=maxiter, mean=biasmean2
          djs_iterstat, biassub2 - biassub1, sigrej=sigrej, maxiter=maxiter, sigma=biasdifsig

          biasdifsig = corrfac * biasdifsig
          
          gainarr[ix,iy] = (flatmean1 + flatmean2 - biasmean1 - biasmean2) / (flatdifsig^2 - biasdifsig^2)
          rdnoisearr[ix,iy] = gainarr[ix,iy] * biasdifsig / sqrt(2.)
          
       endfor

    endfor

; compute the median gain and read noise

    xstart = xskip
    xend = nx-xskip-1L

    gain = median(gainarr[xstart:xend,*])
    rdnoise = median(rdnoisearr[xstart:xend,*])

return
end
