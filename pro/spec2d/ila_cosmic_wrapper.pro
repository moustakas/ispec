;+
; NAME:
;	ILA_COSMIC_WRAPPER
;       
; PURPOSE:
;	Call ILA_COSMIC on an array of images in iSPEC2d format. 
;
; CALLING SEQUENCE:
;       ila_cosmic_wrapper, crlist, datapath=, prefix=, iaxis=, $
;          sigclip=, objlim=, sigfrac=, /gzip, /wfits, _extra=extra
;
; INPUTS:
;	crlist  - FITS list of images to clean
;
; OPTIONAL INPUTS:
;	datapath - path to the data
;	prefix   - word or letter to prepend to the output FITS files
;                  (default 'c')
;       iaxis    - interpolation direction (see DJS_MASKINTERP)
;                  (default 0)
;	sigclip  - LA_COSMIC sigma clipping parameter (scalar or array)
;	objlim   - LA_COSMIC parameter (scalar or array)
;	sigfrac  - LA_COSMIC parameter (scalar or array)
;	extra    - additional inputs to ILA_COSMIC
;
; KEYWORD PARAMETERS:
;       gzip     - write GZIPPED FITS files
;	wfits    - write out each image as an ISPEC format FITS file
;
; OUTPUTS:
;       If WFITS is set then write out cosmic-ray cleaned images.
;
; COMMENTS:
;       Cosmic-ray pixels are assigned a bad pixel mask value
;       according to IMASK_BITS().  Recall that 0 is a good pixel. 
;
;	If SIGCLIP, OBJLIM, or SIGFRAC are passed as arrays, and the
;	number of elements does not equal the number of elements in
;	CRLIST, then each parameter is set to the zeroth value.
;
; EXAMPLE:
;	IDL> ila_cosmic_wrapper, ['im1.fits','im2.fits'], $
;	     sigclip=[4.0,3.0], sigfrac=0.1, /wfits
;
;	IDL> ila_cosmic_wrapper, 'im1.fits', sigclip=3.5, iaxis=1, /wfits
;
;
; PROCEDURES USED:
;	ILA_COSMIC, RD2DSPEC(), ICLEANUP, MWRFITS, CWD(), WRT2DSPEC,
;	SPLOG, SXADDHIST, SXADDPAR, IM_TODAY(), IMASK_BITS(), MEDIAN() 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 December 18, U of A
;       jm03apr10uofa - checked out ISPEC v1.0.0
;       jm03dec07uofa - added GZIP keyword
;       jm03dec25uofa - pass the sky value to ILA_COSMIC
;       jm05jun29uofa - use a format code when writing OBJLIM,
;                       SIGCLIP, and SIGFRAC to the header
;
; Copyright (C) 2001, 2003, 2005, John Moustakas
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

pro ila_cosmic_wrapper, crlist, datapath=datapath, prefix=prefix, iaxis=iaxis, $
  sigclip=sigclip, objlim=objlim, sigfrac=sigfrac, gzip=gzip, wfits=wfits, _extra=extra

    nlist = n_elements(crlist)
    if nlist eq 0L then begin
       print, 'Syntax - ila_cosmic_wrapper, crlist, datapath=, prefix=, iaxis=, $'
       print, '   sigclip=, objlim=, sigfrac=, /gzip, /wfits, _extra=extra'
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(iaxis) eq 0L then iaxis = 0L
    
    outlist = strarr(nlist)

; error checking    
    
    nsigclip = n_elements(sigclip)
    if nsigclip gt 0L then if nsigclip ne nlist then sigclip = sigclip[0]
       
    nobjlim = n_elements(objlim)
    if nobjlim gt 0L then if nobjlim ne nlist then objlim = objlim[0]
    
    nsigfrac = n_elements(sigfrac)
    if nsigfrac gt 0L then if nsigfrac ne nlist then sigfrac = sigfrac[0]

; flag cosmic rays
    
    splog, 'Identifying cosmic rays...'
    for k = 0L, nlist-1L do begin

       scube = rd2dspec(crlist[k],datapath=datapath) ; read the structure cube

       image = scube.image
       sigmap = scube.sigmamap
       sky = scube.sky
       mask = scube.mask
       header = *scube.header

       crimage = image
       
       ncols = scube.naxis1
       nrows = scube.naxis2

       if nsigclip gt 0L then sigc = sigclip[k<(nsigclip-1L)] else sigc = 0L
       if nobjlim gt 0L then objl = objlim[k<(nobjlim-1L)] else objl = 0L
       if nsigfrac gt 0L then sigf = sigfrac[k<(nsigfrac-1L)] else sigf = 0L

; remove the object spectrum (fit along the dispersion axis)

;      imsub = image*0.0
;      for i = 0L, nrows-1L do begin
;
;         fitcoeff = func_fit(findgen(ncols),reform(image[*,i]),3,function_name='flegendre',yfit=yfit)
;         imsub[*,i] = image[*,i] - yfit
;
;      endfor

       ila_cosmic, crimage, outlist=cleanim, masklist=crmask, /zeroindexed, $
         sigclip=sigc, objlim=objl, sigfrac=sigf, skyval=median(sky), _extra=extra
       
; CRMASK 0B for good pixels and 1B for cosmic ray pixels
       
       splog, 'Image '+string(k+1,format='(I0)')+'/'+string(nlist,format='(I0)')+': '+$
         string(long(total(crmask)),format='(I0)')+' bad pixels.'
       
; linearly interpolate over image cosmic ray pixels in the wavelength 
; direction.  assign the appropriate bit number (8) to the bad pixel
; mask

       outimage = djs_maskinterp(image,crmask,iaxis=iaxis)
       outsigmap = sigmap
;      outsigmap = djs_maskinterp(sigmap,crmask,iaxis=0L) ; <-- wrong!

       tempmask = crmask*fix(0)
       wcrmask = where(crmask eq 1B,nwcrmask)
       if nwcrmask ne 0L then tempmask[wcrmask] = imask_bits(/crpix)
;      tempmask[where(crmask eq 1B)] = fix(8)
       outmask = mask + tempmask ; 0 is good

;      outmask = ((mask eq 0B) + crmask) eq 0B ; 1B is good

; update the header with the rejection parameters and the input files

       sxaddpar, header, 'ICOSMIC', crlist[k], ' input to ILA_COSMIC_WRAPPER', before='HISTORY'
       sxaddpar, header, 'OBJLIM', objl, ' cosmic ray rejection parameter', $
         before='HISTORY', format='(F12.2)'
       sxaddpar, header, 'SIGCLIP', sigc, ' cosmic ray rejection parameter', $
         before='HISTORY', format='(F12.2)'
       sxaddpar, header, 'SIGFRAC', sigf, ' cosmic ray rejection parameter', $
         before='HISTORY', format='(F12.2)'
       sxaddhist, "'Cosmic rays cleaned "+im_today()+"'", header ; history note

       if keyword_set(wfits) then begin ; write out

          if keyword_set(prefix) then outname = prefix+crlist[k] else outname = 'c'+crlist[k]
          outname = repstr(outname,'.gz','')

          splog, 'Writing '+datapath+outname+'.'
          wrt2dspec, outname, float(outimage), float(outsigmap), outmask, $
            skyimage=float(scube.sky), header, datapath=datapath, gzip=gzip
          
       endif

       icleanup, scube ; clean up memory
       
    endfor 

return
end
