;+
; NAME:
;       WRITE_BINNED_STANDARDS
;
; PURPOSE:
;       Sum the CALSPEC and SPEC50CAL standard-star spectra in
;       BINSIZE-Angstrom bandpasses.  
;
; CALLING SEQUENCE:
;       write_binned_standards, binsize=, /doplot, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       binsize - bin size (default 10 Angstroms)
;
; KEYWORD PARAMETERS:
;       doplot - visualize the binned spectra
;       write  - write the resultant binned spectra to the default
;                path 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       REPSTR(), MRDFITS(), SPLOG, FLAM_UNITS()
;
; COMMENTS:
;       Bandpass wavelengths are defined at the *middle* of the
;       wavelength interval given by BINSIZE.
;
;       The maximum wavelength is fixed at 1.2 microns, even though
;       the calspec data extend much redder.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Apr 14, U of A - written, based on previous
;          versions
;       jm05jun01uofa - better formatting code on output; use 48.59 as 
;                       the AB magnitude zero-point; better
;                       speed-of-light
;       jm07jun20nyu  - Calspec database updated; AB magnitude
;                       zero-point updated to 48.6
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

pro write_binned_standards, binsize=binsize, doplot=doplot, write=write

    if (n_elements(binsize) eq 0L) then binsize = 10.0

    rootpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/data')
    pathlist = rootpath+['spec50cal','calspec']+'/SPECTRA'
    rootlist = ['spec50cal','calspec']
    npath = n_elements(pathlist)

    light = 2.99792458D18       ; speed of light [A/s]

    if keyword_set(write) then doplot = 0L

    for ipath = 0L, npath-1L do begin

       if keyword_set(write) then dfpsplot, rootpath+rootlist[ipath]+'.ps', /landscape

       path = pathlist[ipath]
       root = rootlist[ipath]

       pushd, path

       if strmatch(rootlist[ipath],'*calspec*') then begin
          flist = file_search('*.fits',count=fcount)
          outlist = root+'_'+repstr(flist,'.fits','.dat')
       endif else begin
          flist = file_search('*.tab',count=fcount)
          outlist = root+'_'+repstr(flist,'.tab','.dat')
       endelse
          
       finerpix = 10.0          ; oversampling factor
       for i = 0L, fcount-1L do begin

          if strmatch(rootlist[ipath],'*calspec*') then begin

             spec = mrdfits(flist[i],1,h,/silent)
             stdwave = spec.wavelength ; [Angstrom]
             stdflux = spec.flux       ; [erg/s/cm^2/A]

             minwave = min(stdwave) > 1000.0
             maxwave = max(stdwave) < 12000.0 ; NOTE!

          endif else begin

             readcol, flist[i], stdwave, abmag, format='F,F', /silent
             stdflux = 10.0^(-0.4*(abmag+48.6))*light/stdwave/stdwave ; [erg/s/cm^2/A]
             
             minwave = min(stdwave) > 3300.0
             maxwave = max(stdwave) < 12000.0

          endelse
          
;         print, flist[i], minmax(stdwave), format='(A30,2x,F10.1,2x,F10.1)'
          
; compute the bandpasses

          bandstart = minwave + (100 - (minwave mod 100))
          bandend = maxwave - (maxwave mod 100)

          nbins = long((bandend-bandstart)/binsize)
          bandwave = findgen(nbins)*binsize+bandstart ;+binsize/2.0

          starflux = fltarr(nbins)
;         print, flist[i], minmax(stdwave), minmax(bandwave), format='(A30,2x,2F10.1,2x,2F10.1)'

          for j = 0L, nbins-1L do begin
             
             wavearray = bandstart + j*binsize - binsize/2.0 + findgen(binsize*finerpix)/finerpix + 1.0/finerpix/2.0
             binflux = interpolate(stdflux,findex(stdwave,wavearray))
;            binflux = interpol(stdflux,stdwave,wavearray)
             starflux[j] = total(binflux)/binsize/finerpix

          endfor

          starmag = -2.5 * alog10(bandwave * bandwave * starflux / light) - 48.6 ; AB magnitude

          if keyword_set(write) then begin
             
             splog, 'Writing '+outlist[i]+'.'
             openw, lun, '../'+outlist[i], /get_lun
             for k = 0L, nbins-1L do printf, lun, bandwave[k], starmag[k], binsize, $
               format='(F7.1,3x,F6.3,3x,F4.1)'
             free_lun, lun

          endif

          if keyword_set(doplot) or keyword_set(write) then begin
             
             plot, bandwave, starflux, ps=4, xsty=3, ysty=3, xtitle='Wavelength [\AA]', $
               ytitle='Flux ('+flam_units()+')', title=flist[i], charsize=2.0, charthick=2.0, $
               xthick=2.0, ythick=2.0, xmargin=[12,3], xrange=[3000,7000]
             oplot, stdwave, stdflux
             if keyword_set(doplot) then cc = get_kbrd(1)

          endif
          
       endfor       

       popd

       if keyword_set(write) then dfpsclose
       
    endfor 
       
return
end
