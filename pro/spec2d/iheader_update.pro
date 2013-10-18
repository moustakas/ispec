;+
; NAME:
;	IHEADER_UPDATE()
;
; PURPOSE:
;	Change and update the data headers with the following:
;	airmass at mid-exposure, zenith distance, Julian date,
;	parallactic angle, UT time at mid-exposure, sidereal time, and
;	[optionally] the detector gain and read noise.  Will also
;	precess the (RA,DEC) coordinates to J2000.0 epoch, and add an
;	EQUINOX keyword.
;
; CALLING SEQUENCE:
;       header = iheader_update(oldheader,gain=,rdnoise=,pscale=)
;
; INPUTS:
;	oldheader - header to update
;
; OPTIONAL INPUTS:
;       pscale    - spatial plate scale [arcsec/pixel]
;	gain      - detector gain [e/ADU]
;	rdnoise   - detector read noise [e]
;
; OUTPUTS:
;	header    - updated header
;
; COMMENTS:
;	Computes the effective airmass according to the IRAF
;	SETAIRMASS formulae.
;
; PROCEDURES USED:
;	SXPAR(), PRECESS, IM_HMS2DEC(), SXADDPAR, OBSERVATORY, SXADDHIST,
;	SXDELPAR , CT2LST, ISPEC_VERSION(), DJS_MEAN()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 18, U of A
;       jm02nov06uofa - compute LST from the UT time rather than grab
;                       it from the header
;       jm03apr10uofa - checked out ISPEC v1.0.0
;       jm05jun17uofa - handle both EPOCH and EQUINOX keywords; added
;                       MJD-OBS keyword
;       jm05jun30uofa - this routine now removes *ALL* IRAF header
;                       keywords 
;
; Copyright (C) 2001-2003, 2005, John Moustakas
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

function iheader_update, oldheader, gain=gain, rdnoise=rdnoise, pscale=pscale

    header = oldheader
    header = header[where(strcompress(oldheader,/remove) ne '')] ; take out blanks

    sxaddhist, "'ISPEC "+ispec_version()+"'", header
    sxaddhist, "'Header keywords updated "+im_today()+"'", header

    if (n_elements(gain) ne 0L) then sxaddpar, header, 'GAIN', gain
    if (n_elements(rdnoise) ne 0L) then sxaddpar, header, 'RDNOISE', rdnoise

    if (n_elements(pscale) eq 0L) then begin
       splog, 'WARNING: No plate scale information given.' 
       pscale = 1.0 ; plate scale ["/pixel]
    endif

    ra = sxpar(header,'RA',count=count) 
    if (count ne 0L) then begin ; only update objects

       epoch = sxpar(header,'EPOCH',count=nepoch)
       equinox = sxpar(header,'EQUINOX',count=nequinox)

; precess the coordinates to the 2000.0 equinox and update the header
; with the new RA, DEC, and EPOCH

       old_ra = 15.0*im_hms2dec(sxpar(header,'RA')) ; [degrees]
       old_dec = im_hms2dec(sxpar(header,'DEC'))    ; [degrees]

       new_epoch = 2000.0
       new_ra = old_ra
       new_dec = old_dec

       precess, new_ra, new_dec, epoch, new_epoch

       sxaddpar, header, 'RA', strjoin(strsplit(im_dec2hms(new_ra/15.0D),' ',/extract),':'), $
         ' right ascension [HMS]'
       sxaddpar, header, 'DEC', strjoin(strsplit(im_dec2hms(new_dec),' ',/extract),':'), $
         ' declination [DMS]', after='RA'
       sxaddpar, header, 'EPOCH', new_epoch, ' coordinate epoch', format='(F6.1)', after='DEC'
       sxaddpar, header, 'EQUINOX', new_epoch, ' coordinate epoch', format='(F6.1)', after='EPOCH'

; UT time at the middle of the exposure
       
       exptime = sxpar(header,'EXPTIME') ; [s]
       sxaddpar, header, 'EXPTIME', float(exptime), ' exposure time [seconds]', format='(F12.3)'

       ut = sxpar(header,'UT')
       sxaddpar, header, 'UT', strjoin(strsplit(im_dec2hms(im_hms2dec(ut)),' ',/extract),':'), $
         ' universal time'
       
       ut_decimal = im_hms2dec(ut)
       utmiddle = strjoin(strsplit(im_dec2hms(ut_decimal+0.5*exptime/3600.0D),' ',/extract),':')
       sxaddpar, header, 'UTMIDDLE', utmiddle, ' universal time at mid-exposure ', after='UT'

; add the Julian date to the header

       date = sxpar(header,'DATE-OBS',count=ndate)
       sxaddpar, header, 'DATE-OBS', date, ' date of observation (yyyy-mm-dd)'
       
       date = strsplit(date,'[-,T]',/extract)

       time = strsplit(utmiddle,':',/extract)
       jd = julday(date[1],date[2],date[0],time[0],time[1],time[2])

       sxaddpar, header, 'JD', jd, ' Julian Date at mid-exposure', after='UTMIDDLE'
       sxaddpar, header, 'MJD-OBS', jd-2400000.5D, ' Modified Julian Date at mid-exposure', after='JD'

; compute the LST time and the airmass at the middle of the exposure
; and update the header; in order to do this calculation we have to
; precess the RA and DEC of the objects to the year of the
; observation 

       obsname = strcompress(sxpar(header,'OBSERVAT'),/remove)
       observatory, obsname, obs

       today_ra = new_ra
       today_dec = new_dec
       
       precess, today_ra, today_dec, new_epoch, float(date[0]) ;, /print
       
       ct2lst, lst, (360.0-obs.longitude), junk, jd ; convert to LST
       lst = strjoin(strsplit(string(im_dec2hms(lst)),' ',/extract),':')
       sxaddpar, header, 'ST', lst, ' sidereal time', after='UTMIDDLE'
;      lst = sxpar(header,'ST')
       lstarr = im_hms2dec(lst)+[0.0,0.5,1.0]*exptime/3600.0 ; LST at the beginning, middle, & end of exposure

       ha = lstarr*15.0D - today_ra ; hour angle [degrees]

       scale = 750.0D ; atmospheric scale height [Allen 1973, p. 125, 133]

       cos_zd = sin(obs.latitude*!dtor) * sin(today_dec*!dtor) + $   ; zd = zenith distance
         cos(obs.latitude*!dtor) * cos(today_dec*!dtor) * cos(ha*!dtor)
       x = scale * cos_zd

       airmass = sqrt(x^2.0D + 2.0*scale + 1.0) - x
       airmass_eff = (airmass[0] + 4.0*airmass[1] + airmass[2])/6.0D

; remove all IRAF header keywords!!

       iraf = where(strmatch(header,'*IRAF*',/fold) eq 1B,niraf,comp=notiraf)
       header = header[notiraf]
;      if (sxpar(header,'DATE') ne 0L) then sxdelpar, header, 'DATE' 
;      if (sxpar(header,'IRAF-TLM') ne 0L) then sxdelpar, header, 'IRAF-TLM'

       sxaddpar, header, 'AIRMASS', float(airmass_eff), ' airmass at mid-exposure', $         
         before='HISTORY', format='(F12.4)'

       meanzd = djs_mean(acos(cos_zd)*!radeg)
       sxaddpar, header, 'ZD', float(meanzd), ' zenith distance at mid-exposure [degrees]', $
         after='AIRMASS', format='(F12.2)'

; add the parallactic angle to the header

       top = sin(ha[0]*!dtor)
       bot = tan(obs.latitude*!dtor)*cos(new_dec*!dtor)-sin(new_dec*!dtor)*cos(ha[0]*!dtor)
       pa = atan(top/bot)*!radeg

       if (pa lt 0.0) then pa = 360D + pa
       if (pa gt 180.0) then pa = pa - 180D

       sxaddpar, header, 'PARANGLE', float(pa), ' parallactic angle at mid-exposure [degrees]', $
         after='AIRMASS', format='(F12.2)'
       
    endif

; insert header information for the spatial plate scale    

    naxis2 = sxpar(header,'NAXIS2')
    
    sxaddpar, header, 'CRVAL2', 0.0, ' position at the reference pixel [arcsec]', before='HISTORY'
    sxaddpar, header, 'CRPIX2', naxis2/2L+1L, ' reference pixel number', before='HISTORY' ; indexed starting at 1
    sxaddpar, header, 'CD2_2', pscale, ' plate scale [arcsec/pixel]', before='HISTORY';, format='(F19.16)'
    sxaddpar, header, 'CDELT2', pscale, ' plate scale [arcsec/pixel]', before='HISTORY'
    sxaddpar, header, 'CTYPE2', 'LINEAR', ' projection type', before='HISTORY'
    
; if a SCAN header keyword exists, calculate the effective exposure
; time (for drift-scanning data)

    scan = sxpar(header,'SCANLEN',count=scancount)
    if (scancount ne 0L) then begin

       aperture = sxpar(header,'APERTURE',count=naperture) ; [arcsec]
       exptime = sxpar(header,'EXPTIME',count=nexptime)    ; [s]

       if (float(scan) gt 0.0) and (naperture ne 0L) and (nexptime ne 0L) then begin
          
          exptime_eff = exptime * aperture / float(scan)

          sxaddpar, header, 'EXPTIME', float(exptime_eff), ' effective exposure time [seconds]', $
            format='(F12.3)'
          sxaddpar, header, 'OEXPTIME', float(exptime), ' original exposure time [seconds]', $
            format='(F12.3)', after='EXPTIME'
;         sxaddhist, "'Original exposure time "+strtrim(string(exptime),2)+" seconds'", header
          sxaddhist, "'Effective exposure time calculated "+im_today()+"'", header

       endif
          
; IHEADER_KEYWORDS makes scanlen a string.  convert it to floating
; point.  do the same for the position angle, POSANGLE, below.
       
       sxaddpar, header, 'SCANLEN', float(scan), ' scan length [arcsec]', format='(F12.2)'
       if (naperture ne 0L) then sxaddpar, header, 'APERTURE', float(aperture), $
         ' slit width [arcsec]', format='(F12.2)'
       
    endif

    posangle = sxpar(header,'POSANGLE',count=pacount)
    if (pacount ne 0L) then sxaddpar, header, 'POSANGLE', float(posangle), $
      ' slit position angle [degrees]', format='(F12.1)'

return, header
end
