;+
; NAME:
;   IRESTORE_SPECLINEFIT()
;
; PURPOSE:
;   Reconstruct the best-fitting absorption- and emission-line
;   spectra determined by ISPECLINEFIT().  Also compute additional
;   other basic/useful info like total stellar masses and mass-
;   and luminosity-weighted ages. 
;
; INPUTS: 
;   specdata - SPECDATA data structure from ISPECLINEFIT() [NGALAXY]
;
; OPTIONAL INPUTS: 
;   templatedir - path name to the spectral synthesis models used in
;                 ISPECLINEFIT() (deafult ${ISPEC_DIR}/templates)
;   instr_vdisp - instrumental velocity dispersion (km/s); the default
;                 is to use the value in the SPECDATA structure, but
;                 by setting this parameter to zero, for example, the
;                 continuum model will be at the resolution of the
;                 templates used (e.g., 3 A FWHM for BC03)
;   linename    - only rebuild the emission line(s) specified by
;                 LINENAME (the default is to rebuild all the lines
;                 that were measured); note that this optional input
;                 is stronger than /STRONG_LINES
;
; KEYWORD PARAMETERS: 
;   nocontinuum    - do not restore the continuum spectrum
;   nolines        - do not restore the emission-line spectrum
;   repair_4959    - force the [OIII] 4959 line to be 1/2.984 times
;                    the strength of [OIII] 5007 (if they were fitted) 
;   negative_lines - restore lines with negative amplitudes (default
;                    is to suppress only restore positive lines); note
;                    that one could alternatively use SNRCUT_LINE
;   strong_lines   - only restore a set of predefined strong (mostly
;                    nebular) emission lines
;   snrcut_line    - only restore lines with S/N>SNRCUT_LINE; default 3
;
; OUTPUTS: 
;   model - [NGALAXY] structure containing the best-fitting continuum
;     plus emission-line model in the observed (redshifted) and rest
;     frame, as well as some other info (see code comments) 
;      .wave          - 
;      .flux          - 
;      .restwave      - 
;      .restflux      -
;      .mass_template -
;      .mass_total    - 
;      .age_mass      - 
;      .age_B         - 
;      .age_V         - 
;      .age_R         - 
;
; OPTIONAL OUTPUTS:
;   outwave      - interpolate the reconstructed model onto this
;                  **rest-frame** wavelength array; the default is to
;                  use the template wavelength vector
;   templateinfo - template informational structure
;
; COMMENTS:
;   Uses a common block to only read the model templates once. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 27, NYU - written based on
;     AGES_RESTORE_BC03()  
;   jm08sep04nyu - added LINENAME optional input
;   jm09aug18ucsd - additional error checking
;
; Copyright (C) 2008-2009, John Moustakas
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

function irestore_speclinefit, specdata, templatedir=templatedir, $
  instr_vdisp=instr_vdisp1, outwave=outwave, linename=linename, $
  templateinfo=templateinfo, nocontinuum=nocontinuum, nolines=nolines, $
  repair_4959=repair_4959, negative_lines=negative_lines, $
  strong_lines=strong_lines, snrcut_line=snrcut_line

    common com_irestore_speclinefit, common_templateinfo, templatewave, $
      templateflux

    light = 2.99792458D5 ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    ngalaxy = n_elements(specdata)
    if (ngalaxy eq 0L) then begin
       doc_library, 'irestore_speclinefit'
       return, -1L
    endif

; call this routine recursively    

    if (ngalaxy gt 1L) then begin
       for igal = 0L, ngalaxy-1L do begin
          print, format='("Rebuilding object ",I0,"/",I0,A10,$)', $
            igal, ngalaxy, string(13b)
          model1 = irestore_speclinefit(specdata[igal],templatedir=templatedir,$
            instr_vdisp=instr_vdisp1,outwave=outwave,linename=linename,$
            templateinfo=templateinfo,nocontinuum=nocontinuum,nolines=nolines,$
            repair_4959=repair_4959,negative_lines=negative_lines,$
            strong_lines=strong_lines,snrcut_line=snrcut_line)
          if (igal eq 0L) then model = replicate(model1,ngalaxy)
          model[igal] = model1
       endfor
       return, model
    endif

; read the templates

    if (not keyword_set(nocontinuum)) then begin
       
       if (n_elements(templatedir) eq 0L) then templatedir = getenv('ISPEC_DIR')+'/templates/'

       if (n_elements(common_templateinfo) eq 0L) then begin
          alltemplates = strtrim(specdata.templatefile,2)
          templates = alltemplates[uniq(alltemplates,sort(alltemplates))]
          ntemplate = n_elements(templates)
          splog, 'Reading info on '+string(ntemplate,format='(I0)')+' templates'
          for it = 0L, ntemplate-1L do begin
             if (file_test(templatedir+templates[it],/regular) eq 0L) then begin
                splog, 'Template file '+templatedir+templates[it]+' not found'
             endif else begin
                templateinfo1 = mrdfits(templatedir+templates[it],1,/silent)
                templateflux1 = mrdfits(templatedir+templates[it],2,templatehdr,/silent)
;               templateres = sxpar(templatehdr,'FWHMRES') ; FWHM resolution [Angstrom] 
                templatewave = make_wave(templatehdr) ; same wavelength vector!
                if (it eq 0L) then begin
                   templateinfo = templateinfo1
                   templateflux = templateflux1
                endif else begin
                   templateinfo = [templateinfo,templateinfo1]
                   templateflux = [[[templateflux]],[[templateflux1]]]
                endelse
             endelse
             common_templateinfo = templateinfo
          endfor
          if (n_elements(templateinfo) eq 0L) then return, -1L
       endif else templateinfo = common_templateinfo

    endif
       
    zabs = specdata.z_abs
    zobj = specdata.z_obj ; emission lines are defined relative to Z_OBJ

; specify the wavelength vector

    if (n_elements(outwave) eq 0L) then restwave = templatewave else $
      restwave = outwave ; OUTWAVE must be in the rest frame!
    wave = restwave*(1.0+zabs)
    
; construct the best-fitting continuum model    

    if keyword_set(nocontinuum) then begin

       modelflux = wave*0.0
       modelrestflux = restwave*0.0

    endif else begin
       
       npix = n_elements(restwave)
       bestindx = specdata.template_bestindx
       
       ntemplate = templateinfo[bestindx].ntemplate
       mpfitcoeff = [specdata.continuum_coeff,specdata.continuum_ebv]

       if (n_elements(instr_vdisp1) ne 0L) then instr_vdisp = instr_vdisp1 else $
         instr_vdisp = specdata.instr_vdisp
       vdisp = specdata.continuum_vdisp
       vdisp_tot = sqrt(vdisp^2.0+instr_vdisp^2.0) ; [km/s]

; redshift and convolve (see IBACKFIT) 

       ztemplatewave = templatewave*(1.0+zabs)
       ztemplateflux = dblarr(npix,ntemplate)
       for ii = 0L, ntemplate-1L do begin
          linterp, alog10(ztemplatewave), templateflux[*,ii,bestindx], $
            alog10(wave), ztemplateflux1
          ztemplateflux[*,ii] = k_smooth(alog10(wave),ztemplateflux1,vdisp_tot)
       endfor

       redcurve = k_lambda(restwave,/charlot) ; rest!
       redcurve = rebin(reform(redcurve,npix,1),npix,ntemplate)

       modelflux = ibackmodel(ztemplateflux,mpfitcoeff,redcurve=redcurve)
       modelrestflux = modelflux*(1.0+zabs)

    endelse
       
; construct the emission-line spectrum

    if keyword_set(nolines) then begin

       linerestflux = restwave*0.0
       lineflux = wave*0.0
       
    endif else begin

       if (n_elements(snrcut_line) eq 0) then snrcut_line = 3.0

       if (n_elements(linename) eq 0L) then linename = specdata.linename
       linename = strupcase(strtrim(linename,2))
       nline = n_elements(linename)

       sigmares = fltarr(nline)
       params = fltarr(3,nline)
       
       for iline = 0L, nline-1L do begin
          fluxtrue = tag_exist(specdata,linename[iline],index=fluxindx)
          wavetrue = tag_exist(specdata,linename[iline]+'_WAVE',index=waveindx)
          sigmatrue = tag_exist(specdata,linename[iline]+'_SIGMA',index=sigmaindx)
          instrtrue = tag_exist(specdata,linename[iline]+'_SIGMA_INSTR',index=instrindx)
          ztrue = tag_exist(specdata,linename[iline]+'_LINEZ',index=zindx)

; the instrumental resolution *could* be, optionally, set to zero or
; to the resolution of the templates
;         sigmares[iline] = (templateres*light/(specdata.(waveindx)))/light/alog(10.0)
          sigmares[iline] = (specdata.(instrindx))[0]/light/alog(10.0) ; instrumental resolution (log-lambda)

          linez = specdata.(zindx)
          if (linez[1] ne -2.0) then begin ; line out of range
             dlinez = linez[0]-zobj ; this offset must match ISPECLINEFIT!
          
             params[1,iline] = alog10(specdata.(waveindx)*(1.0+dlinez))                   ; rest wavelength (log-lambda)
             params[2,iline] = (specdata.(sigmaindx))[0]/alog(10.0)/light                 ; line-width (log-lambda)
             params[0,iline] = (specdata.(fluxindx))[0]/alog(10.0)/10.0^params[1,iline]*$ ; peak amplitude (erg/s/cm2/A)
               ((specdata.(fluxindx))[0]/(specdata.(fluxindx))[1] gt snrcut_line)
          endif
       endfor

; apply a variety of optional quality cuts

       if keyword_set(repair_4959) then begin
          oiii4959 = where(linename eq 'OIII_4959',n4959)
          oiii5007 = where(linename eq 'OIII_5007',n5007)
          if (n4959 ne 0L) and (n5007 ne 0L) then $
            params[0,oiii4959] = params[0,oiii5007]/2.984 
       endif

       if keyword_set(strong_lines) then begin
          strong_linename = ['OII_3727','H_GAMMA','H_BETA','OIII_4959',$
            'OIII_5007','NII_6548','H_ALPHA','NII_6584','SII_6716',$
            'SII_6731']
          strong = where_array(strong_linename,linename)
          if (strong[0] ne -1L) then begin
;            print, linename[strong]
             nline = n_elements(strong)
             params = params[*,strong] 
             sigmares = sigmares[strong]
          endif else begin
             params[0,*] = 0.0  ; suppress
          endelse 
       endif

       if (keyword_set(negative_lines) eq 0) then begin
          neg = where(params[0,*] le 0.0,nneg)
          if (nneg ne 0L) then params[0,neg] = 0.0
       endif
       
; finally build the emission-line spectrum       
       junk = check_math()
       linerestflux = imultigauss(alog10(restwave),params,$
         nline=nline,sigmares=sigmares)
       lineflux = linerestflux/(1.0+zabs)
    endelse

    flux = modelflux + lineflux
    restflux = modelrestflux + linerestflux

;   restwave = templatewave
;   wave = restwave*(1.0+zabs)
;
;   if (n_elements(outwave) ne 0L) then begin
;      full_restwave = restwave & full_wave = wave
;      full_restflux = restflux & full_flux = flux
;
;      restwave = outwave
;      wave = restwave*(1.0+zabs)
;
;      linterp, full_restwave, full_restflux, restwave, restflux, missing=0.0
;      linterp, full_wave, full_flux, wave, flux, missing=0.0
;   endif

; debugging code/AGES
;   ss = read_ages_specfit(specdata.galaxy)
;   ww = where(ss[*,0] gt 0.0)
 
;   djs_plot, ss[ww,0], ss[ww,3], ps=10, xsty=3, ysty=3;, xr=[4700,5100]
;   djs_oplot, templatewave, linerestflux, color='red', ps=10, thick=1
;   cc = get_kbrd(1)

;   djs_plot, ss[ww,0], ss[ww,1], ps=10, xsty=3, ysty=3
;   djs_oplot, templatewave, restflux, color='red'
;   djs_oplot, templatewave, modelrestflux, color='cyan'
;   cc = get_kbrd(1)

; finally compute the stellar masses, the mass-weighted age, and the
; BVR luminosity-weighted ages

    if keyword_set(nocontinuum) then begin

       mass_template = -999.0
       mass_total = -999.0
       age_mass = -999.0
       age_B = -999.0
       age_V = -999.0
       age_R = -999.0
       
    endif else begin

       mass_template = specdata.continuum_coeff*4.0D*!dpi*$ ; stellar mass in each template [M_sun]
         (lumdist(specdata.z_abs,H0=70.0,/silent)*3.085678D24)^2.0 
       mass_total = total(mass_template,/double) ; total stellar mass *within the fiber* [M_sun]

       ages = templateinfo.template_age/1D9 ; template ages [yr]
       ml_B = templateinfo.template_ml_B ; B-band M/L
       ml_V = templateinfo.template_ml_V ; V-band M/L
       ml_R = templateinfo.template_ml_R ; R-band M/L
       
       age_mass = tsum(ages,ages*mass_template)/tsum(ages,mass_template)
       age_B = tsum(ages,ages*mass_template/ml_B)/tsum(ages,mass_template/ml_B)
       age_V = tsum(ages,ages*mass_template/ml_V)/tsum(ages,mass_template/ml_V)
       age_R = tsum(ages,ages*mass_template/ml_R)/tsum(ages,mass_template/ml_R)

    endelse
       
; pack everything into a structure and return

    model = {$
      wave:          wave, $                 ; [A]
      flux:          flux, $                 ; [erg/s/cm2/A]
      restwave:      restwave, $             ; [A]
      restflux:      restflux, $             ; [erg/s/cm2/A]
      mass_template: float(mass_template), $ ; [M_sun]
      mass_total:    float(mass_total), $    ; [M_sun]
      age_mass:      float(age_mass), $      ; [Gyr]
      age_B:         float(age_B), $         ; [Gyr]
      age_V:         float(age_V), $         ; [Gyr]
      age_R:         float(age_R)}           ; [Gyr]

return, model
end
