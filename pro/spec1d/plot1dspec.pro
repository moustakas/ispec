;+
; NAME:
;	PLOT1DSPEC
;
; PURPOSE:
;	Plot one-dimensional extracted spectra.
;
; CALLING SEQUENCE:
;       plot1dspec, specname, datapath=, nsmooth=, outpath=, psname=, $
;          header=, prefix=, /overplot, /left, /right, /all, /normalize, $
;          /telluric, /rem_continuum, /skyclean, /gzip, /postscript, $
;          _extra=extra
;
; INPUTS:
;	specname - name of the multi-extension spectrum
;
; OPTIONAL INPUTS:
;	datapath - path to the spectrum (default PWD)
;	outpath  - 
;	nsmooth  - boxcar smooth with NSMOOTH width
;	psname   - name of the postscript plot (default to OBJECT) 
;	extra    - keywords for SPLOT
;
; KEYWORD PARAMETERS:
;	postscript - generate postscript output
;	all        - plot all FITS files in the current directory, or
;                    in DATAPATH
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;	SPLOT, SOPLOT, SMOOTH, MRDFITS(), TEXTOIDL, PLOTFAVES, CWD(),
;	DFPSPLOT, DFPSCLOSE
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 19, U of A
;       jm04jan15uofa - better error checking; added TELLURIC and
;                       REM_CONTINUUM keywords
;       jm04mar18uofa - syntactical improvements; use DFPS routines
;                       over PS_ routines
;       jm04nov18uofa - added GZIP keyword
;
; Copyright (C) 2001, 2004, John Moustakas
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

pro plotinfo, specname, header, xtitle=xtitle, ytitle=ytitle, $
  legendbox=legendbox, psname=psname, _extra=extra, all=all

    galaxy_str  = strtrim(sxpar(header,'GALAXY',count=gcount),2)
    object_str  = strtrim(sxpar(header,'OBJECT',count=ocount),2)
    scanlen_str = strtrim(sxpar(header,'SCANLEN',count=scount),2)
    extract_str = strtrim(sxpar(header,'EXTRACT',count=ecount),2)

    if (gcount eq 0L) then title = strupcase(specname) else title = strupcase(galaxy_str)
    if (ocount eq 0L) then object = '' else object = ' ('+strupcase(object_str)+')'
    if (scount eq 0L) then scanlen = '' else scanlen = string(scanlen_str,format='(G0.0)')
    if (ecount eq 0L) then extract = '' else extract = string(extract_str,format='(G0.0)')

    if (scount ne 0L) and (ecount ne 0L) then $
      aperture = extract+'" x '+scanlen+'" Aperture' else $
      aperture = ''

;   legendbox = [specname,title]
;   legendbox = title+object
;   legendbox = [title+object,aperture+object]

;   fluxcor = sxpar(header,'FLUXCOR')
;   if fluxcor eq 0L then ytitle = 'Counts' else $
;     ytitle = textoidl('f_{\lambda} (erg s^{-1} cm^{-2} \AA^{-1})')

    psname = repstr(repstr(repstr(specname,'.fits',''),'.ms',''),'.gz','')+'.ps'

return
end

pro plot1dspec, specname, datapath=datapath, nsmooth=nsmooth, outpath=outpath, $
  psname=psname, header=header, prefix=prefix, overplot=overplot, left=left, right=right, $
  all=all, normalize=normalize, telluric=telluric, rem_continuum=rem_continuum, $
  skyclean=skyclean, gzip=gzip, postscript=postscript, _extra=extra, flagspec=flagspec, $
  scale=scale, yrange=yrange

    nspec = n_elements(specname)
    if (nspec eq 0L) and (n_elements(all) eq 0L) then begin
       print, 'Syntax - plot1dspec, specname, datapath=, nsmooth=, outpath=, $'
       print, '   psname=, header=, prefix=, /overplot, /left, /right, /all, $'
       print, '   /normalize, /telluric, /rem_continuum, /skyclean, /gzip, $'
       print, '   /postscript, _extra=extra'
       return
    endif

    if not keyword_set(datapath) then datapath = cwd()
    if not keyword_set(outpath) then outpath = datapath

; read all the FITS files in the current directory
    
    if (n_elements(all) ne 0L) then $
      specname = file_basename(file_search(datapath+'*.fits*',count=nspec))
    
    if keyword_set(overplot) and xregistered('splot') then begin

       scube = rd1dspec(specname[0],datapath=datapath)
       soplot, scube.wave, scube.spec, ps=10, color='green'
       icleanup, scube, /main
       
    endif

    skywaves = [5577.339,5889.950,6300.32] ; sky wavelengths

    if n_elements(normalize) eq 1L then scale = 1.0 else if $
      (n_elements(scale) eq 0L) then scale = 1E17

    if keyword_set(flagspec) then openw, lun1, 'remove_candidates.txt', /get_lun
    
    for i = 0L, nspec-1L do begin
    
       i = i > 0L

       scube = rd1dspec(specname[i],/silent,datapath=datapath)
       header = scube.header

       if keyword_set(normalize) then begin
          scube.spec = im_normalize(scube.spec,scube.wave,const=normconst,_extra=extra)
          scube.sigspec = scube.sigspec / normconst
       endif

       
; interpolate over "bad" sky pixels and telluric absorption

       if keyword_set(skyclean) then begin
          scube.spec = iskymask(scube.spec,scube.sigspec,scube.wave,mask=scube.mask)
          skymask = scube.mask ge 1B
       endif
       
       if keyword_set(telluric) then begin
          tmask = (telluric_mask(scube.wave) eq 0B)
          if keyword_set(skyclean) then skymask = (skymask + tmask) ge 1B else skymask = tmask
       endif

       if (keyword_set(skyclean)) or (keyword_set(telluric)) then $
         scube.spec = djs_maskinterp(scube.spec,skymask,scube.wave)
       
; smooth the spectrum

       if keyword_set(nsmooth) then begin
          scube.spec = smooth(scube.spec,nsmooth)
          scube.sigspec = smooth(scube.sigspec,nsmooth)
          scube.sky = smooth(scube.sky,nsmooth)
       endif

; grab plot information

       plotinfo, specname[i], header, xtitle=xtitle, ytitle=ytitle, $
         legendbox=legendbox, psname=psname

       if (n_elements(prefix) eq 0L) then prefix = ''
       psname = prefix+psname

; S/N statistics

       goodsig = where(scube.sigspec gt 0.0,ngoodsig,ncomp=nbadsig)
       if (ngoodsig ne 0L) then begin
       
          snr = scube.spec[goodsig]/scube.sigspec[goodsig]
          snrmean = mean(snr)
          snrmedian = median(snr)
          snrsig = stddev(snr)

          snrstr = 'S/N = '+strtrim(string(snrmedian,format='(F12.1)'),2)
;         snrstr = 'S/N = '+string(snrmedian,format='(F5.1)')+' ('+string(snrmean,format='(F5.1)')+' +/- '+$
;           string(snrsig,format='(F5.1)')+')'

       endif else snrstr = ''
       
;      legendbox = [legendbox,snrstr]
       
; remove the continuum

       if keyword_set(rem_continuum) then begin

          cc = im_fitcontinuum(scube.wave,scube.spec,scube.sigspec,method=2,/mask,nocflux=nocflux,/silent)
          scube.spec = nocflux
          
       endif
       
       if keyword_set(postscript) then begin
          dfpsplot, outpath+psname, /landscape
          postthick = 5.0
       endif else postthick = 2.0
          
       if n_elements(normalize) eq 1L then $
         ytitle = textoidl('Normalized f_{\lambda}') else $
         ytitle = textoidl('f_{\lambda} (10^{-'+strmid(scale,1,/reverse)+'} '+flam_units()+')')
;      ytitle = 'Relative Flux'
       xtitle = 'Wavelength (\AA)'

       goodpix = where((scube.wave lt 5500.0) or (scube.wave gt 5640.0)) ; avoid bad sky pixels
       yrange = scale*[min(scube.spec[goodpix]),1.1*max(scube.spec[goodpix])]
;      if (n_elements(yrange) eq 0L) then yrange = scale*[min(scube.spec[goodpix]),1.1*max(scube.spec[goodpix])]

       if (nspec eq 1L) and (not keyword_set(postscript)) then begin ; SPLOT cannot be blocked
          splot, scube.wave, scale*scube.spec, ps=10, xtitle=xtitle, $
            ytitle=ytitle, charsize=2.0, charthick=postthick, xsty=3, ysty=3, $
            yrange=yrange, xmargin=[12,3], _extra=extra, thick=postthick
;         soplot, scube.wave, scube.sigspec, ps=10, line=2
       endif else begin
          if (i eq 0L) and (not keyword_set(postscript)) then window, 0, xs=850, ys=600;xs=450, ys=300
          djs_plot, scube.wave, scale*scube.spec, ps=10, xtitle=xtitle, $
            ytitle=ytitle, charsize=2.0, charthick=postthick, xsty=3, ysty=3, $
            yrange=yrange, xmargin=[12,3], _extra=extra, xthick=postthick, ythick=postthick, $
            thick=postthick-2
;         djs_oplot, scube.wave, scube.sigspec, ps=10, line=2
       endelse
          
;      legend, legendbox, left=left, right=right, /top, box=0, $
;        charthick=postthick, charsize=1.5

       if keyword_set(postscript) then begin
          dfpsclose
          if keyword_set(gzip) then spawn, ['gzip -f '+outpath+psname], /sh
       endif

       if (nspec gt 1L) and (i le nspec-2L) and (not keyword_set(postscript)) then begin
          prompt:

          print, 'Object '+strn(i,length=3)+' '+specname[i]+' [Options: b,g,o,r,s,q]'
;         print, 'Object '+strn(i,length=3)+' '+legendbox[0]+' [Options: b,g,o,r,s,q]'
          cc = strupcase(get_kbrd(1))
          case strlowcase(strcompress(cc,/remove)) of
             'm': if keyword_set(flagspec) then printf, lun1, specname[i]
             'b': i = (i-2L) ; back
             'g': begin      ; goto 
                number = ''
                read, number, prompt='Goto spectrum number (0-'+strn(nspec-1L)+'): '
                number = 0 > long(number-1L) < (nspec-2L)
                i = number
             end
             'r': i = i-1L ; redraw
             'o': begin    ; overplot another spectrum
                number = ''
                read, number, prompt='Overplot spectrum number (0-'+strn(nspec-1L)+'): '
                number = 0 > long(number) < (nspec-1L)
                j = number                
                print, 'Overplotting '+specname[j]+'.'
                ospec = rd1dspec(specname[j],/silent,datapath=datapath)
                if keyword_set(normalize) then begin
                   ospec.spec = im_normalize(ospec.spec,ospec.wave,const=normconst,_extra=extra)
                   ospec.sigspec = ospec.sigspec / normconst
                endif
;               norm = interpol(scube.spec,scube.wave,5500.0)
;               djs_oplot, ospec.wave, norm*ospec.spec/interpol(ospec.spec,ospec.wave,5500.0), $
;                 ps=10, color='green', thick=0.5
                djs_oplot, ospec.wave, scale*ospec.spec, ps=10, color='green', thick=0.5
                goto, prompt
             end
             's': begin ; overplot the sky spectrum
                djs_oplot, scube.wave, scale*scube.sky, ps=10, color='red', thick=0.5
;               norm = interpol(scube.spec,scube.wave,5500.0)
;               djs_oplot, scube.wave, scale*norm*scube.sky/interpol(scube.sky,scube.wave,5500.0), $
;                 ps=10, color='red', thick=0.5
                goto, prompt
             end
             'q': begin
                if keyword_set(flagspec) then free_lun, lun1
                return          ; quit
             end
             else: 
          endcase
       endif
        
; clean up memory

       icleanup, scube
    
   endfor 

   if keyword_set(flagspec) then free_lun, lun1
   
return
end
