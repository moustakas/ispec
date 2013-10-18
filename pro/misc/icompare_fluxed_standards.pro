;+
; NAME:
;       ICOMPARE_FLUXED_STANDARDS
;
; PURPOSE:
;       Compare flux-calibrated standard stars to the published
;       standard spectral energy distributions.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;       searchrad    - standard star search radius [arcsec] (if
;                      coordinates are poor)
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Nov 15, U of A - written
;       jm05jan17uofa - generate a single-page postscript file 
;-

function starsearch, ra, dec, stdpath=stdpath, stdfile=stdfile, star=star, _extra=extra

    rootpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards')
    readcol, rootpath+'SUPPORTED', path, format='A', /silent
;   path = ['calspec','spec50cal','oke1990','ctionewcal','bessell1999']
    npath = n_elements(path)
    indx = lindgen(npath)
    
    w = where(strmatch(path,'*'+stdpath+'*',/fold_case) eq 1B)
    remove, w, indx

    path = path[indx]
    npath = n_elements(path)

; search in each standard star path

    for k = 0L, npath-1L do begin

       starinfo = mrdfits(rootpath+path[k]+'/'+'table_info.fits',1,/silent)
       
       ntot = im_djs_angle_match(ra,dec,starinfo.ra,starinfo.dec,mindx=mindx,$
         mdist=mdist,_extra=extra)

       if (ntot eq 1L) then begin
          print, 'found in '+path[k]+'.'
          stdfile = path[k]+'/'+strn(starinfo[mindx].file)
          star = strcompress(strupcase(starinfo[mindx].star),/remove)
          return, mindx
       endif
       
    endfor
    
return, mindx
end

pro icompare_fluxed_standards, starlist, datapath=datapath, stdpath=stdpath, $
  searchrad=searchrad, suffix=suffix, psname=psname, postscript=postscript

    nstar = n_elements(starlist)
    
    if (nstar eq 0L) then begin
       print, 'Syntax - icompare_fluxed_standards'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(stdpath) eq 0L) then stdpath = 'calspec'
    if (n_elements(searchrad) eq 0L) then searchrad = 200.0
    if (n_elements(suffix) eq 0L) then suffix = '' else suffix = '_'+suffix
    if (n_elements(psname) eq 0L) then psname = 'qaplot_fluxed_standards'+suffix+'.ps'

; forage relevant header information
    
    info = iforage(datapath+starlist)

; initialize some standard-star parameters    
    
    rootpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards')
    splog, 'Stellar template path '+rootpath+stdpath+'.'
    if (file_test(rootpath+stdpath,/directory) eq 0L) then begin
       splog, 'Input STDPATH does not exist . . . setting to calspec.'
       stdpath = 'calspec'
    endif
    starinfo = mrdfits(rootpath+stdpath+'/table_info.fits',1,/silent)

    if keyword_set(postscript) then begin

       dfpsplot, datapath+psname, /color, /landscape
       postthick = 5.0
       
    endif else begin

       window, xsize=750, ysize=500
       postthick = 2.0

    endelse

    pagemaker, nx=2<nstar, ny=round(nstar/2.0), position=position, /normal, $
      xspace=1.0, xmargin=[1.0,0.1], ymargin=[0.2,0.9], yspace=1.0

; loop on each object    
    
    for istar = 0L, nstar-1L do begin

; read the one-dimensional spectrum 

       object = strtrim(info[istar].object,2)
       
       print, format='("Reading ",A0," . . . ",$)', info[istar].file
       spec = rd1dspec(starlist[istar],datapath=datapath,/silent)
       wave = spec.wave
       flux = spec.spec

       ra = 15.0*hms2dec(info[istar].ra)
       dec = hms2dec(info[istar].dec)
       
; search for the current standard star and read it.  search radius is
; SEARCHRAD

       ntot = im_djs_angle_match(ra,dec,starinfo.ra,starinfo.dec,$
         dtheta=searchrad/3600.0,units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)

       case long(ntot) of
          0L: begin
             print, format='(" . . . no standard star match in ", A0, " . . . ",$)', stdpath
             mindx = starsearch(ra,dec,dtheta=searchrad/3600.0,$
               units='degrees',stdpath=stdpath,stdfile=stdfile,star=star)
             if (mindx[0] eq -1L) then begin
                splog, 'Star '+info[istar].object+' '+info[istar].file+' not found.'
                return
             endif
          end
          1L: begin
             print, ' . . . found '+strn(starinfo[mindx].star)+'.'
             stdfile = stdpath+'/'+strn(starinfo[mindx].file)
          end
          2L: begin
             splog, 'Multiple stars found for '+info[istar].object+' '+info[istar].file+'.'
             return
          end
       endcase

; all the standard star data files have the same format (see
; UPDATE_STANDARDS_DATABASE)

       readcol, rootpath+stdfile, bandwave, stdmag, binsize, /silent
       stdflux = 10D^(-0.4*(stdmag+48.59D)) * 2.99793D18 / bandwave / bandwave ; [erg/s/cm2/A]

       if (istar gt 0L) then noerase = 1L

       scale = 1E13
       djs_plot, wave, scale*flux, ps=10, xsty=3, ysty=3, xrange=xrange, $
         yrange=yrange, charsize=1.3, charthick=postthick, thick=2, $
         xtitle='Wavelength [\AA]', ytitle='f_{\lambda} [10^{-13} '+flam_units()+']', $
         xthick=postthick, ythick=postthick, position=position[*,istar], $
         noerase=noerase
       djs_oplot, bandwave, scale*stdflux, ps=4, color='red', syms=2.0

       legend, object, /right, /top, charsize=1.3, charthick=postthick, box=0

       if (not keyword_set(postscript)) and (nstar gt 1L) then cc = get_kbrd(1)
       
    endfor

    title = 'Fluxed Standard Stars - '+repstr(suffix,'_','')

    if n_elements(title) ne 0L then $
      xyouts, position[2,0]+(position[0,1]-position[2,0])/2.0, charthick=postthick, $
      position[3,0]*1.03, title, /normal, charsize=1.3, align=0.5
    
    if keyword_set(postscript) then dfpsclose
    
return
end
    
