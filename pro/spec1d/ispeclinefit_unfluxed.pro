;+
; NAME:
;   ISPECLINEFIT_UNFLUXED()
;
; PURPOSE:
;   Model the emission-line spectra of galaxies and QSOs based on
;   *unfluxed* spectra.  
;
; INPUTS:
;   wave    - wavelength array (observed frame, A) [NPIX,NGAL]
;   flux    - corresponding flux vector (erg/s/cm2/A) [NPIX,NGAL]  
;   invvar  - corresponding inverse variance array; masked pixels 
;             should have INVVAR=0 [NPIX,NGAL] 
;   zobj    - redshift [NGAL]
;   specres - FWHM instrumental resolution, which we assume does not
;             vary with wavelength (Angstrom)
;
; OPTIONAL INPUTS:
;   sigmax       - maximum intrinsic line-width for ILINEFIT(); can be
;                  a scalar or an [NGAL] vector (default 500 km/s); we 
;                  include this parameter here to allow SIGMAX to be
;                  different for each object (e.g., galaxies vs QSOs)
;   galaxy       - optional galaxy name for the output structure and QA
;                  plot [NGAL]
;   outpath      - output data path name (default current directory)
;   suffix       - suffix for all the output files (default '')
;   linefile     - emission line file; see READ_LINEPARS() (default
;                  ${ISPEC_DIR}+'etc/elinelist.dat')   
;   extra        - keywords and parameters for IBACKFIT(), ILINEFIT(),
;                  and SPECTRAL_INDICES() 
;
; KEYWORD PARAMETERS:
;   nologfile - do not generate a log file
;   clobber   - overwrite existing files
;   doplot    - plot everything on the screen, suppressing the
;               postscript output 
;   silent    - suppress extensive SPLOGing
;
; OUTPUTS:
;   specdata - line-fitting results (fluxes, EWs, widths, etc.)
;
; OPTIONAL OUTPUTS:
;    specdatafile - name of the output SPECDATA file
;    specfit - fitting results for the last object, useful for
;              debugging ([6,NPIX,NGALAXY] array)
;      specfit[0,*,igal] = rest wavelength
;      specfit[1,*,igal] = rest flux
;      specfit[2,*,igal] = best-fitting emission-line spectrum 
;      specfit[3,*,igal] = smooth (residual) continuum correction 
;      specfit[4,*,igal] = rest inverse variance spectrum
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Sep 01, NYU - based on earlier code and on the
;     latest version of ISPECLINEFIT()
;
; Copyright (C) 2008, John Moustakas
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

function ispeclinefit_unfluxed, wave1, flux1, invvar1, zobj=zobj, $
  specres=specres, sigmax=sigmax1, galaxy=galaxy, outpath=outpath, $
  suffix=suffix, linefile=linefile, specfit=specfit, specdatafile=specdatafile, $
  _extra=extra, nologfile=nologfile, doplot=doplot, clobber=clobber, $
  silent=silent

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    
    if (n_params() ne 3L) then begin
       doc_library, 'ispeclinefit'
       return, -1 
    endif

    dims = size(wave1,/dim)
    if (size(wave1,/n_dim) eq 1L) then ngalaxy = 1L else $
      ngalaxy = dims[1]

    if (total(size(wave1,/dim) ne size(flux1,/dim)) ne 0.0) or $
      (total(size(wave1,/dim) ne size(invvar1,/dim)) ne 0.0) then begin
       splog, 'Dimensions of WAVE, FLUX, and INVVAR must match'
       return, -1
    endif
    
; check for ZOBJ and SPECRES
    
    nzobj = n_elements(zobj)
    case nzobj of
       0L: begin
          splog, 'ZOBJ must be defined'
          return, -1L
       end
       else: begin
          if (nzobj ne ngalaxy) then begin
             splog, 'Dimensions of ZOBJ must match WAVE'
             return, -1L
          endif
       end 
    endcase

    if (n_elements(specres) eq 0L) then begin
       splog, 'SPECRES must be defined'
       return, -1L
    endif

; handle SIGMAX and GALAXY

    nsigmax = n_elements(sigmax1)
    case nsigmax of
       0L: sigmax = replicate(500.0,ngalaxy) ; [km/s]
       1L: sigmax = replicate(sigmax1,ngalaxy)
       else: begin
          if (nsigmax ne ngalaxy) then begin
             splog, 'Dimensions of SIGMAX must match the number of objects'
             return, -1L
          endif
          sigmax = sigmax1
       end  
    endcase

    if (n_elements(galaxy) eq 0L) then begin
       digits = string(ceil(alog10(ngalaxy)),format='(I0)')
       galaxy = 'Galaxy_'+string(findgen(ngalaxy),format='(I'+$
         digits+'.'+digits+')')
    endif
    
; define path and filename defaults; could check whether the output
; files exist
    
    if (n_elements(outpath) eq 0L) then begin
       spawn, 'pwd', outpath
       outpath = outpath[0]+'/'
    endif

    if (n_elements(suffix) eq 0L) then suffix = ''

    mjdstr = get_mjdstr()
    if (suffix ne '') then begin
       logfile = outpath+mjdstr+'_'+suffix+'.log'
       specfitfile = outpath+mjdstr+'_'+suffix+'_specfit.fits'
       specdatafile = outpath+mjdstr+'_'+suffix+'_specdata.fits'
       psfile = outpath+mjdstr+'_'+suffix+'_specfit.ps'
    endif else begin
       logfile = outpath+mjdstr+'.log'
       specfitfile = outpath+mjdstr+'_specfit.fits'
       specdatafile = outpath+mjdstr+'_specdata.fits'
       psfile = outpath+mjdstr+'_specfit.ps'
    endelse

    if file_test(specdatafile,/regular) and $
      (not keyword_set(clobber)) then begin
       splog, 'Output file name '+specdatafile+' exists; try /CLOBBER'
       return, -1L
    endif
    
    if (not keyword_set(nologfile)) then begin
       splog, filename=logfile
       splog, 'Log file '+logfile+' opened '+systime()
    endif

    splog, 'IDL version: ' + string(!version,format='(99(A," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: '+uname[0]
    splog, 'Output datapath is '+outpath

; read the emission line parameters and constraints

    if (n_elements(linefile) eq 0L) then linefile = filepath('',root_dir=$
      getenv('ISPEC_DIR'),subdirectory='etc')+'elinelist.dat'
    splog, 'Reading emission-line parameters '+linefile
    linepars = read_linepars(linefile=linefile)
    if size(linepars,/type) ne 8L then begin
       splog, 'Problem reading '+linefile
       return, -1L
    endif

; open the postscript file    

    if keyword_set(doplot) then begin
       im_window, 0, xratio=0.7, yratio=0.5
    endif else begin
;      splog, 'Opening plot file '+file_basename(psfile)
;      ps_start, psfile, /landscape, font=1, charsize=1.5, default_thickness=4
;      dfpsplot, psfile, /landscape, /color
       im_plotconfig, 1, psfile=psfile, /landscape
    endelse

; loop on each object
    
    stime0 = systime(1)
;   for igal = 0L, 3 do begin
    for igal = 0L, ngalaxy-1 do begin

       splog, format='("Fitting object #",I0,"/",I0)', igal+1, ngalaxy

       wave = wave1[*,igal]*1.0D ; note double!
       flux = flux1[*,igal]*1.0D
       invvar = invvar1[*,igal]*1.0D
       in_zobj = zobj[igal]

       splog, galaxy[igal]+', z = '+strtrim(string(in_zobj,format='(F12.6)'),2)

; model the continuum as a smooth function 

       notzero = where(invvar gt 0)
       medsnr = djs_median(flux[notzero]*sqrt(invvar[notzero])) > 0.0 ; note!

       emask = iemission_mask(wave,z=in_zobj,vdisp=100.0,$
         /sky,/nebular,/telluric,/qso)
       good = where(emask ne 0.0,ngood)
       if (ngood ne 0L) then begin
          smooth1 = medsmooth(flux[good],151)
          linterp, wave[good], smooth1, wave, smooth1
          continuum = smooth(smooth1,51,/edge_truncate)
       endif
       smooth_continuum = continuum*0.0 ; for symmetry with ISPECLINEFIT()
       
; now fit the nebular emission lines; define all the rest-frame
; quantities we will need later

       restwave = wave/(1.0+in_zobj)
       restflux = (flux - smooth_continuum)*(1.0+in_zobj)
       restcontinuum = continuum*(1.0+in_zobj)
       restinvvar = invvar/(1.0+in_zobj)^2.0

       instr_lineres = light*(specres/(1.0+in_zobj))/fwhm2sig/linepars.wave ; [km/s]

       if (not keyword_set(silent)) then splog, 'Measuring the emission-line strengths'
       linefit = ilinefit(restwave,restflux-restcontinuum,restinvvar,$
         linepars,zguess=0.0,instr_lineres=instr_lineres,_extra=extra,$
         sigmax=sigmax[igal],speclinefit=speclinefit,linefit_chi2=linefit_chi2,$
         linefit_niter=linefit_niter,linefit_status=linefit_status,silent=silent)
       linefit.linez = (linefit.linez+in_zobj)*(linefit.linez ne 0.0)

; measure upper limits and EWs; could optionally combine the [OII]
; doublet here

       if (not keyword_set(silent)) then splog, 'Computing upper limits and EWs'
       ulimit = where(linefit.linearea_err gt -2.0,nulimit)
       if (nulimit ne 0L) then begin
          ulinefit = linefit[ulimit]
          ulinefit = iupper_limits(restwave,restflux-speclinefit,ulinefit,$
            snrcut=1.0,glosigma=5.0,ghisigma=20.0,/telluric,debug=debug)
          linefit[ulimit].linecontlevel     = ulinefit.linecontlevel
          linefit[ulimit].linecontlevel_err = ulinefit.linecontlevel_err
          linefit[ulimit].linelimit         = ulinefit.linearea
       endif

       linefit = fill_ew(temporary(linefit))

       if (not keyword_set(silent)) then struct_print, $
         struct_trimtags(linefit,select=$
;        ['*name','*area*','*sigma*','*limit*','*npix*'])
         ['*name','*area*','*sigma*','*linez*','*chi2*'])

; parse and pack everything into a structure
       
       linefit = struct_addtags({linefit_chi2: float(linefit_chi2), $
        linefit_niter: fix(linefit_niter), linefit_status: $
         fix(linefit_status)},parse_ilinefit(linefit))

; pack everything into the SPECDATA structure; also make SPECFIT 

       specfit1 = float([[wave/(1.0+in_zobj)],[flux*(1.0+in_zobj)],$
         [continuum*(1.0+in_zobj)],[speclinefit],$
         [smooth_continuum*(1.0+in_zobj)],[invvar/(1.0+in_zobj)^2.0]])
       if (igal eq 0L) then specfit = specfit1 else $
         specfit = [[[temporary(specfit)]],[[specfit1]]]
       
       if (igal eq 0L) then begin
          linefit1 = im_empty_structure(struct_trimtags(linefit,$
            except=['*BOX*']))
          
          specdata1 = {$
            galaxy:                  '', $
            z_obj:                  0.0, $
            specres:            specres, $
            minwave:                0.0, $
            maxwave:                0.0, $
            z_abs:                  0.0, $
            continuum_snr:       medsnr}
          specdata = struct_addtags(specdata1,linefit1)
          specdata = replicate(specdata,ngalaxy)
       endif

       specdata[igal].galaxy = galaxy[igal]
       specdata[igal].z_obj = in_zobj
       specdata[igal].z_abs = in_zobj
       specdata[igal].minwave = min(wave[where(invvar gt 0)])
       specdata[igal].maxwave = max(wave[where(invvar gt 0)])
       
       junk = specdata[igal] & struct_assign, linefit, $
         junk, /nozero & specdata[igal] = junk

; continuously overwrite SPECDATAFILE in case the code crashes
          
       if (not keyword_set(silent)) then splog, 'Updating '+file_basename(specdatafile)
       mwrfits, reform(specdata), specdatafile, /create 

       if (not keyword_set(silent)) then splog, 'Updating '+file_basename(specfitfile)
       mwrfits, specfit, specfitfile, /create 

; finally build the QA plot

       qaplot_ispeclinefit, specdata[igal], specfit[*,*,igal], $
         postscript=(not keyword_set(doplot)), /unfluxed

    endfor
    
    splog, format='("Total time for ISPECLINEFIT = ",G0," '+$
      'minutes.")', (systime(1)-stime0)/60.0
    splog, /close

; compress SPECDATAFILE, close the postscript file, and return     
    
    spawn, 'gzip -f '+specdatafile, /sh
    spawn, 'gzip -f '+specfitfile, /sh
    im_plotconfig, psfile=psfile, /pdf, /psclose, /landscape

;   ps_end, /pdf, /delete_ps
;   dfpsclose
;   spawn, 'gzip -f '+psfile, /sh
;   spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
;   spawn, '/bin/rm -f '+psfile

return, reform(specdata)
end
