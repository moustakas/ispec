;+
; NAME:
;   ISPECLINEFIT()
;
; PURPOSE:
;   Model the emission- and absorption-line spectra of galaxies. 
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
;   vdisp        - stellar velocity dispersion; can be a scalar or an 
;                  [NGAL] vector (default 150 km/s) 
;   sigmax       - maximum intrinsic line-width for ILINEFIT(); can be
;                  a scalar or an [NGAL] vector (default 500 km/s); we 
;                  include this parameter here to allow SIGMAX to be
;                  different for each object (e.g., galaxies vs QSOs)
;   vmaxtweak    - maximum update to the input absorption-line
;                  redshift for  IBACKFIT(); can be a scalar or an
;                  [NGAL] vector (default 500 km/s); we  include this
;                  parameter here to allow it to be different for each
;                  object 
;   vlinemaxtweak- maximum shift in the emission-line center for
;                  ILINEFIT(); can be a scalar or an [NGAL] vector
;                  (default 500 km/s); we  include this parameter here
;                  to allow it to be different for each object
;   zsnrcross    - minimum required S/N to improve the input redshift
;                  in IBACKFIT(); can be a scalar or an [NGAL] vector
;                  (default 2.0); include this parameter here to allow
;                  ZSNRCROSS to be different for each object (e.g.,
;                  galaxies vs QSOs) 
;   galaxy       - optional galaxy name for the output structure and QA
;                  plot [NGAL]
;   outpath      - output data path name (default current directory)
;   suffix       - suffix for all the output files (default '')
;   templatefile - template file name(s) (see WRITE_BC03_TEMPLATES);
;                  pass a vector of file names to fit different
;                  metallicities (default ${ISPEC_DIR}/templates/
;                  BC03_Z02_salpeter_templates.fits)
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
;      specfit[2,*,igal] = best-fitting continuum spectrum
;      specfit[3,*,igal] = best-fitting emission-line spectrum 
;      specfit[4,*,igal] = smooth (residual) continuum correction 
;      specfit[5,*,igal] = rest inverse variance spectrum
;
; COMMENTS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Aug 14, NYU - totally rewritten version of an
;     earlier code of the same name
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

function ispeclinefit, wave1, flux1, invvar1, zobj=zobj, specres=specres, $
  vdisp=vdisp1, sigmax=sigmax1, vmaxtweak=vmaxtweak1, vlinemaxtweak=vlinemaxtweak1, $
  zsnrcross=zsnrcross1, galaxy=galaxy, outpath=outpath, suffix=suffix, $
  templatefile=templatefile, linefile=linefile, specfit=specfit, $
  specdatafile=specdatafile, _extra=extra, nologfile=nologfile, doplot=doplot, $
  clobber=clobber, silent=silent

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

; handle VDISP, SIGMAX, VMAXTWEAK, VLINEMAXTWEAK, ZSNRCROSS, and
; GALAXY 

    nvdisp = n_elements(vdisp1)
    case nvdisp of
       0L: vdisp = replicate(150.0,ngalaxy) ; [km/s]
       1L: vdisp = replicate(vdisp1,ngalaxy)
       else: begin
          if (nvdisp ne ngalaxy) then begin
             splog, 'Dimensions of VDISP must match the number of objects'
             return, -1L
          endif
          vdisp = vdisp1
       end 
    endcase

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

    nvmaxtweak = n_elements(vmaxtweak1)
    case nvmaxtweak of
       0L: vmaxtweak = replicate(500.0,ngalaxy) ; [km/s]
       1L: vmaxtweak = replicate(vmaxtweak1,ngalaxy)
       else: begin
          if (nvmaxtweak ne ngalaxy) then begin
             splog, 'Dimensions of VMAXTWEAK must match the number of objects'
             return, -1L
          endif
          vmaxtweak = vmaxtweak1
       end  
    endcase

    nvlinemaxtweak = n_elements(vlinemaxtweak1)
    case nvlinemaxtweak of
       0L: vlinemaxtweak = replicate(500.0,ngalaxy) ; [km/s]
       1L: vlinemaxtweak = replicate(vlinemaxtweak1,ngalaxy)
       else: begin
          if (nvlinemaxtweak ne ngalaxy) then begin
             splog, 'Dimensions of VLINEMAXTWEAK must match the number of objects'
             return, -1L
          endif
          vlinemaxtweak = vlinemaxtweak1
       end  
    endcase

    nzsnrcross = n_elements(zsnrcross1)
    case nzsnrcross of
       0L: zsnrcross = replicate(2.0,ngalaxy)
       1L: zsnrcross = replicate(zsnrcross1,ngalaxy)
       else: begin
          if (nzsnrcross ne ngalaxy) then begin
             splog, 'Dimensions of ZSNRCROSS must match the number of objects'
             return, -1L
          endif
          zsnrcross = zsnrcross1
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

; define and read the templates

    imf = 'salpeter' ; 'chabrier'
    templatedir = filepath('',root_dir=getenv('ISPEC_DIR'),$
      subdirectory='templates')
    if (n_elements(templatefile) eq 0L) then begin
       templatefile = templatedir+'BC03_Z02_'+imf+'_templates.fits'  ; solar metallicity
;      templatefile = templatedir+'BC03_Z004_'+imf+'_templates.fits' ; LMC metallicity 
;      templatefile = templatedir+'BC03_Z05_'+imf+'_templates.fits'  ; twice-solar metallicity 
;      templatefile = templatedir+['BC03_Z004_'+imf+'_templates.fits', $
;        'BC03_Z02_'+imf+'_templates.fits','BC03_Z05_'+imf+'_templates.fits']
    endif
    ntemplatefile = n_elements(templatefile)

    for ii = 0L, ntemplatefile-1L do begin

       if file_test(templatefile[ii],/regular) eq 0L then begin
          splog, 'Template file '+templatefile[ii]+' not found'
          return, -1L
       endif

       splog, 'Reading template file '+templatefile[ii]
       templateinfo1 = mrdfits(templatefile[ii],1,/silent)           ; info structure
       templateflux1 = mrdfits(templatefile[ii],2,templatehead,/silent) ; flux

       if (ii eq 0L) then begin
          templatewave = make_wave(templatehead)      
          templateres = sxpar(templatehead,'FWHMRES') ; FWHM resolution [Angstrom] 
          templateinfo = templateinfo1
          templateflux = templateflux1
       endif else begin
          templateinfo = [[templateinfo],[templateinfo1]]
          templateflux = [[[templateflux]],[[templateflux1]]]
       endelse
    endfor

    templateinfo = reform(templateinfo)
    templateflux = reform(templateflux)

    splog, 'Fitting with '+string(templateinfo[0].ntemplate,$
      format='(I0)')+' templates'

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
       splog, 'Opening plot file '+file_basename(psfile)
       dfpsplot, psfile, /landscape, /color
    endelse

; loop on each object
    
    stime0 = systime(1)
    for igal = 0L, ngalaxy-1L do begin

       splog, format='("Fitting object ",I0,"/",I0)', igal+1, ngalaxy

       wave = wave1[*,igal]*1.0D ; note double!
       flux = flux1[*,igal]*1.0D
       invvar = invvar1[*,igal]*1.0D
       in_zobj = zobj[igal]
       in_vdisp = vdisp[igal]

       splog, galaxy[igal]+', z = '+strtrim(string(in_zobj,format='(F12.6)'),2)+$
         ', vdisp = '+strtrim(string(in_vdisp,format='(F10.1)'),2)+' km/s'

; fit the stellar continuum
       backfit = ibackfit(wave,flux,invvar,templatewave,templateflux,$
         zobj=in_zobj,vdisp=in_vdisp,templateres=templateres,specres=specres,$
         vmaxtweak=vmaxtweak[igal],zsnrcross=zsnrcross[igal],silent=silent,$
         _extra=extra)

; median-smooth the residuals; do this by median-smoothing the
; spectrum free of emission lines and other wacky regions, interpolate
; onto the full wavelength array, and then box-car smooth; it works
; nicely! 

       residuals = flux-backfit.continuum
       smooth_continuum = residuals*0.0

       emask = iemission_mask(wave,z=in_zobj,vdisp=100.0,$
         /sky,/nebular,/telluric,/qso)
       good = where((emask ne 0.0),ngood)
       if (ngood ne 0L) then begin
          smooth1 = medsmooth(residuals[good],151)
          linterp, wave[good], smooth1, wave, smooth1
          smooth_continuum = smooth(smooth1,51,/edge_truncate)
       endif
       
;      smooth_continuum = smooth(medsmooth(residuals,$
;        250),100,/nan,/edge_truncate)
;      med = im_medxbin(wave,residuals,200,weight=emask*1.0)
       
; now fit the nebular emission lines; define all the rest-frame
; quantities we will need later; shift the emission lines to IN_ZOBJ,
; in case Z_ABS was spurious

       z_abs = backfit.z_abs
       restwave = wave/(1.0+z_abs)
       restflux = (flux - smooth_continuum)*(1.0+z_abs)
       restcontinuum = backfit.continuum*(1.0+z_abs)
       restinvvar = invvar/(1.0+z_abs)^2.0

       instr_lineres = light*(specres/(1.0+z_abs))/fwhm2sig/linepars.wave ; [km/s]

       if (not keyword_set(silent)) then splog, 'Measuring the emission-line strengths'
       linefit = ilinefit(restwave,restflux-restcontinuum,restinvvar,$
         linepars,zguess=z_abs-in_zobj,instr_lineres=instr_lineres,_extra=extra,$
         sigmax=sigmax[igal],vlinemaxtweak=vlinemaxtweak[igal],speclinefit=speclinefit,$
         linefit_chi2=linefit_chi2,linefit_niter=linefit_niter,$
         linefit_status=linefit_status,silent=silent)
;      linefit.linez = (linefit.linez+z_abs)*(linefit.linez ne 0.0)
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

; measure spectral indices; when subtracting the model emission-line
; spectrum don't allow (spurious) negative lines

       if (not keyword_set(silent)) then splog, 'Measuring spectral indices'
       good = where(restinvvar gt 0.0)

       indices = spectral_indices(restwave[good]*1.0D,1D0*restflux[good]-(1D0*speclinefit[good]>0.0),$ ; note!
         ferr=1.0D/sqrt(restinvvar[good]),debug=debug,_extra=extra,/silent)
       raw_indices = spectral_indices(restwave[good]*1.0D,restflux[good]*1.0D,$
         ferr=1.0D/sqrt(restinvvar[good]),debug=debug,_extra=extra,/silent)
       model_indices = spectral_indices(restwave[good]*1.0D,restcontinuum[good]*1.0D,$
         ferr=1.0D/sqrt(restinvvar[good]),debug=debug,_extra=extra,/silent)

       newindices = [indices.indices,strtrim(model_indices.indices,2)+'_model',$
         strtrim(raw_indices.indices,2)+'_raw']
       
       indices = struct_addtags({indices: newindices},struct_addtags(struct_trimtags(indices,$
         except='INDICES'),struct_trimtags(im_struct_trimtags(raw_indices,select=$
         tag_names(raw_indices),newtags=tag_names(raw_indices)+'_raw'),except='INDICES_RAW')))
       indices = struct_addtags({indices: newindices},struct_addtags(struct_trimtags(indices,$
         except='INDICES'),struct_trimtags(im_struct_trimtags(model_indices,select=$
         tag_names(model_indices),newtags=tag_names(model_indices)+'_model'),except='INDICES_MODEL')))

; pack everything into the SPECDATA structure; also make SPECFIT 

       specfit1 = float([[wave/(1.0+z_abs)],[flux*(1.0+z_abs)],$
         [backfit.continuum*(1.0+z_abs)],[speclinefit],$
         [smooth_continuum*(1.0+z_abs)],[invvar/(1.0+z_abs)^2.0]])
       if (igal eq 0L) then specfit = specfit1 else $
         specfit = [[[temporary(specfit)]],[[specfit1]]]
       
       if (igal eq 0L) then begin
          backfit1 = im_empty_structure(struct_trimtags(backfit,$
            except=['CONTINUUM']));,'TEMPLATE_BESTINDX']))
          linefit1 = im_empty_structure(struct_trimtags(linefit,$
            except=['*BOX*']))
          indices1 = im_empty_structure(indices)
          
          specdata1 = {$
            galaxy:                  '', $
            z_obj:                  0.0, $
            specres:            specres, $
            minwave:                0.0, $
            maxwave:                0.0, $
            templatefile:  file_basename(templatefile)}

          specdata = struct_addtags(specdata1,backfit1)
          specdata = struct_addtags(temporary(specdata),linefit1)
          specdata = struct_addtags(temporary(specdata),indices1)
          specdata = replicate(specdata,ngalaxy)
       endif

       specdata[igal].galaxy = galaxy[igal]
       specdata[igal].z_obj = in_zobj
       specdata[igal].minwave = min(wave[where(invvar gt 0)])
       specdata[igal].maxwave = max(wave[where(invvar gt 0)])
       
       junk = specdata[igal] & struct_assign, backfit, $
         junk, /nozero & specdata[igal] = junk
       junk = specdata[igal] & struct_assign, linefit, $
         junk, /nozero & specdata[igal] = junk
       junk = specdata[igal] & struct_assign, indices, $
         junk, /nozero & specdata[igal] = junk

; continuously overwrite SPECDATAFILE in case the code crashes
          
       if (not keyword_set(silent)) then splog, 'Updating '+file_basename(specdatafile)
       mwrfits, reform(specdata), specdatafile, /create 

       if (not keyword_set(silent)) then splog, 'Updating '+file_basename(specfitfile)
       mwrfits, specfit, specfitfile, /create 
       
; finally build the QA plot
       
       qaplot_ispeclinefit, specdata[igal], specfit[*,*,igal], $
         templateinfo[specdata[igal].template_bestindx], $
         postscript=(not keyword_set(doplot))

    endfor
    
    splog, format='("Total time for ISPECLINEFIT = ",G0," '+$
      'minutes.")', (systime(1)-stime0)/60.0
    splog, /close

; compress SPECDATAFILE, close the postscript file, and return     
    
    spawn, 'gzip -f '+specdatafile, /sh
    spawn, 'gzip -f '+specfitfile, /sh

    dfpsclose
    spawn, 'gzip -f '+psfile, /sh
;   spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
;   spawn, '/bin/rm -f '+psfile

return, reform(specdata)
end
