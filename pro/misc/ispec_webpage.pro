;+
; NAME:
;	ISPEC_WEBPAGE
;
; PURPOSE:
;	Generate web page visualizations for iSPEC2d.
;
; CALLING SEQUENCE:
;       ispec_webpage, logbase, rootname=, weblist=, red_path=, $
;          html_path=, /html_only, /cleanpng, /mainmask, /nopsconvert, $
;          /silent
;
; INPUTS:
;       logbase - base filename for the webpages; a LOGBASE.html file
;                 is created in the HTML_PATH directory and all the
;                 output of this routine is placed in the
;                 HTML_PATH/LOGBASE subdirectory, which will be
;                 created if it does not exist
;
; OPTIONAL INPUTS:
;	rootname   - root prefix of the reduced 2D images (default
;                    'f')
;       weblist    - as an alternative to ROOTNAME, you can specify
;                    the list images to include in the web page using
;                    this keyword; for example, you may not want the
;                    standard stars on the webpage (overwrites
;                    ROOTNAME) 
;	red_path   - path name to the 2D FITS spectra (default CWD)
;	html_path  - path where html and png files should be written 
;                    (default './WWW/'); this directory is created if
;                    it does not exist
;
; KEYWORD PARAMETERS:
;	html_only  - only create html files, do not create png images;
;                    useful for quickly updating html files
;       cleanpng   - delete all *.png files in HTML_PATH/LOGBASE 
;       mainmask   - show the bad pixel mask on the main web page
;                    (useful for identifying spurious cosmic-ray
;                    features, etc.)
;       nopsconvert - do not convert the PS files to PNG files
;
; OUTPUTS:
;       The main HTML page includes a single image for each ;
;       observation, along with some basic observing information.
;       This image links to a dedicated HTML file showing the
;       reduction history of each object.  
;
; PROCEDURES USED:
;       CWD(), SPLOG, REPSTR(), IFORAGE(), REMOVE, ISPEC_CREATE_PNG,
;       HEADFITS(), SXPAR(), CMSET_OP(), NICEPRINT
;
; INTERNAL SUPPORT ROUTINES:
;       FILE_HISTORY()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas 2003 July 15, U of A, written based in *large*
;          part on Karl Gordon's KGWEBPAGE
;       jm04jul15uofa - remove one-dimensional spectra from the file
;                       list 
;       jm04sep27uofa - bug fix: the unix command RM chokes when the
;                       argument list is too long
;       jm05jun16uofa - updated documentation; bug fix when organizing
;                       the PS files; streamlined calling of
;                       ISPEC_CREATE_PNG
;       jm05jun20uofa - handle the new sequence of reductions for
;                       iSPEC2d v2.0 (e.g., case when sky subtraction
;                       occurs after cosmic-ray rejection) and
;                       improved postscript file display
;       jm05jul20uofa - added NOPSCONVERT keyword
;       jm05jul28uofa - when files are missing from WEBLIST, list them
;
; Copyright (C) 2003-2005, John Moustakas
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

function file_history, file
; NB.  if sky-subtraction took place after ICRCOMBINE then the
; ccd-processed images will not be recognized properly

    if strmatch(file,'*.gz*') then gzip = '.gz' else gzip = ''
    
    header = headfits(file)
    history = {file: file, header: header, ncomb:  fix(sxpar(header,'INCOMB'))}
    
; calibration history header keywords

    calibkeys = ['IBIAS','IFLAT','ISFLAT','WMAPNAME','TRACENAM','IEXTNAME','SENSNAME']
    for k = 0L, n_elements(calibkeys)-1L do begin
       entry = sxpar(header,calibkeys[k],count=count)
       if (count eq 0L) then $
         history = create_struct(history,calibkeys[k],'') else $
         history = create_struct(history,calibkeys[k],strtrim(entry,2))
    endfor

; the following keywords have to be treated individually to get the
; file history right

    iraw = sxpar(header,'IRAW',count=irawcount)
    icosmic = sxpar(header,'ICOSMIC',count=icosmiccount)
    iskysb = sxpar(header,'ISKYSB',count=iskysbcount)
    iwave = sxpar(header,'IWAVE',count=iwavecount)
    itrace = sxpar(header,'ITRACE',count=itracecount)
    iflux = sxpar(header,'IFLUX',count=ifluxcount)

; the case of NCOMB>0 is treated below
    
    if (history.ncomb eq 0L) then begin
       if (irawcount eq 0L) then begin
          history = create_struct(history,'IRAW','') 
          message, 'This should never happen.'
       endif else begin
          history = create_struct(history,'IRAW',iraw)
          icomb = 'r'+strtrim(iraw,2)+gzip
          history = create_struct(history,'ICOMB',icomb)
       endelse
    endif
       
    if (icosmiccount eq 0L) then $
      history = create_struct(history,'ICOSMIC','') else $
      history = create_struct(history,'ICOSMIC',icosmic)

    if (iskysbcount eq 0L) then $
      history = create_struct(history,'ISKYSB','') else $
      history = create_struct(history,'ISKYSB',iskysb)

    if (iwavecount eq 0L) then $
      history = create_struct(history,'IWAVE','') else $
      history = create_struct(history,'IWAVE',iwave)

    if (itracecount eq 0L) then $
      history = create_struct(history,'ITRACE','') else $
      history = create_struct(history,'ITRACE',itrace)

    if (ifluxcount eq 0L) then $
      history = create_struct(history,'IFLUX','') else $
      history = create_struct(history,'IFLUX',iflux)

; write down the proper file history and comments

    filehist = ''
    comments = ''
    
    if (irawcount eq 1) then begin
       filehist = [filehist,iraw]
       comments = [comments,'raw image']
       filehist = [filehist,icomb]
       comments = [comments,'ccd-processed']
    endif
    
    if (icosmiccount eq 1) then begin
       filehist = [filehist,'c'+icosmic]
       comments = [comments,'cosmic-ray rejected']
    endif
    
    if (iskysbcount eq 1) then begin
       filehist = [filehist,'s'+iskysb]
       comments = [comments,'sky-subtracted']
    endif

    if (iwavecount eq 1) and (itracecount eq 0) and (ifluxcount eq 0) then begin
       filehist = [filehist,'w'+iwave] ; <-- bug!
       comments = [comments,'wavelength calibrated']
    endif

    if (iwavecount eq 0) and (itracecount eq 0) and (ifluxcount eq 1) then begin
       filehist = [filehist,'f'+iflux] ; <-- bug!
       comments = [comments,'flux calibrated']
    endif

    if (iwavecount eq 0) and (itracecount eq 1) and (ifluxcount eq 0) then begin
       filehist = [filehist,'d'+itrace]
       comments = [comments,'distortion corrected']
    endif

    if (iwavecount eq 1) and (itracecount eq 0) and (ifluxcount eq 1) then begin
       filehist = [filehist,'fw'+iwave]
       comments = [comments,'flux calibrated, wavelength calibrated']
    endif

    if (iwavecount eq 1) and (itracecount eq 1) and (ifluxcount eq 1) then begin
       filehist = [filehist,'fdw'+iwave]
       comments = [comments,'flux calibrated, distortion corrected, wavelength calibrated']
    endif

    filehist = filehist[1L:n_elements(filehist)-1L]
    comments = comments[1L:n_elements(comments)-1L]
    
; the following keywords can appear NCOMB times, so treat them
; separately 
    
    keys = ['IRAW','ISKYSB','ICOMB'] ; the order of these keywords is important!

    rfilehist = ''
    rcomments = ''

    if (history.ncomb gt 0L) then begin

       for i = 0L, history.ncomb-1L do begin

          for k = 0L, n_elements(keys)-1L do begin
             
             loopkey = keys[k]+string(i+1,format='(I2.2)') ; <-- assumes fewer than 99 images
             entry = sxpar(header,loopkey,count=count)

             if (count eq 0L) then $
               history = create_struct(history,loopkey,'') else begin

                history = create_struct(history,loopkey,strtrim(entry,2))

                case k of
                   0L: begin
                      rfilehist = [rfilehist,strtrim(entry,2)]
                      rcomments = [rcomments,'raw image']
                   end
                   1L: begin
                      rfilehist = [rfilehist,['','s']+strtrim(entry,2)]
                      rcomments = [rcomments,'ccd-processed','sky subtracted']
                   end
                   2L: begin
                      rfilehist = [rfilehist,''+strtrim(entry,2)]
                      rcomments = [rcomments,'ccd-processed']
                   end
                   else: 
                endcase

             endelse

          endfor

       endfor
    
       rfilehist = rfilehist[1L:n_elements(rfilehist)-1L]
       rcomments = rcomments[1L:n_elements(rcomments)-1L]

       history = create_struct(history,'FILEHIST',reverse([rfilehist,filehist]),$
         'COMMENTS',reverse([rcomments,comments]))

    endif else begin

       history = create_struct(history,'FILEHIST',reverse(filehist),$
         'COMMENTS',reverse(comments))

    endelse

return, history
end

pro ispec_webpage, logbase, rootname=rootname, weblist=weblist, red_path=red_path, $
  html_path=html_path, html_only=html_only, png_only=png_only, cleanpng=cleanpng, $
  mainmask=mainmask, nopsconvert=nopsconvert, silent=silent

    nlogbase = n_elements(logbase)
    if (nlogbase ne 1L) then begin
       doc_library, 'ispec_webpage'
       return
    endif
    
    if (n_elements(rootname) eq 0L) then rootname = 'f'
    if (n_elements(red_path) eq 0L) then red_path = cwd()
    if (n_elements(html_path) eq 0L) then html_path = './WWW/'

    if file_test(red_path,/directory) eq 0L then begin
       splog, 'RED_PATH '+red_path+' does not exist.'
       return
    endif

    if file_test(html_path,/directory) eq 0L then begin
       if not keyword_set(silent) then splog, 'HTML_PATH '+html_path+' does not exist - creating it.'
       spawn, ['mkdir -p '+html_path], /sh
    endif

    if file_test(html_path+logbase,/directory) eq 0L then begin
       if not keyword_set(silent) then splog, 'Path '+html_path+logbase+' does not exist - creating it.'
       spawn, ['mkdir -p '+html_path+logbase], /sh
    endif else if not keyword_set(silent) then splog, 'Writing web pages to path '+html_path+logbase+'.'

    if n_elements(weblist) eq 0L then begin
    
       pushd, red_path
       files = file_search(rootname+'*.fits*',count=n_files)
       popd

       if (n_files eq 0L) then begin
          splog, 'No FITS files found in RED_PATH '+red_path+' with ROOTNAME '+rootname+'.'
          return
       endif

    endif else begin

       files = weblist
       tempfiles = file_search(red_path+weblist,count=n_files)
       if (n_files ne n_elements(weblist)) then begin
          splog, 'The following files in WEBLIST do not exist:'
          missingfiles = cmset_op(files,'AND',/not2,file_basename(tempfiles))
          niceprint, missingfiles
          return
       endif

    endelse
    
; clean up PNG, PS, and HTML files from calling this routine
; previously; *note* I've run into an error of "Argument list too
; long" for /bin/rm, so below you see the work-around

    if keyword_set(cleanpng) then begin
       if not keyword_set(silent) then splog, 'Deleting all PNG and PS files from '+html_path+logbase+'.'
       pushd, html_path+logbase
       spawn, ['/bin/ls -f *.png*  | xargs /bin/rm -f '], /sh
       spawn, ['/bin/ls -f *.ps*   | xargs /bin/rm -f '], /sh
       spawn, ['/bin/ls -f *.html* | xargs /bin/rm -f '], /sh
       popd
    endif

    pushd, red_path & psfiles = file_search('*.ps',count=n_psfiles) & popd
    if (n_psfiles eq 0L) then if not keyword_set(silent) then $
      splog, 'No postscript files found in RED_PATH '+red_path+'.'

; exclude "FITSLIST" files

    if (n_psfiles ne 0L) then begin
       junk = where(strmatch(psfiles,'*fitslist*',/fold) eq 1B,njunk,comp=keep,ncomp=nkeep)
       if (nkeep ne 0L) then begin
          psfiles = psfiles[keep]
          n_psfiles = nkeep
       endif
    endif
    
    pngfiles = repstr(psfiles,'.ps','.png')
    if keyword_set(html_only) then nopsconvert = 1L

    if (keyword_set(nopsconvert) eq 0L) then begin
       for i = 0L, n_psfiles-1L do begin
          if not keyword_set(silent) then splog, 'Converting '+psfiles[i]+' to '+pngfiles[i]+'.'
          spawn, ['cp -pf '+red_path+psfiles[i]+' '+html_path+logbase+'/'], /sh
          spawn, ['convert '+red_path+psfiles[i]+' '+html_path+logbase+'/'+pngfiles[i]], /sh
       endfor
    endif

; forage header information from the spectra and sort by Julian date

    if not keyword_set(silent) then splog, 'Foraging header information for '+$
      string(n_files,format='(I0)')+' files.'
    forage = iforage(red_path+files) ; <-- HERE!!!!

    keep = where(forage.naxis eq 2L,nkeep)
    if (nkeep eq 0L) then begin
       splog, 'No two-dimensional images found!'
       return
    endif 

    if (nkeep ne n_files) then splog, 'Removing '+string(n_files-nkeep,format='(I0)')+$
      ' one-dimensional images from the file list.'

    n_files = nkeep
    files = files[keep]
    forage = forage[keep]

    srt = sort(forage.jd)
    files = files[srt]
    forage = forage[srt]
    
    fitsfiles = files
    for k = 0L, n_files-1L do files[k] = strmid(fitsfiles[k],0,strpos(fitsfiles[k],'.fits'))

; generate the postscript HTML page

    if (n_elements(n_psfiles) ne 0L) and (not keyword_set(png_only)) then begin
       
       npscols = 3L
       xwidth = string(fix(100/float(npscols)/3.0),format='(I0)')
       iwidth = string(fix(100/float(npscols)),format='(I0)')

; attempt to sort the PSFILES into iSPEC order

       lastindx = lindgen(n_psfiles)
       psbias = where(strmatch(psfiles,'qaplot*bias*') eq 1B,npsbias)
       psillum = where(strmatch(psfiles,'qaplot*illum*') eq 1B,npsillum)
       psdome = where(strmatch(psfiles,'qaplot*domeflat*') eq 1B,npsdome)
       pssky = where(strmatch(psfiles,'qaplot*skyflat*') eq 1B,npssky)
       psresp = where(strmatch(psfiles,'qaplot*response*') eq 1B,npsresp)
       pslamp = where(strmatch(psfiles,'qaplot*arc*') eq 1B,npslamp)
       pssens = where(strmatch(psfiles,'qaplot*sens*') eq 1B,npssens)
       pstell = where(strmatch(psfiles,'qaplot*telluric*') eq 1B,npstell)

       newpsfiles = ''
       newpsindx = -1L
       if (npsbias ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[psbias]]
          newpsindx = [newpsindx,psbias]
       endif
       if (npsillum ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[psillum]]
          newpsindx = [newpsindx,psillum]
       endif
       if (npsdome ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[psdome]]
          newpsindx = [newpsindx,psdome]
       endif
       if (npssky ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[pssky]]
          newpsindx = [newpsindx,pssky]
       endif
       if (npsresp ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[psresp]]
          newpsindx = [newpsindx,psresp]
       endif
       if (npslamp ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[pslamp]]
          newpsindx = [newpsindx,pslamp]
       endif
       if (npssens ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[pssens]]
          newpsindx = [newpsindx,pssens]
       endif
       if (npstell ne 0L) then begin
          newpsfiles = [newpsfiles,psfiles[pstell]]
          newpsindx = [newpsindx,pstell]
       endif

       newpsindx = newpsindx[1L:n_elements(newpsindx)-1L]
       
; subscript any additional files

       if (n_elements(newpsfiles) lt n_elements(psfiles)) then begin

          remove, newpsindx, lastindx
          newpsfiles = [newpsfiles,psfiles[lastindx]]
          
       endif

       psfiles = newpsfiles[1L:n_elements(newpsfiles)-1L]
       n_psfiles = n_elements(psfiles)

       rootpsfiles = repstr(psfiles,'.ps','')
       pngfiles = repstr(psfiles,'.ps','.png')
       nlistcols = 3L
       
       if not keyword_set(silent) then splog, 'Writing '+html_path+logbase+'_ps.html'
       
       openw, lun1, html_path+logbase+'_ps.html', /get_lun
       printf, lun1, '<html>'
       printf, lun1, '<head>'
       printf, lun1, '<title>iSPEC Postscript Output for '+logbase+'</title>'
       printf, lun1, '<style type="text/css">'
       printf, lun1, '<!--'
       printf, lun1, 'body {'
       printf, lun1, '   color: black;'
       printf, lun1, '   background-color: white;'
       printf, lun1, '   text-align: center;'
       printf, lun1, '}'
       printf, lun1, '-->'
       printf, lun1, '</style>'
       printf, lun1, '</head>'
       printf, lun1, '<body>'
       printf, lun1, '<h1>iSPEC2d Postscript Output for '+logbase+'</h1>'
       printf, lun1, '</br>'

; generate the menu       
       
       printf, lun1, '<table border="1" align="center" frame="box" width="100%">'
       printf, lun1, '<tbody>'
       for i = 0L, ceil(n_psfiles/float(nlistcols))-1L do begin

          printf, lun1, '<tr width="100%">'
          for j = 0L, nlistcols-1L do begin
             indx = j + i*nlistcols
             if (indx le n_psfiles-1L) then begin
                printf, lun1, '<td width="1%">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')+'</td>'
                printf, lun1, '<td width="10%"><a href="'+logbase+'/'+$
                  psfiles[indx]+'">'+rootpsfiles[indx]+'</a></td>'
             endif
          endfor
          printf, lun1, '</tr>'

       endfor
       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '<br />'

       printf, lun1, '<table align="center" cellpadding="10" border="2" cellspacing="10" width="100%">'
       printf, lun1, '<tbody>'

       for i = 0L, ceil(n_psfiles/float(npscols))-1L do begin

; images          
          
          printf, lun1, '<tr width="100%">'
          for j = 0L, npscols-1L do begin
             indx = j + i*npscols
             if (indx le n_psfiles-1L) then begin
                printf, lun1, '<td width='+iwidth+'%><a href="'+logbase+'/'+$
                  pngfiles[indx]+'"><img width="100%" src="'+$
                  logbase+'/'+pngfiles[indx]+'"></a></td>'
             endif
          endfor
          printf, lun1, '</tr>'

; image captions          
          
          printf, lun1, '<tr width="100%">'
          for j = 0L, npscols-1L do begin
             indx = j + i*npscols
             if (indx le n_psfiles-1L) then begin
                printf, lun1, '<td class="center">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')+' <a href="'+logbase+'/'+psfiles[indx]+$
                  '">'+rootpsfiles[indx]+'</a></td>'
             endif
          endfor
          printf, lun1, '</tr>'
          
       endfor 

       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '</body>'
       printf, lun1, '</html>'
       free_lun, lun1

    endif

; generate the spectral HTML pages

    if (not keyword_set(png_only)) then begin
       
       if not keyword_set(silent) then splog, 'Writing '+html_path+logbase+'.html'
       openw, unit1, html_path+logbase+'.html',/get_lun
       printf, unit1, '<html>'
       printf, unit1, '<head>'
       printf, unit1, '<title>Visual Log for '+logbase+'</title>'
       printf, unit1, '<style type="text/css">'
       printf, unit1, '<!--'
       printf, unit1, 'body {'
       printf, unit1, '   color: black;'
       printf, unit1, '   background-color: white;'
       printf, unit1, '   text-align: center;'
       printf, unit1, '}'
       printf, unit1, '-->'
       printf, unit1, '</style>'
       printf, unit1, '</head>'
       printf, unit1, '<body>'
       printf, unit1, '<h1>Visual Log for '+logbase+'</h1>'
       printf, unit1, '<h3><a href="'+logbase+'_ps.html">iSPEC2d Postscript Output</a></h3>'
       printf, unit1, '<table align="center" cellpadding="2" border="2" cellspacing="0">'
       printf, unit1, '<tbody>'

    endif
       
    n_cols = 0
    max_cols = 1

    for i = 0L, n_files-1L do begin

       if file_test(red_path+fitsfiles[i]) then begin

          if not keyword_set(silent) then splog, 'Working on ' + files[i]+'.'

; determine the file history then create the individual HTML file

          history = file_history(fitsfiles[i])

          if (not keyword_set(png_only)) then begin
             
             openw, unit2, html_path+logbase+'/'+files[i]+'.html', /get_lun
             printf, unit2, '<html>'
             printf, unit2, '<head>'
             printf, unit2, '<title>Log page for '+files[i]+'</title>'
             printf, unit2, '<style type="text/css">'
             printf, unit2, '<!--'
             printf, unit2, 'body {'
             printf, unit2, '   color: black;'
             printf, unit2, '   background-color: white;'
             printf, unit2, '   text-align: center;'
             printf, unit2, '   font-family: "Times New Roman";'
             printf, unit2, '}'
             printf, unit2, 'th {text-align: center;}'
             printf, unit2, 'td {text-align: center;}'
             printf, unit2, 'td.left {text-align: left;}'
             printf, unit2, '.heavy {'
             printf, unit2, '   font-weight: bold;'
             printf, unit2, '   margin-bottom: 1ex;'
             printf, unit2, '   font-size: 110%;'
             printf, unit2, '}'
             printf, unit2, '-->'
             printf, unit2, '</style>'
             printf, unit2, '</head>'
             printf, unit2, '<body>'
             printf, unit2, '<h1>Log page for '+files[i]+'</h1>'
             printf, unit2, '<p>'
             printf, unit2, '<table align="center" cellpadding="2" border="1" cellspacing="0">'
             printf, unit2, '<tbody>'
             printf, unit2, '<tr>'
             printf, unit2, '<th>object</th><th>date</th><th>time</th><th>airmass</th><th>exptime</th>'
             printf, unit2, '</tr>'
             printf, unit2, '<tr>'
             object = forage[i].object & if strcompress(object,/remove) eq '' then object = ' - '
             printf, unit2, '<td>'+object+'</td>'
             printf, unit2, '<td>'+forage[i].date+'</td>'
             printf, unit2, '<td>'+forage[i].ut+'</td>'
             printf, unit2, '<td>'+string(forage[i].airmass,format='(F5.2)')+'</td>'
             printf, unit2, '<td>'+string(forage[i].exptime,format='(F10.2)')+'</td>'
             printf, unit2, '</tr>'
             printf, unit2, '</tbody>'
             printf, unit2, '</table>'
             printf, unit2, '</p>'

; setup sub files

             printf, unit2, '<p>'
             printf, unit2, '<table cellpadding="4" border="1" cellspacing="0">'
             printf, unit2, '<tbody>'

          endif
             
          subfiles = history.filehist
          nsubfiles = n_elements(subfiles)

          for j = 0L, nsubfiles-1L do begin

             rootfile = strmid(subfiles[j],0,strpos(subfiles[j],'.fits'))
             comment = history.comments[j]
             
             if file_test(red_path+subfiles[j],/regular) then begin

                png_file = strmid(subfiles[j],0,strpos(subfiles[j],'.fits'))+'.png'
                png_file_tb = repstr(png_file,'.png','_tb.png')

                png_mask_file = repstr(png_file,'.png','_mask.png')
                png_mask_file_tb = repstr(png_file_tb,'_tb.png','_mask_tb.png')

                if (j eq 0L) then begin

                   main_png_file = png_file
                   main_png_file_tb = png_file_tb

                endif
                
; generate PNG images of the images and bad pixel masks

                if not keyword_set(html_only) then begin

                   print, '   Writing '+html_path+logbase+'/'+png_file+'.'

; spectrum plus thumbnail                   

                   ispec_create_png, red_path+subfiles[j], png_file=html_path+logbase+'/'+png_file, $
                     png_thumbfile=html_path+logbase+'/'+png_file_tb, file_ext=0

; mask plus thumbnail                   

                   if (strmatch(comment,'*raw*',/fold) eq 0B) then begin ; the raw image has no mask

                      ispec_create_png, red_path+subfiles[j], png_file=html_path+logbase+'/'+png_mask_file, $
                        png_thumbfile=html_path+logbase+'/'+png_mask_file_tb, file_ext=2, /mask

                   endif
                      
                endif

                if (not keyword_set(png_only)) then begin

                   printf, unit2,'<tr><td class="left"><span class="heavy">'+rootfile+' ('+comment+')</span><br /><br />'
                   printf, unit2,'<a href="'+png_file+'"><img width="100%" src="'+png_file_tb+'"></a><br /><br />'
                   if (strmatch(comment,'*raw*',/fold) eq 0B) then $
                     printf, unit2,'<a href="'+png_mask_file+'"><img width="100%" src="'+png_mask_file_tb+'"></a>'
                   printf, unit2, '</td></tr>'

                endif
                   
             endif 

          endfor

          if (not keyword_set(png_only)) then begin
             
             printf, unit2, '</tbody>'
             printf, unit2, '</table>'
             printf, unit2, '</p>'
             printf, unit2, '</body>'
             printf, unit2, '</html>'
             free_lun, unit2

; add this observation to the log

             if (n_cols EQ 0) then printf, unit1, '<tr>'

             printf, unit1, '<td>'
             printf, unit1, '<a href="'+logbase+'/'+files[i]+'.html">' + files[i] + '</a>'
             printf, unit1, '<br />'+object
             printf, unit1, '<br />'+forage[i].date
             printf, unit1, '<br />'+forage[i].ut
             printf, unit1, '<br /> AM = '+string(forage[i].airmass,format='(F5.2)')
             printf, unit1, '<br />'+string(forage[i].exptime,format='(F10.2)')+' secs'
             printf, unit1, '</td>'
             printf, unit1, '<td>'
             if keyword_set(mainmask) then $
               printf, unit1, '<a href="'+logbase+'/'+files[i]+'.html"> <img width="100%" src="'+logbase+'/'+$
               repstr(main_png_file_tb,'_tb.png','_mask_tb.png')+'"></a>' else $
               printf, unit1, '<a href="'+logbase+'/'+files[i]+'.html"> <img width="100%" src="'+logbase+'/'+$
                 main_png_file_tb+'"></a>'
             printf, unit1, '</td>'
             n_cols = n_cols + 1
             if (n_cols EQ max_cols) then begin
                printf, unit1, '</tr>'
                n_cols = 0
             endif

          endif
             
       endif 

;      print, main_png_file+'/'+main_png_file_tb
       
    endfor

    if (not keyword_set(png_only)) then begin
       
       if (n_cols LT max_cols) then printf, unit1, '</tr>'

       printf, unit1, '</tbody>'
       printf, unit1, '</table>'
       printf, unit1, '</body>'
       printf, unit1, '</html>'
       free_lun, unit1

    endif
       
return
end
