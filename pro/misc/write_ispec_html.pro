;+
; NAME:
;       WRITE_ISPEC_HTML
;
; PURPOSE:
;       Generate the iSPEC web pages.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
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
;       J. Moustakas, 2004 Feb 24, U of A
;-

pro write_ispec_html

    webpath = getenv('HOME')+'/public_html/ispec/'

    pushd, webpath
    
    splog, 'Generating the iSPEC home page.'
    openw, lun, 'index.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://cerebus.as.arizona.edu/~ioannis">'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="http://cerebus.as.arizona.edu/~ioannis/research/ispec/" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="http://cerebus.as.arizona.edu/~ioannis/styles/ispec.css" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>iSPEC: IDL Spectroscopic Data Reduction Software</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>iSPEC</h1>'
    printf, lun, '<h2>IDL Spectroscopic Data Reduction Software</h2>'
    printf, lun, '<br />'
;   printf, lun, '<a id="contents"><h2 class="noindent">table of contents</h2></a>'
    printf, lun, '<div class="contents">'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="#introduction">introduction</a></li>'
    printf, lun, '<li><a href="#download">download & install</a></li>'
    printf, lun, '<li><a href="#documentation">documentation</a></li>'
    printf, lun, '<li><a href="#examples">examples</a></li>'
    printf, lun, '</ul>'
    printf, lun, '</div>'
    printf, lun, '<br />'
    printf, lun, '<div class="section">'
    printf, lun, '<h1 class="center"><a id="introduction">introduction</a></h1>'
    printf, lun, '<p>iSPEC is a long-slit data reduction and analysis software '+$
      'package written in IDL. The purpose of iSPEC is to reduce long-slit '+$
      'spectroscopic observations to produce fully calibrated two- and one-dimensional spectra. '+$
      ''
;     'iSPEC can be used both in IDL command-line mode for 100% control over each step in the data '+$
;     'reduction process or in batch mode for fast, flexible, and repeatable data reduction. '+$
;     'One particular strength of iSPEC is that errors and bad pixels are propagated through '+$
;     'every stage of the reduction, enabling a quantitative assessment of each pixel in the '+$
;     'final one-dimensional spectrum.  iSPEC is still in the development phase. Although the IDL '+$
;     'code itself has been extensively documented, a minimal amount of other documentation yet '+$
;     'exists. "Recipe" documents will follow soon! My thanks go to Karl Gordon and Andy Marble '+$
;     'for testing my code extensively on their own spectroscopic
;     observations.'
    printf, lun, '</p>'
    printf, lun, '</div>'
    printf, lun, '<hr />'
    printf, lun, '<br />'

    printf, lun, '<div class="section">'

    printf, lun, '<h1 class="center"><a id="download">download & install</a></h1>'

    printf, lun, '<p>In developing iSPEC I relied on a number of indispensable'
    printf, lun, 'astronomical IDL libraries that must be installed in your local IDL'
    printf, lun, 'path before iSPEC can be used.  I do not include these libraries'
    printf, lun, 'transparently with iSPEC because most users will already have, for'
    printf, lun, 'example, the Goddard IDL astronomical library.</p>'

    printf, lun, '<p>Below are links to gzipped tar files containing all the IDL'
    printf, lun, 'routines you will need to successfully use iSPEC.  You can download'
    printf, lun, 'the libraries directly here or from the originating website.'
    printf, lun, 'Installation instructions are given below the table.  Eventually I'
    printf, lun, 'will be providing an anonymous CVS server to all these libraries.</p>'

    printf, lun, '<table> '
    printf, lun, '<tr>'
    printf, lun, '<td><a href="http://idlastro.gsfc.nasa.gov">Goddard IDL Library</a></td>'
    printf, lun, '<td><a href="goddard.tar.gz">[goddard.tar.gz (~1.1 MB)]</a></td>'
    printf, lun, '</tr> '
    printf, lun, '<tr>'
    printf, lun, '<td><a href="http://spectro.princeton.edu/idlspec2d_install.html">'
    printf, lun, "Dave Schlegel's idlutils \& idlspec2d</a></td>"
    printf, lun, '<td><a href="sdsspro.tar.gz">[sdsspro.tar.gz (~18.6 MB)]</a></td>'
    printf, lun, '</tr> '
    printf, lun, '<tr>'
    printf, lun, '<td><a href="http://cow.physics.wisc.edu/~craigm/idl/idl.html">Craig Markwardt IDL Programs</a></td>'
    printf, lun, '<td><a href="markwardt.tar.gz">[markwardt.tar.gz (~459 kB)]</a></td>'
    printf, lun, '</tr> '
    printf, lun, '<tr>'
    printf, lun, '<td><a href="http://cerebus.as.arizona.edu/~ioannis/research/ispec/ispec.html">iSPEC</a></td>'
    printf, lun, '<td><a href="ispec.tar.gz">[ispec.tar.gz (~12 MB)]</a></td>'
    printf, lun, '</tr> '
    printf, lun, '</table>'
    printf, lun, '<br />'

    printf, lun, '<p>After downloading the tar files choose where you want to unpack'
    printf, lun, 'them.  For example, all of my IDL libraries are in ${HOME}/idl.  Once'
    printf, lun, 'that is done read the <a href="INSTALL">INSTALL</a> file to learn how'
    printf, lun, 'to set up all the required shell and IDL environment variables.'
    printf, lun, 'Now go reduce some data!</p>'

    printf, lun, '</div>'
    printf, lun, '<hr />'
    printf, lun, '<br />'
    printf, lun, '<div class="section">'
    printf, lun, '<h1 class="center"><a id="documentation">documentation</a></h1>'
    printf, lun, '<p><a href="documentation.html">documentation</p>'
    printf, lun, '</div>'
    printf, lun, '<hr />'
    printf, lun, '<br />'
    printf, lun, '<div class="section">'
    printf, lun, '<h1 class="center"><a id="examples">examples</a></h1>'
    printf, lun, '<p><a href="examples.html">examples</p>'
    printf, lun, '</div>'
    printf, lun, '<hr />'
    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
      'J. Moustakas</a> with questions or comments.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd
    
stop    

return
end
