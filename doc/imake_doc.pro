pro imake_doc, outdir=outdir
; build HTML documentation for iSPEC.  an output directory can be
; specified with OUTDIR
    
    idir = getenv('ISPEC_DIR')

;   spawn, 'find '+idir+' -name "*.pro"  -print', prolist
    spawn, 'find '+idir+'/idl/analysis -name "*.pro"  -print', prolist1
    spawn, 'find '+idir+'/idl/misc -name "*.pro"  -print', prolist2
    spawn, 'find '+idir+'/idl/preproc -name "*.pro"  -print', prolist3
    spawn, 'find '+idir+'/idl/spec1d -name "*.pro"  -print', prolist4
    spawn, 'find '+idir+'/idl/spec2d -name "*.pro"  -print', prolist5

    prolist = [prolist1,prolist2,prolist3,prolist4,prolist5]

    if n_elements(outdir) eq 0L then outdir = idir+'/doc/'
    
    make_html_help, prolist, outdir+'ispec_doc.html', $
      /strict, /link_files, title='IDL Help for ISPEC', /verbose

return
end
