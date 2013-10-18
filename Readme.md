J. Moustakas, 2003 July 7, U of A

Assume that iSPEC and its associated libraries have been downloaded
into your ${HOME}/idl directory.  Define the following environment
variables in your ~/.cshrc file:

	setenv IMPRO_DIR     ${HOME}/idl/impro
	setenv ISPEC_DIR     ${HOME}/idl/ispec
	setenv IDLUTILS_DIR  ${HOME}/idl/idlutils
	setenv IDLSPEC2D_DIR ${HOME}/idl/idlspec2d
	setenv PATH ${IDLUTILS_DIR}/bin:${IDLSPEC2D_DIR}/bin:${PATH}

Now make sure that ${HOME}/idl is in your local IDL path.  In your
~/.idlenv file add the line

	setenv IDL_PATH         $IDL_PATH{:}+${HOME}/idl

To create a local HTML version of the documentation of ISPEC,
including links to the programs, go into ${ISPEC_DIR}/doc, start IDL,
and type 'imake_doc'.  This routine will create a file called
'ispec_doc.html' in the same directory which can then be viewed with a
browser.
