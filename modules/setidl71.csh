#/bin/tcsh
#source "/home/jprchlik/.cshrc"
#71 needed for CMS2 
setenv IDL_VER "71"

############################################################################
# re-assign idl specific content based on the desired version
if ( -d "/usr/local/etc/profile.d" ) then
  # re-source the profile.d environment for idl-only
  if ( -r "/usr/local/etc/profile.d/49-exelis-idl.csh" ) then
    source /usr/local/etc/profile.d/49-exelis-idl.csh
  endif
else
  # re-source the environment for idl for non-profile.d machines
  ssxg_setup_exelis_idl
endif
############################################################################
############################################################################

#Set up the name of a file with your initial IDL commands.
setenv IDL_STARTUP .idl_startup
if ( ! -e ${IDL_STARTUP} ) then
   touch ${IDL_STARTUP}
endif

############################################################################
#SSWIDL settings
#Set instrument and packages to be used after  $SSW_INSTR:
#setenv SSW_INSTR    "GEN TRACE SXT CDS CHIANTI SUMER BINARIES SOT XRT EIS EIT MDI SOT AIA IRIS LASCO SECCHI"

# You can find the names of the packages by typing the following at the
# sswidl prompt:
# sswidl> ssw_packages
# You can find the names of the instruments by typing the following at the
# sswidl prompt:
# sswidl> print, ssw_instruments()
############################################################################
#Setup SOHO paths for image_tool (in sswidl)
if ( ! -d "/usr/local/etc/profile.d" ) then
  ssxg_setup_image_tool_soho
endif
setenv IDL_PATH /usr/local/itt/idl/idl${IDL_VER}/lib 

setenv IDL_DLM_PATH "/home/jprchlik/personaladditions/code/idl/cdawlib/source:<IDL_DEFAULT>"

idl

