Here lies the code for converting 1D solar wind ion spectra measured by DSCOVR into fit parameters.
Each heading contains the name to a module that is an important piece in the DSCOVR data pipeline.



adv_kp.csh
===========
C-shell wrapper script for the DSCOVR datapipeline. It is useful because this script allows the DSCOVR data pipeline to be run as a cron job.
The script compiles the required libraries so the pipeline can run automatically an calls the IDL program calc_through_errors.


calc_through_errors.pro
=======================
This code calls dsc_advanced_kp (from dsc_advanced_kp_fit_v2.pro) over a given time range. This call is done in a way that days with errors
are skipped, but the skipping will fail if the Wind data lags behind the DSCOVR pipeline. When this happens email Mike Stevens and tell 
him you need more Wind data processed to run the DSCOVR pipeline.

```IDL
.r calc_through_error
calc_through_errors,curv
```

calc_through_errors,curv,clobber=clobber,reflag=reflag,start=start

curv : curv is the only required argument which is the current working verision of the
program as an integer
clobber : clobber keyword refits the proton spectra for days already computed by kp_fit
reflag  : reflag keyword tells whether to reproduce the cdf files and their associated flags
start   : start is the start time in Julian days Assumes June 4, 2016 if no day given (use double(julday(M,D,YYYY)))


dsc_advanced_kp_fit_v2.pro
==========================
This is the primary IDL code for extracting the DSCOVR solar wind parameters from the ion spectra.
The naming convention for the file name is a bit off because the main program is dsc_advanced_kp_fit;
however, the program also contains a dummy program call dsc_advanced_kp_fit_v2. The purpose of the 
dummy program is that it allows the full file to compile when using the resolve routine funciton.


dsc_advanced_kp_fit
------------------
dsc_advanced_kp_fit is the primary wrapper which allows all other programs to do their dirty work.
Use the following code in IDL (not SSWIDL) to run the program.

```IDL
.r dsc_advanced_kp_fit_v2
dsc_advanced_kp_fit,2018,1,-1,/save 
```

Above I specified the year (2018), doy (1), version (-1), and the save keyword;
however, dsc_advanced_kp_fit may take many more keywords. The usage for dsc_advanced_kp_fit is    
 dsc_advanced_kp_fit, year, doy, version, show=show, hold=hold, $     
                         save=save, rezero=rezero, wind_show=wind_show$     
                         averaging_length = averaging_length, $     
                         flybacks=flybacks, $     
                         neg_offset=neg_offset, verbose=verbose, $     
                         clobber = clobber           

Now here are the descriptions for each variable and keyword:    
year : The year of the ion spectra you want to reduce in YYYY format. Value must be an integer.    
doy  : The day in the above year you want to reduce the spectrum for in D format. Value must be an integer.   
version : The version of the cdf files you want to create. Value must be an integer (I use negative values for testing). 
 You should increment this value whenever you deliver a new distribution 
or you make modifications to the pipeline that you will pass along to the NASA cdf archive. NASA cares about this a lot and you 
should be careful not to overwrite version of a file you already distributed to the archive. People will not be happy.    
show : Use to show the ionspectrum and fit over the course of the day. This slows down the pipeline, but it useful if you
are curious about parameters returned by the pipeline.    
hold : Deprecated, does nothing.    
save : Saves the analysis to a cdf file and plots. By default the program will not save the analysis. If you flip the
save keyword the program will not overwrite a previously existing cdf file unless you also set the clobber keyword. 
Setting the save keyword performs a couple checks before deciding to save. First, it checks if the out cdf exists 
and if not go ahead and save normally. Then it checks if the idl save file exists if it does just restore the save
file and save the cdf and plot. Finally, it checks if the clobber keyword is set and if it run the routine and overwrite
any existing cdf files.    
rezero : This keyword is set by default and removes a base line from the DSCOVR ion spectrum, which is likely due to a 
poor ground. Set rezero = 0 to remove base line subtraction (but do not do that).    
wind_show : Overplots the Wind ion spectrum on the DSCOVR ion spectrum. Only works if show keyword already set.    
averaging_length : A keyword that sets the a number of spectra to use before storing the parameters from the DSCOVR spectra. Value must be an integer.
The default value is 15, which in the current observing mode corresponds to 1 minute.    
flybacks : Deprecated, for pre-June 2016 DSCOVR data to zero out fly back noise.    
neg_offset : Maximum Value to use in each channel to find the offset (Default = 0). As the name implies, this should be a negative
number or a postive addition to the spectra.    
verbose : Prints more information about the analysis.    
clobber : If using save, using clobber will over write any existing cdf file with the given year,doy, and version. You should only
ever set clobber if you are working on a version that is NOT distributed to the NASA archive.     



compare_wind_dscovr.pro
=======================
Compares WIND and DSCOVR data for a given day and plot puts plots with 
WIND and DSCOVR overplotted (_ontop_)
and DSCOVR-WIND (_delta_) in ../out_plots/. This program is automatically called by dsc_advanced_kp_fit when using the save keyword.
Use these plots to verify DSCOVR data by eye. To run the program as a stand alone use the following syntax:    

```IDL
.r compare_wind_dscovr
compare_wind_dscovr,2018,1,-1
```

Below is the full set of inputs for the program:    

year :  The year of the parameters you want to plot in YYYY format. Value must be an integer.  
doy  :  The day in the above year where you want to compare the DSCOVR parameters to Wind visually in D format. Value must be an integer. 
version : The version of the cdf files you want to create. Value must be an integer (I use negative values for testing). 
span :  span in days to compare WIND and DSCOVR (default = 1)
filefmt : filefmt is the format of the idl save file you need to restore (default = '("dscovr_h1_fc_",I4,I02,I02,"_v",I02)')
archive : archive is the location of the cdf files archive (default = '/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected')
ffmt : ffmt is the file format of the output png file with plasma dscovr plotted ontop of wind (default =  '("../out_plots/compare_wind_dsc_ontop_",I4,"_",I03,".png")')
dfmt : dfmt is the file format of the output png file with dscovr wind plasma difference plotted (default = '("../out_plots/compare_wind_dsc_delta_",I4,"_",I03,".png")')
mfmt : mfmt is the file format of the output png file with magnetic field dscovr plotted ontop of wind (default = '("../out_plots/compare_wind_dsc_ontop_mag_",I4,"_",I03,".png")')


idl_dscovr_to_cdf.pro
=====================
Coverts idl save files to CDF format and stores the cdf files at orchive in the call.

If not ran from dsc_advanced_kp_fit_v2 then the follow contingencies are required following:    
```IDL
@'/crater/utilities/idl/mike/idlstartup'
@compile_cdaweb
@compile_IDLmakecdf
```

Normally to program runs for dsc_advanced_kp_fit, so the contigencies are not required. The call for the program is as follow:    
```IDL
.r idl_dscovr_to_cdf
idl_dscovr_to_cdf,2018,1,-1
```

Below is the full call for the program:    
idl_dscovr_to_cdf,year,doy,version,archive=archive,skeleton=skeleton,filefmt=filefmt, $     
                  outfmt=outfmt,orchive=orchive,outdom=outdom    

Below is the full set of inputs for the program:    
year :  The year of the parameters you want to covert to cdf in YYYY format. Value must be an integer.       
doy  :  The day in the above year where you want to convert the idl sav file to a cdf in D format. Value must be an integer.      
version  : The version of the cdf files you want to create. Value must be an integer (I use negative values for testing).      
archive  : archive is the location of the idl sav file archive (default ='/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected')     
orchive  : orchive is the location of the output cdf file directory (default ='/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected')
skeleton : skeleton is the full path to the cdf skeleton (default ='../skeleton/dscovr_h1_fc_0000000_v01.cdf')    
filefmt  : filefmt is the file format of the files in the idl save directory (default ='("dsc_fc_advkp_1minute_corrected_",I4,"_",I03,".idl")')      
outfmt   : outfmt is the file format of the output cdf files (default   = '("dscovr_h1_fc_",I4,I02,I02,"_v",I02,".cdf")')      

