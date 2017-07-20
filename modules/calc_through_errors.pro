;----------------------------------------------------------------
;
;USAGE
;calc_through_errors,clobber=clobber,reflag=reflag,start=start
;
;COMMENTS
;clobber keyword refits the proton spectra for days already computed by kp_fit
;reflag keyword tells whether to reproduce the cdf files and their associated flags
;start is the start time in Julian days Assumes June 4, 2016 if no day given
;
;----------------------------------------------------------------
pro calc_through_errors,clobber=clobber,reflag=reflag,start=start
if keyword_set(clobber) then clobber = 1 else clobber = 0
if keyword_set(reflag) then reflag = 1 else reflag = 0
if keyword_set(start) then start = start else start = DOUBLE(JULDAY(06,04,2016)) ; Start of good DSCOVR observations
set_plot,'Z'

;resolve/compile dsc_advanced_kp_fit program
resolve_routine,'dsc_advanced_kp_fit_v2',/compile_full_file

curv = 5 ; Current working version
efmt = '("***Skipping ",I4,"/",I03," due to missing DSCOVR data***")';format print day errors
dfmt = '("Analyzing ",I4,"/",I03)';format print day
mess = '#############################################################' ;mark error prints
;start = DOUBLE(JULDAY(01,01,2017)) ; Start of good DSCOVR observations
;start = DOUBLE(JULDAY(05,04,2017)) ; Start of good DSCOVR observations
ender = DOUBLE(systime(/JULIAN,/UTC)-14) ; loop until 14 days before today
recom = DOUBLE(systime(/JULIAN,/UTC)-14-45) ; automatically recompute all spectra in the last 45 days



;Loop over all DSCOVR days
for i=start,ender-1 do begin
;get year of observation
    caldat,i,mon,day,year
;Get doy of observation
    doy = fix(i-JULDAY(01,01,year))+1
;Check if day current day is within auto recal regime unless reflag is set
    if ((i gt recom) and (reflag eq 0)) then clobber=1 else clobber=0

;Run dsc_advanced_kp_fit 
    dstr = string([year,doy],format=dfmt) 
    print,mess
    print,dstr

    case 1 of 
        ;if clobber set recompute spectral fitting and wind adjustment to spectra
        (clobber eq 1): dsc_advanced_kp_fit,year,doy,curv,/save,/clobber 
        ;if clobber not set and reflag not set run fitting on days without idl save files (skip previously computed days)
        ((clobber eq 0) and (reflag eq 0)):  dsc_advanced_kp_fit,year,doy,curv,/save
        ((clobber eq 0) and (reflag eq 1)):  begin
            idl_dscovr_to_cdf,year,doy,curv
            compare_wind_dscovr,year,doy,curv
        end



    endcase
    catch,Error_status

;Print the doy and year if anything goes wrong but
;continue looping through other days
    if Error_status ne 0 then begin
        estr = string([year,doy],format=efmt)
        print,mess
        print,estr
        print,mess
        catch,/cancel
    endif

    print,mess
endfor

set_plot,'X'

end
