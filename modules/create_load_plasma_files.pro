
;===================================================
;
;
;
;Uses
;    IF not ran from dsc_advanced_kp_fit_v2 then the follow contingencies are required
;    following 
;    @'/crater/utilities/idl/mike/idlstartup' 
;    @compile_cdaweb
;    @compile_IDLmakecdf
;
;ROUTINE for creating IDL sav file for load plasma which integrates plasma and magnetic info
;
;USAGE
;create_load_plasma_files,year,doy,marchive=marchive,parchive=parchive,tarchive=tarchive,mfilefmt=mfilefmt,pfilefmt=pfilefmt,tfilefmt=tfilefmt,outpath=outpath,outfmt=outfmt
;
;===================================================


;===================================================
; 
; Create IDL structure containing DSCOVR plasma and magnetic field information
;
;USAGE
;plasmag_struct
;
;
;===================================================

;pro plasmag_struct
;
;tmp = {plasmag_dat, $
;     
;
;
;end


;===================================================
;
;Return newest version of file
;
;USAGE
;new_file = new_ver(archive,file)
;
;===================================================
function new_ver,archive,file


new_file = file_search(archive+file)
if n_elements(size(new_file)) gt 3 then new_file = new_file[n_elements(new_file)-1]


return,new_file
end

;===================================================
;MAIN loop
;===================================================
pro create_load_plasma_files,year,doy,marchive=marchive,parchive=parchive,tarchive=tarchive,mfilefmt=mfilefmt,pfilefmt=pfilefmt,tfilefmt=tfilefmt,outpath=outpath,outfmt=outfmt

if keyword_set(parchive) then parchive=parchive else parchive='/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected';set output cdf location
parchive = parchive+'/'

if keyword_set(marchive) then marchive=marchive else marchive='/crater/observatories/dscovr/plasmag/mag';set output cdf location
marchive = marchive+'/'

if keyword_set(tarchive) then tarchive=tarchive else tarchive='/crater/observatories/dscovr/traj/pre_or';set output cdf location
tarchive = tarchive+'/'

if keyword_set(pfilefmt) then pfilefmt=pfilefmt else pfilefmt = '("dscovr_h1_fc_",I4,I02,I02,"_v*.cdf")' 
if keyword_set(mfilefmt) then mfilefmt=mfilefmt else mfilefmt = '(I4,"/dscovr_h0_mag_",I4,I02,I02,"_v*.cdf")'
if keyword_set(tfilefmt) then tfilefmt=tfilefmt else tfilefmt = '(I4,"/dscovr_orbit_pre_",I4,I02,I02,"_v*.cdf")'


if keyword_set(outpath) then outpath = outpath else outpath = "/crater/observatories/dscovr/plasmag/merged/plsmag"
outpath = outpath+"/"
if keyword_set(outfmt) then outfmt = outfmt else outfmt = '("dsc.plsmag.",I4,"_",I03,".",I03,".idl")'



;GET JULDAY FROM DOY AND YEAR

obs_day = double(JULDAY(1,1,year,0,0,0))
obs_day = double(obs_day+doy-1) ;convert day of year into Julian day

;OUTPUT observed day in format in YYYYMMDD format
CALDAT,obs_day,omon,oday,oyear

pcdf = string([oyear,omon,oday],format=pfilefmt)
mcdf = string([oyear,oyear,omon,oday],format=mfilefmt)
tcdf = string([oyear,oyear,omon,oday],format=tfilefmt)


magf = new_ver(marchive,mcdf)
plsf = new_ver(parchive,pcdf)
orbf = new_ver(tarchive,tcdf)

;read in cdf files for given day
print,plsf
mag = read_mycdf("",magf,/all)
pls = read_mycdf("",plsf,/all)
orb = read_mycdf("",orbf,/all)


;Set up arrays for plasma. 
vx = pls.v_gse.dat[0,*]
vy = pls.v_gse.dat[1,*]
vz = pls.v_gse.dat[2,*]
wd = pls.THERMAL_SPD.dat
np = pls.np.dat
dq = pls.dqf.dat


;find flagged values
rvx = where((vx le -9999.0) or (dq lt 0))
rvy = where((vy le -9999.0) or (dq lt 0))
rvz = where((vz le -9999.0) or (dq lt 0))
rwd = where((wd le -9999.0) or (dq lt 0))
rnp = where((np le -9999.0) or (dq lt 0))

;Replace -1.e30 with -9999.0
if n_elements(size(rvx)) gt 3 then vx[rvx] = -9999.0
if n_elements(size(rvy)) gt 3 then vy[rvy] = -9999.0
if n_elements(size(rvz)) gt 3 then vz[rvz] = -9999.0
if n_elements(size(rwd)) gt 3 then wd[rwd] = -9999.0
if n_elements(size(rnp)) gt 3 then np[rnp] = -9999.0

jdp = CDF_EPOCH_TOJULDAYS(pls.epoch.dat)
jdm = CDF_EPOCH_TOJULDAYS(mag.epoch1.dat)
jdo = CDF_EPOCH_TOJULDAYS(orb.epoch.dat)

;Interpolate mag. data
bx = interpol(mag.b1gse.dat[0,*], jdm ,jdp)
by = interpol(mag.b1gse.dat[1,*], jdm, jdp)
bz = interpol(mag.b1gse.dat[2,*], jdm, jdp)
bm = interpol(mag.b1f1.dat, jdm, jdp)

;Interpolate pos. data
RE = 6371; Earth radius in km, for conversion of trajectories
x = interpol(orb.gse_pos.dat[0,*], jdo, jdp)/RE 
y = interpol(orb.gse_pos.dat[1,*], jdo, jdp)/RE
z = interpol(orb.gse_pos.dat[2,*], jdo, jdp)/RE

;convert epoch into doy fraction
ddoy = jdp - JULDAY(1,1,year,0,0,0)+1.

DOY_DSC   = ddoy
BX_DSC    = bx
BY_DSC    = by
BZ_DSC    = bz
B_DSC     = bm
VX_DSC    = vx
VY_DSC    = vy
VZ_DSC    = vz
V_DSC     = sqrt(vx^2.+vy^2.+vz^2.)
EW_DSC    =  fltarr(n_elements(bx))-9999.0 
NS_DSC    =  fltarr(n_elements(bx))-9999.0 
W_DSC     =  wd
N_DSC     =  np
X_DSC     =  x 
Y_DSC     =  y 
Z_DSC     =  z 
HEOVH_DSC =  fltarr(n_elements(bx))-9999.0 ;He/H


save,DOY_DSC,BX_DSC,BY_DSC,BZ_DSC,B_DSC,VX_DSC,VY_DSC,VZ_DSC,V_DSC,EW_DSC,NS_DSC,EW_DSC,NS_DSC,W_DSC,N_DSC,X_DSC,Y_DSC,Z_DSC,HEOVH_DSC,filename=outpath+string([year,doy,doy+1],format=outfmt)





end
