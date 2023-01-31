; use NCL to calculate atmosheric zonal stream function
; Yiling: modified from Eli's
;
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"		; functions required to
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"           	; plot.  include before
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"         	; begin
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
begin

PI=3.1415926

; read in data
;filename="mH.tier1.CN.011.6hrly.cam2.h0.050115_FSDS_PRECT_PS_Q_QREFHT_SOLIN_T_TREFHT_TS_U_V.nc"
filename="mH.tier1.CN.002.6hrly.cam2.h0.120115_FSDS_PRECT_PS_Q_QREFHT_SOLIN_T_TREFHT_TS_U_V.nc"
;filename="cesmpifv1vdpCMIP6.002.6hrly.cam2.h0.470115_FSDS_PRECT_PS_Q_QREFHT_SOLIN_T_TREFHT_TS_U_V.nc"
f=addfile(filename,"r");
lat=f->lat
lon=f->lon
time=f->time
;pres_=f->lev
; vinth2p requires the lev_p to be expressed in mb [hPa]
lev_p           = (/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000/)         
lev_p!0         = "lev_p"                  ; variable/dim name 
lev_p&lev_p     =  lev_p                   ; create coordinate variable
lev_p@long_name = "pressure"               ; attach some attributes
lev_p@units     = "hPa"
lev_p@positive  = "down"

hyam            = f->hyam         ; read from file
hybm            = f->hybm
P0mb            = f->P0*0.01
ps              = f->PS           ; these are in Pa
u               = f->U
v               = f->V
  
dv    = new ( dimsizes(u), typeof(u))  	; divergence
ud    = new ( dimsizes(dv),typeof(dv))	; zonal divergent wind 
vd    = new ( dimsizes(dv),typeof(dv))	; meridional divergent wind
  
uv2dvg   (u,v,dv)          ; u,v ==> divergence
dv2uvg   (dv,ud,vd)        ; dv  ==> divergent wind components
; convert to pressure levels
up   = vinth2p(ud,hyam,hybm,lev_p,ps,1,P0mb,1,False)
indlat=ind(lat.gt.-5 .and. lat.lt.5) ; tropical from 5S to 5N
lat_=lat(indlat)
ntim=dimsizes(time)
nlev=dimsizes(lev_p)
nlon=dimsizes(lon)
up_=up(:,:,indlat,:)
ps_=ps(:,indlat,:)
; calculate meridional zonal streamfunction
uvmsg = 1e+36
dsizes_x = dimsizes(up_)
nlat     = dsizes_x(2)
up_T    = new ( (/ntim,nlev,nlon,nlat /), float, uvmsg )  ; transpose up
ps_T    = new ( (/ntim,nlon,nlat /), float, uvmsg )  ; transpose ps
up_T1   = new ( (/ntim,nlev,nlon,nlat /), double, uvmsg )  ; transpose up_T*coslat
lon_fake= lon*0  ; fake lon for zonal_mpsi with cos(lon)=1
do i=0,ntim-1,1
do j=0,nlev-1,1
up_T(i,j,:,:)=transpose(up_(i,j,:,:))
end do
ps_T(i,:,:)=transpose(ps_(i,:,:))
end do
do i=0,nlat-1,1
up_T1(:,:,:,i)=up_T(:,:,:,i)*cos(lat_(i)/180*PI)
end do
;;; Calculate zonal circulation
mzpsi= zonal_mpsi(up_T1,lon_fake,lev_p*100.,ps_T)
mzpsi!0="time"
mzpsi&time=time
mzpsi!1="lev"
mzpsi&lev=lev_p
mzpsi&lev@units="mb"
mzpsi!2="lon"
mzpsi&lon=lon
mzpsi@units="kg/s"
mzpsi@long_name="zonal Streamfunction"
mzpsi=(/mzpsi/1e11/)
mzpsi@units="svp"

;***********************************
;output netcdf data
system("rm -f mH.tier1.CN.002.6hrly.cam2.h0.120115_mzpsi.nc")
ncdf= addfile("mH.tier1.CN.002.6hrly.cam2.h0.120115_mzpsi.nc" ,"c")
ncdf->lon=lon
ncdf->time=time
up!0    = "time"
up!1    = "lev"
up!2    = "lat"
up!3    = "lon"
up&time= time
up&lat= lat
up&lev= lev_p
up&lon= lon
up@long_name= "Pressure Level U"
up@units="m/s"
ncdf->mzpsi=mzpsi
ncdf->U_lev=up
end
