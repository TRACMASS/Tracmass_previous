!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
 
SUBROUTINE readfields

IMPLICIT none

#include "/Applications/Utilities/netcdf-3.6.0-p1/include/netcdf.inc"
!#include "/sw/include/netcdf.inc"
!#include "netcdf.inc"
#include "../param.h"

INTEGER IMU,JMV,IT
PARAMETER(IMU=IMT-1,JMV=JMT-1,IT=10)

common/vel/u(IMT,0:JMAX,KM,NST),v(IMT,0:JMAX,KM,NST),hs(IMT,JMAX,NST),w(0:KM),ff
real u,v,hs
real*8 w,ff

common /coord/dx,dy,deg,stlon1,stlat1,csu(JMT),cst(JMT),zw(0:KM)
REAL*8        dx,dy,deg,stlon1,stlat1,csu,cst,zw

common /grid/dxdy(IMT,JMT),dztb(IMT,JMT,KD),dz(KM),rmin,dr,tmin,dtemp,smin,dsalt,kmt(IMT,JMT)
REAL*8 dxdy,dztb,dz,rmin,dr,tmin,dtemp,smin,dsalt
INTEGER kmt

#ifdef tempsalt
common/dens/ tem(IMT,JMT,KM,2),sal(IMT,JMT,KM,2),rho(IMT,JMT,KM,2)
real tem,sal,rho
real tempb(KM),saltb(KM),rhob(KM)
integer kmm
#endif

common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
integer    ints,intstart,intend,intrun,intspin,intstep,intmin,intmax

common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=27) :: directory

CHARACTER :: dirdata*23,dates(62)*17
CHARACTER(LEN=63) :: rfil
!data dirdata/'/Volumes/sjo5/data/tun/'/
 data dirdata/'/Users/doos/data/abcde/'/

integer i,im,j,jm,k,kk,nread,ndates,ints2

REAL LAT_U(IMU,JMT),LON_U(IMU,JMT)
REAL LAT_V(IMt,JMv),LON_V(IMt,JMv)  
REAL LAT_r(IMt,JMT),LON_r(IMT,JMT)  

INTEGER MAX_IT,VAR_ID
REAL ZVR(IMT,JMV,KM),ZUR(IMU,JMT,KM),ZRR(IMT,JMT,KM),TIME(It)

INTEGER MASK_U(IMU,JMT),MASK_V(IMT,JMV),MASK_RHO(IMT,JMT)      	

integer ncid !output ID index of netCDF file
integer ierr !error 
integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
integer varid,varidx,varidy,varidz,varidt !output ID index of variable
integer startA(2),startB(3),startC(4) !input index vector of position to start reading
integer countA(2),countB(3),countC(4) !input lengths of 'volume' to be retrieved
integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
integer p, x1, y1, z1, t1 !?
    
real fieldh(IMT,JMT,1),fieldx(IMU,JMT,KM,1),fieldy(IMT,JMV,KM,1),fieldr(IMT,JMT,KM,1)
real dzu(IMT,JMT,KM),dzv(IMT,JMT,KM),dxv,dyu(JMT),pi

logical around

save nread,rfil,ndates,dzu,dzv,dyu,dxv

!_____________ swap between datasets ___________________________________

hs(:,:,1)=hs(:,:,2)
u(:,:,:,1)=u(:,:,:,2)
v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
tem(:,:,:,1)=tem(:,:,:,2)
sal(:,:,:,1)=sal(:,:,:,2)
rho(:,:,:,1)=rho(:,:,:,2)
#endif

!____________________________ initialise ___________________________
if(ints.eq.intstart) then
 pi=2.*asin(1.)
 call coordinat
 hs=0.
 u=0.
 v=0.
#ifdef tempsalt
 tem=0.
 sal=0.
 rho=0.
#endif
 ndates=0

!______________________________ Read ROMS grid horizontal ________________________
ierr=NF_OPEN(dirdata//'test_grd2.nc',NF_NOWRITE,ncid)
print *,dirdata//'test_grd2.nc'
if(ierr.ne.0) stop 3001

! Lire les longitudes:
ierr=nf_inq_varid(NCID,'lon_u',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LON_U)
ierr=nf_inq_varid(NCID,'lon_v',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LON_V)      
ierr=nf_inq_varid(NCID,'lon_rho',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LON_r)     
! Lire les latitudes:
ierr=nf_inq_varid(NCID,'lat_u',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LAT_U)
ierr=nf_inq_varid(NCID,'lat_v',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LAT_V)         
ierr=nf_inq_varid(NCID,'lat_rho',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,LAT_r) 

do j=2,jmt
 dyu(j)=Lat_v(1,j)-Lat_v(1,j-1)
! print *,j,Lat_v(1,j),dyu(j)
enddo
dyu(1)=dyu(2)

dxv=Lon_v(2,1)-Lon_v(1,1)

 print *,dxv
  
do i=1,imt
 do j=1,jmt
  dxdy(i,j)=dxv*cos(pi/180.*lat_v(i,j))*dyu(j)*deg**2 
 enddo
enddo

! Lire les masks:
ierr=nf_inq_varid(NCID,'mask_u',VAR_ID)
ierr=nf_get_var_int(NCID,VAR_ID,MASK_U)
ierr=nf_inq_varid(NCID,'mask_rho',VAR_ID)
ierr=nf_get_var_int(NCID,VAR_ID,MASK_RHO)      
ierr=nf_inq_varid(NCID,'mask_v',VAR_ID)
ierr=nf_get_var_int(NCID,VAR_ID,MASK_V)
! lire zu_r
ierr=nf_inq_varid(NCID,'zu_r',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,ZUR)	
! lire zv_r
ierr=nf_inq_varid(NCID,'zv_r',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,ZVR)
! lire zr_r
ierr=nf_inq_varid(NCID,'zr_r',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,ZRR)
! Close the file up
ierr=NF_CLOSE(ncid)
if(ierr.ne.0) stop 3040

do k=1,KM
 kk=KM+1-k

 do j=1,JMT
  do i=1,IMU
   dzu(i+1,j,kk)=zur(i,j,k)*dyu(j)*deg
  enddo
  dzu(1,j,kk)=dzu(2,j,k)
 enddo

 do i=1,IMT
  do j=1,JMV
   dzv(i,j+1,kk)=zvr(i,j,k)*dxv*cos(pi/180.*lat_v(i,j))*deg
  enddo
  dzv(i,1,kk)=dzv(i,1,k)
 enddo

 do i=1,IMT
  do j=1,JMT
   dztb(i,j,kk)=zrr(i,j,k)
  enddo
 enddo

enddo


endif !____________________________________ endif initstart _____________________________

!print *,'nread1=',nread

if(mod(ints,10).eq.0) then

ints2=ints
666  continue
if(ints2.gt.intmax .or. ints2.lt.intmin) then
 ints2=ints2-intend+intstart-intstep
 goto 666
endif

if(mod(ints,360).eq.0) then
print *,'psi written for ints=',ints,ints2
call writepsi
endif

! ndates=ndates+1
 ndates=ints2/90+1
 rfil=dirdata//'roms_hisM2_0100.nc'
! rfilv=dirdata//'KAB042j_5d_'//dates(ndates)//'_grid_V.nc'
 print *,'rfil=',ints,ints2,ndates,ints2/90+1,dates(ndates)
 inquire(file=rfil,exist=around)
 if(.not.around) stop 4556
endif
ierr=NF_OPEN(rfil,NF_NOWRITE,ncid)
if(ierr.ne.0) stop 3750

nread=mod(ints,10)+1
print *,'nread=',nread,ints

! Lire le temps:
ierr=nf_inq_varid(NCID,'ocean_time',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,TIME)

print *,'ocean_time',TIME
!print *,((time(n)-time(1))/60.,n=1,it)

! Lire zeta:
ierr=nf_inq_varid(NCID,'zeta',VAR_ID)

startB(1)=1
startB(2)=1
startB(3)=nread
countB(1)=IMT
countB(2)=JMT
countB(3)=1
ierr=NF_GET_VARA_REAL(ncid,varid,startB,countB,fieldh)
!ierr=NF_CLOSE(ncid)

!ierr=nf_get_var_real(NCID,VAR_ID,ZETA)
print *,NCID,'zeta',VAR_ID
      
! Lire U
ierr=nf_inq_varid(NCID,'u',VAR_ID)
print *,NCID,'u',VAR_ID
startC(1)=1
startC(2)=1
startC(3)=1
startC(4)=nread
countC(1)=IMU
countC(2)=JMT
countC(3)=KM
countC(4)=1
ierr=NF_GET_VARA_REAL(ncid,varid,startC,countC,fieldx)
print *,NCID,'u',VAR_ID,ierr
if(ierr.ne.0) stop 3798


! Lire V
ierr=nf_inq_varid(NCID,'v',VAR_ID)
ierr=nf_get_var_real(NCID,VAR_ID,V)
print *,NCID,'v',VAR_ID



ierr=NF_INQ_VARID(ncid,'vozocrtx',varid) ! the main data fields
if(ierr.ne.0) stop 3762


!print *,'read fields'
!print *,(fieldx(i,100,1,1),i=1,IMT/3)

ierr=NF_CLOSE(ncid)



do j=1,JMT
 do i=1,IMT
  hs(i,j,2)=0.01*fieldh(i,j,1)
 enddo
enddo

!print *,(ssh(i,250),i=1,100)

do j=1,JMT
 jm=j-1
 if(jm.eq.0) jm=1
 do i=1,IMT
  im=i-1
  if(im.eq.0) im=IMT
  do k=1,KM-1
   kk=KM+1-k
   u(i,j,k,2)=fieldx(i,j,kk,1)*dzu(i,j,k)
   v(i,j,k,2)=fieldy(i,j,kk,1)*dzv(i,j,k)
   rho(i,j,k,2)=fieldr(i,j,kk,1)
  enddo
   u(i,j,KM,2)=fieldx(i,j,1,1)*(dzu(i,j,KM)+0.5*(hs(i,j,2)+hs(im,j,2)))
   v(i,j,KM,2)=fieldy(i,j,1,1)*(dzv(i,j,KM)+0.5*(hs(i,j,2)+hs(i,jm,2)))
   rho(i,j,KM,2)=fieldr(i,j,1,1)
!   if(i.eq.210 .and. j.eq.252) print *,'density=',rho(i,j,KM,2)
 enddo
enddo


return
end subroutine readfields

!_______________________________________________________________________
