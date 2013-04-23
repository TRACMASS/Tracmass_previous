SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel

  USE mod_time
  USE mod_grid
  USE mod_name
  
  
#ifdef tempsalt
USE mod_dens
#endif
IMPLICIT none

CHARACTER :: dates(62)*17
!CHARACTER(LEN=65) :: rfilu,rfilv,rfilh,rfilr

!REAL :: e1v(IMT+2,jmt),e1t(IMT+2,jmt),e2u(IMT+2,jmt),e2t(IMT+2,JMT)

INTEGER, SAVE :: nread,ndates
CHARACTER(LEN=65), SAVE :: rfilu,rfilv,rfilh,rfilr

integer i,ip,j,jp,k,kk,ints2,im

integer ncid !output ID index of netCDF file
integer ierr !error 
integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
integer varid,varidx,varidy,varidz,varidt !output ID index of variable
integer startA(1),startB(4),startC(2) !input index vector of position to start reading
integer countA(1),countB(4),countC(2) !input lengths of 'volume' to be retrieved
integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
integer p, x1, y1, z1, t1 !?

REAL*4, DIMENSION(KM) :: valsz
!REAL*4, ALLOCATABLE, DIMENSION(:) :: ssh
REAL*4, DIMENSION(IMT+2,JMT) :: ssh
REAL*4, DIMENSION(IMT+2,JMT,KM,1) :: fieldx,fieldy,fieldr
    
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: e1v,e1t,e2u,e2t
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dzu,dzv,dzt

logical around

 data dates/'19900702_19901001','19901001_19901231',                             & ! 1990
'19910101_19910402','19910402_19910702','19910702_19911001','19911001_19911231', & ! 1991
'19920101_19920402','19920402_19920702','19920702_19921001','19921001_19921231', & ! 1992
'19930101_19930402','19930402_19930702','19930702_19931001','19931001_19931231', & ! 1993
'19940101_19940402','19940402_19940702','19940702_19941001','19941001_19941231', & ! 1994
'19950101_19950402','19950402_19950702','19950702_19951001','19951001_19951231', & ! 1995
'19960101_19960402','19960402_19960702','19960702_19961001','19961001_19961231', & ! 1996
'19970101_19970402','19970402_19970702','19970702_19971001','19971001_19971231', & ! 1997
'19980101_19980402','19980402_19980702','19980702_19981001','19981001_19981231', & ! 1998
'19990101_19990402','19990402_19990702','19990702_19991001','19991001_19991231', & ! 1999
'20000101_20000402','20000402_20000702','20000702_20001001','20001001_20001231', & ! 2000
'20010101_20010402','20010402_20010702','20010702_20011001','20011001_20011231', & ! 2001
'20020101_20020402','20020402_20020702','20020702_20021001','20021001_20021231', & ! 2002
'20030101_20030402','20030402_20030702','20030702_20031001','20031001_20031231', & ! 2003
'20040101_20040402','20040402_20040702','20040702_20041001','20041001_20041231', & ! 2004
'20050101_20050402','20050402_20050702','20050702_20051001','20051001_20051231'/   ! 2005

!  if ( .NOT. ALLOCATED(valsz) ) then
!     allocate ( valsz(KM) )
!     allocate ( ssh(IMT+2,JMT) )
!     allocate ( fieldx(IMT+2,JMT,KM,1),fieldy(IMT+2,JMT,KM,1),fieldr(IMT+2,JMT,KM,1) )
!  end if

 if ( .not. allocated (e1v) ) then
   allocate ( e1v(IMT+2,JMT),e1t(IMT+2,JMT),e2u(IMT+2,JMT),e2t(IMT+2,JMT) )
   allocate ( dzu(IMT+2,JMT,KM),dzv(IMT+2,JMT,KM),dzt(IMT+2,JMT,KM) )
 end if

!_____________ swap between datasets ___________________________________

    hs(:,:,1)=hs(:,:,2)
    uflux(:,:,:,1)=uflux(:,:,:,2)
    vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1)=wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)=tem(:,:,:,2)
    sal(:,:,:,1)=sal(:,:,:,2)
    rho(:,:,:,1)=rho(:,:,:,2)
#endif

!_______________________ update the time counting ________________________________________
 iday=iday+5
 if(iday.gt.idmax(imon,iyear)) then
  iday=iday-idmax(imon,iyear)
  imon=imon+1
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
! if kan skrivas här om man vill börja om från iyear0
  endif
 endif
ntime=10000*iyear+100*imon+iday

!____________________________ initialise ___________________________
!print *,'ints=',ints,intstart
if(ints.eq.intstart) then
! call coordinat
 hs=0.
 uflux=0.
 vflux=0.
#ifdef tempsalt
 tem=0.
 sal=0.
 rho=0.
#endif
 ndates=0

!______________________________ Read ORCA grid horizontal ________________________
ierr=NF90_OPEN(trim(inDataDir)//'topo/mesh_hgr.nc',NF90_NOWRITE,ncid)
!ierr=NF_OPEN('/Users/doos/data/orc12/topo/mesh_hgr.nc',NF_NOWRITE,ncid)
if(ierr.ne.0) stop 3001

ierr=NF90_INQ_VARID(ncid,'nav_lev',varidz)
if(ierr.ne.0) stop 3087
startA(1)=1
countA(1)=km
ierr=NF90_GET_VAR(ncid,varidz,valsz,startA,countA)
if(ierr.ne.0) stop 3100
zw(0)=0.d0
!open(77,file='/Users/doos/data/orc12/topo/depthorca5')
do k=1,km
 zw(k)=dble(valsz(k))
 kk=km+1-k
 dz(kk)=zw(k)-zw(k-1) 
! print *,k,zw(k),kk,dz(kk)
!write(77,*) dz(kk)
enddo

startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=1
countB(1)=IMT+2
countB(2)=jmt
countB(3)=1
countB(4)=1
ierr=NF90_INQ_VARID(ncid,'e1t',varid)
if(ierr.ne.0) stop 3002
ierr=NF90_GET_VAR(ncid,varid,e1t,startB,countB)
if(ierr.ne.0) stop 3005
!print *,(e1t(i,jmt-2),i=300,310)

ierr=NF90_INQ_VARID(ncid,'e2u',varid)
if(ierr.ne.0) stop 3010
ierr=NF90_GET_VAR(ncid,varid,e2u,startB,countB)
if(ierr.ne.0) stop 3015
!do j=1,jmt
!print *,j,e2u(1,j)
!enddo
!stop 4956

ierr=NF90_INQ_VARID(ncid,'e2t',varid)
if(ierr.ne.0) stop 3020
ierr=NF90_GET_VAR(ncid,varid,e2t,startB,countB)
if(ierr.ne.0) stop 3025
!print *,(e2t(i,jmt-2),i=300,310)

ierr=NF90_INQ_VARID(ncid,'e1v',varid)
if(ierr.ne.0) stop 3030
ierr=NF90_GET_VAR(ncid,varid,e1v,startB,countB)
if(ierr.ne.0) stop 3035
!print *,(e1v(i,jmt-2),i=300,310)
!do j=1,jmt
!print *,j,e1v(1,j)/deg
!enddo
!stop 4568

! Close the file up
ierr=NF90_CLOSE(ncid)
if(ierr.ne.0) stop 3040

do i=1,IMT
 do j=1,jmt
  dxdy(i,j)=e1t(i,j)*e2t(i,j)
 enddo
enddo

!______________________________ Read ORCA grid horizontal ________________________
ierr=NF90_OPEN(trim(inDataDir)//'topo/mesh_zgr.nc',NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 4001
! dzt
ierr=NF90_INQ_VARID(ncid,'e3t_ps',varid)
if(ierr.ne.0) stop 4004
startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=1
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
if(ierr.ne.0) stop 4002
!do i=200,250
!print *,i,(fieldx(i,1,k,1),k=16,35)
!enddo
!print *,(dz(km+1-k),k=1,km)
dzt(:,:,:)=fieldx(:,:,:,1)
!print *,'dzt',(dzt(200,200,k),k=1,km)

! dzu
ierr=NF90_INQ_VARID(ncid,'e3u_ps',varid)
if(ierr.ne.0) stop 4014
startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=1
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
if(ierr.ne.0) stop 4012
dzu(:,:,:)=fieldx(:,:,:,1)
!print *,'dzu',(dzu(200,200,k),k=1,km)

! dzv
ierr=NF90_INQ_VARID(ncid,'e3v_ps',varid)
if(ierr.ne.0) stop 4024
startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=1
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
if(ierr.ne.0) stop 4032
dzv(:,:,:)=fieldx(:,:,:,1)
!print *,'dzv',(dzv(200,200,k),k=1,km)

ierr=NF90_CLOSE(ncid)
if(ierr.ne.0) stop 4031

!stop 4596

endif !____________________________________ endif initstart _____________________________

if(mod(ints,18).eq.1) then

ints2=ints
666  continue
if(ints2.gt.intmax .or. ints2.lt.intmin) then
 ints2=ints2-intend+intstart-intstep
 goto 666
endif

 if(iday0.eq.2-10 .and. imon0.eq.7 .and. iyear0.eq.1990) then
   ndates=ints2/18+1
 elseif(iday0.eq.1-10 .and. imon0.eq.10 .and. iyear0.eq.1992) then
  ndates=ints2/18+10
 else
  print *,iyear0,imon0,iday0,ndates
  print *,iyear,imon,iday,ndates
  stop 9567
 endif

 rfilu=trim(inDataDir)//'gcm/KAB042j_5d_'//dates(ndates)//'_grid_U.nc'
 rfilv=trim(inDataDir)//'gcm/KAB042j_5d_'//dates(ndates)//'_grid_V.nc'
 rfilr=trim(inDataDir)//'gcm/KAB042j_5d_'//dates(ndates)//'_sigma.nc'
 rfilh=trim(inDataDir)//'gcm/KAB042j_5d_'//dates(ndates)//'_SSH.nc'
 print *,'rfilu=',rfilu,ints,ints2,ndates,dates(ndates)
 inquire(file=rfilu,exist=around)
 if(.not.around) stop 4556
 inquire(file=rfilv,exist=around)
 if(.not.around) stop 4557
 inquire(file=rfilr,exist=around)
 if(.not.around) stop 4558
 inquire(file=rfilh,exist=around)
 if(.not.around) stop 4559
endif
!nread=mod(ints/5,18)+1
nread=mod(ints,18)+1

! zonal velocity
ierr=NF90_OPEN(rfilu,NF90_NOWRITE,ncid)
!print *,'rfilu=',rfilu
if(ierr.ne.0) stop 3750


ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) ! the main data fields
if(ierr.ne.0) stop 3762

startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=nread
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
if(ierr.ne.0) stop 3798

!print *,'read fields'
!print *,(fieldx(i,100,1,1),i=1,IMT+2/3)

ierr=NF90_CLOSE(ncid)


! meridional velocity

ierr=NF90_OPEN(rfilv,NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3751

ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! the main data fields
if(ierr.ne.0) stop 3763

startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=nread
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldy,startB,countB)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

! density

ierr=NF90_OPEN(rfilr,NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3758

ierr=NF90_INQ_VARID(ncid,'sigma',varid) ! the main data fields
if(ierr.ne.0) stop 3767

startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=nread
countB(1)=IMT+2
countB(2)=jmt
countB(3)=km
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,fieldr,startB,countB)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

! sea surface height

ierr=NF90_OPEN(rfilh,NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3759

ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) stop 3763

startB(1)=1
startB(2)=1
startB(3)=1
startB(4)=nread
countB(1)=IMT+2
countB(2)=jmt
countB(3)=1
countB(4)=1
ierr=NF90_GET_VAR(ncid,varid,ssh,startB,countB)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

do j=1,jmt
 do i=1,IMT
  hs(i,j,2)=0.01*ssh(i,j)
 enddo
enddo

!print *,(ssh(i,250),i=1,100)
!do j=1,jmt
!print *,j,ssh(1,j),ssh(2,j),ssh(IMT+2-1,j),ssh(IMT+2,j)
!enddo
!stop 4568

do j=1,jmt
 jp=j+1
 if(jp.gt.jmt) jp=jmt
 do i=1,IMT
  ip=i+1
  if(ip.eq.IMT+1) ip=1
  do k=1,km-1
   kk=km+1-k
   uflux(i,j,k,2)=fieldx(i,j,kk,1)*e2u(i,j)*dzu(i,j,kk)
   vflux(i,j,k,2)=fieldy(i,j,kk,1)*e1v(i,j)*dzv(i,j,kk)
   rho(i,j,k,2)=fieldr(i,j,kk,1)
  enddo
   uflux(i,j,km,2)=fieldx(i,j,1,1)*e2u(i,j)*(dzu(i,j,1)+0.5*(hs(i,j,2)+hs(ip,j,2)))
   vflux(i,j,km,2)=fieldy(i,j,1,1)*e1v(i,j)*(dzv(i,j,1)+0.5*(hs(i,j,2)+hs(i,jp,2)))
   rho(i,j,km,2)=fieldr(i,j,1,1)
!   if(i.eq.210 .and. j.eq.252) print *,'density=',rho(i,j,km,2)
 enddo
enddo

if(ints.eq.intstart) then
do j=1,jmt
 do i=1,IMT
  do k=1,km
   kk=km+1-k
   if(fieldr(i,j,k,1).ne.0.) kmt(i,j)=k
   if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
  enddo
  if(kmt(i,j).ne.0) then
   dztb(i,j,1)=dzt(i,j,kmt(i,j))
  else
   dztb(i,j,1)=0.
  endif
 enddo
enddo

! testa kmt

!stop 32456

!print *,dz
!print *,'eastern indian'
!do j=260,200,-1
!print 679,j,(kmt(i,j),i=1,50)
!enddo
!print *,'western pacific'
!do j=260,200,-1
!print 679,j,(kmt(i,j),i=50,250)
!enddo
!print *,'eastern pacific'
!do j=JMT,250,-1
!print 679,j,(kmt(i,j),i=251,450)
!enddo
!print *,'atlantic'
!do j=260,250,-1
!print 679,j,(kmt(i,j),i=451,641)
!enddo
!print *,'western indian'
!do j=JMT,250,-1
!print 679,j,(kmt(i,j),i=641,imt)
!enddo


679 format(i3,1x,400i1)
!stop 4966

!do j=292,289,-1
!print 679,j,(kmt(i,j),i=139,143)
!679 format(400i3)
!enddo

endif

do j=2,jmt
 do i=1,IMT
  im=i-1
  if(im.eq.0) im=IMT
  do k=1,km
   if(uflux(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'u',i,j,k,kmt(i,j),uflux(i,j,k,2)
   if(vflux(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'v',i,j,k,kmt(i,j),vflux(i,j,k,2)
   if(uflux(im,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'u',im,j,k,kmt(i,j),uflux(im,j,k,2)
   if(vflux(i ,j-1,k,2).ne.0..and. kmt(i,j).eq.0) print *,'v',i,j-1,k,kmt(i,j),vflux(i,j-1,k,2)
  enddo
 enddo
enddo


return
end subroutine readfields

!_______________________________________________________________________
