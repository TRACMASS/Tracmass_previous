SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  use mod_seed

#ifdef tempsalt
  USE mod_dens
#endif
  IMPLICIT none
  
  INTEGER, PARAMETER :: NTID=73
#ifdef orca1
  INTEGER, PARAMETER :: IJKMAX2=254
#else
  INTEGER, PARAMETER :: IJKMAX2=281
#endif

  ! = Loop variables
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm,ipp,jmm,jpp
  INTEGER, SAVE                                :: n  

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
!  CHARACTER (len=50)                         :: varName
!  INTEGER                                   :: start1d, count1d
 ! INTEGER, DIMENSION(2)                     :: start2d, count2d
 ! INTEGER, DIMENSION(3)                     :: start3d, count3d
  
    ! = Variables for filename generation
!  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix !, dstamp
!  INTEGER                                    :: intpart1 ,intpart2 ,subYr
!  INTEGER                                    :: filePos ,fileJD ,subYrD
!  INTEGER                                    :: yr1 ,mn1 ,dy1
 ! INTEGER                                    :: yr2 ,mn2 ,dy2
  

  REAL*4, ALLOCATABLE, DIMENSION(:,:)         :: ssh
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp
  REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: temp2d_doub
  
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: botbox
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:)   :: e1t,e2t
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: kmu,kmv
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: ntimask
  REAL*4 :: temp2d_simp(IMT+2,JMT)
  INTEGER itemp(IMT+2,JMT)
  INTEGER :: varid, ncid
  REAL*4 dd,dmult

LOGICAL around

!REAL*8 temp3d_doub(IMT+2,JMT,KM) ! to be used when botbox is written

  
  alloCondGrid: if ( .not. allocated (botbox) ) then
     allocate (  botbox(IMT,JMT,3) )
     allocate ( kmu(IMT,JMT)    ,kmv(IMT,JMT) )
     allocate ( ntimask(NTID,IJKMAX) )
  end if alloCondGrid

  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt),temp3d_simp(IMT+2,JMT,KM), temp2d_doub(IMT+2,JMT)   )
     allocate ( e1t(IMT+2,JMT) , e2t(IMT+2,JMT) )
  end if alloCondUVW
    
  start1d  = [ 1]
  count1d  = [km]
  start2d  = [  1   ,   1,    1,    1]
  count2d  = [imt+2 , jmt,    1,    1]
  start3d  = [  1   ,   1,    1,    1]
  count3d  = [imt+2 , jmt,   km,    1]
  
    ! === swap between datasets ===
  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif

  ! === initialise ===
  initFieldcond: if(ints.eq.intstart) then
     ! call coordinat
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     
gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'

ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3751
ierr=NF90_INQ_VARID(ncid,'e1t',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e1t(i,:)=temp2d_doub(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2u',varid) ! dyu in meters
if(ierr.ne.0) stop 3764
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  dyu(i,:)=temp2d_doub(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2t',varid) ! the main data fields
if(ierr.ne.0) stop 3765
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e2t(i,:)=temp2d_doub(i,:)
  dxdy(i,:) = e1t(i,:) * e2t(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! the main data fields
if(ierr.ne.0) stop 3766
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  dxv(i,:)=temp2d_doub(i,:)
enddo
ierr=NF90_CLOSE(ncid)

! The bathymetry
gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
if(ierr.ne.0) stop 3767
ierr=NF90_GET_VAR(ncid,varid,itemp,start2d,count2d)

do i=1,IMT
  kmt(i,:)=itemp(i,:)
enddo

kmu=0 ; kmv=0

do j=1,jmt
jp=j+1
if(jp.eq.jmt+1) jp=jmt  ! should be north fold instead
 do i=1,imt
   kmu(i,j)=min(itemp(i,j),itemp(i+1,j))
   kmv(i,j)=min(itemp(i,j),itemp(i,jp))
 enddo
enddo

!botbox=0.
!! layer thickness at u points
!ierr=NF90_INQ_VARID(ncid,'e3u',varid) 
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmu(i,j).ne.0) botbox(i,j,1)=temp3d_doub(i,j,kmu(i,j))
!enddo
!enddo
!
!! layer thickness at v points
!ierr=NF90_INQ_VARID(ncid,'e3v',varid)
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmv(i,j).ne.0) botbox(i,j,2)=temp3d_doub(i,j,kmv(i,j))
!enddo
!enddo
!
!! layer thickness at T points
!ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmt(i,j).ne.0) botbox(i,j,3)=temp3d_doub(i,j,kmt(i,j))
!enddo
!enddo


open(21,file=trim(inDataDir)//'topo/botbox',form='unformatted')
!write(21) botbox
read(21) botbox
close(21)


!ierr=NF90_INQ_VARID(ncid,'nav_lon',varid) ! londitude mesh
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!open(21,file=trim(inDataDir)//'topo/long',form='unformatted')
!write(21) temp2d_simp
!close(21)
!
!ierr=NF90_INQ_VARID(ncid,'nav_lat',varid) ! londitude mesh
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!open(21,file=trim(inDataDir)//'topo/lat',form='unformatted')
!write(21) temp2d_simp
!close(21)


ierr=NF90_CLOSE(ncid)



!open(21,file=trim(inDataDir)//'topo/kmt',form='unformatted')
!write(21)kmt
!close(21)
!stop 4985


! Depth coordinates
zw(0)    = 0.d0
zw(1:km) = get1DfieldNC(trim(gridFile) ,'e3t_0')
do k=1,km
 kk=km+1-k
 dz(kk)=zw(k)
 zw(k)=zw(k)+zw(k-1)
! print *,k,zw(k),kk,dz(kk)
end do

! Bottom box 
do j=1,JMT
 do i=1,IMT
  if(kmt(i,j).ne.0) then
   dztb(i,j,1)=botbox(i,j,3)
  else
   dztb(i,j,1)=0.
  endif
 enddo
enddo

#ifdef initxyt
! Time for individual start positions
if(ijkmax.ne.IJKMAX2) then
 print *,ijkmax,IJKMAX2
 stop 4396
endif
open(84,file=trim(inDataDir)//'topo/masktime005b',form='unformatted')
 read(84) ntimask
close(84)
#endif

iday=startDay-5
imon=startMon
iyear=startYear
     
endif initFieldcond

! date update
  iday=iday+5
  if(iday.gt.idmax(imon,1999)) then
     iday=iday-idmax(imon,1999)
     imon=imon+1
     if(imon.eq.13) then
        imon=1
        iyear=iyear+1
        ! if kan skrivas här om man vill börja om från iyear0
     endif
  endif
!iyear=2000 ! quick and dirty fix for years
ntime=10000*iyear+100*imon+iday
n=n+1 ! quick and dirty fix to be generalised

! file names
#ifdef orca1
 dataprefix='xxxx/ORCA1-N202_xxxxxxxx'
write(dataprefix(17:24),'(i8)') ntime
write(dataprefix(1:4),'(i4)') iyear
#else
 dataprefix='xxxx/ORCA025-N112_xxxxxxxx'
write(dataprefix(19:26),'(i8)') ntime
write(dataprefix(1:4),'(i4)') iyear
#endif

  fieldFile = trim(inDataDir)//trim(dataprefix)
    
! temp, salt and ssh
!temp2d_simp = get2DfieldNC(trim(fieldFile)//'d05T.nc' ,'sossheig')
!print *,trim(fieldFile)//'d05T.nc'

ierr=NF90_OPEN(trim(fieldFile)//'d05T.nc',NF90_NOWRITE,ncid)
ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) stop 3768
ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

do i=1,IMT+1
do j=1,JMT
  hs(i,j,2)=temp2d_simp(i,j)
enddo
enddo

do i=4,IMT+1
hs(i,jmt+1,2) =hs(IMT+4-i,jmt-3,2)  !  north fold 
enddo
#ifdef tempsalt 
temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05T.nc' ,'votemper')

#endif     

dmult=1.  ! amplification of the velocity amplitude by simple multiplication

! u velocity
!temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05U.nc' ,'vozocrtx')
gridFile = trim(fieldFile)//'d05U.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i+1,j,2))
   if(k.eq.kmu(i,j)) dd = botbox(i,j,1)
   uflux(i,j,kk,2)=temp3d_simp(i,j,k) * dyu(i,j) * dd * dmult
  enddo
 enddo
enddo

! v velocity
!temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05V.nc' ,'vomecrty')
gridFile = trim(fieldFile)//'d05V.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! kmt field
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i,j+1,2))
   if(k.eq.kmv(i,j)) dd = botbox(i,j,2)
   vflux(i,j,kk,2)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
  enddo
 enddo
enddo


#ifdef coarse

temp3d_simp=0.
do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   temp3d_simp(i,j,kk)=uflux(i,j,kk,2)
  enddo
 enddo
enddo

uflux(:,:,:,2)=0.
do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   uflux(i,j,kk,2)=temp3d_simp(i,j,kk)
  enddo
 enddo
enddo

! 2x2 grid cells into one for U
!do i=2,IMT,2
! im=i-1
! ip=i+1
! if(im.lt.1  ) im=im+IMT
! if(ip.gt.IMT) ip=ip-IMT
! do j=1,JMT,2
!  jm=j-1
!  jp=j+1
!  if(jm.lt.1  ) jm=1
!  if(jp.gt.JMT) jp=JMT
!  do k=1,KM
!   temp3d_simp(im,j ,k)=0.5*( uflux(im,j,k,2) + uflux(im,jp,k,2) )    ! mean on western wall
!   temp3d_simp(im,jp,k)=temp3d_simp(im,j,k)				     	      ! mean on western wall
!   temp3d_simp(ip,j ,k)=0.5*( uflux(ip,j,k,2) + uflux(ip,jp,k,2) )    ! mean on eastern wall
!   temp3d_simp(ip,jp,k)=temp3d_simp(ip,j,k)						      ! mean on eastern wall
!   temp3d_simp(i ,j ,k)=0.5*(temp3d_simp(im,j,k)+temp3d_simp(ip,j,k)) ! linear interpolation in the middle
!   temp3d_simp(i ,jp,k)=temp3d_simp(i,j,k)						      ! linear interpolation in the middle
!  enddo
! enddo
!enddo



! 3x3 grid cells into one for U
!do i=3,IMT,3
! im =i-1
! imm=i-2
! ip =i+1
! if(im .lt.1  ) im =im +IMT
! if(imm.lt.1  ) imm=imm+IMT
! if(ip .gt.IMT) ip =ip -IMT
! do j=1,JMT,3
!  jm =j-1
!  jp =j+1
!  if(jm .lt.1  ) jm =1
!  if(jp .gt.JMT) jp =JMT
!  do k=1,KM
!   temp3d_simp(imm,j ,k)=(uflux(imm,jm,k,2)+uflux(imm,j,k,2)+uflux(imm,jp,k,2))/3. ! western
!   temp3d_simp(imm,jm,k)=temp3d_simp(imm,j,k)							  
!   temp3d_simp(imm,jp,k)=temp3d_simp(imm,j,k)							 
!					  
!   temp3d_simp(ip ,j ,k)=(uflux(ip ,jm,k,2)+uflux(ip ,j,k,2)+uflux(ip ,jp,k,2))/3. ! eastern
!   temp3d_simp(ip ,jm,k)=temp3d_simp(ip ,j,k)							  
!   temp3d_simp(ip ,jp,k)=temp3d_simp(ip ,j,k)							
!
!   temp3d_simp(im ,j ,k)=2./3.*temp3d_simp(imm,j,k)+1./3.*temp3d_simp(ip,j,k) ! linear interpolation in the middle
!   temp3d_simp(im ,jm,k)=temp3d_simp(i,j,k)							  
!   temp3d_simp(im ,jp,k)=temp3d_simp(i,j,k)							  
!
!   temp3d_simp(ip ,j ,k)=1./3.*temp3d_simp(imm,j,k)+2./3.*temp3d_simp(ip,j,k) ! linear interpolation in the middle
!   temp3d_simp(ip ,jm,k)=temp3d_simp(ip,j,k)							  
!   temp3d_simp(ip ,jp,k)=temp3d_simp(ip,j,k)							  
!  enddo
! enddo
!enddo


! 4x4 grid cells into one for U
do i=3,IMT,4
 im =i-1
 imm=i-2
 ip =i+1
 ipp=i+2
 if(im .lt.1  ) im =im +IMT
 if(imm.lt.1  ) imm=imm+IMT
 if(ip .gt.IMT) ip =ip -IMT
 if(ipp.gt.IMT) ipp=ipp-IMT
 do j=2,JMT,4
  jm =j-1
  jp =j+1
  jpp=j+2
  if(jm .lt.1  ) jm =1
  if(jp .gt.JMT) jp =JMT
  if(jpp.gt.JMT) jpp=JMT
  do k=1,KM
   temp3d_simp(imm,j  ,k)=0.25*(uflux(imm,jm,k,2)+uflux(imm,j,k,2)+uflux(imm,jp,k,2)+uflux(imm,jpp,k,2)) ! western
   temp3d_simp(imm,jm ,k)=temp3d_simp(imm,j,k)							  
   temp3d_simp(imm,jp ,k)=temp3d_simp(imm,j,k)							 
   temp3d_simp(imm,jpp,k)=temp3d_simp(imm,j,k)		
					  
   temp3d_simp(ipp,j  ,k)=0.25*(uflux(ipp,jm,k,2)+uflux(ipp,j,k,2)+uflux(ipp,jp,k,2)+uflux(ipp,jpp,k,2)) ! eastern
   temp3d_simp(ipp,jm ,k)=temp3d_simp(ipp,j,k)							  
   temp3d_simp(ipp,jp ,k)=temp3d_simp(ipp,j,k)							
   temp3d_simp(ipp,jpp,k)=temp3d_simp(ipp,j,k)							    

   temp3d_simp(im ,j  ,k)=0.75*temp3d_simp(imm,j,k)+0.25*temp3d_simp(ipp,j,k) ! linear interpolation in the middle
   temp3d_simp(im ,jm ,k)=temp3d_simp(im,j,k)							  
   temp3d_simp(im ,jp ,k)=temp3d_simp(im,j,k)							  
   temp3d_simp(im ,jpp,k)=temp3d_simp(im,j,k)	

   temp3d_simp(i  ,j  ,k)=0.5*temp3d_simp(imm,j,k)+0.5*temp3d_simp(ipp,j,k) ! linear interpolation in the middle
   temp3d_simp(i  ,jm ,k)=temp3d_simp(i,j,k)							  
   temp3d_simp(i  ,jp ,k)=temp3d_simp(i,j,k)							  
   temp3d_simp(i  ,jpp,k)=temp3d_simp(i,j,k)	

   temp3d_simp(ip ,j  ,k)=0.25*temp3d_simp(imm,j,k)+0.75*temp3d_simp(ipp,j,k) ! linear interpolation in the middle
   temp3d_simp(ip ,jm ,k)=temp3d_simp(ip,j,k)							  
   temp3d_simp(ip ,jp ,k)=temp3d_simp(ip,j,k)							  
   temp3d_simp(ip ,jpp,k)=temp3d_simp(ip,j,k)	
  enddo
 enddo
enddo

! put back into velocity field array
do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   uflux(i,j,kk,2)=temp3d_simp(i,j,kk) 
  enddo
 enddo
enddo

!___________ v


temp3d_simp=0.
do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   temp3d_simp(i,j,kk)=vflux(i,j,kk,2)
  enddo
 enddo
enddo

vflux(:,:,:,2)=0.
do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   vflux(i,j,kk,2)=temp3d_simp(i,j,kk)
  enddo
 enddo
enddo

! 2x2 grid cells into one for V
!do i=1,IMT-1,2
! im=i-1
! ip=i+1
! if(im.lt.1  ) im=im+IMT
! if(ip.gt.IMT) ip=ip-IMT
! do j=2,JMT-1,2
!  jm=j-1
!  jp=j+1
!  if(jm.lt.1  ) jm=1
!  if(jp.gt.JMT) jp=JMT
!  do k=1,KM
!   temp3d_simp(i ,jm,k)=0.5*( vflux(i,jm,k,2) + vflux(ip,jm,k,2) )     ! mean on southern wall
!   temp3d_simp(ip,jm,k)=temp3d_simp(i,jm,k)						       ! mean on southern wall
!   temp3d_simp(i ,jp,k)=0.5*( vflux(i,jp,k,2) + vflux(ip,jp,k,2) )     ! mean on northern wall
!   temp3d_simp(ip,jp,k)=temp3d_simp(i,jp,k)						       ! mean on northern wall
!   temp3d_simp(i ,j ,k)=0.5*(temp3d_simp(i,jm,k)+temp3d_simp(i,jp,k))  ! linear interpolation in the middle
!   temp3d_simp(ip,j ,k)=temp3d_simp(i,j,k)							   ! linear interpolation in the middle
!  enddo
! enddo
!enddo


! 3x3 grid cells into one for V
!do i=2,IMT-1,3
! im =i-1
! ip =i+1
! if(im .lt.1  ) im =im +IMT
! if(ip .gt.IMT) ip =ip -IMT
! if(ipp.gt.IMT) ipp=ipp-IMT
! do j=3,JMT-1,3
!  jmm=j-2
!  jm =j-1
!  jp =j+1
!  if(jmm.lt.1  ) jmm=1
!  if(jm .lt.1  ) jm =1
!  if(jp .gt.JMT) jp =JMT
!  do k=1,KM
!   temp3d_simp(i  ,jmm,k)=1./3.*(vflux(im,jmm,k,2)+vflux(i,jmm,k,2)+vflux(ip,jmm,k,2)) ! southern
!   temp3d_simp(im ,jmm,k)=temp3d_simp(i,jmm,k)							  
!   temp3d_simp(ip ,jmm,k)=temp3d_simp(i,jmm,k)							 
!					  
!   temp3d_simp(i  ,jp ,k)=1./3.*(vflux(im,jp,k,2)+vflux(i,jp,k,2)+vflux(ip,jp,k,2)) ! northern
!   temp3d_simp(im ,jp ,k)=temp3d_simp(i,jp ,k)							  
!   temp3d_simp(ip ,jp ,k)=temp3d_simp(i,jp ,k)							
!
!   temp3d_simp(i  ,jm ,k)=2./3.*temp3d_simp(i,jmm,k)+1./3.*temp3d_simp(i,jp,k) ! linear interpolation in the middle
!   temp3d_simp(im ,jm ,k)=temp3d_simp(i,jm,k)							  
!   temp3d_simp(ip ,jm ,k)=temp3d_simp(i,jm,k)							  
!
!   temp3d_simp(i  ,j  ,k)=1./3.*temp3d_simp(i,jmm,k)+2./3.*temp3d_simp(i,jp,k) ! linear interpolation in the middle
!   temp3d_simp(im ,j  ,k)=temp3d_simp(i,j,k)							  
!   temp3d_simp(ip ,j  ,k)=temp3d_simp(i,j,k)							  
!					  
!  enddo
! enddo
!enddo

! 4x4 grid cells into one for V
do i=2,IMT-1,4
 im =i-1
 imm=i-2
 ip =i+1
 ipp=i+2
 if(im .lt.1  ) im =im +IMT
 if(imm.lt.1  ) imm=imm+IMT
 if(ip .gt.IMT) ip =ip -IMT
 if(ipp.gt.IMT) ipp=ipp-IMT
 do j=3,JMT-2,4
  jm =j-1
  jmm=j-2
  jp =j+1
  jpp=j+2
  if(jm .lt.1  ) jm =1
  if(jmm.lt.1  ) jmm=1
  if(jp .gt.JMT) jp =JMT
  if(jpp.gt.JMT) jpp=JMT
  do k=1,KM
   temp3d_simp(i  ,jmm,k)=0.25*(vflux(im,jmm,k,2)+vflux(i,jmm,k,2)+vflux(ip,jmm,k,2)+vflux(ipp,jmm,k,2)) ! southern
   temp3d_simp(im ,jmm,k)=temp3d_simp(i,jmm,k)							  
   temp3d_simp(ip ,jmm,k)=temp3d_simp(i,jmm,k)							 
   temp3d_simp(ipp,jmm,k)=temp3d_simp(i,jmm,k)		
					  
   temp3d_simp(i  ,jpp,k)=0.25*(vflux(im,jpp,k,2)+vflux(i,jpp,k,2)+vflux(ip,jpp,k,2)+vflux(ipp,jpp,k,2)) ! northern
   temp3d_simp(im ,jpp,k)=temp3d_simp(i,jpp,k)							  
   temp3d_simp(ip ,jpp,k)=temp3d_simp(i,jpp,k)							
   temp3d_simp(ipp,jpp,k)=temp3d_simp(i,jpp,k)							    

   temp3d_simp(i  ,jm ,k)=0.75*temp3d_simp(i,jmm,k)+0.25*temp3d_simp(i,jpp,k) ! linear interpolation in the middle
   temp3d_simp(im ,jm ,k)=temp3d_simp(i,jm,k)							  
   temp3d_simp(ip ,jm ,k)=temp3d_simp(i,jm,k)							  
   temp3d_simp(ipp,jm ,k)=temp3d_simp(i,jm,k)	

   temp3d_simp(i  ,j  ,k)=0.5*temp3d_simp(i,jmm,k)+0.5*temp3d_simp(i,jpp,k) ! linear interpolation in the middle
   temp3d_simp(im ,j  ,k)=temp3d_simp(i,j,k)							  
   temp3d_simp(ip ,j  ,k)=temp3d_simp(i,j,k)							  
   temp3d_simp(ipp,j  ,k)=temp3d_simp(i,j,k)	

   temp3d_simp(i  ,jp ,k)=0.25*temp3d_simp(i,jmm,k)+0.75*temp3d_simp(i,jpp,k) ! linear interpolation in the middle
   temp3d_simp(im ,jp ,k)=temp3d_simp(i,jp,k)							  
   temp3d_simp(ip ,jp ,k)=temp3d_simp(i,jp,k)							  
   temp3d_simp(ipp,jp ,k)=temp3d_simp(i,jp,k)	
  enddo
 enddo
enddo

! put back into velocity field array
do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   vflux(i,j,kk,2)=temp3d_simp(i,j,kk) 
  enddo
 enddo
enddo


#endif


#ifdef initxyt
ijkst(:,5)=ntimask(n,:)
#endif


!stop 4967
!uflux=0.0001 ; vflux=0.0001  ! special case with no verlocities

!do j=1,jmt
! do i=1,imt
!  do k=1,km
!   kk=km+1-k
!   if(kk.gt.kmu(i,j)) uflux(i,j,k,2)=0.
!   if(kk.gt.kmv(i,j)) vflux(i,j,k,2)=0.
!  enddo
! enddo
!enddo

     deallocate ( ssh , temp3d_simp, temp2d_doub, e1t , e2t)


  return
end subroutine readfields



