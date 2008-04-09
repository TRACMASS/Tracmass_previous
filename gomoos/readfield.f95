SUBROUTINE readfields
  
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none

  CHARACTER :: dirdata*10
  CHARACTER(len=63), SAVE :: ncFile
  data dirdata/'/data/GOM/'/
  !data dirdata/'/Users/doos/data/abcde/'/

  integer i,im,j,jm,k,kk,ints2

  INTEGER, SAVE :: nread,ndates
  integer ncid !output ID index of netCDF file
  integer ierr !error 
  integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
  integer varid,varidx,varidy,varidz,varidt !output ID index of variable
  integer startA(1),startB(4),startC(2) !input index vector of position to start reading
  integer countA(1),countB(4),countC(2) !input lengths of 'volume' to be retrieved
  integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
  integer p, x1, y1, z1, t1 !?
  integer intpart1,intpart2
  
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)     :: depth, ssh
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)     :: mask
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)     :: u0mask,u2mask,v0mask,v2mask
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:)     :: e1v,e1t,e2u,e2t,ang
  REAL,          ALLOCATABLE, DIMENSION(:)       :: valsz,fort
  REAL,          ALLOCATABLE, DIMENSION(:,:)     :: NCfieldkmt 
  REAL,          ALLOCATABLE, DIMENSION(:,:,:,:) :: NCfieldx,NCfieldy,NCfieldr 
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: dzu,dzv,dzt
  
  REAL :: x_offset,y_offset,temp_offset,salt_offset, ssh_offset
  REAL :: x_scale,y_scale,temp_scale,salt_scale, ssh_scale
  
  logical around
  
  real, dimension(22) :: dS = (/ 0.008, 0.008, 0.017, 0.033, 0.067, &
       0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, &
       0.067, 0.067, 0.067, 0.033, 0.017, 0.008, 0.008, 0.000 /)
  
  integer year,month,day
  character*16 :: fileName
  
  real :: pi=3.1416E0
  real :: rad,  radi
 
  if ( .NOT. ALLOCATED(ang) ) then
     allocate ( ang(imt,jmt),mask(imt,jmt),depth(imt,jmt) )
     allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km),dzt(imt,jmt,km) )
     allocate ( u0mask(imt,jmt),u2mask(imt,jmt) )
     allocate ( v0mask(imt,jmt),v2mask(imt,jmt) )
  end if
  allocate ( valsz(km),fort(imt*jmt*2) )
  allocate ( NCfieldx(1,imt,jmt,km),NCfieldy(1,imt,jmt,km) )
  allocate ( NCfieldr(1,imt,jmt,km),NCfieldkmt(imt,jmt) )
  allocate ( ssh(imt,jmt) )

  !_____________ swap between datasets ___________________________________
  u(:,:,:,1)=u(:,:,:,2)
  v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif
  
  !____________________________ initialise ___________________________
  if(ints.eq.intstart) then
     call coordinat
     u=0.
     v=0.
#ifdef tempsalt
     tem=0.
     sal=0.
     rho=0.
#endif
     ndates=0
     
     rad=pi/180.
     radi=180./pi
     

     ! === Read  and setup horizontal Grid===
     ierr=nf90_OPEN('/data/GOM/grid.cdf',nf90_NOWRITE,ncid)
     if(ierr.ne.0) stop 3001

     startB(1)=1
     startB(2)=1
     startB(3)=1
     startB(4)=1
     countB(1)=imt
     countB(2)=jmt
     countB(3)=1
     countB(4)=1

     OPEN(41,FILE='/data/GOM/fort.41')
     READ(41,*) fort

     e2t=reshape(fort(1:imt*jmt),(/ imt, jmt/) )
     e2u=e2t

     e1t=reshape(fort(imt*jmt+1:imt*jmt*2),(/ imt, jmt/) )
     e1v=e1t

     dxdy=e1t*e2t

     !=====    Velocity angles      ===== 
     ierr=nf90_INQ_VARID(ncid,'ang',varid)
     if(ierr.ne.0) stop 3002
     ierr=NF90_GET_VAR (ncid,varid,ang,startB,countB)
     if(ierr.ne.0) stop 3005
     ang=rad*ang

     !=========================
     !===== Set up depths =====
     !=========================

     ierr=nf90_OPEN('/data/GOM/grid.cdf',nf90_NOWRITE,ncid)
     if(ierr.ne.0) stop 4001

     startB(1)=1
     startB(2)=1
     startB(3)=1
     startB(4)=1
     countB(1)=imt
     countB(2)=jmt
     countB(3)=km
     countB(4)=1

     ! === Load the bathymetry ===
     ierr=nf90_INQ_VARID(ncid,'depth',varid)
     if(ierr.ne.0) stop 4002

     ierr=NF90_GET_VAR (ncid,varid,NCfieldkmt)
     if(ierr.ne.0) stop 4003

!     fieldkmt=100    !kuk even depth

     depth=int(floor(NCfieldkmt))
     mask=1
     do i=1,imt
        do j=1,jmt
           if (depth(i,j) .lt. 0) then
              depth(i,j)=0
              mask(i,j)=0
           end if
        end do
     end do
     kmt=depth

     ierr=nf90_CLOSE(ncid)
     if(ierr.ne.0) stop 4031
     !===========================================
     !stop 4596

     !============================
     !===== Read in uv-masks =====
     !============================
     
     ierr=nf90_OPEN('/data/GOM/uvmasks.cdf',nf90_NOWRITE,ncid)
     if(ierr.ne.0) stop 5001
     
     startB(1)=1
     startB(2)=1
     startB(3)=1
     startB(4)=1
     countB(1)=imt
     countB(2)=jmt
     countB(3)=km
     countB(4)=1
     
     ! === Load the masks ===
     ierr=nf90_INQ_VARID(ncid,'u0mask',varid)
     if(ierr.ne.0) stop 5002
     ierr=NF90_GET_VAR (ncid,varid,u0mask)
     if(ierr.ne.0) stop 5003
     u0mask=-(u0mask-1)*9999

     ierr=nf90_INQ_VARID(ncid,'u2mask',varid)
     if(ierr.ne.0) stop 5012
     ierr=NF90_GET_VAR (ncid,varid,u2mask)
     if(ierr.ne.0) stop 5013
     u2mask=(u2mask-1)*9999
     
     ierr=nf90_INQ_VARID(ncid,'v0mask',varid)
     if(ierr.ne.0) stop 5022
     ierr=NF90_GET_VAR (ncid,varid,v0mask)
     if(ierr.ne.0) stop 5023
     v0mask=-(v0mask-1)*9999
     
     ierr=nf90_INQ_VARID(ncid,'v2mask',varid)
     if(ierr.ne.0) stop 5032
     ierr=NF90_GET_VAR (ncid,varid,v2mask)
     if(ierr.ne.0) stop 5033
     v2mask=(v2mask-1)*9999
 
  endif !=== endif initstart ===

  intpart1=mod(ints+1,8)
  if (intpart1 .eq. 0) then
     intpart1=8
  endif
  intpart2=floor((ints)/8.)+1

  x_scale     = 9.1556e-5
  y_scale     = 9.1556e-5
  temp_scale  = 0.0005340739
  salt_scale  = 0.0006103702
  ssh_scale   = 0.000244

  x_offset    =-3.3845e-19
  y_offset    =-5.4561e-18
  temp_offset = 12.5
  salt_offset = 20.0
  ssh_offset  = 0

ndates=intpart2

call gdate (2453005+1+ndates, year,month,day)
fileName='unh.20000000.cdf'
write(fileName(5:12),'(i4i2.2i2.2)') year,month,day
ncFile=dirdata//''//fileName//''
print *,ncFile,intpart1

inquire(file=ncFile,exist=around)
if(.not.around) stop 4556

nread=mod(ints/5,18)+1

  startB(1)=1
  startB(2)=1
  startB(3)=1
  startB(4)=intpart1
  countB(1)=imt
  countB(2)=jmt
  countB(3)=km
  countB(4)=1

 ierr=nf90_OPEN(ncFile,nf90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3750

  ! === zonal velocity ===
  ierr=nf90_INQ_VARID(ncid,'u',varid) ! the main data fields
  if(ierr.ne.0) stop 3762
  ierr=NF90_GET_VAR (ncid,varid,NCfieldx,startB,countB)
  if(ierr.ne.0) then  
     print *,'ierr',ierr,'ncid',ncid,'varid',varid,'startB',startB,'countB',countB
     stop 3798
  end if
  NCfieldx=NCfieldx*x_scale+x_offset

  ! === meridional velocity ===
  ierr=nf90_INQ_VARID(ncid,'v',varid) ! the main data NCfields
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR (ncid,varid,NCfieldy,startB,countB)
  if(ierr.ne.0) stop 3799
  NCfieldy=NCfieldy*y_scale+y_offset

  do k=1,km
     NCfieldx(1,:,:,k)=NCfieldx(1,:,:,k)*cos(ang)+NCfieldy(1,:,:,k)*sin(ang)
     NCfieldy(1,:,:,k)=-NCfieldx(1,:,:,k)*sin(ang)+NCfieldy(1,:,:,k)*cos(ang)
  end do

  if (intstep .le. 0) then
     NCfieldx=-NCfieldx
     NCfieldy=-NCfieldy
  end if
  
  u0: do k=1,km
     NCfieldx(1,:,:,k)=min(NCfieldx(1,:,:,k),real(u0mask))
  end do u0
  u2: do k=1,km
     NCfieldx(1,:,:,k)=max(NCfieldx(1,:,:,k),real(u2mask))
  end do u2
  v0: do k=1,km
     NCfieldy(1,:,:,k)=min(NCfieldy(1,:,:,k),real(v0mask))
  end do v0
  v2: do k=1,km
     NCfieldy(1,:,:,k)=max(NCfieldy(1,:,:,k),real(v2mask))
  end do v2

!  NCfieldy=0

  ! === Salinity ===
  ierr=nf90_INQ_VARID(ncid,'salt',varid) ! the main data fields
  if(ierr.ne.0) stop 3767
  ierr=NF90_GET_VAR (ncid,varid,NCfieldr,startB,countB)
  if(ierr.ne.0) stop 3799
  NCfieldr=NCfieldr*salt_scale+salt_offset
  sal(:,:,:,2)=NCfieldr(1,:,:,:)

  ! === Sea Surface height ===
  countB(1)=imt
  countB(2)=jmt
  countB(3)=1
  countB(4)=1  
  ierr=nf90_INQ_VARID(ncid,'elev',varid) ! the main data fields
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR (ncid,varid,ssh,startB,countB)
  if(ierr.ne.0) stop 3799
  
  ssh=ssh *ssh_scale+ssh_offset

  do  k=1,km
     dzt(:,:,km-k+1)=dS(k)*(depth(:,:)+ssh(:,:))
  end do
  
  do i=1,imt-1
     dzv(i,:,:)=0.5*(dzt(i,:,:)+dzt(i+1,:,:))
  enddo
  dzv(imt,:,:)=dzv(imt-1,:,:)
  do j=1,jmt-1
     dzu(:,j,:)=0.5*(dzt(:,j,:)+dzt(:,j+1,:))
  enddo
  dzu(:,jmt,:)=dzu(:,jmt-1,:)

 ierr=nf90_CLOSE(ncid)
   
 do k=1,km
    kk=km-k+1
    u(:,1:80,k,2)  =NCfieldx(1,:,:,kk)*e2u(:,:)*dzu(:,:,k)*mask
    v(:,1:80,k,2)  =NCfieldy(1,:,:,kk)*e1v(:,:)*dzv(:,:,k)*mask
    rho(:,1:80,k,2)=NCfieldr(1,:,:,k)
 enddo

 !Sensitivity run
! u=u*e10.0
! v=v*10.0

 
end subroutine readfields


!=================================================


subroutine gdate (jd, year,month,day)
  !
  !---computes the gregorian calendar date (year,month,day)
  !   given the julian date (jd).
  !
  integer jd,year,month,day,i,j,k
  l= jd+68569
  n= 4*l/146097
  l= l-(146097*n+3)/4
  i= 4000*(l+1)/1461001
  l= l-1461*i/4+31
  j= 80*l/2447
  k= l-2447*j/80
  l= j/11
  j= j+2-12*l
  i= 100*(n-49)+i+l
  
  year= i
  month= j
  day= k
  
  return
end subroutine gdate
!_______________________________________________________________________
