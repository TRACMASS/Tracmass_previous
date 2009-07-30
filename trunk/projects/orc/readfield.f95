SUBROUTINE readfields

!  USE io_ezcdf
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  

#ifdef tempsalt
  USE mod_dens
#endif
  IMPLICIT none
! INCLUDE 'netcdf.inc'
! INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'
! include "/Applications/Utilities/netcdf-3.6.2_old/include/netcdf.inc"

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos

!  integer p, x1, y1, z1, t1 !?

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName
!  INTEGER                                   :: start1d, count1d
 ! INTEGER, DIMENSION(2)                     :: start2d, count2d
 ! INTEGER, DIMENSION(3)                     :: start3d, count3d
  
    ! = Variables for filename generation
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2 ,subYr
  INTEGER                                    :: filePos ,fileJD ,subYrD
  INTEGER                                    :: yr1 ,mn1 ,dy1
  INTEGER                                    :: yr2 ,mn2 ,dy2
  



  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: ssh
!  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof
  
  REAL*4, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,e1t,e2u,e2t
!  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzt

!NET CDF STUFF
!INTEGER ierr !error 
INTEGER dimidx,dimidy,dimidz,dimidt !output ID index of dimension
INTEGER varid,varidx,varidy,varidz,varidt !output ID index of variable
INTEGER startA(1),startB(4),startC(2) !input index vector of position to start reading
INTEGER countA(1),countB(4),countC(2) !input lengths of 'volume' to be retrieved
INTEGER leny,lenz,lent !output Length of dimension
!INTEGER ncid, id_x, id_y,id_idta, id_jdta, jpidta, jpjdta  
    INTEGER                                 :: ncid

real*8 tempxy(IMT+2,JMT),tempxyz(IMT+2,JMT,KM)
REAL temp3d_simp(IMT+2,JMT,KM),temp2d_simp(IMT+2,JMT)
INTEGER itemp(IMT+2,JMT)

  logical around

!  print *,'read boerjar'
  
  alloCondGrid: if ( .not. allocated (e1v) ) then
     allocate ( e1v(IMT+2,JMT)    ,e1t(IMT+2,JMT) )
     allocate ( e2u(IMT+2,JMT)    ,e2t(IMT+2,JMT) )
!     allocate ( dzt(IMT+2,JMT,KM) )
  end if alloCondGrid
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
!     allocate ( rhof(IMT+2,JMT,KM) )
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
ierr=NF90_GET_VAR(ncid,varid,tempxy,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e1t(i,:)=tempxy(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2u',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,tempxy,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e2u(i,:)=tempxy(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2t',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,tempxy,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e2t(i,:)=tempxy(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,tempxy,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e1v(i,:)=tempxy(i,:)
enddo
ierr=NF90_CLOSE(ncid)

dxdy = e1t * e2t
   

! The bathymetry
gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,itemp,start2d,count2d)
ierr=NF90_CLOSE(ncid)
do i=1,IMT
  kmt(i,:)=itemp(i,:)
enddo

!open(21,file=trim(inDataDir)//'topo/kmt',form='unformatted')
!write(21)kmt
!close(21)
!stop 4985

! Depth coordinates
zw(0)    = 0.d0
zw(1:km) = get1DfieldNC(trim(gridFile) ,'e3w_0')
do k=1,km
 kk=km+1-k
 dz(kk)=zw(k)
 zw(k)=zw(k)+zw(k-1)
! print *,k,zw(k),kk,dz(kk)
end do

! Bottom box not implemented yet
do j=1,JMT
 do i=1,IMT
!  do k=1,KM
!   kk=km+1-k
!   if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
!  enddo
  if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dz(kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        enddo
     enddo
     
iday=startDay-5
imon=startMon
iyear=startYear
     
endif initFieldcond

! date update
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
! file names
dataprefix='/9999/ORCA025-N112_20000105'
write(dataprefix(20:27),'(i8)') ntime
write(dataprefix(2:5),'(i4)') iyear
  fieldFile = trim(inDataDir)//trim(dataprefix)
  


! temp, salt and ssh
temp2d_simp = get2DfieldNC(trim(fieldFile)//'d05T.nc' ,'sossheig')
do i=1,IMT+1
  hs(i,:,2)=temp2d_simp(i,:)
enddo
hs(:,jmt+1,2) =hs(:,1,2)
#ifdef tempsalt 
temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05T.nc' ,'votemper')
#endif

! u velocity
temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05U.nc' ,'vozocrtx')
do i=1,IMT
 do k=1,km-1
  kk=km+1-k
  uflux(i,:,k ,2)=temp3d_simp(i,:,kk) * e2u(i,:) * dz(kk)
 enddo
  uflux(i,:,km,2)=temp3d_simp(i,:,1 ) * e2u(i,:) * (dz(km) + 0.5*(hs(i,:,2) + hs(i+1,:,2)))
enddo
! v velocity
temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05V.nc' ,'vomecrty')
do i=1,IMT
 do k=1,km-1
  kk=km+1-k
  vflux(i,:,k ,2)=temp3d_simp(i,:,kk) * e1v(i,:) * dz(kk)
 enddo
  vflux(i,:,km,2)=temp3d_simp(i,:,1 ) * e2u(i,:) * (dz(km) + 0.5*(hs(i,:,2) + hs(i,:+1,2)))
enddo
  
!  print *,'read slut'

  return
end subroutine readfields



