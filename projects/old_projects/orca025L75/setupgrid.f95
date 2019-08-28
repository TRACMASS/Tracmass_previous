SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------



  ! === Init local variables for the subroutine ===
  INTEGER                                     :: i ,j ,k, ip, jp, im, jm
  INTEGER                                     :: kk, ii
  REAL*8                                      :: dp,dd
  REAL*8,  ALLOCATABLE, DIMENSION(:,:)        :: temp2d_doub
  REAL*8,  ALLOCATABLE, DIMENSION(:,:,:)      :: temp3d_doub
  
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t !,rhom
  CHARACTER (len=200)                         :: gridFile
!  REAL*4									  :: long(IMT,JMT),lat(IMT,JMT)

  ! === Start and count mask for reading netCDF files ===
  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count2D  = [         imt,        jmt ,  1 , 1 ]
  map2D    = [          1 ,          2 ,  3 , 4 ]  
  start3D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count3D  = [         imt,        jmt , KM , 1 ]
  map3D    = [          1 ,          2 ,  3 , 4 ] 
  
  call coordinat


  ! === Open mesh file ===
  gridFile = trim(inDataDir)//'topo/mesh_mask_ORCA025.L75_nov2010.nc'
  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3751

  ! === Read dx and dy for T points ===
  ! === Compute area of grid box    ===
  allocate( temp2d_doub(IMT,JMT) )
  allocate ( e1t(IMT,JMT) , e2t(IMT,JMT) )
  
  ierr=NF90_INQ_VARID(ncid,'e1t',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
    e1t(i,:)=temp2d_doub(i,:)
  enddo

  ierr=NF90_INQ_VARID(ncid,'e2t',varid) 
  if(ierr.ne.0) stop 3765
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
    e2t(i,:)=temp2d_doub(i,:)
    dxdy(i,:) = e1t(i,:) * e2t(i,:)
  enddo

  ! === Read dy for u points ===
  ierr=NF90_INQ_VARID(ncid,'e2u',varid) 
  if(ierr.ne.0) stop 3764
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
    dyu(i,:)=temp2d_doub(i,:)
  enddo

  ! === Read dx for v points ===
  ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! dx for v points
  if(ierr.ne.0) stop 3766
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
    dxv(i,:)=temp2d_doub(i,:)
  enddo
  
  ! extras for analysis to be commented out
! === Read dy for u points ===
!  ierr=NF90_INQ_VARID(ncid,'glamu',varid) 
!  if(ierr.ne.0) stop 3764
!  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
!  if(ierr.ne.0) stop 3799
!  do i=1,IMT
!  do  j=1,JMT
!    long(i,j)=temp2d_doub(i,j)
!    if(long(i,j).lt.0.) long(i,j)=long(i,j)+360.
!  enddo
!  enddo
!  ierr=NF90_INQ_VARID(ncid,'gphiv',varid) 
!  if(ierr.ne.0) stop 3764
!  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
!  if(ierr.ne.0) stop 3799
!  do i=1,IMT
!  do  j=1,JMT
!    lat(i,j)=temp2d_doub(i,j)
!  enddo
!  enddo
!do j=JMT,1,-1
! write (*,"(i4,2x,99f6.2)") j,(lat(i,j),i=1,IMT,200)
!enddo
    
  deallocate( temp2d_doub, e1t, e2t )
  
!#ifdef orca025
!  ! === In orca025 the vertical mesh is ===
!  ! === in a separate file              ===
!  ierr=NF90_CLOSE(ncid)
!  gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
!  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!#endif


  ! === Read and compute depth coordinates ===
  zw(0) = 0.d0
  ierr=NF90_INQ_VARID(ncid ,'e3t_0',varid)
  if(ierr.ne.0) stop 3777
  ierr=NF90_GET_VAR(ncid ,varid ,zw(1:KM) ,start1d ,count1d)
  if(ierr.ne.0) stop 3778
  do k=1,km
    kk=km+1-k
    dz(kk)=zw(k)
    zw(k)=zw(k)+zw(k-1)
!    print *,k,kk,dz(kk),zw(k)
  end do
!stop 936
  ! === Read kmt - number of vertical T points ===
  ierr=NF90_INQ_VARID(ncid,'mbathy',varid)
  if(ierr.ne.0) stop 3767
  ierr=NF90_GET_VAR(ncid,varid,kmt,start2d,count2d)
  
  ! === Read and compute the depth of the bottom T point ===
  allocate (  botbox(IMT,JMT,3) )
  allocate( temp3d_doub(IMT,JMT,KM) )
  
  botbox=0.
  ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
  do i=1,IMT
    do j=1,JMT
      if(kmt(i,j).ne.0) then
        botbox(i,j,3)=temp3d_doub(i,j,kmt(i,j))
        if(botbox(i,j,3).eq.0.) then
          print *,i,j,kmt(i,j),botbox(i,j,3),temp3d_doub(i,j,kmt(i,j))
          botbox(i,j,3)=dz(km+1-kmt(i,j))
        endif
      endif
    enddo
  enddo

  ! === Read and compute the depth of the bottom u and v point ===
  allocate ( kmu(IMT,JMT)    ,kmv(IMT,JMT) )
  kmu=0 ; kmv=0

  do j=1,JMT
    jp=j+1
    if(jp.eq.jmt+1) jp=jmt
    do i=1,imt
      ip=i+1
      if(ip.eq.IMT+1) ip=1
      ! Number of depth levels for u and v
      kmu(i,j)=min(kmt(i,j),kmt(ip,j),KM)
      kmv(i,j)=min(kmt(i,j),kmt(i,jp),KM)
      ! Depth of bottom box for u
      if(kmu(i,j).ne.0) then
        if(kmu(i,j) .eq. kmt(i,j) .and. kmu(i,j) .eq. kmt(ip,j)) then
          botbox(i,j,1)=min(dztb(i,j,1),dztb(ip,j,1))
        elseif(kmu(i,j) .eq. kmt(i,j)) then
          botbox(i,j,1)=dztb(i,j,1)
        elseif(kmu(i,j) .eq. kmt(ip,j)) then
          botbox(i,j,1)=dztb(ip,j,1)
        else
          print*,kmu(i,j),kmt(i,j),kmt(ip,j)
        endif
      else
        botbox(i,j,1) = 0.d0
      endif
      ! Depth of bottom box for v
      if(kmv(i,j).ne.0) then
        if(kmv(i,j) == kmt(i,j) .and. kmv(i,j) == kmt(i,jp)) then
          botbox(i,j,2) = min(dztb(i,j,1),dztb(i,jp,1))
        elseif(kmv(i,j) == kmt(i,j)) then
          botbox(i,j,2) = dztb(i,j,1)
        elseif(kmv(i,j) == kmt(i,jp)) then
          botbox(i,j,2) = dztb(i,jp,1)
        else
          print*,kmv(i,j),kmt(i,j),kmt(i,jp)
        endif
      else
        botbox(i,j,2) = 0.d0
      endif
    enddo
  enddo

  !  north fold 
  do i=4,IMT
    ii=IMT+4-i
    kmv(i,JMT)=kmv(ii,JMT-3)
    botbox(i,JMT,2)=botbox(ii,JMT-3,2)
  enddo

  ierr=NF90_CLOSE(ncid)


  ! Bottom box 
  do j=1,JMT
    do i=1,IMT
      if(kmt(i,j).ne.0) then
        dztb(i,j,1)=botbox(i+subGridImin-1,j+subGridJmin-1,3)
        if(botbox(i+subGridImin-1,j+subGridJmin-1,3).eq.0.) then
          print *,i,j,kmt(i,j),botbox(i+subGridImin-1,j+subGridJmin-1,3)
          stop 4957
        endif
      else
        dztb(i,j,1)=0.
      endif
    enddo
  enddo
  detta ligger helt galet!

! Initialise time
!  currSec = startSec 
!  currMin = startMin 
!  currHour = startHour 
!  currDay = startDay 
!  currMon = startMon 
!  currYear = startYear 
  
  
!open(21,file=trim(inDataDir)//'topo/longlat',form='unformatted')
!write(21) long
!write(21) lat
!write(21) kmt
!write(21) kmu
!write(21) kmv
!close(21)

!kmt(940:1200,499)=1
!do j=JMT,1,-1
! write (*,"(i4,999i1)") j,(kmt(i,j),i=800,1300,2)
!enddo
!stop 3956


end SUBROUTINE setupgrid
