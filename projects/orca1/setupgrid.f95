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
!  INTEGER, ALLOCATABLE, DIMENSION(:,:)        :: temp2d_int
  REAL*8                                      :: dp,dd
!  REAL*8,  ALLOCATABLE, DIMENSION(:,:)        :: temp2d_doub
!  REAL*8,  ALLOCATABLE, DIMENSION(:,:,:)      :: temp3d_doub
!  REAL*4									  :: long(IMT,JMT),lat(IMT,JMT)
  
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t !,rhom
  CHARACTER (len=200)                         :: gridFile

  ALLOCATE ( e1t(IMT,JMT) , e2t(IMT,JMT) )
  ALLOCATE ( dztb(imt,jmt,km),dzu(imt,jmt,km),dzv(imt,jmt,km) )
  ALLOCATE ( kmu(IMT,JMT), kmv(IMT,JMT) )
  map2D    = [3, 4,  1, 1 ]
  map3D    = [2, 3,  4, 1 ]
  ncTpos = 1

  ! === Open mesh file ===
!  gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'
  if(KM==42) then
 !  gridFile = trim(inDataDir)//'topo/mesh_mask_ecearth_42l.nc'
   gridFile = '/Users/doos/data/ecearth/topo/mesh_mask_ORCA1_ecearth2_42l.nc'
  elseif(KM==64) then
  ! gridFile = trim(inDataDir)//'topo/mesh_mask_orca1l64.nc'
  endif
  e1t  = get2DfieldNC(gridFile, 'e1t')
  e2t  = get2DfieldNC(gridFile, 'e2t')
  dxdy = e1t * e2t

  dyu  = get2DfieldNC(gridFile, 'e2u')
  dxv  = get2DfieldNC(gridFile, 'e1v')
  dx   = dxv(imt/2, jmt/2)
  dy   = dyu(imt/2, jmt/2)

!  gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
  dz = get1DfieldNC(gridFile, 'e3t_0')
  do k=1,km
    kk=km+1-k
    dz(kk)=zlev(k)
    zlev(k)=zlev(k)+zlev(k-1)
  end do
  
  kmt = get2DfieldNC(gridFile, 'mbathy')
  
!  print *,(kmt(i,JMT/3),i=800,810)
  
  kmu=0 ; kmv=0
  do j=1,jmt-1
     jp=j+1
     if(jp.eq.jmt+1) jp=jmt
     do i=1,imt
        ip=i+1
        if(ip.eq.IMT+1) ip=1
        kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
        kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
     enddo
  enddo

  do i=4, imt
     ii = imt + 4 - i
     kmv(i,jmt) = kmv(ii,jmt-3)
  enddo

  dztb = get3DfieldNC(gridFile, 'e3t')
  dzu  = get3DfieldNC(gridFile, 'e3u')
  dzv  = get3DfieldNC(gridFile, 'e3v')
  
  do k=1,km
     where (k > kmt)
        dztb(:,:,k) = 0
     end where
     where (k > kmu)
        dzu(:,:,k) = 0
     end where
     where (k > kmv)
        dzv(:,:,k) = 0
    end where
  enddo 
  
  
!  do j=1,jmt-1
!     jp=j+1
!     if(jp.eq.jmt+1) jp=jmt
!     do i=1,imt
!        ip=i+1
!        if(ip.eq.IMT+1) ip=1
!        do k=1,KM
!         dzu(i,j,k)=min(dztb(i,j,k), dztb(ip,j,k))
!         dzv(i,j,k)=min(dztb(i,j,k), dztb(i,jp,k))
!        enddo
!     enddo
!  enddo


! do j=2500,2300,-1
!  print "(400i1)",(kmu(i,j),i=3520,3680)
! enddo
! stop 3895


!  print *,'depth'
!  print *,(depth(i,JMT/3),i=800,810)
!  print *,'dztb'
!  print *,(dztb(i,JMT/3,kmt(i,JMT/3)),i=800,810)
!  Print *,'kolla en punkt'
!  print *,kmt(800,JMT/3)
!  print *, dztb(800,JMT/3,:)
!  print *, depth(800,JMT/3)

    
 ! dztb(:,:,2) = dztb(:,:,1) ! ????? what is this????

 !print *,dztb(1,1000,1:5)

  
  
!  gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'
!  print *,trim(gridFile)
!  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!  if(ierr.ne.0) stop 3751

  ! === Read dx and dy for T points ===
  ! === Compute area of grid box    ===
  !allocate( temp2d_doub(IMT,JMT) )

  
  !ierr=NF90_INQ_VARID(ncid,'e1t',varid) 
  !if(ierr.ne.0) stop 3763
  !ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  !if(ierr.ne.0) stop 3799
  !do i=1,IMT
  !  e1t(i,:)=temp2d_doub(i,:)
  !enddo

  !ierr=NF90_INQ_VARID(ncid,'e2t',varid) 
  !if(ierr.ne.0) stop 3765
  !ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  !if(ierr.ne.0) stop 3799
  !do i=1,IMT
  !  e2t(i,:)=temp2d_doub(i,:)
  !  dxdy(i,:) = e1t(i,:) * e2t(i,:)
  !enddo

  ! === Read dy for u points ===
  !ierr=NF90_INQ_VARID(ncid,'e2u',varid) 
  !if(ierr.ne.0) stop 3764
  !ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  !if(ierr.ne.0) stop 3799
  !do i=1,IMT
  !  dyu(i,:)=temp2d_doub(i,:)
  !enddo

  ! === Read dx for v points ===
  !ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! dx for v points
  !if(ierr.ne.0) stop 3766
  !ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  !if(ierr.ne.0) stop 3799
  !do i=1,IMT
  !  dxv(i,:)=temp2d_doub(i,:)
  !enddo
  
    
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
!  print *,(lat(IMT/2,j),j=1,JMT)
  
  
  !deallocate( temp2d_doub, e1t, e2t )
  
  
!  ! === In orca025 the vertical mesh is ===
!  ! === in a separate file              ===
!  ierr=NF90_CLOSE(ncid)
!  gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
!  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)


  ! === Read and compute depth coordinates ===
!  zw(0) = 0.d0
!  ierr=NF90_INQ_VARID(ncid ,'e3t_0',varid)
!  if(ierr.ne.0) stop 3777
!  ierr=NF90_GET_VAR(ncid ,varid ,zw(1:KM) ,start1d ,count1d)
!  if(ierr.ne.0) stop 3778
!  do k=1,km
!    kk=km+1-k
!    dz(kk)=zw(k)
!    zw(k)=zw(k)+zw(k-1)
!    print *,k,kk,dz(kk),zw(k)
!  end do

!  ALLOCATE( temp2d_int(IMT,JMT), kmu(IMT,JMT), kmv(IMT,JMT) )

  ! === Bathymetry ===
  !if(ierr.ne.0) stop 5751
  !ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
  !if(ierr.ne.0) stop 3767
  !ierr=NF90_GET_VAR(ncid,varid,temp2d_int,start2d,count2d)
  !do i=1,IMT
  !    kmt(i,:)=temp2d_int(i,:)
  !enddo
  
  !kmu=0 ; kmv=0
  !do j=1,jmt
  !    jp=j+1
  !    if(jp.eq.jmt+1) jp=jmt
  !    do i=1,imt
  !        ip=i+1
  !        if(ip.eq.IMT+1) ip=1
  !        kmu(i,j)=min(temp2d_int(i,j),temp2d_int(ip,j),KM)
  !        kmv(i,j)=min(temp2d_int(i,j),temp2d_int(i,jp),KM)
  !    enddo
  !enddo
  
!  north fold 
!  do i=4,IMT
!    ii=IMT+4-i
!    kmv(i,JMT)=kmv(ii,JMT-3)
!  enddo


!  DEALLOCATE( temp2d_int ) 
!  ! === Read and compute the depth of the bottom T point ===
!  ALLOCATE( temp3d_doub(IMT+2,JMT,KM) )
  
  ! === dz at T points ===
!  ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
!  if(ierr.ne.0) stop 3763
!  ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!  do i=1,IMT
!      do j=1,JMT
!          if(kmt(i,j).ne.0) dztb(i,j,1)=temp3d_doub(i,j,kmt(i,j))
!      enddo
!  enddo
  
!  DEALLOCATE( temp3d_doub )
!  ierr=NF90_CLOSE(ncid)


  ! dz is independent of x,y except at bottom
!  do j=1,JMT
!      do i=1,IMT
!          if(kmt(i,j).eq.0) then
!              dztb(i,j,1)=0.
!          endif
!      enddo
!  enddo
  
!  dztb(:,:,2)=dztb(:,:,1)  ! this and the third "time" dimension is probably unnecessary

  
!open(21,file='/Users/doos/Dropbox/data/ecearth/topo/kmt_ecearth2_42.bin',form='unformatted')
!!write(21) long
!!write(21) lat
!!write(21) e3t_0
!write(21) kmt
!!write(21) kmu
!!write(21) kmv
!!write(21) dxdy
!!write(21) botbox
!close(21)
!!print *,dxdy
!stop 4596



  
end SUBROUTINE setupgrid
