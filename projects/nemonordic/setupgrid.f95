SUBROUTINE setupgrid
  
   USE mod_precdef
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
   !    ===  Set up the grid for ORCA0083 configuration ===
   ! =============================================================
   ! Subroutine for defining the grid of the ORCA0083 config. 
   ! Run once before the loop starts.
   ! -------------------------------------------------------------
   ! The following arrays will be populated:
   !
   !  dxdy - Horizontal area of cells (T points)
   !  dz   - Thickness of standard level (T point) 
   !  dzt  - Time-invariant thickness of level (T point)
   !  dzu  - Time-invariant thickness of level (U point) not needed when u is in m3/s
   !  dzv  - Time-invariant thickness of level (V point)
   !  kmt  - Number of levels from surface to seafloor (T point)
   !  kmu  - Number of levels from surface to seafloor (U point)
   !  kmv  - Number of levels from surface to seafloor (V point)
   !
   ! -------------------------------------------------------------
    
   ! === Init local variables for the subroutine ===
   INTEGER                                      :: i ,j ,k, n, kk, ii, &
   &                                               ip, jp, im, jm !! Loop indices
   REAL(DP), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t        !! dx, dy [m]
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)    :: tmp4D
   CHARACTER (len=200)                          :: gridFile 
   INTEGER, PARAMETER :: IMTG=619,JMTG=523,ILAG=430,JLAG=274
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)    :: temp2d_doub
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: temp3d_doub
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)    :: temp2d_int

!   map2D    = [3, 4,  1, 1 ]
!   map3D    = [2, 3,  4, 1 ]
   ncTpos   = 1
   
  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [  1 ,   1,  1 , 1 ]
  count2D  = [imtg,jmtg,  1 , 1 ]
  map2D    = [  1 ,  2 ,  3 , 4 ]  
  start3D  = [1   ,1   ,  1 , 1 ]
  count3D  = [imtg,jmtg, KM , 1 ]
  map3D    = [1   ,   2,  3 , 4 ] 
   ! --- Read dx, dy at T points --- 
   !
   allocate (temp2d_doub(IMTG, JMTG)  )
!   allocate (dzt0(imt,jmt,km) )

 gridFile = trim(inDataDir)//'topo/mesh_mask_baltix_56_levels.nc'
  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3751
  ierr=NF90_INQ_VARID(ncid,'e1t',varid)
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799

allocate ( e1t(imt,jmt) , e2t(imt,jmt) )
do j=1,JMT
 do i=1,IMT
  e1t(i,j)=temp2d_doub(i+ILAG,j+JLAG) ! [m]
 enddo
enddo

  ierr=NF90_INQ_VARID(ncid,'e2t',varid)
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799

do j=1,JMT
 do i=1,IMT
  e2t(i,j)=temp2d_doub(i+ILAG,j+JLAG) ! [m]
 enddo
enddo

   
   dxdy(1:imt,1:jmt) = e1t(1:imt,1:jmt) * e2t(1:imt,1:jmt)
   deallocate ( e1t, e2t )
     
   !
   ! --- Read dy at U points and dx at V points --- 
   !
!   dyu  = get2DfieldNC(gridFile, 'e2u')
!   dxv  = get2DfieldNC(gridFile, 'e1v')
!   dx   = dxv(imt/2, jmt/2)
!   dy   = dyu(imt/2, jmt/2)
   
   !
   ! Read dz at T points without considering 
   ! bottom partial cells and variable volume  
   !
!   allocate (zlev(KM)  )

   zlev = get1DfieldNC(gridFile, 'e3t_0')
   do k=1,km
      kk=km+1-k
      dz(kk)=zlev(k)
      zlev(k)=zlev(k)+zlev(k-1)
   end do
   
!   print *,dz,zlev
   
   
   !
   ! Read number of valid levels at U, V, T points
   ! as 2D array
   !
!   kmt = get2DfieldNC(gridFile, 'mbathy') 
  ALLOCATE( temp2d_int(IMTG,JMTG), kmu(IMT,JMT), kmv(IMT,JMT) )
   
  if(ierr.ne.0) stop 5751
  ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
  if(ierr.ne.0) stop 3767
  ierr=NF90_GET_VAR(ncid,varid,temp2d_int,start2d,count2d)
  
do j=1,JMT
 do i=1,IMT
  kmt(i,j)=temp2d_int(i+ILAG,j+JLAG)
 enddo
enddo

!open(78,file=TRIM(inDataDir)//'topo/topo.bin',form='unformatted')
!write(78) kmt
!close(78)



  DEALLOCATE( temp2d_int ) 
   
   kmu=0 ; kmv=0
   do j=1,jmt-1
      jp=j+1
      do i=1,imt-1
         ip=i+1
         kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
         kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
      enddo
   enddo

   !
   ! Read layer thickness at U, V, T points 
   ! without considering variable volume.
   !
   ! SSH variability is accounted for in readfield each time step
   !
   ALLOCATE( temp3d_doub(IMTG,JMTG,KM),dzt0(imt,jmt,km) )
   dzt0=0.
     ! === dz at T points ===
  ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
   do k=1,km
  do i=1,IMT
      do j=1,JMT
          if(k<=kmt(i,j)) dzt0(i,j,km+1-k)=temp3d_doub(i+ILAG,j+JLAG,k)
!          if(K==km) print *,i,k,temp3d_doub(i+ILAG,j+JLAG,km+1-k)
      enddo
  enddo
  enddo
  
!  print *,(temp3d_doub(imtg/2,j,1),j=1,jmtg)
!  print *,(dzt0(imt/2,j,km),j=1,jmt)
!  stop 3957
  
  DEALLOCATE( temp3d_doub )

   ! Ensure thickness is zero in invalid points
   !
!      do k=1,km
!         where (k > kmt(1:imt,1:jmt))
!            dzt0(:,:,k) = 0
!         end where
!      enddo 

   
end SUBROUTINE setupgrid
