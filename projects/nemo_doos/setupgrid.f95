SUBROUTINE setupgrid
   ! =============================================================
   ! 
   ! Purpose
   ! -------
   !
   ! Set up the NEMO mesh, i.e. longitudes, latitudes, cell size etc. 
   ! 
   ! Method
   ! ------
   !
   ! Read information about the horizontal and vertical mesh from 
   ! NEMO output or input files. 
   ! Four files are needed: 
   !   coordFile - lon,lat
   !   hgridFile - e1t,e2t,e1u,e2u etc
   !   zgridFile - dz (1-D), dzt, dzu, dzv (3D)
   !   bathyFile - bathymetry
   ! In NEMO, a file named mesh_mask usually contains all of the above
   ! except bathymetry. 
   !
   ! The following arrays will be populated:
   !
   !  dxdy - Horizontal area of cells (T points)
   !  dz   - Thickness of standard level (T point) 
   !  dzt  - Time-invariant thickness of level (T point)
   !  dzu  - Time-invariant thickness of level (U point)
   !  dzv  - Time-invariant thickness of level (V point)
   !  kmt  - Number of levels from surface to seafloor (T point)
   !  kmu  - Number of levels from surface to seafloor (U point)
   !  kmv  - Number of levels from surface to seafloor (V point)
   !
   ! History:
   ! --------
   !
   ! 02.2019:  J.Kjellsson unifies code for all NEMO configurations
   !
   !
   ! -------------------------------------------------------------
   
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
    
   ! === Init local variables for the subroutine ===
   INTEGER                                      :: i ,j ,k, n, kk, ii, &
   &                                               ip, jp, im, jm !! Loop indices
   REAL(DP), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t        !! dx, dy [m]
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)    :: tmp4d
   CHARACTER (len=200)                          :: gridFile 
   
   !
   ! Print some settings
   !
   print*,' Read 3D dz:                   ',read3Ddz
   print*,' Vertical grid is upside down: ',gridIsUpsideDown
   print*,' Grid data directory:          ',trim(topoDataDir)
   print*,' Coordinates file:             ',trim(coordFile)
   print*,' Hor. grid file:               ',trim(hgridFile)
   print*,' Vert. grid file:              ',trim(zgridFile)
   print*,' Bathy. file:                  ',trim(bathyFile)
   
   !
   ! Set up which positions to read in the netCDF files
   !
   map2D    = [3, 4,  1, 1 ]
   map3D    = [2, 3,  4, 1 ]
   ncTpos   = 1
         
   !
   ! --- Read dx, dy at T points --- 
   !
   allocate ( e1t(imt,jmt) , e2t(imt,jmt) )
   e1t  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dx_name)
   e2t  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dy_name)
   dxdy(1:imt,1:jmt) = e1t(1:imt,1:jmt) * e2t(1:imt,1:jmt)
   deallocate ( e1t, e2t )
  
   !
   ! --- Read dy at U points and dx at V points --- 
   !
   dyu  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dyu_name)
   dxv  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dxv_name)
   dx   = dxv(imt/2, jmt/2)
   dy   = dyu(imt/2, jmt/2)
   
   !
   ! Read dz at T points without considering 
   ! bottom partial cells and variable volume  
   !
   dz = get1DfieldNC(trim(topoDataDir)//trim(zgridFile), dz_1D_name)
   if (gridIsUpsideDown) then
      do k=1,km
         kk=km+1-k
         dz(kk)=zlev(k)
         zlev(k)=zlev(k)+zlev(k-1)
      end do
   else
      do k=1,km
         dz(k)=zlev(k)
         zlev(k)=zlev(k)+zlev(k-1)
      end do
   end if
   
   !
   ! Read number of valid levels at U, V, T points
   ! as 2D array
   !
   kmt = get2DfieldNC(trim(topoDataDir)//trim(bathyFile), kBathy_name)
   allocate ( kmu(imt,jmt), kmv(imt,jmt) )
   
   kmu=0 ; kmv=0
   do j=1,jmt
      jp=j+1
      if(jp == jmt+1) jp=jmt
      do i=1,imt
         ip=i+1
         if(ip == imt+1) ip=1
         kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
         kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
      enddo
   enddo
   
   ! Land-sea mask at north fold 
   do i=4, imt
      ii = imt + 4 - i
      kmv(i,jmt) = kmv(ii,jmt-3)
   enddo
   
   !
   ! Read layer thickness at U, V, T points 
   ! without considering variable volume.
   !
   ! SSH variability is accounted for in readfield each time step
   !
   allocate ( dzu(imt,jmt,km,2),dzv(imt,jmt,km,2), dzt0(imt,jmt,km) )
   
   if (read3Ddz) then
      
      print*,'  Reading 3D dz for u,v,t points ' 
      if (gridIsUpsideDown) then
         allocate ( tmp4d(imt,jmt,km,1) )
         tmp4d(:,:,:,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzt_3D_name)
         dzt0(:,:,km:1:-1)  = tmp4d(:,:,1:km,1)
         
         tmp4d(:,:,:,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzu_3D_name)
         dzu(:,:,km:1:-1,1) = tmp4d(:,:,1:km,1)
         
         tmp4d(:,:,:,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzv_3D_name)
         dzv(:,:,km:1:-1,1) = tmp4d(:,:,1:km,1)
         deallocate( tmp4d )
      else
         dzt0(:,:,1:km)  = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzt_3D_name)
         dzu(:,:,1:km,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzu_3D_name)
         dzv(:,:,1:km,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzv_3D_name)
      end if
      
   else
      
      print*,'  Set dz horizontally constant for u,v,t points ' 
      print*,'  i.e. no partial steps '
      do j=1,jmt
         do i=1,imt
            dzt0(i,j,1:km)  = dz(1:km)
            dzu(i,j,1:km,1) = dz(1:km)
            dzv(i,j,1:km,1) = dz(1:km)
         end do
      end do
      
   end if
   
   !i = 844
   !j = 1410
   !print*,'dzu ',dzu(i,j,:,1)
   !print*,'dzt ',dzt0(i,j,:)
   !print*,'dzv ',dzv(i,j,:,1)
   
      
   !
   ! Ensure thickness is zero in invalid points
   !
   do n=1,2
      do k=1,km
         kk = km-k+1
         where (kk > kmt(1:imt,1:jmt))
            dzt0(:,:,k) = 0
         end where
         where (kk > kmu(1:imt,1:jmt))
            dzu(:,:,k,n) = 0
         end where
         where (kk > kmv(1:imt,1:jmt))
            dzv(:,:,k,n) = 0
         end where
      enddo 
   end do
   
   !do k=1,km
   !print*,'k, kk, kmu, dzu ',k,kk,kmu(i,j), dzu(i,j,k,1)
   !end do
   
   !i = 844
   !j = 1410
   !print*,'dzu2 ',dzu(i,j,:,1)   
   
end SUBROUTINE setupgrid
