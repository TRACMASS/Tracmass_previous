SUBROUTINE readfields
   !!
   !!
   !!  Read velocities from satellite altimetry
   !!
   !!
  
   USE mod_precdef
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_traj
   USE mod_getfile
   use mod_seed
   use mod_tempsalt
   
#ifdef tempsalt
   USE mod_dens
   USE mod_stat
#endif
   
   IMPLICIT none
  
   ! = Loop variables
   INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
   INTEGER                                      :: kbot,ktop
   INTEGER, SAVE                                :: ntempus=0,ntempusb=0,nread
   ! = Variables used for getfield procedures
   CHARACTER (len=200)                          :: fieldFile
   ! = Variables for filename generation
   CHARACTER (len=200)                          :: dataprefix
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)        :: ke
   REAL(PP) dd,hu,hv,uint,vint,zint,hh,h0
  
#ifdef tempsalt
   REAL(PP), ALLOCATABLE, DIMENSION(:)     :: rhozvec, depthzvec, latvec
   REAL(PP), ALLOCATABLE, DIMENSION(:)     :: tmpzvec, salzvec
#endif

   LOGICAL around

!---------------------------------------------------------------

   alloCondUVW: if(.not. allocated (ke)) then
      allocate ( ke(imt,jmt,km))
#ifdef tempsalt
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
#endif
   endif alloCondUVW
   
 
   call datasetswap
   
   call updateClock
   
   
   dataprefix = 'YYYY/dt_global_allsat_phy_l4_YYYYMMDD_YYYYMMDD.nc'
   write(dataprefix(1:4),'(i4.4)')   currYear
   write(dataprefix(30:33),'(i4.4)') currYear
   write(dataprefix(34:35),'(i2.2)') currMon
   write(dataprefix(36:37),'(i2.2)') currDay
   write(dataprefix(39:42),'(i4.4)')   2018
   write(dataprefix(43:44),'(i2.2)')   1
   write(dataprefix(45:46),'(i2.2)')   15
   

   fieldFile = trim(physDataDir)//trim(dataprefix)
   
   do k=1,km
      dzt(:,:,km,nsp) = 1. 
   end do
   
   ncTpos = 1
   
      
   map2d = [3, 4, 1, 2]
   !! Read u, v 
   uvel(1:imt,1:jmt,1) = get2DfieldNC(trim(fieldFile), ueul_name)
   vvel(1:imt,1:jmt,1) = get2DfieldNC(trim(fieldFile), veul_name)
   
! Comment this out next time
!   do j = 1, jmt
!    do i = 1, imt
!     if(uvel(i,j,1)> -214748364) then 
!      kmt(i,j)=1.
!     end if
!    end do
!   end do

   
   !! Land points are indicated by FillValue = -2147483647
   !! Find abnormal velocities and set them to zero to indicate land instead
   where (uvel(1:imt,1:jmt,1) <= -2147483647)
      uvel(1:imt,1:jmt,1) = 0.
      vvel(1:imt,1:jmt,1) = 0.
!      kmt(1:imt,1:jmt)    = 0
!      kmu(1:imt,1:jmt)    = 0
!      kmv(1:imt,1:jmt)    = 0
   end where
   

   
   !! Multiply by scale factor
   uvel(1:imt,1:jmt,1) = uvel(1:imt,1:jmt,1) * 0.0001  ![m/s]
   vvel(1:imt,1:jmt,1) = vvel(1:imt,1:jmt,1) * 0.0001  ![m/s]
   
   
   ke(:,:,:) = uvel(1:imt,1:jmt,:) * uvel(1:imt,1:jmt,:) + &
             & vvel(1:imt,1:jmt,:) * vvel(1:imt,1:jmt,:)

   !! uvel, vvel come on an A grid, so we need to interpolate to 
   !! staggered C grid
   do i = 1 , imt
      ip = i+1
      if (i == imt) ip=1
      do j = 1, jmt
         jp = j+1
         if (j == jmt) jp=jmt
         do k = 1, km
            uflux(i,j,k,nsp) = 0.5 * (uvel(i,j,k) + uvel(ip,j,k)) * dyu(i,j) * dzt(i,j,k,nsp) 
            vflux(i,j,k,nsp) = 0.5 * (vvel(i,j,k) + vvel(i,jp,k)) * dxv(i,j) * dzt(i,j,k,nsp) 
         enddo
      enddo
   enddo
   
   !! Zero velocity at coastlines
   do j = 1, jmt
    do i = 1, imt
     if(kmu(i,j) ==  0) then 
      uflux(i,j,:,:) = 0.
     end if
     if (kmv(i,j) == 0) then 
      vflux(i,j,:,:) = 0.
     end if
    end do
   end do
   

! Store only for a first separate run using the entire time series   
! Comment this out next time
!   open(unit=111,file='/Users/doos/data/aviso/topo/kmt.bin',form='unformatted')
!   write(111) kmt
!   close(111) 
      
   !! Zero meridional flux at j=0 and j=jmt
   vflux(:,0  ,:,:) = 0.
   vflux(:,jmt,:,:) = 0.
   
   
   return
   
end subroutine readfields



