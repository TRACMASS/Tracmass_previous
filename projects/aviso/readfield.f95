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
   REAL(P4) dd,hu,hv,uint,vint,zint,hh,h0
  
#ifdef tempsalt
   REAL(P4), ALLOCATABLE, DIMENSION(:)     :: rhozvec, depthzvec, latvec
   REAL(P4), ALLOCATABLE, DIMENSION(:)     :: tmpzvec, salzvec
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
   
   dataprefix = 'YYYY/dt_global_allsat_madt_uv_YYYYMMDD_YYYYMMDD.nc'
   write(dataprefix(1:4),'(i4.4)')   currYear
   write(dataprefix(31:34),'(i4.4)') currYear
   write(dataprefix(35:36),'(i2.2)') currMon
   write(dataprefix(37:38),'(i2.2)') currDay
   write(dataprefix(40:43),'(i4.4)')   2014
   write(dataprefix(44:45),'(i2.2)')   1
   write(dataprefix(46:47),'(i2.2)')   6
   
   fieldFile = trim(inDataDir)//'uv/'//trim(dataprefix)
   do k=1,km
      dzt(:,:,km,nsp) = 1. 
   end do
   
   ncTpos = 1
   
   kmt(:,:) = 1
   kmu(:,:) = 1
   kmv(:,:) = 1
      
   map2d = [3, 4, 1, 2]
   !! Read u, v 
   uvel(1:imt,1:jmt,1) = get2DfieldNC(trim(fieldFile), 'u')
   vvel(1:imt,1:jmt,1) = get2DfieldNC(trim(fieldFile), 'v')
   
   !! Land points are indicated by FillValue = -2147483647
   !! Find abnormal velocities and set them to zero to indicate land instead
   where (uvel(1:imt,1:jmt,1) <= -2147483647)
      uvel(1:imt,1:jmt,1) = 0.
      vvel(1:imt,1:jmt,1) = 0.
      kmt(1:imt,1:jmt)    = 0
      kmu(1:imt,1:jmt)    = 0
      kmv(1:imt,1:jmt)    = 0
   end where
      
   !! Multiply by scale factor
   uvel(1:imt,1:jmt,1) = uvel(1:imt,1:jmt,1) * 0.0001  ![m/s]
   vvel(1:imt,1:jmt,1) = vvel(1:imt,1:jmt,1) * 0.0001  ![m/s]
   
   print*,fieldFile
   
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
   do j = 2, jmt-1
      do i = 1, imt
         ip = i+1
         if (i == imt) ip=1
         im = i-1
         if (i == 1) im=imt
         if (kmt(i,j) == 1 .and. kmt(ip,j) == 0) then !western boundary
            uflux(i,j,1,:) = 0.
            kmu(i,j) = 0
         end if
         if (kmt(i,j) == 1 .and. kmt(im,j) == 0) then !eastern boundary
            uflux(im,j,1,:) = 0.
            kmu(im,j) = 0
         end if
         if (kmt(i,j) == 1 .and. kmt(i,j-1) == 0) then !southern boundary
            vflux(i,j-1,1,:) = 0.
            kmv(i,j-1) = 0
         end if
         if (kmt(i,j) == 1 .and. kmt(i,j+1) == 0) then !northern boundary
            vflux(i,j,1,:) = 0.
            kmv(i,j) = 0
         end if
      end do
   end do
   
   
   
   !open(unit=111,file='/Users/joakim/Downloads/kmt.bin',form='unformatted')
   !write(111) kmt
   !close(111) 
   
   !open(unit=111,file='/Users/joakim/Downloads/kmu.bin',form='unformatted')
   !write(111) kmu
   !close(111) 
   
   !open(unit=111,file='/Users/joakim/Downloads/kmv.bin',form='unformatted')
   !write(111) kmv
   !close(111) 
   
   !! Zero meridional flux at j=0 and j=jmt
   vflux(:,0  ,:,:) = 0.
   vflux(:,jmt,:,:) = 0.
   
   do i=1,IMT
   do j=1,JMT
   do k=1,KM
   if(kmv(i,j) == 0 .and. vflux(i,j,km+1-k,nsp)/=0. ) then
      print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
      stop 4966
   endif
   if(kmu(i,j) == 0 .and. uflux(i,j,km+1-k,nsp)/=0. ) then
      print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j)
      stop 4967
   endif
   enddo
   enddo
   enddo
   
   return
   
end subroutine readfields



