module mod_tracer_functions

   use mod_precdef
   use mod_log

   use mod_time
   use mod_grid, only     : imt, jmt, km, subGridJmin, nsp   
   use mod_tempsalt, only : tracers3D, tracers2D, n3Dtracers, n2Dtracers

   use mod_dens
   use mod_stat

contains
   
   subroutine calculate_3Dtracer(itracer)
      
      implicit none
      
      integer, intent(in) :: itracer
      integer             :: i,j,k,jt
      integer             :: item, isal, ihum
      real(PP), allocatable, DIMENSION(:)   :: rhozvec, depthzvec, latvec,   &
                                             & tmpzvec, salzvec,             &
                                             & prevec , humvec , temvec, dsevec
      
      if (tracers3D(itracer)%desc == 'sigma' ) then
         
         item = 0
         isal = 0
         do jt=1,n3Dtracers
            if (tracers3D(jt)%desc == 'temperature') item=jt
            if (tracers3D(jt)%desc == 'pr_salinity') isal=jt
         end do            
         if (item == 0) then 
            print*,' Could not find temperature, which is needed for density! '
            stop
         end if
         if (isal == 0) then
            print*,' Could not find practical salinity which is needed for density! '
            stop
         end if
         
         allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))         
         
         ! Calculate potential density  
         depthzvec = 0.
         do j=1,jmt
            ! wait, here 1/12 degree step is hard coded
            ! thats not even true for ORCA12 
            latvec=-80+1./12.*float(j+subGridJmin-1)
            do i=1,imt
               tmpzvec = tracers3D(item)%data(i,j,:,nsp)
               salzvec = tracers3D(isal)%data(i,j,:,nsp)
               call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)               
               tracers3D(itracer)%data(i,j,:,nsp) = rhozvec - 1000.
            end do
         end do
         
         deallocate( tmpzvec,salzvec,rhozvec,depthzvec,latvec )
      
      end if
      
      ! 
      ! Atmosphere thermodynamic quantities
      ! 
      if (tracers3D(itracer)%desc == 'DSE' ) then
         
         ! This function is just a placeholder for how you would
         ! calculate dry static energy. 
         
         item = 0
         ihum = 0
         
         do jt=1,n3Dtracers
            if (tracers3D(jt)%desc == 'temperature') item=jt
            if (tracers3D(jt)%desc == 'sp_humidity') ihum=jt
         end do            
         if (item == 0) then 
            print*,' Could not find temperature, which is needed for DSE! '
            stop
         end if
         if (ihum == 0) then
            print*,' Could not find specific humidity which is needed for DSE! '
            stop
         end if
         
         allocate ( prevec(km), humvec(km), temvec(km), dsevec(km) )
         
         ! Calculate dry static energy
         do j=1,jmt
            do i=1,imt
               !prevec = a + b * ps(i,j)
               temvec = tracers3D(item)%data(i,j,:,nsp)
               humvec = tracers3D(ihum)%data(i,j,:,nsp)
               ! do the calculations 
               ! dse = cp * T + g * Z 
               tracers3D(itracer)%data(i,j,:,nsp) = dsevec(:)
            end do
         end do
         
         deallocate( tmpzvec,salzvec,rhozvec,depthzvec,latvec )
      
      end if
      
   end subroutine calculate_3Dtracer

end module mod_tracer_functions
