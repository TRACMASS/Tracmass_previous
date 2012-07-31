module mod_pos
  USE mod_param
  USE mod_grid
  USE mod_vel
  USE mod_loopvars
  USE mod_time
#ifdef streamxy
  USE mod_streamxy
#endif
#ifdef streamv
  USE mod_streamv
#endif
#ifdef streamr
  USE mod_streamr
#endif
#ifdef stream_thermohaline
  USE mod_stream_thermohaline
#endif
  USE mod_psi
    
  IMPLICIT none
  
contains
  
  subroutine  pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)
    
    INTEGER                                    :: mra,mta,msa
    INTEGER                                    :: mrb,mtb,msb
    REAL                                       :: uu
    INTEGER                                    :: ia, iam, ja, ka,k
    INTEGER                                    :: ib, jb, kb
    REAL                                       :: temp,salt,dens
    REAL*8                                     :: dza,dzb, zz
    REAL*8, INTENT(IN)                         :: x0, y0, z0
    REAL*8, INTENT(OUT)                        :: x1, y1, z1
        
    ! === calculate the new positions ===
    ! === of the trajectory           ===    
    scrivi=.false.
    if(ds==dse) then ! eastward grid-cell exit 
       scrivi=.false.
       uu=(rbg*uflux(ia,ja,ka,nsp)+rb*uflux(ia ,ja,ka,nsm))*ff
       if(uu.gt.0.d0) then
          ib=ia+1
          if(ib.gt.IMT) ib=ib-IMT 
       endif
       x1=dble(ia)
#if defined timeanalyt
       call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr) 
       call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr)
#endif /*timeanalyt*/

#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
       mrb=int((dens-rmin)/dr)+1
       if(mrb.lt.1 ) mrb=1
       if(mrb.gt.MR) mrb=MR
#if defined streamts 
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
#endif 
#endif 
#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 

       call savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,1,1,real(subvol*ff))
        
    elseif(ds==dsw) then ! westward grid-cell exit
       scrivi=.false.
       uu=(rbg*uflux(iam,ja,ka,nsp)+rb*uflux(iam,ja,ka,nsm))*ff
       if(uu.lt.0.d0) then
          ib=iam
       endif
       x1=dble(iam)
#if defined timeanalyt
       call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr) ! meridional position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vertical position
#endif
!       scrivi=.true.      
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
       mrb=int((dens-rmin)/dr)+1
       if(mrb.lt.1 ) mrb=1
       if(mrb.gt.MR) mrb=MR
#if defined streamts 
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
#endif 
#endif 
#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(iam,ja,ka,mrb,mta,mtb,msa,msb,1,-1,real(subvol*ff))

    elseif(ds==dsn) then ! northward grid-cell exit
       
       scrivi=.false.
       uu=(rbg*vflux(ia,ja,ka,nsp)+rb*vflux(ia,ja,ka,nsm))*ff
       if(uu.gt.0.d0) then
          jb=ja+1
       endif
       y1=dble(ja)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr) ! zonal position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vertical position
#endif
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
       mrb=int((dens-rmin)/dr)+1
       if(mrb.lt.1 ) mrb=1
       if(mrb.gt.MR) mrb=MR
#if defined streamts 
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
#endif 
#endif 
#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,2,1,real(subvol*ff))

    elseif(ds==dss) then ! southward grid-cell exit
       
       scrivi=.false.
       uu=(rbg*vflux(ia,ja-1,ka,nsp)+rb*vflux(ia,ja-1,ka,nsm))*ff
       if(uu.lt.0.d0) then
          jb=ja-1
#ifndef ifs 
        if(jb==0) stop 34578
#endif
       endif
       y1=dble(ja-1)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr) ! zonal position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vertical position
#endif
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
       mrb=int((dens-rmin)/dr)+1
       if(mrb.lt.1 ) mrb=1
       if(mrb.gt.MR) mrb=MR
#if defined streamts 
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
#endif 
#endif 
#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(ia,ja-1,ka,mrb,mta,mtb,msa,msb,2,-1,real(subvol*ff))
       
    elseif(ds==dsu) then ! upward grid-cell exit
       scrivi=.false.
       call vertvel(rb,ia,iam,ja,ka)
#ifdef full_wflux
       uu=wflux(ia,ja,ka,nsm)
#else
       uu=rbg*wflux(ka,nsp)+rb*wflux(ka,nsm)
#endif
       if(uu.gt.0.d0) then
          kb=ka+1
       endif
       z1=dble(ka)
       if(kb==KM+1) then    ! prevent "evaporation" and put particle from the surface
          kb=KM           
          z1=dble(KM)-0.5d0 ! to the middle of the surface layer
       endif
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr)
       call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr)
#endif

#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,3,1,real(subvol*ff))

    elseif(ds==dsd) then ! downward grid-cell exit
       scrivi=.false.
       call vertvel(rb,ia,iam,ja,ka)
       
#ifdef full_wflux
       if(wflux(ia,ja,ka-1,nsm).lt.0.d0) kb=ka-1
#else
       if(rbg*wflux(ka-1,nsp)+rb*wflux(ka-1,nsm).lt.0.d0) kb=ka-1
#endif              
       z1=dble(ka-1)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr)
       call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr)
#endif
#ifdef sediment
       if(kb==KM-kmt(ia,ja)) then
          nsed=nsed+1
          nrj(ntrac,6)=2
          trj(ntrac,1)=x1
          trj(ntrac,2)=y1
          trj(ntrac,3)=z1
          trj(ntrac,4)=tt
          trj(ntrac,5)=subvol
          nrj(ntrac,1)=ib
          nrj(ntrac,2)=jb
          nrj(ntrac,3)=ka
          nrj(ntrac,4)=niter
          nrj(ntrac,5)=idint(ts)
          nrj(ntrac,7)=1
          !call writedata(13)
          !cycle ntracLoop
       endif
#endif
#if defined stream_thermohaline
! calculate the layers of temperature and salinity for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(ia,ja,ka-1,mrb,mta,mtb,msa,msb,3,-1,real(subvol*ff))

    elseif( ds==dsc .or. ds==dsmin) then  
       ! shortest time is the time-steping 
       scrivi=.true.
#ifdef timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds,rr)
       call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds,rr)
#else           
       ! If there is no spatial solution, 
       ! which should correspond to a convergence zone
       if(dse==UNDEF .and. dsw==UNDEF .and. dsn==UNDEF .and. & 
          dss==UNDEF .and. dsu==UNDEF .and. dsd==UNDEF ) then
          
          ! move if atmosphere, freeze if ocean
          ib=ia ; jb=ja ; kb=ka
!          print *,'convergence for ',ib,jb,kb,x0,y0,z0
!#ifdef ifs
!          call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr) ! zonal crossing 
!          call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr) ! merid. crossing 
!          call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vert. crossing 
!#else
!          x1=x0 ; y1=y0 ; z1=z0 
!          print *,ib,jb,kb,x1,y1,z1
!#endif  
          ! If there is at least one spatial solution 
          ! but the shortest cross time is the time step
       endif
!       else
          call pos_orgn(1,ia,ja,ka,x0,x1,ds,rr) ! zonal crossing 
          call pos_orgn(2,ia,ja,ka,y0,y1,ds,rr) ! merid. crossing 
          call pos_orgn(3,ia,ja,ka,z0,z1,ds,rr) ! vert. crossing 
!       endif
#endif
    endif
    
!    
!! This is just a try and needs to be implemented and tested thoroughly
!! It will neet do be implemented for all varbottombox
!#ifdef baltix
!! depth conversion for bottom box
!#ifdef varbottombox
!    if(  (ds==dse .or. ds==dsw .or. ds==dsn .or. ds==dss)  .and.  &
!         (ka==KM+1-kmt(ia,ja) .or. kb==KM+1-kmt(ib,jb))             ) then
!#ifdef zgrid3Dt 
!        dza=dz(ka)
!        dzb=dz(kb)
!        if(ka==KM+1-kmt(ia,ja)) dza=dztb(ia,ja,1)
!        if(kb==KM+1-kmt(ib,jb)) dzb=dztb(ib,jb,1)
!#elif  zgrid3D
!        dza=dzt(ia,ja,ka)
!        dzb=dzt(ib,jb,kb)
!#else
! stop 4967
!#endif /*zgrid3Dt*/
!       if(dza.ne.dzb) then
!		zz=dble(int(z1))+1.d0 - (1.d0-z1+dble(int(z1)))*dza/dzb 
!	
!		if( zz.le.dble(int(z1)) .or. zz.ge.dble(int(z1)+1) ) then
!	 		print *, 'fel',zz,z1,dble(int(z1))
!	 		print *,'z0,z1=',z0,z1
!	 		print *,'nytt z1=',zz
!	 		print *, 'dz',dza,dzb,dza/dzb
!	 		print *, 'ds',ds,dse,dsw,dsn,dss
!	 		print *, 'kmt',kmt(ia,ja),kmt(ib,jb),ka,kb
!	 		print *, 'distance from top of a box in m',(1.-(z1-int(z1)))*dza
!	 		print *, 'distance from top of b box in m',(1.-(zz-int(zz)))*dzb
!	 		print *, 'ska vara mellan 0 och 1',(1.-z1+dble(int(z1)))*dza/dzb
!	 		print *, 'ia,ib,ja,jb=',ia,ib,ja,jb
!	 		print *, 'x0,x1,y0,y1=',x0,x1,y0,y1
!	 		
!		 	stop 4956
!		endif
!	
!		z1=zz
!	endif
!	endif
!#endif /*varbottombox*/
!#endif /*baltix*/

    
  end subroutine pos
  

end module mod_pos
