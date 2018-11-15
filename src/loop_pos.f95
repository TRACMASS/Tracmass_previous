module mod_pos
  USE mod_precdef
  USE mod_param
  USE mod_grid
  USE mod_vel
  USE mod_loopvars
  USE mod_time, only: intrpb, intrpbg
  USE mod_streamfunctions
  USE mod_psi
  USE mod_tempsalt
  USE mod_traj, only: ntrac

  IMPLICIT none
  
contains
  
  subroutine  pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)
    
    INTEGER                                    :: mra,mta,msa
    INTEGER                                    :: mrb,mtb,msb
    REAL                                       :: uu
    INTEGER                                    :: ia, iam, ja, ka,k
    INTEGER                                    :: ib, jb, kb
    REAL                                       :: temp,salt,dens
    REAL(DP)                                   :: dza,dzb, zz
    REAL(DP), INTENT(IN)                       :: x0, y0, z0
    REAL(DP), INTENT(OUT)                      :: x1, y1, z1
      
    ! === calculate the new positions ===
    ! === of the trajectory           ===
    scrivi=.false.
    if(ds==dse) then ! eastward grid-cell exit 
       scrivi=.false.
       uu=(intrpbg*uflux(ia,ja,ka,nsp)+intrpb*uflux(ia ,ja,ka,nsm))*ff
       if(uu.gt.0.d0) then
          ib=ia+1
          if(ib.gt.IMT) ib=ib-IMT 
       endif
       x1=dble(ia)
#if defined timeanalyt
       call pos_time(2,ia,ja,ka,y0,y1)
       call pos_time(3,ia,ja,ka,z0,z1)
#else
       call pos_orgn(2,ia,ja,ka,y0,y1,ds) 
       call pos_orgn(3,ia,ja,ka,z0,z1,ds)
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
       !PRINT *, ia, ja, ka, temp, salt
#endif 

       call savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,1,1,real(subvol*ff))
        
    elseif(ds==dsw) then ! westward grid-cell exit
       scrivi=.false.
       uu=(intrpbg*uflux(iam,ja,ka,nsp)+intrpb*uflux(iam,ja,ka,nsm))*ff
       if(uu.lt.0.d0) then
          ib=iam
       endif
       x1=dble(iam)
#if defined timeanalyt
       call pos_time(2,ia,ja,ka,y0,y1)
       call pos_time(3,ia,ja,ka,z0,z1)
#else
       call pos_orgn(2,ia,ja,ka,y0,y1,ds) ! meridional position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds) ! vertical position
#endif
!       scrivi=.true.      
#if defined streamr 
       call interp2(ib,jb,kb,temp,salt,dens)
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       !PRINT *, '2', ia, ja, ka, temp, salt
#endif 
       call savepsi(iam,ja,ka,mrb,mta,mtb,msa,msb,1,-1,real(subvol*ff))

    elseif(ds==dsn) then ! northward grid-cell exit
       scrivi=.false.
       uu=(intrpbg*vflux(ia,ja,ka,nsp)+intrpb*vflux(ia,ja,ka,nsm))*ff
       if(uu.gt.0.d0) then
          jb=ja+1
       endif
       y1=dble(ja)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1)
       call pos_time(3,ia,ja,ka,z0,z1)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds) ! zonal position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds) ! vertical position
#endif
#if defined streamr 
       call interp2(ib,jb,kb,temp,salt,dens)
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       uu=(intrpbg*vflux(ia,ja-1,ka,nsp)+intrpb*vflux(ia,ja-1,ka,nsm))*ff
       if(uu.lt.0.d0) then
          jb=ja-1
#ifndef atmospheric 
        if(jb==0) stop 34578
#endif
       endif
       y1=dble(ja-1)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1)
       call pos_time(3,ia,ja,ka,z0,z1)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds) ! zonal position
       call pos_orgn(3,ia,ja,ka,z0,z1,ds) ! vertical position
#endif
#if defined streamr 
       call interp2(ib,jb,kb,temp,salt,dens)
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
       call vertvel(ia,iam,ja,ka)
#if defined explicit_w || full_wflux
       uu=wflux(ia,ja,ka,nsm)
#else
       uu=intrpbg*wflux(ka,nsp)+intrpb*wflux(ka,nsm)
#endif
       if(uu.gt.0.d0) then
          kb=ka+1
       endif
       z1=dble(ka)
   
       if(kb==KM+1) then    ! prevent particles to cross the sea surface 
          kb=KM             
#if defined hydro 
          z1=dble(KM)       ! put them exactly at the surface for hydro
#else
          z1=dble(KM)-0.5d0 ! ! put them in the midle of the surface layer 
#endif
       endif
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1)
       call pos_time(2,ia,ja,ka,y0,y1)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds)
       call pos_orgn(2,ia,ja,ka,y0,y1,ds)
#endif
! SARA
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
        !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
        !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
        !call interp(ia,ja,ka,x1,y1,z1,temp,salt,dens,nsm)
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
       call vertvel(ia, iam, ja, ka)
     
#if defined explicit_w || full_wflux
       if(wflux(ia,ja,ka-1,nsm).lt.0.d0) kb=ka-1
#else
       if(intrpbg*wflux(ka-1,nsp)+intrpb*wflux(ka-1,nsm).lt.0.d0) kb=ka-1
#endif              
       z1=dble(ka-1)
#if defined timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1)
       call pos_time(2,ia,ja,ka,y0,y1)
#else
       call pos_orgn(1,ia,ja,ka,x0,x1,ds)
       call pos_orgn(2,ia,ja,ka,y0,y1,ds)
#endif
#ifdef sediment
       if(kb==KM-kmt(ia,ja)) then
          nsed=nsed+1
          nrj(ntrac,6)=2
          trj(1,ntrac)=x1
          trj(2,ntrac)=y1
          trj(3,ntrac)=z1
          trj(4,ntrac)=tt
          trj(5,ntrac)=subvol
          nrj(1,ntrac)=ib
          nrj(2,ntrac)=jb
          nrj(3,ntrac)=ka
          nrj(4,ntrac)=niter
          nrj(5,ntrac)=idint(ts)
          nrj(7,ntrac)=1
          !call writedata(13)
          !cycle ntracLoop
       endif
#endif

! SARA
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1)
       call interp2(ib,jb,kb,temp,salt,dens)
        !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
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
! calculate the layers of temp and salt for both a-box and b-box
       call interp2(ib,jb,kb,temp,salt,dens)
        !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,nsm)
       mtb=int((temp-tmin)/dtemp)+1
       if(mtb.lt.1 ) mtb=1
       if(mtb.gt.MR) mtb=MR
       msb=int((salt-smin)/dsalt)+1
       if(msb.lt.1 ) msb=1
       if(msb.gt.MR) msb=MR
       call interp2(ia,ja,ka,temp,salt,dens)
 !       call interp(ia,ja,ka,x1,y1,z1,temp,salt,dens,nsm)
       mta=(temp-tmin)/dtemp+1
       if(mta.lt.1 ) mta=1
       if(mta.gt.MR) mta=MR
       msa=(salt-smin)/dsalt+1
       if(msa.lt.1 ) msa=1
       if(msa.gt.MR) msa=MR
#endif 
       call savepsi(ia,ja,ka-1,mrb,mta,mtb,msa,msb,3,-1,real(subvol*ff))

    elseif( ds==dsc .or. ds==dsmin) then ! shortest time is the time-steping 
       scrivi=.true.
#ifdef timeanalyt
       call pos_time(1,ia,ja,ka,x0,x1)
       call pos_time(2,ia,ja,ka,y0,y1)
       call pos_time(3,ia,ja,ka,z0,z1)
       if(dble(int(x0))/=x0 .and. abs(int(x0)-int(x1))>=1.) then
        print *,'pangx=',x0,x1,y0,y1,z0,z1
!        print *,'ds',dse,dsw,dsn,dss,dsu,dsd,dsc,dsmin
        x1=x0 
     !   stop 3854
       endif
       if(dble(int(y0))/=y0 .and. abs(int(y0)-int(y1))>=1.) then
        print *,'pangy=',y0,y1
    !    stop 3854
         y1=y0  
       endif
       if(dble(int(z0))/=z0 .and. abs(int(z0)-int(z1))>=1.) then
        print *,'pangz=',ntrac,z0,z1
         z1=z0
        stop 3854
       endif
!       print *,'x0x1',x0,x1,y0,y1,z0,z1
#else           
       ! If there is no spatial solution, i.e. a convergence zone
       if(dse==UNDEF .and. dsw==UNDEF .and. dsn==UNDEF .and. & 
          dss==UNDEF .and. dsu==UNDEF .and. dsd==UNDEF ) then       
          ib=ia ; jb=ja ; kb=ka
       endif
!       else
          call pos_orgn(1,ia,ja,ka,x0,x1,ds) ! zonal crossing 
          call pos_orgn(2,ia,ja,ka,y0,y1,ds) ! merid. crossing 
          call pos_orgn(3,ia,ja,ka,z0,z1,ds) ! vert. crossing 
!       endif
#endif
    endif
    
  end subroutine pos
  

end module mod_pos
