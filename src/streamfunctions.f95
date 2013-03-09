
module mod_pos
  USE mod_param
  USE mod_grid
  USE mod_vel
  USE mod_loopvars
  USE mod_time
  USE mod_streamfunctions
  USE mod_psi
    
  IMPLICIT none
  
contains

subroutine pos_stream_funcs_1()
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,1)
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

     end subroutine pos_stream_funcs_1


     subroutine pos_stream_funcs_2()
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,1)
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
     end subroutine pos_stream_funcs_2


     subroutine pos_stream_funcs_3()
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,1)
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
     end subroutine pos_stream_funcs_3


     subroutine pos_stream_funcs_4()
#if defined streamr 
!       call interp(ib,jb,kb,x1,y1,z1,1)
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
     end subroutine pos_stream_funcs_4


     subroutine pos_stream_funcs_5()

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

     end subroutine pos_stream_funcs_5


     subroutine pos_stream_funcs_6()

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
     end subroutine pos_stream_funcs_6
