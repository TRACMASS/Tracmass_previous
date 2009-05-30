c23456789012345678901234567890123456789012345678901234567890123456789012

c_______________________________________________________________________
 
      subroutine readfields

      implicit none

#include "netcdf.inc"
#include "param.h"

      common/vel/u(imt,jmt,km,2),v(imt,jmt,km,2),w(0:km),ff
      real u,v
      real*8  w,ff

      common /coord/ dx,dy,yu,deg,csu,cst,dz,stlon,stlat
      real*8 dx,dy(jmt),yu(jmt),deg,cst(jmt),csu(jmt),dz(km),stlon,stlat

#ifdef tempsalt
      common/dens/ tem(imt,jmt,km,2),
     &             sal(imt,jmt,km,2),rho(imt,jmt,km,2)
      real tem,sal,rho
      real tempb(km),saltb(km),rhob(km)
      integer kmm
#endif

c      common/bathy/kmt(imt,jmt)
      integer kmt(imt,jmt)

      integer i,j,k,kk,im,ints,jnts,intstart,istat
#ifdef bolus
      real field(imt,jmt,km,6)
#else
      real field(imt,jmt,km,4)
#endif
c      character*46 datanamn
c      character*39 datanamn
      character*43 datanamn
c      data datanamn/'/net/navier/hav2/doos/trajec/paleo/data/eo2fc5'/
c      data datanamn/'/net/navier/hav2/doos/trajec/paleo/data/eo3fc5'/
       data datanamn/'/net/knudsen/export/home1/doos/paleo/modern'/
c      data datanamn/'/net/navier/zee2/doos/data/paleo/eo3fc5'/
c      data datanamn/'/net/navier/hav2/doos/trajec/paleo/data/modern'/

      real ax,ay
      integer id,idu,idv,idt,ids ! NetCDF identifiers
      integer start(4),count(4)
      save idu,idv,idt,ids 
#ifdef bolus
      integer iduiso,idviso
      save iduiso,idviso
#endif

      do j=1,jmt
       do i=1,imt
        do k=1,km
         u(i,j,k,1)=u(i,j,k,2)
         v(i,j,k,1)=v(i,j,k,2)
#ifdef tempsalt 
         tem(i,j,k,1)=tem(i,j,k,2)
         sal(i,j,k,1)=sal(i,j,k,2)
         rho(i,j,k,1)=rho(i,j,k,2)
#endif
        enddo
       enddo
      enddo
c      if(ff.eq.1.) then
       ints=jnts
c       if(ints.eq.nsmax+1) stop 2567
c      elseif(ff.eq.-1.) then
c       ints=nsmax-jnts+1
c       if(ints.eq.0) stop 2568
c      endif

c open the NetCDF datasets
      if (ints.eq.intstart) then
	print *,datanamn//'.u.nc'
       istat = nf_open(datanamn//'.u.nc' ,nf_nowrite, idu)
       istat = nf_open(datanamn//'.v.nc' ,nf_nowrite, idv)
       istat = nf_open(datanamn//'.t.nc' ,nf_nowrite, idt)
       istat = nf_open(datanamn//'.s.nc' ,nf_nowrite, ids)
#ifdef bolus
       istat = nf_open(datanamn//'.uiso.nc' ,nf_nowrite, iduiso)
       istat = nf_open(datanamn//'.viso.nc' ,nf_nowrite, idviso)
#endif

c        call readtopo

      endif


c read data
      start(1)=2 
      start(2)=1 
      start(3)=1
c      start(4)=1117-ints 
c      start(4)=817-ints 
      start(4)=ints 
      count(1)=imt 
      count(2)=jmt 
      count(3)=km
      count(4)=1 
      istat = nf_inq_varid(idu,"U",id)
      istat = nf_get_vara_real(idu,id,start,count,field(1,1,1,1))
      istat = nf_inq_varid(idv,"V",id)
      istat = nf_get_vara_real(idv,id,start,count,field(1,1,1,2))
      istat = nf_inq_varid(idt,"T",id)
      istat = nf_get_vara_real(idt,id,start,count,field(1,1,1,3))
      istat = nf_inq_varid(ids,"S",id)
      istat = nf_get_vara_real(ids,id,start,count,field(1,1,1,4))
#ifdef bolus
      istat = nf_inq_varid(iduiso,"UISO",id)
      istat = nf_get_vara_real(iduiso,id,start,count,field(1,1,1,5))
      istat = nf_inq_varid(idviso,"VISO",id)
      istat = nf_get_vara_real(idviso,id,start,count,field(1,1,1,6))
#endif
c      print *,'ints=',ints,' yr=',float(ints)/12.
c      if(mod(ints,12).eq.0) print *,'ints=',ints,' yr=',float(ints)/12.
c write kmt
c      do j=1,jmt
c       do i=1,imt
c       do k=1,km
c         if(field(i,j,k,3).lt.100.) kmt(i,j)=k
c        enddo
c       enddo
c      enddo
c      print *,kmt
c      open(86,file='../data/kmt.modern',form='unformatted')
c      write(86) kmt
c      close(86)
c      stop 2678

c transports
      do k=1,km
       kk=km+1-k
       do j=2,jmt
        ax=0.5*dy(j)*deg*dz(kk)*0.01
        ay=0.5*dx*deg*dz(kk)*csu(j)*0.01
c        ax=dy(j)*deg*dz(kk)*0.01
c        ay=dx*deg*dz(kk)*csu(j)*0.01
        do i=1,imt
c        if(field(i,j  ,k,2).ne.0.) print *,field(i,j  ,k,2)
         im=i-1
         if(im.eq.0) im=imt
         if(field(i,j  ,k,1).gt.1.e6) field(i,j  ,k,1)=0.
         if(field(i,j-1,k,1).gt.1.e6) field(i,j-1,k,1)=0.
         if(field(i ,j ,k,2).gt.1.e6) field(i ,j ,k,2)=0.
         if(field(im,j ,k,2).gt.1.e6) field(im,j ,k,2)=0.
         if(field(i,j-1,k,3).gt.1.e6) field(i,j-1,k,3)=0.
         if(field(i ,j ,k,3).gt.1.e6) field(i ,j ,k,3)=0.
         if(field(i ,j ,k,4).gt.1.e6) field(i ,j ,k,4)=0.
         if(field(i,j-1,k,4).gt.1.e6) field(i,j-1,k,4)=0.
#ifdef bolus
c         if(field(i,j-1,k,5).gt.1.e6) field(i,j-1,k,5)=0.
         if(field(i,j,k,5).gt.1.e6) field(i,j,k,5)=0.
c         if(field(im,j,k,6).gt.1.e6) field(im,j,k,6)=0.
         if(field(i,j,k,6).gt.1.e6) field(i,j,k,6)=0.
#endif
c         print *,field(i,j,k,1),field(i,j,k,2),
c     &           field(i,j,k,3),field(i,j,k,4)
         u(i,j,kk,2)=(field(i,j,k,1)+field(i,j-1,k,1))*ax
         v(i,j,kk,2)=(field(i,j,k,2)+field(im,j ,k,2))*ay
c         u(i,j,kk,2)=field(i,j,k,1)*ax
c         v(i,j,kk,2)=field(i,j,k,2)*ay
#ifdef bolus
c         u(i,j,kk,2)=u(i,j,kk,2)+field(i,j-1,k,5)*ax*2.
         u(i,j,kk,2)=u(i,j,kk,2)+field(i,j,k,5)*ax*2.
c         v(i,j,kk,2)=v(i,j,kk,2)+field(im,j,k,6)*ay*2.
         v(i,j,kk,2)=v(i,j,kk,2)+field(i,j,k,6)*ay*2.
#endif

#if defined tempsalt
         tem(i,j,kk,2)=field(i,j,k,3)
         sal(i,j,kk,2)=1000.*field(i,j,k,4)+35.
#endif
c         if(v(i,j,kk,2).ne.0.) print *,v(i,j,kk,2)
        enddo
       enddo
      enddo

#if defined tempsalt
c density
      do j=1,jmt
       do i=1,imt
        kmm=0
        do k=1,km
         tempb(k)=field(i,j,k,3)
         saltb(k)=field(i,j,k,4)
         if(tempb(k).ne.0.) kmm=k
        enddo
        if(kmm.ne.0) then
         call stat(tempb,saltb,rhob,kmm)
         do k=1,kmm
          kk=km+1-k
          rho(i,j,kk,2)=rhob(k)
         enddo
        endif
       enddo
      enddo
#endif

      return
      end

c_______________________________________________________________________
