	From: 	  markus.meier@smhi.se
	Subject: 	Re: read rco
	Date: 	tisdag 26 sep 2006 14.44.03 GMT+02:00
	To: 	  doos@misu.su.se, hanna@misu.su.se

Hej Krisofer, hej Hanna!

Försök använda följande program.
Glöm inte "setenv F_UFMTENDIAN big" innan ni kör programmet.
Om ni behöver en jättebra läsprogramm återkommer gärna till mig.

Mvh, Markus


Kristofer Döös wrote:

Hej Markus!

Jag och Hanna håller på att försöka läsa snap-filerna i 2nm_122.
Dom är inte skrivna på samma sätt som de körningar vi haft tidigare.
Har du lite fortran-kod hur dom är skrivna?

K.



#
#
# %Z% SCCS module: %M%  version: %I%
#    Creation date: %G%
.SUFFIXES :
.SUFFIXES : .o .F
OFILES= read_slarvig.o
DOPTNS=
#
#========================================================
#Monolith ifc
#GPP=  $(DOPTNS)
#FC =    ifc  -I /opt/scali/include -O2
#INCLUDE=
#LDR = ifc 
#LIBS =   -L /opt/scali/lib -lfmpi -lmpi -lpthread  -lsvml -lPEPCF90
#========================================================
#Tornado
 GPP=  $(DOPTNS)
 FC =    ifort -O2 -assume byterecl -Nscampi
 INCLUDE=
 LDR = ifort
 LIBS =  -Nscampi
#========================================================
#
.F.o :
	$(FC) $(GPP)   $(INCLUDE)  -c $<

full: $(OFILES)
	$(LDR) $(OFILES)  $(LIBS) -o ./extract_files

# Remove intermediate files
clean:
	rm -f *.o  *.log
c
c     f90 -o extract extract_fluxes_rco2-134.F
c
      program extract

      parameter(imt=320,jmt=362,km=41,nt=4,NLEN=29557,
     &          NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,mstat=13,
     &          NSNAPS=1,nlength=2379,many3=6,many=4,
     &          salt0=17.0,
     &          l_tke=0,
     &          klev0=1,i_halocline=0)
      real*4 zdzz(km),zw(km),dz(km),dzw(0:km),
     &       dxt(imt),dyt(jmt),phi(jmt),phit(jmt),yu(jmt),
     &       cst(jmt),csu(jmt),dsalt(km),
     &       fvn(km),fvs(km),fue(km),fuw(km),fw(0:km),
     &                  hr(imt,jmt),
     &                  a_2d(imt,jmt),
     &                  b_2d(imt,jmt),vmask(4),
     &                  eta(imt,jmt),
     &                  ubt(imt,jmt),
     &                  vbt(imt,jmt),
     &                  temp(imt,jmt,km),
     &                  salt(imt,jmt,km),
     &                  tke(imt,jmt,km),
     &                  epsi(imt,jmt,km),
     &                  u3(imt,jmt,km),
     &                  v3(imt,jmt,km),
     &                  w3(imt,jmt,km),
     &                  tracer(imt,jmt,km,many),
     &                  transp1(km,nlength),
     &                  transp2(km,nlength),
     &                  transp3(km,nlength),
     &                  transp4(km,nlength),
     &                  vol_by2(nlength),vol_by5(nlength),
     &                  transp(many3,nlength),
     &                  s_transp(many3,nlength)
      real*4  ird,ird2,ird3,ird4
      real*4  snap1d(NLEN)

      real*8  ird0,ird20,ird30,ird40
      real*4  rd2d(imt,jmt)
      real*4  rd1d_a(NSNAPS), rd1d_b(NSNAPS)
      integer kmt(imt,jmt),kmu(imt,jmt),kpn(imt,jmt),
     &        ispvar(NSNAPS),isplev(NSNAPS),klev(imt,jmt),
     &        npts(imt,jmt)
      integer itt0,year,month,day,hour,minute,second,
     &imt0,jmt0,km0,nt0,NLEN0,NSNAPS0
      integer idmax(12)
      character tfile*120,infile*120,ofil*20,ifile*23
      logical around
      data vmask/1.e10,2.e10,3.e10,4.e10/
      data idmax/31,28,31,30,31,30,31,31,30,31,30,31/
      data ifile/'box_RCO2_3m.map'/

      tfile=
     &'topo_RCO2.asc'
c     tfile=
c    &'../topo/baltic_2nm/'//tfile
      print*,tfile
      open(unit=2,file=tfile,status='old',form='formatted')
      read(2,'(40i3)') kmt
      close(2)

      do i=1,imt
      do j=1,jmt
        if(kmt(i,j).gt.km) kmt(i,j)=km
      enddo
      enddo

      do i=1,imt-1
      do j=1,jmt-1
        kmu(i,j)=min(kmt(i,j),kmt(i+1,j),
     &               kmt(i,j+1),kmt(i+1,j+1))
      enddo
      enddo

      iunit=10
        open(unit=iunit,file=ifile,status='old',
     &       form='formatted')
        read(iunit,200) kpn
200     format(40i3)
        close(iunit)

      do i=1,imt
      do j=1,jmt
        if(kpn(i,j).eq.0 .and. kmt(i,j).ne.0) then
          if(kpn(i-1,j  ).ne.0) kpn(i,j)=kpn(i-1,j  )
          if(kpn(i  ,j-1).ne.0) kpn(i,j)=kpn(i  ,j-1)
          if(kpn(i+1,j  ).ne.0) kpn(i,j)=kpn(i+1,j  )
          if(kpn(i  ,j+1).ne.0) kpn(i,j)=kpn(i  ,j+1)
          if(kpn(i,j).eq.0) print*,'error in kpn at',i,j
        endif
      end do
      end do

      dxdeg_1 = 6.6666667E-02
      dydeg_1 = 3.3333333E-02
      stlon_1 = 9.016700
      stlat_1 = 53.84163

      radius = 6370.e5
      pi=4.*atan(1.)
      radian = 360./(2.*pi)

      do i=1,imt
        dxt(i)=dxdeg_1*radius/radian
      enddo
      do j=1,jmt
        dyt(j)=dydeg_1*radius/radian
      enddo
      do j=1,jmt
        yu(j)=stlat_1+j*dydeg_1
      enddo
      phi(1) =yu(1)/radian
      phit(1)=phi(1)-0.5*dydeg_1/radian
      do j=1,jmt-1
        phi(j+1) = (yu(1)+j*dydeg_1)/radian
      enddo
      do j=2,jmt
        phit(j) = 0.5*(phi(j-1)+phi(j))
      enddo
      do j=1,jmt
        cst(j)  = cos(phit(j))
        csu(j)  = cos(phi (j))
      enddo

      dz( 1)=   3.E2
      dz( 2)=   3.E2
      dz( 3)=   3.E2
      dz( 4)=   3.E2
      dz( 5)=   3.E2
      dz( 6)=   3.E2
      dz( 7)=   3.E2
      dz( 8)=   3.E2
      dz( 9)=   3.E2
      dz(10)=   3.E2
      dz(11)=   3.E2
      dz(12)=   3.E2
      dz(13)=   3.E2
      dz(14)=   3.007080E2
      dz(15)=   3.063581E2
      dz(16)=   3.175872E2
      dz(17)=   3.342542E2
      dz(18)=   3.561495E2
      dz(19)=   3.829976E2
      dz(20)=   4.144610E2
      dz(21)=   4.501440E2
      dz(22)=   4.895979E2
      dz(23)=   5.323265E2
      dz(24)=   5.777925E2
      dz(25)=   6.254241E2
      dz(26)=   6.746222E2
      dz(27)=   7.247683E2
      dz(28)=   7.752317E2
      dz(29)=   8.253778E2
      dz(30)=   8.745760E2
      dz(31)=   9.222075E2
      dz(32)=   9.676735E2
      dz(33)=   10.10402E2
      dz(34)=   10.49856E2
      dz(35)=   10.85539E2
      dz(36)=   11.17002E2
      dz(37)=   11.43851E2
      dz(38)=   11.65746E2
      dz(39)=   11.82413E2
      dz(40)=   11.93642E2
      dz(41)=   11.99292E2

      dzw(0) = 0.5*dz(1)
      do k = 1,km-1
        dzw(k) = 0.5*(dz(k) + dz(k+1))
      end do
      dzw(km)  = 0.5*dz(km)

      zw(1) = dz(1)
      do k = 2,km
        zw(k) = zw(k-1) + dz(k)
      end do

      do i=1,imt-1
      do j=1,jmt-1
        if(kmu(i,j).ne.0) then
          hr(i,j) = 1./zw(kmu(i,j))
        endif
      enddo
      enddo

      zdzz(1)=0.5*dz(1)
      do k=2,km
        zdzz(k)= zdzz(k-1)+dzw(k-1)
      end do

      do j=1,jmt
      do i=1,imt
        a_2d(i,j)=vmask(1)
        npts(i,j)=0
        if(kmt(i,j).ge.1) then
          a_2d(i,j)=0.
        endif
      enddo
      enddo

      num=0
      iyeara=1999
      iyeare=iyeara
      do iyear=iyeara,iyeare
        if(iyear.eq.iyeara) then
          ima=12
        else
          ima=1
        endif 
        if(iyear.eq.iyeare) then
          ime=12
        else
          ime=12
        endif
      do imon=ima,ime
        idmax0=idmax(imon)
        if(mod(iyear,4).eq.0.and.imon.eq.2) idmax0=29
      do iday=1,idmax0
      do ihour=0,18,6
      ofil='d0000000000.snap1'
      write(ofil(2:5),'(i4)') iyear
      if(imon.lt.10) then
        write(ofil(7:7),'(i1)') imon
      else
        write(ofil(6:7),'(i2)') imon
      endif
      if(iday.lt.10) then
        write(ofil(9:9),'(i1)') iday
      else
        write(ofil(8:9),'(i2)') iday
      endif
      if(ihour.lt.10) then
        write(ofil(11:11),'(i1)') ihour
      else
        write(ofil(10:11),'(i2)') ihour
      endif

      infile=
     &'/nobackup/rossby3/sm_meier/rco2/work_rco2_baltic_2nm_122/'//ofil
      print*,infile
      inquire(file=infile,exist=around)
      if(.not.around) goto 3700
      print*,infile
      open(unit=1,file=infile,status='old',
     &     form='unformatted')
      num=num+1
c     start read header:
      read(1) ird
      itt0 = ird
      read(1) ird
      km0 = ird
      read(1) ird
      nt0 = ird
      read(1) ird
      imt0 = ird
      read(1) ird
      jmt0 = ird
      read(1) ird
      NLEN0 = ird
      read(1) ird
      NSNAPS0 = ird
      read(1) ird
      year = ird
      read(1) ird
      month = ird
      read(1) ird
      day = ird
      read(1) ird
      hour = ird
      read(1) ird
      minute = ird
      read(1) ird
      second = ird
c
      read(1) ird0,ird20,ird30
      dtts = ird0
      totsec = ird20
      snapd = ird30
      read(1) ird0,ird20,ird30,ird40
      dx = ird0
      dy = ird20
      dxdeg = ird30
      dydeg = ird40
      read(1) ird0,ird20
      stlon = ird0
      stlat = ird20
c
      read(1) rd1d_a,rd1d_b
      do i=1,NSNAPS
         ispvar(i) = rd1d_a(i)
         isplev(i) = rd1d_b(i)
      enddo
      read(1) rd2d
      do i=1,imt
         do j=1,jmt
            kmt(i,j) = rd2d(i,j)
         enddo
      enddo
c     end read header:

        read(1,err=2000) ird
        i0=nint(ird)
        read(1,err=2000) snap1d(1:i0)
        ii=0
        do j=1,jmt
        do i=1,imt
          eta(i,j)=vmask(1)
          if(kmt(i,j).ge.1) then
            ii=ii+1
            eta(i,j)=snap1d(ii)
          endif
        enddo
        enddo

        read(1,err=2000) ird
        i0=nint(ird)
        read(1,err=2000) snap1d(1:i0)
        ii=0
        do j=1,jmt
        do i=1,imt
          ubt(i,j)=vmask(1)
          if(kmu(i,j).ge.1) then
            ii=ii+1
            ubt(i,j)=snap1d(ii)/hr(i,j)
          endif
        enddo
        enddo

        read(1,err=2000) ird
        i0=nint(ird)
        read(1,err=2000) snap1d(1:i0)
        ii=0
        do j=1,jmt
        do i=1,imt
          vbt(i,j)=vmask(1)
          if(kmu(i,j).ge.1) then
            ii=ii+1
            vbt(i,j)=snap1d(ii)/hr(i,j)
          endif
        enddo
        enddo

      do m=4,58
        read(1,err=2000) ird
        i0=nint(ird)
        read(1,err=2000) snap1d(1:i0)
c       print*,'2d snapshots',m,i0,snap1d(1)
      enddo

      do k=1,km
        read(1,err=2001) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2001) snap1d(1:i0)
        endif
c       print*,'temperature',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
          temp(i,j,k)=vmask(1)
          if(kmt(i,j).ge.k) then
            ii=ii+1
            temp(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2001  continue

      do k=1,km
        read(1,err=2002) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2002) snap1d(1:i0)
        endif
c       print*,'salinity',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
          salt(i,j,k)=vmask(1)
          if(kmt(i,j).ge.k) then
            ii=ii+1
            salt(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2002  continue

      if(l_tke.eq.1) then

      do k=1,km
        read(1,err=2006) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2006) snap1d(1:i0)
        endif
c       print*,'tke',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
          tke(i,j,k)=vmask(1)
          if(kmt(i,j).ge.k) then
            ii=ii+1
            tke(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2006  continue

      do k=1,km
        read(1,err=2007) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2007) snap1d(1:i0)
        endif
c       print*,'epsi',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
          epsi(i,j,k)=vmask(1)
          if(kmt(i,j).ge.k) then
            ii=ii+1
            epsi(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2007  continue

      endif

      do k=1,km
        read(1,err=2003) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2003) snap1d(1:i0)
        endif
c       print*,'u3',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
c         u3(i,j,k)=vmask(1)
          u3(i,j,k)=0.
          if(kmu(i,j).ge.k) then
            ii=ii+1
            u3(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2003  continue

      do k=1,km
        read(1,err=2004) ird
        i0=nint(ird)
        if(i0.gt.0) then
          read(1,err=2004) snap1d(1:i0)
        endif
c       print*,'v3',k,i0,snap1d(1)
        ii=0
        do j=1,jmt
        do i=1,imt
c         v3(i,j,k)=vmask(1)
          v3(i,j,k)=0.
          if(kmu(i,j).ge.k) then
            ii=ii+1
            v3(i,j,k)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
2004  continue

      do m=1,many
      do k=1,km
          read(1,err=2005) ird
          i0=nint(ird)
          if(i0.gt.0) then
            read(1,err=2005) snap1d(1:i0)
          endif
        ii=0
        do j=1,jmt
        do i=1,imt
          tracer(i,j,k,m)=vmask(1)
          if(kmt(i,j).ge.k) then
            ii=ii+1
            tracer(i,j,k,m)=snap1d(ii)
          endif
        enddo
        enddo
      enddo
 2005 continue
      enddo

      goto 3500
2000  print*,'error 2000'
3500  continue
      close(1)

3700  continue
      enddo
      enddo
      enddo
      enddo
4000  continue

      print*,'extraction ready! number of snapshots:',num

      end
