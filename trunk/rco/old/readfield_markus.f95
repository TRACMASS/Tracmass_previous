!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine readfields

IMPLICIT none
# include "../param.h"

common/vel/u(imt,0:JMAX,km,NST),v(imt,0:JMAX,km,NST),hs(imt,JMAX,NST),w(0:km),ff,rand(6)
real u,v,hs,rand
real*8 w,ff

common /coord/dx,dy,deg,stlon1,stlat1,csu(jmt),cst(jmt),zw(0:km)
REAL*8        dx,dy,deg,stlon1,stlat1,csu,cst,zw

common /grid/dxdy(imt,jmt),dztb(imt,jmt),dz(km),rmin,dr,tmin,dtemp,smin,dsalt,kmt(imt,jmt)
REAL*8 dxdy,dztb,dz,rmin,dr,tmin,dtemp,smin,dsalt
INTEGER kmt

#ifdef tempsalt
common/dens/ tem(imt,jmt,km,2),sal(imt,jmt,km,2),rho(imt,jmt,km,2)
real tem,sal,rho
real tempb(km),saltb(km),rhob(km)
integer kmm
#endif

#ifdef time
common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
integer    ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
#endif

common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=23) :: directory

integer NLEN,dtts,NSNAPS,NSNAPS_2D,NSNAPS_3D,nd,l_tke,nlength

parameter (NLEN=29557,NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,NSNAPS=1,nlength=2379)
parameter (l_tke=0)

integer itt0,year,month,day,hour,minute,second,imt0,jmt0,km0,nt0,NLEN0,NSNAPS0
integer i,j,k,m,inum,num,kz,ii,ints2,kk,i0
integer iyear0,imon0,iday0,ihour0,iyear,imon,iday,ihour
integer kmu(imt,jmt),idmax(12,1980:2007)

real*8  ird0,ird20,ird30,ird40,dxa,dya,stlon,stlat,dxdeg,dydeg

real snapd,totsec

real*4 ird,ird2,ird3,ird4,rd2d(imt,jmt),rd1d_a(NSNAPS), rd1d_b(NSNAPS)
real*4 zdzz(km),dzw(0:km),dxt(imt),dyt(jmt),phi(jmt),phit(jmt),yu(jmt)
real*4  snap1d(NLEN)

real snap2d(imt,jmt)

character ofile*20,infile*48,zfile*123,rfile*39
character*3 a_exp1
character*2 a_exp2
logical around

data idmax/31,29,31,30,31,30,31,31,30,31,30,31, & ! 1980 
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1981
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1982
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1983
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 1984
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1985
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1986
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1987
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 1988
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1989
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1990
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1991
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 1992
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1993
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1994
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1995
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 1996
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1997
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1998
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 1999
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 2000
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 2001
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 2002
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 2003
           31,29,31,30,31,30,31,31,30,31,30,31, & ! 2004
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 2005
           31,28,31,30,31,30,31,31,30,31,30,31, & ! 2006
           31,28,31,30,31,30,31,31,30,31,30,31/   ! 2007


integer ittstart,itt

save kmu,inum,ihour,iday,imon,iyear,dxa,dya

ihour=ihour+6
if(ihour.eq.24) then
 ihour=0
 iday=iday+1
 if(iday.gt.idmax(imon,iyear)) then
  iday=1
  imon=imon+1
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
! if kan skrivas här om man vill börja om från 1981
  endif
 endif
endif

!____________________________ initialise ___________________________
if(ints.eq.intstart) then
hs=0.
u=0.
v=0.
#ifdef tempsalt
tem=0.
sal=0.
rho=0.
#endif

dya=0.005*dy*deg
dxa=0.005*dx*deg

!print *,directory//'kmt/kmt'
!open (unit=86,file=directory//'kmt/kmt',form='unformatted')
!read (86) kmt
!close (86)


! === Create velocity grid ===
!do i=1,imt-1
! do j=1,jmt-1
!  kmu(i,j)=min(kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
! enddo
!enddo

! starting date ...ska flyttas till main.f95 så småningom
ihour0=0
iday0=15
imon0=5
iyear0=1990

ihour=ihour0
iday=iday0
imon=imon0
iyear=iyear0

!inum=1981010100

endif

!_______________________________________________________________________

! swap between datasets

      do j=1,jmt
       do i=1,imt
        hs(i,j,1)=hs(i,j,2)
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

! write psi every year
      if(mod(ints,360).eq.0) then
       print *,'psi written for ints=',ints
       call writepsi
      endif

!=======================================================================
!     ===  Start main loop for extraction 
!=======================================================================

!      iday=-1
!      ittstart=  22176         ! 19800528    
      ints2=ints
  666  continue
      if(ints2.gt.intmax .or. ints2.lt.intmin) then
       ints2=ints2-intend+intstart-intstep
       goto 666
      endif


inum=1000000*iyear+10000*imon+100*iday+ihour
!print 99,ihour,iday,imon,iyear,ints,inum
99 format('Hour: ',i2,' Date: ',i2,'/',i2,' ',i4,' ints=',2i12)
                     
!     === Create filenames for the snap-files to be used ===
ofile='d0000000000.snap1'
write(ofile(2:11),'(i10)') inum

infile=directory//'2nm_122/'//ofile
print *,inum,infile
inquire(file=infile//'.gz',exist=around)
if(.not.around) stop 4555
zfile='gunzip -c '//infile//'.gz > '//directory//'tmp/'//name
!print *,zfile
CALL system(zfile)
rfile=directory//'tmp/'//name
!print *,rfile
inquire(file=rfile,exist=around)
if(.not.around) stop 4556

!     === open snap file ===
open(unit=30,file=rfile,status='old',form='unformatted',err=4000)

! start read header:
read(30) ird
itt0 = ird
read(30) ird
km0 = ird
if(km0.ne.km) stop 2844
read(30) ird
nt0 = ird
read(30) ird
imt0 = ird
if(imt0.ne.imt) stop 2845
read(30) ird
jmt0 = ird
if(jmt0.ne.jmt) stop 2846
read(30) ird
NLEN0 = ird
read(30) ird
NSNAPS0 = ird
read(30) ird
year = ird
read(30) ird
month = ird
read(30) ird
day = ird
read(30) ird
hour = ird
read(30) ird
minute = ird
read(30) ird
second = ird

read(30) ird0,ird20,ird30
dtts = ird0
totsec = ird20
snapd = ird30
read(30) ird0,ird20,ird30,ird40
dx = ird0
dy = ird20
dxdeg = ird30
dydeg = ird40
read(30) ird0,ird20
stlon1 = ird0
stlat1 = ird20

read(30) rd1d_a,rd1d_b
!do i=1,NSNAPS
! ispvar(i) = rd1d_a(i)
! isplev(i) = rd1d_b(i)
!enddo
read(30) rd2d

if(ints.eq.intstart) then

do i=1,imt
 do j=1,jmt
  kmt(i,j) = rd2d(i,j)
 enddo
enddo

!open(21,file=directory//'kmt/kmtny',form='unformatted')
!write(21)kmt
!close(21)

 do i=1,imt-1
  do j=1,jmt-1
   kmu(i,j)=min(kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
!   if(kmu(i,j).ne.0) then
!    hr(i,j) = 1./zmw(kmu(i,j))
!   else
!    hr(i,j) = 0.
!   endif
  enddo
 enddo
!print *,hr

call coordinat

endif

! ssh 
read(30,err=2000) ird
i0=nint(ird)
read(30,err=2000) snap1d(1:i0)
ii=0
do j=1,jmt
 do i=1,imt
  hs(i,j,2)=0.
  if(kmt(i,j).ge.1) then
   ii=ii+1
   hs(i,j,2)=0.01*snap1d(ii)
!   print *,i,j,hs(i,j,2),hs(i,j,1)
  endif
 enddo
enddo

! ubt
read(30,err=2000) ird
i0=nint(ird)
read(30,err=2000) snap1d(1:i0)
!ii=0
!do j=1,jmt
! do i=1,imt
!  ubt(i,j)=vmask(1)
!  if(kmu(i,j).ge.1) then
!   ii=ii+1
!   ubt(i,j)=snap1d(ii)/hr(i,j)
!  endif
! enddo
!enddo

! vbt
read(30,err=2000) ird
i0=nint(ird)
read(30,err=2000) snap1d(1:i0)
!ii=0
!do j=1,jmt
! do i=1,imt
!  vbt(i,j)=vmask(1)
!  if(kmu(i,j).ge.1) then
!   ii=ii+1
!   vbt(i,j)=snap1d(ii)/hr(i,j)
!  endif
! enddo
!enddo

! vad är detta?
do m=4,58
 read(30,err=2000) ird
 i0=nint(ird)
 read(30,err=2000) snap1d(1:i0)
! print*,'2d snapshots',m,i0,snap1d(1)
enddo

! temperature
do k=1,km
 read(30,err=2001) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2001) snap1d(1:i0)
 endif
! print*,'temperature',k,i0,snap1d(1)
#ifdef tempsalt
 kk=km+1-k
 ii=0
 do j=1,jmt
  do i=1,imt
   tem(i,j,kk,2)=0.
   if(kmt(i,j).ge.k) then
    ii=ii+1
    tem(i,j,kk,2)=snap1d(ii)
   endif
  enddo
 enddo
#endif
enddo
2001  continue

do k=1,km
 read(30,err=2002) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2002) snap1d(1:i0)
 endif
! print*,'salinity',k,i0,snap1d(1)
#ifdef tempsalt
 kk=km+1-k
 ii=0
 do j=1,jmt
  do i=1,imt
   sal(i,j,kk,2)=0.
   if(kmt(i,j).ge.k) then
    ii=ii+1
    sal(i,j,kk,2)=snap1d(ii)
   endif
  enddo
 enddo
#endif
enddo
2002  continue

! zonal velocity
do k=1,km
 kk=km+1-k
 read(30,err=2003) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2003) snap1d(1:i0)
 endif
! print*,'u3',k,i0,snap1d(1)
 ii=0
 do j=1,jmt
  do i=1,imt
   snap2d(i,j)=0.
   if(kmu(i,j).ge.k) then
    ii=ii+1
    snap2d(i,j)=snap1d(ii)
   endif
  enddo
 enddo
! u -> transports
 do j=2,jmt
  do i=1,imt-1
    if(kk.ne.km) then
     u(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))*dz(kk)
    else
     u(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))*( dz(kk)+0.5*(hs(i,j,2)+hs(i+1,j,2)) )
    endif
  enddo
 enddo

enddo
2003  continue

! meridional velocity
do k=1,km
 kk=km+1-k
 read(30,err=2004) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2004) snap1d(1:i0)
 endif
!print*,'v3',k,i0,snap1d(1)
 ii=0
 do j=1,jmt
  do i=1,imt
   snap2d(i,j)=0.
   if(kmu(i,j).ge.k) then
    ii=ii+1
    snap2d(i,j)=dxa*csu(j)*snap1d(ii)
   endif
  enddo
 enddo
! v -> transport
 do j=1,jmt-1
  do i=2,imt
   if(kk.ne.km) then
    v(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*dz(kk)
   else
    v(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*( dz(kk)+0.5*(hs(i,j,2)+hs(i,j+1,2)) )
   endif
  enddo
 enddo

enddo

!print *,'uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu',ints
!print *,(u(i,jmt/2,km,2),i=1,imt)
!print *,'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',ints
!print *,(v(i,jmt/2,km,2),i=1,imt)


2004  continue

goto 9000
    
4000  continue
2000  print*,'error 2000'
stop 2000
9000 continue

close(30) 


#ifdef tempsalt
! the density
      do i=1,imt
       do j=1,jmt
        if(kmt(i,j).ne.0) then

         kmm=kmt(i,j)
         do k=1,kmm
          kk=km+1-k
!          if(tem(i,j,kk,2).lt.-3. .or. sal(i,j,kk,2).gt.45.) then
!           print *,i,j,kk,kmm,tem(i,j,kk,2),sal(i,j,kk,2)
!     &                       ,tem(i,j,km,2),sal(i,j,km,2)
!     &                       ,tem(i,j,km-1,2),sal(i,j,km-1,2)
!     &                       ,tem(i,j,km-2,2),sal(i,j,km-2,2)
!           stop 5096
!          endif
          tempb(k)=tem(i,j,kk,2)
          saltb(k)=(sal(i,j,kk,2)-35.)/1000.
         enddo
         call statv(tempb,saltb,rhob,kmm)
         do k=1,kmm
          kk=km+1-k
          rho(i,j,kk,2)=rhob(k)
         enddo

        endif
       enddo
      enddo


      return
      end subroutine readfields
#endif

!_______________________________________________________________________
      
