!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine readfields

IMPLICIT none
# include "../param.h"

common/vel/u(imt,0:JMAX,km,NST),v(imt,0:JMAX,km,NST),hs(imt,JMAX,NST),w(0:km),ff
real u,v,hs
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
CHARACTER(LEN=27) :: directory


      integer nt,ktot,MAX_EVER_SNAPS,NLEN,dtts

      parameter (nt=2,ktot=214,MAX_EVER_SNAPS=(nt+2)*1+3,NLEN=28237)
      real snap2d(imt,jmt),vmask(4),snap1d(NLEN)
      real*4 vmask4(4)
      integer i,j,k,m,inum,iday,num,iyear0,imon0,iday0,kz,ii,ints2,kk
integer kmu(imt,jmt)
integer idmax(12,80:96),ima(80:96),ime(80:96),idmin(12,80:96)
character ofile*20,infile*44,zfile*123,rfile*39
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
                31,29,31,30,31,30,31,31,30,31,30,31/   ! 1996

data idmin/1, 1, 1, 1,26, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & 
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/ &
          ima/5,16*1/,ime/17*12/,  &
          a_exp1/'246'/, a_exp2/'57'/,dtts/150./ 

      integer ittstart,itt

      save kmu,vmask,vmask4

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

open (unit=86,file='/Volumes/sjo5/data/rco/kmt/kmt',form='unformatted')
read (86) kmt
close (86)

! === Create velocity grid ===
do i=1,imt-1
 do j=1,jmt-1
  kmu(i,j)=min(kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
 enddo
enddo

! === prepare for landmask === 
vmask=0.
vmask4=0.
 
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

      iday=-1
!     ittinit =  21024        ! 19800526
      ittstart=  22176         ! 19800528    
!     ittstart=2401056         ! 19910918   
!      ints2=mod(ints-1,2452)+1

      ints2=ints
  666  continue
      if(ints2.gt.intmax .or. ints2.lt.intmin) then
       ints2=ints2-intend+intstart-intstep
       goto 666
      endif


!      itt=(ints2-1)*2*24*nint(3600./dtts)+ittstart
      itt=ints2*24*nint(3600./dtts)+ittstart
!      print *,ints,itt
      do iyear0=80,96
       do imon0=ima(iyear0),ime(iyear0)
        do iday0=idmin(imon0,iyear0),idmax(imon0,iyear0)
         iday=iday+1
         num=iday*24*nint(3600./dtts)+21024
         if(itt.eq.num) goto 7890
        enddo
       enddo
      enddo
 7890 continue

      inum=itt
!      print 99,iday0,imon0,iyear0,itt
 99   format('Date: ',i2,'/',i2,' 19',i2,' itt=',i10)
                     
!     === Create filenames for the snap-files to be used ===
                     
      ofile='d0000000000.snap1'
      if(inum.lt.10) then
       write(ofile(11:11),'(i1)') inum
      else if(inum.lt.100) then
       write(ofile(10:11),'(i2)') inum
      else if(inum.lt.1000) then
       write(ofile(9:11),'(i3)') inum
      else if(inum.lt.10000) then
       write(ofile(8:11),'(i4)') inum
      else if(inum.lt.100000) then
       write(ofile(7:11),'(i5)') inum
      else if(inum.lt.1000000) then
       write(ofile(6:11),'(i6)') inum
      else if(inum.lt.10000000) then
       write(ofile(5:11),'(i7)') inum
      else if(inum.lt.100000000) then
       write(ofile(4:11),'(i8)') inum
      else if(inum.lt.1000000000) then
       write(ofile(3:11),'(i9)') inum
      else
       write(ofile(2:11),'(i10)') inum
      endif

infile='/Volumes/sjo5/data/rco/old/'//ofile
!print *,infile
zfile='uncompress -c '//infile//'.Z > /Volumes/sjo5/data/rco/old/tmp/'//name
!print *,zfile
CALL system(zfile)
rfile='/Volumes/sjo5/data/rco/old/tmp/'//name
!print *,rfile
inquire(file=rfile,exist=around)
if(.not.around) stop 4556

!     === open snap file ===
      open(unit=30,file=rfile,status='old',form='unformatted',err=4000)

! === Read sealevel ===
      k=1
      read (30) snap1d  
      ii=0
      do j=1,jmt
       do i=1,imt
        hs(i,j,2)=vmask4(1)
        if(kmt(i,j).ge.1) then
         ii=ii+1
         hs(i,j,2)=0.01*snap1d(ii)
        endif
       enddo
      enddo
 
! === Read ubt,vbt ===
      do k=2,3
       read (30) snap1d
      enddo

! === Read fluxes ===
                
      do k=4,50                 
       read (30) snap1d                
      enddo

! === Read 3D-fields ===
                
      k=51
      do 5000 kz=1,km
       kk=km+1-kz
                  
! === Read Temp
       k=k+1
       read (30) snap1d
#ifdef tempsalt
       ii=0
       do j=1,jmt
        do i=1,imt
         snap2d(i,j)=vmask4(1)
         if(kmt(i,j).ge.kz) then
          ii=ii+1
!          if(snap1d(ii).lt.-3.) then
!           print *,'tempproblem',i,j,ii,snap1d(ii)
!           stop 4968
!          endif
          snap2d(i,j)=snap1d(ii)
         else
          snap2d(i,j)=0.
         endif
         tem(i,j,kk,2)=snap2d(i,j)
        enddo
       enddo
#endif

! === Read Sal (3D) ===
       k=k+1
       read (30) snap1d  
#ifdef tempsalt
       ii=0
       do j=1,jmt
        do i=1,imt
         snap2d(i,j)=vmask4(1)
         if(kmt(i,j).ge.(kz)) then
          ii=ii+1
!          if(snap1d(ii).gt.45.) then
!           print *,'saltproblem',i,j,ii,snap1d(ii),tem(i,j,kk,2)
!           stop 4969
!          endif
          snap2d(i,j)=snap1d(ii)
         endif
         if(snap2d(i,j).lt.0.) snap2d(i,j)=0.  ! no overshoot neg sal
         sal(i,j,kk,2)=snap2d(i,j)
        enddo
       enddo  
#endif

! === Real u (3D) ===

       k=k+1
       read (30) snap1d
       ii=0
       do j=1,jmt
        do i=1,imt
         snap2d(i,j)=vmask4(1)
         if(kmu(i,j).ge.kz) then
          ii=ii+1
          snap2d(i,j)=snap1d(ii)
         else
          snap2d(i,j)=0.
         endif
! u -> transport
         if(i.ne.imt .and. j.ne.1) then
          if(kk.ne.km) then
           u(i,j,kk,2)=0.005*dy*deg*(snap2d(i,j)+snap2d(i,j-1))*dz(kk)
          else
u(i,j,kk,2)=0.005*dy*deg*(snap2d(i,j)+snap2d(i,j-1))*( dz(kk)+0.5*(hs(i,j,2)+hs(i+1,j,2)) )
!         if(i.eq.imt/2 .and. j.eq.jmt/2  .and. kk.eq.km) print *,u(i,j,kk,1)
          endif
         endif
        enddo
       enddo

! === Read v (3D) ==
      k=k+1
      read (30) snap1d
      ii=0
      do j=1,jmt
       do i=1,imt
        snap2d(i,j)=vmask4(1)
        if(kmu(i,j).ge.(kz)) then
         ii=ii+1
         snap2d(i,j)=snap1d(ii)
        else
         snap2d(i,j)=0.
        endif
! v -> transport
        if(i.ne.1 .and. j.ne.jmt) then
         if(kk.ne.km) then
v(i,j,kk,2)=0.005*dx*deg*csu(j)*(snap2d(i,j)+snap2d(i-1,j))*dz(kk)
         else
v(i,j,kk,2)=0.005*dx*deg*csu(j)*(snap2d(i,j)+snap2d(i-1,j))*( dz(kk)+0.5*(hs(i,j,2)+hs(i,j+1,2)) )
         endif
        endif
       enddo
      enddo
                  
!     === end of k loop u and v ===
5000  continue

        
4000  continue


!        zfile='compress '//infile
!        CALL system(zfile)

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
      
