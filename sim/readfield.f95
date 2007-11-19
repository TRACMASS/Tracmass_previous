!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine readfields

USE mod_param
USE mod_coord
USE mod_time
USE mod_grid
USE mod_name
USE mod_vel
USE mod_dens
IMPLICIT none

#ifdef tempsalt
real tempb(km),saltb(km),rhob(km)
integer kmm
#endif

INTEGER Iend,jend,kend,m,kk,iints,i,j,k,im,jm

REAL :: Time,TempdS,stlat,stlon
REAL :: uu,vv,ss,tt,zz,dhu,dhv
!real ust(imt,jmt,km,2)
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ddz
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depth
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: dxyz
REAL, SAVE :: dxd

CHARACTER :: FileName*128,TName*14,SName*11,VName*11,hour*3,month*2,day*4

LOGICAL around

 if ( .not. allocated (ddz) ) then
     allocate ( ddz(IMT,JMT,KM) )
  end if

 if ( .not. allocated (depth) ) then
     allocate ( depth(IMT,JMT) )
     allocate ( dxyz(KM) )
  end if

ihour=ihour+int(ff)
if(ihour.eq.25) then
 ihour=1
 iday=iday+1
 if(iday.gt.idmax(imon,iyear)) then
  iday=1
  imon=imon+1
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
  endif
 endif
elseif(ihour.eq.0) then
 ihour=24
 iday=iday-1
 if(iday.eq.0) then
  imon=imon-1
  if(imon.eq.0) then
   imon=12
   iyear=iyear-1
  endif
  iday=idmax(imon,iyear)
 endif
endif


ntime=1000000*iyear+10000*imon+100*iday+ihour
!print 99,ihour,iday,imon,iyear,ints,ntime
99 format('Hour: ',i2,' Date: ',i2,'/',i2,' ',i4,' ints=',2i12)

if(mod(ints,30*24).eq.0) then
print *,'psi written for ints=',ints
call writepsi
endif

! swap between time steps
hs(:,:,1)=hs(:,:,2)
u(:,:,:,1)=u(:,:,:,2)
v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
tem(:,:,:,1)=tem(:,:,:,2)
sal(:,:,:,1)=sal(:,:,:,2)
rho(:,:,:,1)=rho(:,:,:,2)
#endif

      open(45,file='../tmp/'//name)
 81   format(i1,i1,i1)
 82   format(i1,i2)
 83   format(i3)
 84   format(i2)
 89   format(3a)

      i=24*(iday-1)+ihour
      rewind(45)
      if(i.ge.100) then
       write(45,83) i
      elseif(i.ge.10) then
       write(45,82) 0,i
      elseif(i.lt.10) then
       write(45,81) 0,0,i
      endif
      rewind(45)
      read (45,89) hour

      i=imon
      rewind(45)
      if(i.ge.10) then
       write(45,84) i
      elseif(i.lt.10) then
       write(45,81) 0,i
      endif
      rewind(45)
      read (45,89) month

#ifdef sim
FileName=directory//'gcm/81'//month//'/APY'//hour//'.fld'
#endif
#ifdef for
FileName=directory//'gcm/88'//month//'/APL'//hour//'.fld'
#endif

!      print 29,ihour,iday,imon,iyear,ints,FileName
 29   format('Time: ',i2,'h ',i2,'/',i2,i5,' ints=',i10,2x,a128)

inquire(file=FileName,exist=around)
if(.not.around) then
 print *,'data missing',FileName
 stop 4556
endif

! initialise
if(ints.eq.intstart) then
!open(unit=28,file=directory//'trint'//name)
! call coordinat
 dxd=dx*deg  ! grid length in m
 zw=0.
 dxyz=dxd*dz   ! grid wall in m2
 hs=0.
 kmt=0
 depth=0.
 ddz=0.
 u=0.
 v=0.
! ust=0.
#ifdef tempsalt
 tem=0.
 sal=0.
 rho=0.
#endif
endif

! read new data
open(unit=50,file=FileName,status='old',err=4000)
read(50,88) Time,i,iend,j,jend,kend,TempdS,stlat,stlon,TName,SName,VName,(zw(km+1-k),k=2,km+1)
88 format(f8.1,5i5,3f8.1,a14,2a11,/,12f6.1,/,12f6.1,/,12f6.1,/,3f6.1)
!print 88,Time,i,iend,j,jend,kend,TempdS,stlat,stlon,TName,SName,VName,(zw(km+1-k),k=2,km+1)
!stop 2345
 100  continue
      read(50,*) i,j,zz

      if(i.ne.999) then
       depth(i,j)=zz
       goto 100
      endif

 200  continue

      read(50,*) i,j,kend,zz
      if(i.ne.999) then
      hs(i,j,2)=zz  ! sea surface height in m
      kmt(i,j)=kend
      do k=1,kend
       read(50,*) uu,vv,tt,ss
       kk=km+1-k         
       u(i,j,kk,2)=uu*0.01  ! convert from cm/s to m/s
       v(i,j,kk,2)=vv*0.01
       tem(i,j,kk,2)=tt
       sal(i,j,kk,2)=ss
!       if(i.eq.200 .and.j.eq.200) print *,kend,uu,vv,tt,ss,zz
      enddo

       goto 200
      endif     

!      print *,'slutlaest'

      close(50)

 1000 continue
	

 4000	continue

! transports
! grid box volumes in m3
      do i=1,imt
       do j=1,jmt
        do k=1,kmt(i,j)
         kk=km+1-k    
         ddz(i,j,kk)=dz(kk)
         if(k.eq.kmt(i,j)) then
!          ddz(i,j,kk)=depth(i,j)-zw(kk+1)
          ddz(i,j,kk)=depth(i,j)-zw(kk)
          dztb(i,j,1)  =depth(i,j)-zw(kk)
!         if(ddz(i,j,kk).le.1.) 
!        print *,'vad?',i,j,kk,ddz(i,j,kk),depth(i,j),zw(kk+1),kmt(i,j)
         endif
         if(k.eq.1) ddz(i,j,kk)=ddz(i,j,kk)+hs(i,j,2)
         if(ddz(i,j,kk).lt.0.) then
            print *,i,j,k,kk,ddz(i,j,kk),dz(kk)
            stop 2345
         endif
        enddo
       enddo
      enddo

!trint=0.
do i=1,imt-1
 do j=1,jmt-1
  do k=km+1-kmt(i,j),km
   u(i,j,k,2)=u(i,j,k,2)*dxd*0.5*(ddz(i,j,k)+ddz(i+1,j,k))
   v(i,j,k,2)=v(i,j,k,2)*dxd*0.5*(ddz(i,j,k)+ddz(i,j+1,k))
!   ust(i,j,k,1)=ust(i,j,k,1)+u(i,j,k,2)
!   ust(i,j,k,2)=ust(i,j,k,2)+v(i,j,k,2)
!   trint=trint+tra(i,j,k)*ddz(i,j,k)
  enddo
 enddo
enddo

!print 98,ints,trint
!write(28,98) ints,trint
!98 format('tr',i5,1x,f11.0)

! integrate tracer


!print *,FileName
!open(21,file='/Volumes/sjo5/data/skb/simp/uv.euler',form='unformatted')
!write(21) ust
!close(21)
!stop 3567
 
deallocate ( ddz )


return
end subroutine readfields


      
