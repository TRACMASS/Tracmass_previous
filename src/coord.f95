!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine coordinat

USE mod_param
USE mod_grid
USE mod_name
USE mod_tempsalt

IMPLICIT none

INTEGER i,j,k,kk,mois(12)
REAL*8 rlatt,rlatu!,radian,radius
REAL :: a



data mois/31,28,31,30,31,30,31,31,30,31,30,31/
!__________________________________________________________________________________________
! Earth constants
!radius = 6371229.d0 ! earth radius in metre
!radian = pi/180.d0
!deg=radius*radian ! ~ 111000 metre
!tday=24.d0 * 3600.d0

! month lengths including leap years
do i=1900,3000
do k=1,12
idmax(k,i)=mois(k)
enddo
if((mod(i,4).eq.0 .and. mod(i,100).ne.0) .or.  mod(i,400).eq.0 ) idmax(2,i)=29
!#if defined orc
! do k=1,12  ! Julian days
!  idmax(k,i)=idmax(k-1,i)+idmax(k,i)
! enddo
!#endif
!print *,i,(idmax(k,i),k=1,12)
enddo

! x,y gridsize in degrees

rmin=19.d0
rmax=28.5d0

#if defined rco || baltix
dx = 1./15.d0
dy = 1./30.d0
!stlon1= 9.0d0+0.25d0*dx
!stlat1=53.8d0+1.25d0*dy
!print *,'coord stlon,stlat=',stlon1,stlat1
rmin=-2.d0
rmax=10.d0
#endif



dr=(rmax-rmin)/dble(MR-1)

tmin=-3. ; tmax=33
smin=33. ; smax=38.
dtemp=(tmax-tmin)/dble(MR-1)
dsalt=(smax-smin)/dble(MR-1)

! cosines relating to corners of grid box: ...u
! cosines relating to center of grid box:   ...t
do j=1,JMT
 rlatu=stlat1+dy*(j-1)+dy
 rlatt=rlatu-0.5d0*dy
 csu(j)=cos(rlatu*radian)
 cst(j)=cos(rlatt*radian)
enddo




!________________________________________________________________________________________
return
end subroutine coordinat

