!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
subroutine diffusion(x1,y1,z1,ib,jb,kb,dt,snew,st0,st1)
USE mod_param

IMPLICIT none
!#include "param.h"
!#include "part.h"
!=========================================================================================
! subroutine which adds a random position to the new trajectory position
! written by Richard Levine and David Webb
!
! x1,y1,z1 : original non-dimensional position of particle 
!            (fractions of a grid box side in the corresponding direction),
!            and is updated in subroutine to include random element      
! ib,jb,kb : original position in integers
! dt       : crossing time of particle through grid-box (in units of seconds)
! snew     : position of particle within time interval [st0,st1] (in units of seconds/m^3)
! [st0,st1]: is current time interval for interpolation (in units of seconds/m^3)
!
! param.h  : contains grid property definitions
! part.h   : contains velocity, property & status array definitions for particles
! ========================================================================================

REAL*8 xd,yd,zd,xe,ye,ze,x2,y2,z2,um,um0,um1,snew,st0,st1
REAL*8 ah,az
common /diffc/   ah,az   
integer irandom
common /randomc/ irandom
               
integer iloop,ic,jc,kc
logical ltest
      
! horizontal & vertical diff coefficients in m^2/s, which may be set as the same as 
! those of the OGCM 
ah = 2.0d2        
az = 1.0d-4      

!  check particle location in grid
if( kb.le.KM .and. ib.gt.1 .and. ib.le.IMT .and. jb.gt.1 .and. jb.le.JMT )then 
 ltest = .true.
 iloop = 0
 do while(ltest)
  iloop = iloop+1
	  
! calculate diffusive component in xd,yd,zd in metres for time-step dt	  
  call randomstep(dt,xd,yd,zd)
  xe = xd/(dx*cst(jb)*deg)
  ye = yd/(dy*deg)
  ze = zd/dz(kb)
	  
! update new position in x2,y2,z2 & ic,jc,kc for checks on particle position
  x2 = x1+xe
  y2 = y1+ye
  z2 = z1+ze
  ic = int(x2)+1
  jc = int(y2)+1
  kc = KM-int(z2)

!  check that column is deep enough and that at the particle
!  level there is at least one adjacent open ocean velocity point.
  
  if( ic.gt.1. and. ic.le.IMT .and. jc.gt.1.and.jc.le.JMT.and.kc.gt.0) then
   if (kc.le.kmt(ic,jc).and.(kc.le.kmu(ic,jc).or.kc.le.kmu(ic-1,jc).or. &
       kc.le.kmu(ic,jc-1).or.kc.le.kmu(ic-1,jc-1))) ltest=.false.
 endif

!  open boundary
  if(ic.eq.1.and.jc.ge.1.and.jc.le.jmt.and.kc.gt.0) then
   if (kc.le.kmt(ic,jc)) ltest = .false.
  endif

!  If iloop is very large, the particle is almost certainly stuck 
!  somewhere. Print error message and stop track
 
  if(iloop.gt.1000.and.ltest)then
   pstatus0(itrack) = 10
   ltest = .false.
!            write(99,*)'particle ',itrack,' stuck'           
  endif
 enddo

!  accept new position and update track
 x1 = x2
 y1 = y2
 z1 = z2
 ib = ic
 jb = jc
 kb = int(z1)+1

!  check that vertical velocities still OK if on boundary
 if(z1.eq.int(z1))then
  um0 = w0(ib,jb,kb-1)
  um1 = w1(ib,jb,kb-1)
  um  = um0 + (um1-um0)*(snew-st0)/(st1-st0)
  if(um.lt.0d0)then
   kb = kb-1
  endif
 endif

!  boundaries should be caught below - this is in case anything
!  is missed (but remember to reset below if track continues)

else
 pstatus0(itrack) = 8
 pstatus(itrack)  = 1
endif

end subroutine diffusion

!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
subroutine randomstep(dt,xd,yd,zd)
implicit none


!  subroutine to calculate the 'diffusion' steps to make in the x, y and z directions 
!  to represent diffusion effects over a timestep.
!
!  Input:
!
!    dt   -  The length of the timestep in seconds
!
!  Parameters:
!
!    ah   -  Horizontal diffusion coefficient (units of m^2/s)
!    az   -  vertical diffusion coefficient (m^2/s)
!
!  Output:
!
!    xd   - step in x direction (units m)
!    yd   - step in y direction (units m)
!    zd   - step in z direction (units m)
!
!  Uses:
!
!    ran1 - random number generator (produces a uniform distribution between 0.0 and 1.0)
!
!========================================================================================

real*8 dt,xd,yd,zd
real*8 a,radius,pi,theta
parameter (pi = 3.141592654d0)
real*8 ran1,qr

real*8  ah,az
integer irandom
common /randomc/ irandom
common /diffc/   ah,az

! horizontal
qr = dble(ran1(irandom))
a  = sqrt(4*dt*ah) 
radius = a*sqrt(-log(dble(1.-qr)))
theta = dble(2.0d0*pi*ran1(irandom))
xd = sin(theta)*radius
yd = cos(theta)*radius

! vertical
qr = dble(ran1(irandom))
a  = sqrt(4*dt*az) 
radius = a*sqrt(-log(dble(1.-qr)))
theta = dble(2.0d0*pi*ran1(irandom))
zd = sin(theta)*radius

end subroutine randomstep
      
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
     
function ran1(idum)
implicit none

!========================================================================================
!
! uniform random number generator
!
! Input/Output:
!
!   idum   - integer seed updated on each call to the routine
!
!  Output
!
!    ran1  - function returned value.  Real*8 random number in the range 0.0 to 1.0
!
!  Based on routine in Numerical Algorithms - itself based on Knuth
!  and K. Park and K.W. Miller, Comm. of the ACM 31 (1988) P.1192.
!
!========================================================================================

real*8  ran1
integer idum
integer IA,IM,IQ,IR,NTAB,NDIV
real*8  AM,EPS,RNMX

parameter(IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32)
parameter(AM=1./IM,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

integer j,k,iv(NTAB),iy
save iv,iy
data iv /NTAB*0/, iy /0/

if(idum.le.0.or.iy.eq.0)then
 idum=max(-idum,1)
 do 11 j=NTAB+8,1,-1
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if(idum.lt.0) idum=idum+IM
  if(j.le.NTAB) iv(j)=idum
 11 enddo
 iy=iv(1)
endif

k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if(idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum

ran1=min(AM*iy,RNMX)

end function ran1
