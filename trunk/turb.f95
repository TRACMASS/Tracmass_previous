!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
#ifdef turb    

subroutine turbu(ia,ja,ka,rr)

! computes the paramterised turbulent velocities u' and v' into upr

USE mod_param
USE mod_vel
USE mod_turb
IMPLICIT none

real*8 uv(6),rr,rg,en
integer ia,ja,ka,im,jm

call random_number(rand)
rand=2.*rand-1. ! Max. amplitude of turb. (varies with the same aplitude as the mean vel)
!rand=1.*rand-0.5

rg=1.d0-rr

im=ia-1
if(im.eq.0) im=IMT
jm=ja-1

! time interpolated velocities
uv(1)=(rg*u(ia,ja,ka,NST)+rr*u(ia,ja,ka,1))*ff ! western u
uv(2)=(rg*u(im,ja,ka,NST)+rr*u(im,ja,ka,1))*ff ! eastern u
uv(3)=(rg*v(ia,ja,ka,NST)+rr*v(ia,ja,ka,1))*ff ! northern v
uv(4)=(rg*v(ia,jm,ka,NST)+rr*v(ia,jm,ka,1))*ff ! southern v

!upr=uv*rand

upr(:,1)=upr(:,2) ! store u' from previous time iteration step
upr(:,2)=uv(:)*rand(:)

 
return
end subroutine turbu
#endif
!_______________________________________________________________________
