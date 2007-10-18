!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine turbu(ia,ja,ka,rr)

! computes the paramterised turbulent velocities  

USE mod_param
USE mod_vel
USE mod_turb
IMPLICIT none

real*8 uv(6),rr,rg,en
integer ia,ja,ka,im,jm

call random_number(rand)
rand=2.*rand-1.
!rand=1.*rand-0.5

rg=1.d0-rr

im=ia-1
if(im.eq.0) im=IMT
jm=ja-1

uv(1)=(rg*u(ia,ja,ka,NST)+rr*u(ia,ja,ka,1))*ff
uv(2)=(rg*u(im,ja,ka,NST)+rr*u(im,ja,ka,1))*ff
uv(3)=(rg*v(ia,ja,ka,NST)+rr*v(ia,ja,ka,1))*ff
uv(4)=(rg*v(ia,jm,ka,NST)+rr*v(ia,jm,ka,1))*ff

!en=0.25*sqrt(uv(1)**2+uv(2)**2+uv(3)**2+uv(4)**2)

upr=uv*rand
!print *,'random number=',rand
!print *,'uv',uv
!print *,'en',en
!print *,'upr',upr

 
return
end subroutine turbu

!_______________________________________________________________________
