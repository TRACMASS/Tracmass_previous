#ifdef timeanalyt 

subroutine cross_time(ijk,ia,ja,ka,r0,sp,sn,ts,tt,dsmin,dxyz,rr)

  ! subroutine to compute time (sp,sn) when trajectory 
  ! crosses face of box (ia,ja,ka) with the time analytical sceme by Vries and Döös (2001)
  ! two crossings are considered for each direction:  
  ! east and west for longitudinal directions, etc.  
  !
  !  Input:
  !
  !  ijk      : Considered direction (i=zonal, 2=meridional, 3=vertical)
  !  ia,ja,ka : Original position in integers
  !  r0       : Original non-dimensional position in the 
  !             ijk-direction of particle 
  !               (fractions of a grid box side in the 
  !                corresponding direction)
  !  ts       : Model dataset time step of trajectory
  !  tt       : Time in seconds of the trajectory relative to the code start 
  !  dsmin    : Timestep in sec divided by Volume of the grid-cell
  !  dxyz     : Volume of current grid-cell
  !  rr       : Time interpolation constant between 0 and 1 
  !
  !
  !  Output:
  !
  !    sp,sn   : crossing time to reach the grid box wall 
  !              (in units of s/m3)
  !
  ! computation of time (ss) at crossing for a>0, a<0 and a=0
  ! using knowledge of possible sign changes of velocities inside gridbox,
  ! the direction (rijk) (east/west or north/south or up/down) is determined
  ! except for boxes containing land points, solutions are traced until 
  ! the end of a season, so sometimes trajectories do not end on an edge! 
  ! depending on velocities there may be no solution in both directions 
  ! ss0(+dsmin): start(end) of season; s=s0: xi=xi0; s=ss: xi=xin
  ! numerical computation of root:
  ! let f(xi) = 0 be the root to be computed giving the time ss
  ! given some value x.neq.xi, it follows that
  ! 0 = f(xi) = (Taylor expansion) f(x) + (xi-x)*f'(x)
  ! with f' being the derivative of f (---> velocity)
  ! this gives the iteration method: xi_(n+1) = xi_n - f(xi_n)/f'(xi_n)
  ! choose initial value xi00 and use direction rijk for boundary condition 
  ! iteration variable xi_n is optimized during calculations for situations
  ! that would degrade convergency
  ! first, the accuracy parameter xxlim is used for minimizing f(x); then,
  ! xxlim is used to minimize the numerical error in the time variable ss 

USE mod_param
USE mod_vel
USE mod_grid
USE mod_turb
IMPLICIT NONE

INTEGER :: iim,loop,iil,ii,ijk,ia,ja,ka
REAL*8  :: uu,um,vv,vm,ss,aa,r0,sp,sn,dsmin,ts,tt,dxyz
REAL*8  :: f0,f1,dzs,dzu1,dzu2,rijk,s0,ss0,rr,rg
REAL*8,  PARAMETER ::  EPS=1.d-10
!_______________________________________________________________________________

sp=UNDEF ; sn=UNDEF

#ifdef twodim  
if(ijk.eq.3) return
#endif 

s0=tt/dxyz
ss0=dble(idint(ts))*tseas/dxyz
rg=1.d0-rr

loop=0 ; rijk=0.d0 ; ss=UNDEF ; f0=0.d0 ; f1=0.d0

if(ijk.eq.1) then
 ii=ia
 iim=ia-1
 iil=iim
 if(iim.eq.0) iil = imt
 uu=uflux(ii ,ja,ka,nsp)
 um=uflux(iil,ja,ka,nsp)
 vv=uflux(ii ,ja,ka,nsm)
 vm=uflux(iil,ja,ka,nsm)
#ifdef turb   
 if(r0.ne.dble(ja)) then
  uu=uu+upr(7,2)  
  vv=vv+upr(1,2)  
 else  ! add u' from previous iterative time step if on box wall
  uu=uu+upr(7,1) 
  vv=vv+upr(1,1)  
 endif
 if(r0.ne.dble(ja-1)) then
  um=um+upr(8,2)
  vm=vm+upr(2,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr(8,1)  
  vm=vm+upr(2,1)  
 endif
#endif
elseif(ijk.eq.2) then
 ii=ja
 iim=ja-1
 iil=iim
 uu=vflux(ia,ii ,ka,nsp)
 um=vflux(ia,iil,ka,nsp)
 vv=vflux(ia,ii ,ka,nsm)
 vm=vflux(ia,iil,ka,nsm)
#ifdef turb   
 if(r0.ne.dble(ii)) then
  vv=vv+upr(3,2)  
  uu=uu+upr(9,2)  
 else  ! add u' from previous iterative time step if on box wall
  vv=vv+upr(3,1)  
  uu=uu+upr(9,1) 
 endif
 if(r0.ne.dble(iim)) then
  um=um+upr( 4,2)
  vm=vm+upr(10,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr( 4,1)  
  vm=vm+upr(10,1)  
 endif
#endif
elseif(ijk.eq.3) then
 ii=ka
 iim=ka-1
 iil=iim
 uu=wflux(ii ,nsp)
 um=wflux(iil,nsp)
 vv=wflux(ii ,nsm)
 vm=wflux(iil,nsm)
#ifdef turb   
 if(r0.ne.dble(ii)) then
  vv=vv+upr( 5,2)  
  uu=uu+upr(11,2)  
 else  ! add u' from previous iterative time step if on box wall
  vv=vv+upr( 5,1)  
  uu=uu+upr(11,1) 
 endif
 if(r0.ne.dble(iim)) then
  um=um+upr( 6,2)
  vm=vm+upr(12,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr( 6,1)  
  vm=vm+upr(12,1)  
 endif
#endif
 if (ka.eq.km) then
  dzs= dz(km)+rg*hs(ia,ja,nsp)+rr*hs(ia,ja,nsm)
  dzu1=dz(km)+hs(ia,ja,nsp)
  dzu2=dz(km)+hs(ia,ja,nsm)
  if(dabs(dzu1).le.eps) stop 4956
  if(dabs(dzu2).le.eps) stop 4957
  f0=dzs/dzu1
  f1=dzs/dzu2
  uu=uu*f0
  um=um*f0
  vv=vv*f1
  vm=vm*f1
  if(dabs(f0).le.eps) stop 4958
  if(dabs(f1).le.eps) stop 4959
  f0=1.d0/f0
  f1=1.d0/f1
 endif 
endif

aa = -(vv-uu-vm+um)
! ----- a > 0
!     if (aa.gt.0.d0) then
     if (aa.gt.EPS) then
      call apos (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)
! ----- a < 0      
     elseif (aa.lt.-EPS) then
!     elseif (aa.lt.0.d0) then
      call amin (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)
! ----- a = 0 
!     elseif (aa.eq.0.d0) then
      else
      if (aa.ne.0.d0) print *,'anil eps=',aa
      call anil (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)
     endif
! translate direction and positions for old to new tracmass
      if (rijk.eq.dble(ii)) then
       sp=ss-s0
       sn=UNDEF
      elseif(rijk.eq.dble(iim)) then
       sn=ss-s0
       sp=UNDEF
      else
       sp=UNDEF
       sn=UNDEF
      endif
      if(sp.eq.0.d0) sp=UNDEF
      if(sn.eq.0.d0) sn=UNDEF
return
end subroutine cross_time
!_______________________________________________________________________

subroutine pos_time(ijk,ia,ja,ka,r0,r1,ts,tt,dsmin,dxyz,ss0,ds,rr)
  ! subroutine to compute  the new position (r0) of the trajectory 
  !
  !  Input:
  !
  !  ijk      : Considered direction (i=zonal, 2=meridional, 3=vertical)
  !  ia,ja,ka : Original position in integers
  !  r0       : Original non-dimensional position in the 
  !             ijk-direction of particle 
  !               (fractions of a grid box side in the 
  !                corresponding direction)
  !  ts       : Model dataset time step of trajectory
  !  tt       : Time in seconds of the trajectory relative to the code start 
  !  dsmin    : Timestep in sec divided by Volume of the grid-cell
  !  dxyz     : Volume of current grid-cell
  !  ss0      : dble(idint(ts))*tseas/dxyz
  !  ds       : Min time to face crossing calc in cross
  !  rr       : Time interpolation constant between 0 and 1 
  !
  !
  !  Output:
  !
  !  r1       : New non-dimensional position in the 
  !             ijk-direction of particle 
  !               (fractions of a grid box side in the 
  !                corresponding direction)
  !

USE mod_param
USE mod_vel
USE mod_grid
USE mod_turb
IMPLICIT NONE

REAL*8,  PARAMETER ::  xilim=3.0d0,xxlim=1.d-7
INTEGER :: ijk,ia,ja,ka,iim,iil,ii
REAL*8  :: uu,um,vv,vm,xi,xi0,const,ga,erf0,aa,bb,daw0,r0,r1,dsmin,dxyz,tt,ts,rijk
REAL*8  :: s15aff,dawson,s15adf,s15aef,errfun
REAL*8  :: f0,f1,dzs,dzu1,dzu2,s0,ss,ss0,ds,rr,rg
REAL*8,  PARAMETER ::  EPS=1.d-10

#ifdef twodim  
if(ijk.eq.3) then
r1=r0
return
endif
#endif 


s0=tt/dxyz-ds
ss=ts*tseas/dxyz
!rr=1.d0-rg
rg=1.d0-rr
f0=0.d0 ; f1=0.d0


if(ijk.eq.1) then
 ii=ia
 iim=ia-1
 iil=iim
 if(iim.eq.0) iil = imt
 uu=uflux(ii ,ja,ka,nsp)
 um=uflux(iil,ja,ka,nsp)
 vv=uflux(ii ,ja,ka,nsm)
 vm=uflux(iil,ja,ka,nsm)
#ifdef turb   
 if(r0.ne.dble(ja)) then
  uu=uu+upr(7,2)  
  vv=vv+upr(1,2)  
 else  ! add u' from previous iterative time step if on box wall
  uu=uu+upr(7,1) 
  vv=vv+upr(1,1)  
 endif
 if(r0.ne.dble(ja-1)) then
  um=um+upr(8,2)
  vm=vm+upr(2,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr(8,1)  
  vm=vm+upr(2,1)  
 endif
#endif
elseif(ijk.eq.2) then
 ii=ja
 iim=ja-1
 iil=iim
 uu=vflux(ia,ii ,ka,nsp)
 um=vflux(ia,iil,ka,nsp)
 vv=vflux(ia,ii ,ka,nsm)
 vm=vflux(ia,iil,ka,nsm)
#ifdef turb   
 if(r0.ne.dble(ii)) then
  vv=vv+upr(3,2)  
  uu=uu+upr(9,2)  
 else  ! add u' from previous iterative time step if on box wall
  vv=vv+upr(3,1)  
  uu=uu+upr(9,1) 
 endif
 if(r0.ne.dble(iim)) then
  um=um+upr( 4,2)
  vm=vm+upr(10,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr( 4,1)  
  vm=vm+upr(10,1)  
 endif
#endif
elseif(ijk.eq.3) then
 ii=ka
 iim=ka-1
 iil=iim
 uu=wflux(ii ,nsp)
 um=wflux(iil,nsp)
 vv=wflux(ii ,nsm)
 vm=wflux(iil,nsm)
#ifdef turb   
 if(r0.ne.dble(ii)) then
  vv=vv+upr( 5,2)  
  uu=uu+upr(11,2)  
 else  ! add u' from previous iterative time step if on box wall
  vv=vv+upr( 5,1)  
  uu=uu+upr(11,1) 
 endif
 if(r0.ne.dble(iim)) then
  um=um+upr( 6,2)
  vm=vm+upr(12,2)
 else  ! add u' from previous iterative time step if on box wall
  um=um+upr( 6,1)  
  vm=vm+upr(12,1)  
 endif
#endif
 if (ka.eq.km) then
  dzs= dz(km)+rg*hs(ia,ja,nsp)+rr*hs(ia,ja,nsm)
  dzu1=dz(km)+hs(ia,ja,nsp)
  dzu2=dz(km)+hs(ia,ja,nsm  )
  if(dabs(dzu1).le.eps) stop 4966
  if(dabs(dzu2).le.eps) stop 4967
  f0=dzs/dzu1
  f1=dzs/dzu2
  uu=uu*f0
  um=um*f0
  vv=vv*f1
  vm=vm*f1
  if(dabs(f0).le.eps) stop 4968
  if(dabs(f1).le.eps) stop 4969
  f0=1.d0/f0
  f1=1.d0/f1
 endif 
endif

if(dabs(dsmin).le.eps) stop 4981
aa = -(vv-uu-vm+um)/dsmin
bb = um-uu


if(aa.ne.0.d0) then
if(f0.eq.0.d0) then
 ga = -dble(iim) + (vm-um)/(vv-uu-vm+um)
 const = (um*vv-uu*vm)/(vv-uu-vm+um)
else
 ga = -dble(iim) + (f1*vm-f0*um)/(vv-uu-vm+um)
 const = (f0*um*(vv-vm)+f1*vm*(um-uu))/(vv-uu-vm+um)
endif
endif

if (aa.gt.0.d0) then
 const=2.d0*const/dsqrt(2.d0*aa)
 xi0=(bb+aa*(s0-ss0))/dsqrt(2.d0*aa)
 xi =(bb+aa*(ss-ss0))/dsqrt(2.d0*aa)
 daw0 = s15aff(xi0) !      daw0 = s15aff(xi0,ifail)
 r1 = dawson (const,daw0,r0+ga,xi0,xi) -ga
elseif (aa.lt.0.d0) then
 const=const*sqrt(pi/(-2.d0*aa))
 xi0=(bb+aa*(s0-ss0))/dsqrt(-2.d0*aa)
 xi =(bb+aa*(ss-ss0))/dsqrt(-2.d0*aa)
 if (xi0.gt.xilim) then ! complementary error function
  erf0 = s15adf(xi0) !       erf0 = s15adf(xi0, ifail)
 elseif(xi0.lt.-xilim) then
  erf0 = - s15adf(-xi0) !       erf0 = - s15adf(-xi0, ifail)
 else ! error function
  erf0 = - s15aef(xi0)     !       erf0 = - s15aef(xi0, ifail)
 endif
 r1 = errfun (const,erf0,r0+ga,xi0,xi) -ga 
elseif (aa.eq.0.d0) then
 if (f0.eq.0.d0) then
  f0=1.d0
  f1=1.d0
 endif 
 aa = um-uu
 bb =(f0*um-f1*vm)/dsmin
 xi = ss-s0
 if (aa.eq.0.d0) then
  r1 = r0 - xi*(-f0*um+0.5d0*bb*(s0-ss0+ss-ss0))
 else
  ga = (r0-dble(iim)) + (bb*(s0-ss0-1.d0/aa)-f0*um)/aa
  r1 = r0 + ga*(dexp(-aa*xi) - 1.d0) -bb*xi/aa
 endif 
endif

if (r1-dble(ii) .gt.0.d0 .and. r1-dble(ii) .lt.xxlim) r1=dble(ii )-xxlim
if (dble(iim)-r1.gt.0.d0 .and. dble(iim)-r1.lt.xxlim) r1=dble(iim)+xxlim


return
end subroutine pos_time
!_______________________________________________________________________

subroutine apos (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)
! computation of time (ss) at crossing for a>0 

USE mod_param
USE mod_vel
IMPLICIT NONE

REAL*8,  PARAMETER ::  EPS=1.d-10
REAL*8,  PARAMETER ::  xxlim=1.d-7
INTEGER :: ii,iim,iconfig,i,loop
REAL*8  :: uu,um,vv,vm,xi0,xin,ssii,ssiim,xerr,xf1,xi,xi00,xia,xib,xibf
REAL*8  :: aa,bb,ga,const,daw0,s15aff,xf,xf2,xiaf,dawson
REAL*8  :: rijk,r0,ss,f0,f1,s0,ss0,dsmin    

aa = -(vv-uu-vm+um)/dsmin
bb = um-uu
xi0=(bb+aa*(s0-ss0))/dsqrt(2.d0*aa)

! 1) ii or iim ---> land point 
if (uu.eq.0.d0 .and. vv.eq.0.d0) then
 xin =-dlog(ii-r0) 
 if (xi0.lt.0.d0 .and. xi0*xi0.ge.xin) then
  rijk = dble(iim)
  xin=-dsqrt(xi0*xi0-xin)
  ss= 2.d0*(xin-xi0)/dsqrt(2.d0*aa) + s0
  return
 else
  return ! no solution
 endif
elseif (um.eq.0.d0 .and. vm.eq.0.d0 ) then
 xin =-dlog(r0-iim)
 if (xi0.lt.0.d0 .and. xi0*xi0.ge.xin) then
  rijk = dble(ii)
  xin=-dsqrt(xi0*xi0-xin)
  ss= 2.d0*(xin-xi0)/dsqrt(2.d0*aa) + s0
  return
 else
! no solution
  return
 endif
endif

if (f0.eq.0.d0) then
 ga = -dble(iim) + (vm-um)/(vv-uu-vm+um)
 const = (um*vv-uu*vm)/(vv-uu-vm+um)
 ssii = uu/(uu-vv)
 ssiim= um/(um-vm)
else
 ga = -dble(iim) + (f1*vm-f0*um)/(vv-uu-vm+um)
 const = (f0*um*(vv-vm)+f1*vm*(um-uu))/(vv-uu-vm+um)
 ssii = (uu+um*(f0-1.d0))/(uu-vv+um*(f0-1.d0)+vm*(1.d0-f1))
 ssiim= f0*um/(f0*um-f1*vm)
endif
const = 2.d0*const/dsqrt(2.d0*aa)
daw0 = s15aff(xi0)     !  daw0 = s15aff(xi0, ifail)
! 'velocity' at xi0
xf1 = const -(xi0+xi0)*(r0+ga)    

! 2) ii direction (four velocity configurations at edge)
if (ssii.le.0.d0 .or. ssii.ge.1.d0) then
 if (uu.gt.0.d0 .or. vv.gt.0.d0) then
! 2a) configuration: +++
  iconfig = 1
  xi00=(vm-vv)/dsqrt(2.d0*aa)
 else
! 2b) configuration: --- no solution, try rijk=iim
!print *,'configuration: --- no solution, try rijk=iim',ssii,uu,vv
  rijk=dble(iim)
  goto 3000
 endif 
else
 if (vv.gt.0.d0) then
! 2c) configuration: -- 0 ++
  iconfig = 3   
  xi00=(vm-vv)/dsqrt(2.d0*aa)
 else
! 2d) configuration: ++ 0 --
  if (s0-ss0.ge.ssii*dsmin) goto 3000
  iconfig = 4
  xi00=(bb+aa*ssii*dsmin)/dsqrt(2.d0*aa) 
 endif 
endif

xf=dawson(const,daw0,r0+ga,xi0,xi00)-ga-dble(ii)

if (xf.lt.0.d0) goto 3000
rijk=dble(ii)
xib =xi00
xibf=xf
 if (xf1.gt.0.d0) then
  if (xf1*(xi00-xi0).gt.-(r0-rijk)) then
   xi00=xi0
   xf=r0-rijk
  endif
 else
  if (ssiim.gt.0.d0 .and. ssiim.lt.1.d0 .and.vm.gt.0.d0) then
!  check crossing at 3c) configuration: -- 0 ++
   xf2=(bb+aa*ssiim*dsmin)/dsqrt(2.d0*aa)
   xf2=dawson(const,daw0,r0+ga,xi0,xf2)-ga-dble(iim)
   if (xf2.lt.0.d0) goto 3000
  endif                
  if (iconfig.eq.4) then
   xi00=(bb+aa*ssiim*dsmin)/dsqrt(2.d0*aa)
   xf=dawson(const,daw0,r0+ga,xi0,xi00)-ga-dble(ii)
  endif
 endif     
 goto 1000

3000  continue


! 3) iim direction (four velocity configurations at edge)
     if (ssiim.le.0.d0 .or. ssiim.ge.1.d0) then
      if (um.gt.0.d0 .or. vm.gt.0.d0) then
! 3a) configuration: +++ no solution
       rijk=dble(iim)
       return
      else
! 3b) configuration: ---
       iconfig = 6
       xi00=(vm-vv)/dsqrt(2.d0*aa)
      endif 
     else
      if (vm.gt.0.d0) then
! 3c) configuration: -- 0 ++
       if (s0-ss0.ge.ssiim*dsmin) return
       iconfig = 7   
       xi00=(bb+aa*ssiim*dsmin)/dsqrt(2.d0*aa) 
      else
! 3d) configuration: ++ 0 --
       iconfig = 8
       xi00=(vm-vv)/dsqrt(2.d0*aa) 
      endif 
     endif
     xf=dawson(const,daw0,r0+ga,xi0,xi00)-ga-dble(iim)
     if (xf.gt.0.d0) return
     rijk=dble(iim)
     xib =xi00
     xibf=xf
     if (xf1.lt.0.d0) then
      if (xf1*(xi00-xi0).lt.-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     elseif (iconfig.eq.7) then
      xi00=(bb+aa*ssii*dsmin)/dsqrt(2.d0*aa)
      xf=dawson(const,daw0,r0+ga,xi0,xi00)-ga-dble(iim)
     endif      

1000  continue
! calculation of root
     if (loop.eq.100) then    ! why 100 ?
      xi = -0.05d0 + xi0
      do i=1,100              ! why 100 ?
       xi=xi+0.05d0
       xf2= dawson (const,daw0,r0+ga,xi0,xi)
       xf = xf2-ga-rijk
       xf1= const-(xi+xi)*xf2 
!       print *, ' loop', xi, xf, xf1
      enddo
      loop = 1000
     else

     xi = xi00
     xf1= const-(xi+xi)*(xf+rijk+ga)
     xia = xi0
     xiaf= r0-rijk 
     xerr = xxlim
     ssii = xxlim
100   continue
     loop=loop+1
     if(loop.eq.80) then
!      print *, ' loop stop', xi0, xi00, xi, xf, ntrac, loop
      if (dabs(xf).gt.100*xxlim) loop=1000
!      print *, 'loopnr', abs(xf),xxlim,100*xxlim, loop
!       if (xerr.eq.UNDEF) loop=1000
      goto 200
     endif
     if(dabs(xf1).le.eps) stop 4992
     xin = xi - xf/xf1
     if ( dabs(xf).gt.100.d0 .or.(xf*xibf.gt.0.d0 .and.(xf1*xibf.lt.0.d0 .or.xin.lt.xia)) ) then
      xib = xi
      xibf= xf
      xi = 0.5d0*(xi+xia)
     elseif (xf*xibf.lt.0.d0 .and. (xf1*xibf.lt.0.d0 .or.xin.gt.xib) ) then
      xia = xi
      xiaf= xf
      xi = 0.5d0*(xi+xib)
     else
 !     ic=3   What is this????????? Should it be iconfig=3 ?????????????????????
      xi = xin 
     endif

     xf2= dawson (const,daw0,r0+ga,xi0,xi)
     xf = xf2-ga-rijk
     xf1= const-(xi+xi)*xf2 
!     if(loop.eq.1000) print *, ' loop', xi, xf, xf1
     if (dabs(xf).ge.xerr) then
      goto 100
     elseif ( dabs(xf).gt.xxlim*dabs(xf1*(xi-xi0)).and. dabs(xf).lt.ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
      xerr = UNDEF
      ssii = dabs(xf)
      ssiim= xi
      goto 100
     endif
     if (xerr.eq.UNDEF) xi = ssiim

     endif
200   continue
     ss= 2.d0*(xi-xi0)/dsqrt(2.d0*aa) + s0
!     nstat(iconfig,1,1)=nstat(iconfig,1,1)+1       ! removed statistics
!     nstat(iconfig,2,1)=nstat(iconfig,2,1)+loop
     if (xi.ne.xi0) then
      xin=dabs( xf/((xi-xi0)*xf1) )
!      accu(1)=accu(1)+1.0
!      accu(2)=accu(2)+xin
!      if (xin.gt.xxlim.and.ntrac.le.100) print *,ntrac,xin,xf1,ss-s0
     endif
     if (loop.eq.1000.or.ss-s0.lt.0.d0) then
      if (ss-s0.lt.0.d0) print *,' ss-s0 is negative '
!      print *, '+ time cross =', ss, s0,loop
!      print *, aa,uu,um,vv,vm
!      print *, iconfig,r0,rijk
!       print *, xi0,xi00,xib,xf1,xibf,xi,xf
!      print *,ii,iim,dabs(xf),xxlim,dabs(xf1*(xi-xi0)),ssii,xerr
!      STOP
     endif

return
end subroutine apos
!_______________________________________________________________________

subroutine amin (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)

! computation of time (ss) at crossing for a<0 

USE mod_param
USE mod_vel
IMPLICIT NONE

REAL*8,  PARAMETER ::  EPS=1.d-10
REAL*8,  PARAMETER ::  xilim=3.0d0,xxlim=1.d-7
INTEGER :: ii,iim,iconfig,i,loop
REAL*8  :: uu,um,vv,vm,xib,xia,xi00,xi0,xf1,xf,xf2,xin,xibf,xi,xerr,xiaf
REAL*8  :: aa,bb,ga,s15aef,s15adf,erf0,const,errfun
REAL*8  :: rijk,r0,ss,f0,f1,s0,ss0,ssii,ssiim,dsmin   

aa = -(vv-uu-vm+um)/dsmin
bb = um-uu
xi0=(bb+aa*(s0-ss0))/dsqrt(-2.d0*aa)


! 1) ii or iim ---> land point 
     if (uu.eq.0.d0 .and. vv.eq.0.d0) then
      rijk = dble(iim)
      xin =-dlog(ii-r0) 
      xin =-dsqrt(xi0*xi0+xin)
      ss= 2.d0*(xi0-xin)/dsqrt(-2.d0*aa) + s0
      return
     elseif (um.eq.0.d0 .and. vm.eq.0.d0 ) then
      rijk = dble(ii)
      xin =-dlog(r0-iim)
      xin =-dsqrt(xi0*xi0+xin)
      ss= 2.d0*(xi0-xin)/dsqrt(-2.d0*aa) + s0
      return
     endif 

     if (f0.eq.0.d0) then
      ga = -dble(iim) + (vm-um)/(vv-uu-vm+um)
      const = (um*vv-uu*vm)/(vv-uu-vm+um)
     if(dabs(uu-vv).le.eps) stop 4993
      ssii = uu/(uu-vv)
     if(dabs(um-vm).le.eps) stop 4994
      ssiim= um/(um-vm)
     else
      ga = -dble(iim) + (f1*vm-f0*um)/(vv-uu-vm+um)
      const = (f0*um*(vv-vm)+f1*vm*(um-uu))/(vv-uu-vm+um)
     if(dabs(uu-vv+um*(f0-1.d0)+vm*(1.-f1)).le.eps) stop 4995
      ssii = (uu+um*(f0-1.d0))/(uu-vv+um*(f0-1.d0)+vm*(1.-f1))
     if(dabs(f0*um-f1*vm).le.eps) stop 4996
      ssiim= f0*um/(f0*um-f1*vm)
     endif
     const = const*sqrt(pi/(-2.d0*aa))
     if (xi0.gt.xilim) then
! complementary error function
!      erf0 = s15adf(xi0, ifail)
      erf0 = s15adf(xi0)
     elseif(xi0.lt.-xilim) then
!      erf0 = -s15adf(-xi0, ifail)
      erf0 = -s15adf(-xi0)
     else
! error function
!      erf0 = -s15aef(xi0, ifail)
!print *,'erf0',xi0
!stop 5906
      erf0 = -s15aef(xi0)
     endif
! minus-velocity at xi0
     xf1 = xi0*(r0+ga)-(const/dsqrt(pi)) 

! 2) ii direction (four velocity configurations at edge)
     if (ssii.le.0.d0 .or. ssii.ge.1.d0) then
      if (uu.gt.0.d0 .or. vv.gt.0.d0) then
! 2a) configuration: +++
       iconfig = 1
       xi00=(vm-vv)/dsqrt(-2.d0*aa)
      else
! 2b) configuration: --- no solution, try rijk=iim
       rijk=dble(iim)
       goto 4000
      endif 
     else
      if (vv.gt.0.d0) then
! 2c) configuration: -- 0 ++
       iconfig = 3   
       xi00=(vm-vv)/dsqrt(-2.d0*aa)
      else
! 2d) configuration: ++ 0 --
       if (s0-ss0.ge.ssii*dsmin) goto 4000
       iconfig = 4
       xi00=(bb+aa*ssii*dsmin)/dsqrt(-2.d0*aa) 
      endif 
     endif
     xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(ii)
     if (xf.lt.0.d0) goto 4000
     rijk=dble(ii)
     xib =xi00
     xibf=xf
     if (xf1.lt.0.d0) then
      if (xf1*(xi00-xi0).gt.-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     else
      if (ssiim.gt.0.d0 .and.ssiim.lt.1.d0 .and.vm.gt.0.d0) then
! check crossing at 3c) configuration: -- 0 ++
       xf2=(bb+aa*ssiim*dsmin)/dsqrt(-2.d0*aa)
       xf2=errfun(const,erf0,r0+ga,xi0,xf2)-ga-dble(iim)
       if (xf2.lt.0.d0) goto 4000
      endif
      if (iconfig.eq.4) then
       xi00=(bb+aa*ssiim*dsmin)/dsqrt(-2.d0*aa)
       xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(ii)
      endif
     endif 
     goto 2000

4000  continue
! 3) iim direction (four velocity configurations at edge)
     if (ssiim.le.0.d0 .or. ssiim.ge.1.d0) then
      if (um.gt.0.d0 .or. vm.gt.0.d0) then
! 3a) configuration: +++ no solution
!     print *,'3a no solution'
       return
      else
! 3b) configuration: ---
       iconfig = 6
       xi00=(vm-vv)/dsqrt(-2.d0*aa)
      endif 
     else
      if (vm.gt.0.d0) then
! 3c) configuration: -- 0 ++
       if (s0-ss0.ge.ssiim*dsmin) return
       iconfig = 7   
       xi00=(bb+aa*ssiim*dsmin)/dsqrt(-2.d0*aa) 
      else
! 3d) configuration: ++ 0 --
       iconfig = 8
       xi00=(vm-vv)/dsqrt(-2.d0*aa) 
      endif 
     endif
     xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(iim)
     if (xf.gt.0.d0) return
     rijk=dble(iim)
     xib =xi00
     xibf=xf
     if (xf1.gt.0.d0) then
      if (xf1*(xi00-xi0).lt.-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     elseif (iconfig.eq.7) then   
      xi00=(bb+aa*ssii*dsmin)/dsqrt(-2.d0*aa)
      xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(iim)
     endif

2000  continue
! calculation of root
     if (loop.eq.100) then
      xi = 0.005d0 + xi0
      do i=1,200
       xi=xi-0.005d0
       xf2= errfun (const,erf0,r0+ga,xi0,xi)
       xf = xf2-ga-rijk
       xf1= xi*xf2-(const/dsqrt(pi)) 
       print *, ' loop', xi, xf, xf1
      enddo
      loop = 1000
     else

     xi = xi00
     xf1= xi*(xf+rijk+ga)-(const/dsqrt(pi))
     xia = xi0
     xiaf= r0-rijk
     xerr = xxlim
     ssii = xxlim      
200   continue
     loop=loop+1
     if(loop.eq.80) then
!      print *, ' loop stop', xi0, xi00, xi, xf
      if (dabs(xf).gt.100*xxlim) loop=1000
!       if (xerr.eq.UNDEF) loop=1000
      goto 300
     endif
     xf1 = xf1+xf1
     if(xf1.eq.0.) stop 59785
     xin = xi - xf/xf1
     if ( dabs(xf).gt.100.d0 .or.(xf*xibf.gt.0.d0 .and. (xf1*xibf.gt.0.d0 .or.xin.gt.xia)) )then
      xib = xi
      xibf= xf
      xi = 0.5d0*(xi+xia)
     elseif (xf*xibf.lt.0.d0 .and.  (xf1*xibf.gt.0.d0 .or.xin.lt.xib) ) then
      xia = xi
      xiaf= xf
      xi = 0.5d0*(xi+xib)
     else
      xi = xin      
     endif
!     stop 49567
     xf2= errfun (const,erf0,r0+ga,xi0,xi)
     xf1= xi*xf2-(const/dsqrt(pi))
     xf = xf2-ga-rijk
!      print *, ' loop', xi, xf, xf1
     if (dabs(xf).ge.xerr) then
      goto 200
     elseif ( dabs(xf).gt.xxlim*dabs(xf1*(xi-xi0)).and. dabs(xf).lt.ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
      xerr = UNDEF
      ssii = dabs(xf)
      ssiim= xi
      goto 200
     endif
     if (xerr.eq.UNDEF) xi = ssiim

     endif
300   continue
     ss= 2.d0*(xi0-xi)/dsqrt(-2.d0*aa) + s0 
!     nstat(iconfig,1,2)=nstat(iconfig,1,2)+1
!     nstat(iconfig,2,2)=nstat(iconfig,2,2)+loop
     if (xi.ne.xi0) then
      xin=dabs( xf/((xi-xi0)*xf1) )
!      accu(1)=accu(1)+1.0
!      accu(2)=accu(2)+xin
!      if (xin.gt.xxlim.and.ntrac.le.100) print *,ntrac,xin,xf1,ss-s0
     endif
     if (loop.eq.1000.or.ss-s0.lt.0.d0) then
      if (ss-s0.lt.0.d0) print *, ' ss-s0 is negative '
!      print *, '- time cross =', ss, s0, loop
!      print *, aa,uu,um,vv,vm
!      print *, iconfig,r0,rijk
!       print *, xi0,xi00,xib,xf1,xibf,xi,xf
!      print *, bb,ga,const,erf0
!      STOP
     endif

return
end subroutine amin
!_______________________________________________________________________

subroutine anil (ii,iim,r0,rijk,s0,ss,ss0,uu,um,vv,vm,f0,f1,loop,dsmin)

! computation of time (ss) at crossing for aa=0 

USE mod_param
USE mod_vel
IMPLICIT NONE

REAL*8,  PARAMETER ::  EPS=1.d-10
REAL*8,  PARAMETER :: xxlim=1.d-7
INTEGER :: ii,iim
REAL*8  :: aa,bb,ga,dsmin,uu,um,vv,vm,ss,ssii,ssiim
REAL*8  :: xerr,xf,xf1,xf2,xi,xi0,xi00,xia,xiaf,xib,xibf,xin
REAL*8  :: r0,rijk,s0,ss0,f0,f1
INTEGER :: loop,i,iconfig

if (f0.eq.0.d0) then
 f0=1.d0
 f1=1.d0
endif

aa = um-uu
bb =(f0*um-f1*vm)/dsmin 

if (aa.eq.0.d0) then
   if (bb.eq.0.d0) then
     ga = f0*um
     if (ga.gt.0.d0) then
       rijk=dble(ii)
       ss=s0+(rijk-r0)/ga
     elseif (ga.lt.0.) then
       rijk=dble(iim)
       ss=s0+(rijk-r0)/ga
     endif
     return
   else
     if(dabs(bb).le.eps) stop 4537
     ga = -(f0*um+aa*dble(iim))/bb -ss0
     if (s0+ga.ge.0.d0) then
       if (bb.gt.0.d0) then
         rijk=dble(iim)
       else
         rijk=dble(ii)
       endif
       ss=-ga+dsqrt( (s0+ga)**2-2.d0*(rijk-r0)/bb )
     else
       if (bb.gt.0.d0) then
         rijk=dble(ii)
       else
         rijk=dble(iim)
       endif         
       if ( (s0+ga)**2 .ge. 2.d0*(rijk-r0)/bb ) then
         ss=-ga-dsqrt( (s0+ga)**2-2.d0*(rijk-r0)/bb )
       else
         if (bb.gt.0.d0) then
           rijk=dble(iim)
         else
           rijk=dble(ii)
         endif
         ss=-ga+dsqrt( (s0+ga)**2-2.d0*(rijk-r0)/bb )
       endif
     endif         
   endif
   return
else
  if (bb.eq.0.d0) then  ! b=0: stationary velocity fields !!!
    if (uu.gt.0.d0) then
      xi=1.d0+(r0-dble(ii))*(uu-um)/(uu-um*(f0-1.d0))        
      if (xi.gt.0.) then
        rijk=dble(ii)
        if (um.ne.uu) then
          ss= s0 - 1.d0/(uu-um)*dlog(xi)
        else
          ss= s0 + (dble(ii)-r0)/(uu-um*(f0-1.)) 
        endif
        if (ss.le.0.d0) ss=UNDEF
        return
      endif
    endif
    if (um.lt.0.d0) then
      xi=1.d0+(r0-dble(iim))*(uu-um)/(f0*um)        
      if (xi.gt.0.d0) then
        rijk=dble(iim)
        if (um.ne.uu) then
          ss= s0 - 1.d0/(uu-um)*dlog(xi)
        else
          ss= s0 + (dble(iim)-r0)/(f0*um)
        endif
        if (ss.le.0.d0) ss=UNDEF
      endif
    endif
!   ssii = dsmin+dsmin
!   ssiim= dsmin+dsmin
        return
       else
        ssii = (f0*um-aa)/bb
        ssiim= f0*um/bb
       endif
       ga = (r0-dble(iim)) + (bb*(s0-ss0-1.d0/aa)-f0*um )/aa
       xi0 = 0.d0
! ''velocity'' at xi0
       xf1 = f0*um-bb*(s0-ss0)-aa*(r0-dble(iim)) 

! 2) ii direction (four velocity configurations at edge)
       if (ssii.le.0.d0 .or. ssii.ge.dsmin) then
        if (uu.gt.0.d0 .or. vv.gt.0.d0) then
! 2a) configuration: +++
         iconfig = 1
         xi00=ss0+dsmin-s0
        else
! 2b) configuration: --- no solution, try rijk=iim
         rijk=dble(iim)
         goto 3000
        endif 
       else
        if (vv.gt.0.d0) then
! 2c) configuration: -- 0 ++
         iconfig = 3   
         xi00=ss0+dsmin-s0
        else
! 2d) configuration: ++ 0 --
         if (s0-ss0.ge.ssii) goto 3000
         iconfig = 4
         xi00=ss0+ssii-s0 
        endif 
       endif
       xf=(r0-dble(ii))+ga*(dexp(-aa*xi00)-1.d0)-bb*xi00/aa
       if (xf.lt.0.) goto 3000
       rijk=dble(ii)
       xib =xi00
       xibf=xf
       if (xf1.gt.0.d0) then
        if (xf1*(xi00-xi0).gt.-(r0-rijk)) then
         xi00=xi0
         xf=r0-rijk
!          if (b.ne.0.) goto 1000
        endif
       else
        if (ssiim.gt.0. .or. ssiim.lt.dsmin .and.vm.gt.0.d0) then
! check crossing at 3c) configuration: -- 0 ++
         xf2=ss0+ssiim-s0
         xf2=(r0-dble(iim))+ga*(dexp(-aa*xf2)-1.d0)-bb*xf2/aa
         if (xf2.lt.0.d0) goto 3000
        endif                
       endif   
       goto 1000

3000    continue
! 3) iim direction (four velocity configurations at edge)
       if (ssiim.le.0.d0 .or. ssiim.ge.dsmin) then
        if (um.gt.0.d0 .or. vm.gt.0.d0) then
! 3a) configuration: +++ no solution
         return
        else
! 3b) configuration: ---
         iconfig = 6
         xi00=ss0+dsmin-s0
        endif 
       else
        if (vm.gt.0.) then
! 3c) configuration: -- 0 ++
         if (s0-ss0.ge.ssiim) return
         iconfig = 7   
         xi00=ss0+ssiim-s0 
        else
! 3d) configuration: ++ 0 --
         iconfig = 8
         xi00=ss0+dsmin-s0 
        endif 
       endif
       xf=(r0-dble(iim))+ga*(dexp(-aa*xi00)-1.d0)-bb*xi00/aa
       if (xf.gt.0.d0) then
!        rijk=dble(iim)  ! nyinlagt
!        print *,'nytt prov'
!        stop 4968
        return
       endif
       rijk=dble(iim)
       xib =xi00
       xibf=xf
       if (xf1.lt.0.d0) then
        if (xf1*(xi00-xi0).lt.-(r0-rijk)) then
         xi00=xi0
         xf=r0-rijk
        endif
       endif

! determining root and crossing time ss        
1000    continue
       xi00 = dabs(aa)*xi00
       xib  = dabs(aa)*xib
       if (loop.eq.100) then
        xi=xi0-0.05d0
        do i=1,100
         xi = xi+0.05d0
         xf = (r0-rijk)+ga*(dexp(-sign(1.d0,aa)*xi)-1.d0)-bb*xi/(aa*dabs(aa)) 
         xf1= (-bb/aa-aa*ga*dexp(-sign(1.d0,aa)*xi))/dabs(aa)         
!         print *, ' loop', xi, xf, xf1
        enddo
        loop=1000
       else

       xi = xi00 
       xf1= (-bb/aa-aa*ga*dexp(-sign(1.d0,aa)*xi))/dabs(aa) 
       xia = xi0
       xiaf= r0-rijk
       xerr = xxlim
       ssii = xxlim
100     continue
!        print *, ' loop', xi, xf, xf1
       loop=loop+1
       if(loop.eq.90) then
!        print *, ' loop stop', xi0, xi00, xi, xf
        if (dabs(xf).gt.100.d0*xxlim) loop=1000
!         if (xerr.eq.UNDEF) loop=1000
        goto 400
       endif
       xin = xi - xf/xf1
     if ( dabs(xf).gt.100.d0 .or.(xf*xibf.gt.0.d0 .and. (xf1*xibf.lt.0.d0 .or.xin.lt.xia)) ) then
        xib = xi
        xibf= xf
        xi = 0.5d0*(xi+xia)
       elseif (xf*xibf.lt.0. .and. (xf1*xibf.lt.0.d0 .or.xin.gt.xib) ) then
        xia = xi
        xiaf= xf
        xi = 0.5d0*(xi+xib)
       else
        xi = xin      
       endif
       xf = (r0-rijk)+ga*(dexp(-sign(1.d0,aa)*xi)-1.d0)-bb*xi/(aa*dabs(aa)) 
       xf1= (-bb/aa-aa*ga*dexp(-sign(1.d0,aa)*xi) )/dabs(aa) 
       if (dabs(xf).ge.xerr) then
        goto 100
       elseif ( dabs(xf).gt.xxlim*dabs(xf1*(xi-xi0)).and. dabs(xf).lt.ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
        xerr = UNDEF
        ssii = dabs(xf)
        ssiim= xi
        goto 100
       endif
       if (xerr.eq.UNDEF) xi = ssiim

       endif
400     continue
       ss= s0 + xi/dabs(aa)
       if (xi.ne.0.) then
        xin=dabs( xf/((xi-xi0)*xf1) )
!        accu(1)=accu(1)+1.0
!        accu(2)=accu(2)+xin
!        if (xin.gt.xxlim.and.ntrac.le.10) print *,ntrac,xin,xf1,ss-s0
       endif
       if (loop.eq.1000.or.ss-s0.lt.0.d0) then
        if (ss-s0.lt.0.d0) print *,' ss-s0 is negative '
        print *, '0 time cross =', ss, s0,loop
         print *, aa,uu,um,vv,vm
         print *, iconfig,r0,rijk,xi0,xi00,xib,xf1,xibf,xi,xf
!        STOP
       endif

      endif 

return
end subroutine anil
!_______________________________________________________________________

function dawson (const,daw0,r0,xi0,xi)

 IMPLICIT NONE
 REAL*8  :: dawson,const,daw0,r0,xi0,xi,daw,s15aff

 daw=s15aff(xi)

 dawson=r0*dexp(xi0**2-xi**2)+const*(daw-dexp(xi0**2-xi**2)*daw0)

 return

end function dawson
!_______________________________________________________________________

function errfun (const,erf0,r0,xi0,xi)

USE mod_param
 IMPLICIT NONE
 REAL*8,  PARAMETER :: xilim=3.0d0
 REAL*8  :: errfun,const,erf0,r0,xi0,xi,erf,s15aef,s15adf,hh0,hh
 INTEGER :: i

if (xi0.gt.xilim) then ! complementary error function
 erf  = s15adf(xi)  
elseif(xi0.lt.-xilim) then
 erf  = -s15adf(-xi)
else
! if(xi.ne.0.) stop 3956
 erf  = -s15aef(xi) ! error function
endif

if (dabs(xi).gt.25.d0) then 
 hh0=1.d0
 hh=1.d0
 do i=1,10
  hh0=-0.5d0*hh0*dble(i+i-1)/(xi0*xi0)
  hh=hh+hh0
 enddo
 erf=dexp(xi**2-xi0**2)*hh/xi0
 hh0=1.d0
 hh=1.d0
 do i=1,10
  hh0=-0.5d0*hh0*dble(i+i-1)/(xi*xi)
  hh=hh+hh0
 enddo       
 erf = 1.d0/dsqrt(pi)*((hh/xi)-erf) 
else
 erf = dexp(xi**2)*(erf-erf0)
endif 
errfun = r0*dexp(xi**2-xi0**2) + const*erf 

return

end function errfun


!______________________________________________________________________________
!function dawson(x)
function s15aff(x)

! Returns Dawson's integral for any real x.
! From http://imf.ing.ucv.ve/_1numerical_recipe/Fortran77/f6-10.pdf

USE mod_param
IMPLICIT NONE
INTEGER, PARAMETER ::  NMAX=6  ! Denna ska kollas och testas med högre värden
REAL*8,  PARAMETER ::  H=0.4d0,A1=2.d0/3.d0,A2=0.4d0,A3=2.d0/7.d0

REAL*8 :: s15aff,x,dd1,dd2,e1,e2,sum,x2,xp,xx,c(NMAX),pisqin
INTEGER i,init,n0

SAVE init,c
DATA init/0/ !  Flag is 0 if we need to initialize, else 1.

init=0
pisqin=1.d0/dsqrt(pi)

if(init.eq.0) then
 init=1
 do i=1,NMAX
  c(i)=dexp(-((2.d0*dble(i)-1.d0)*H)**2)
 enddo 
endif

if(dabs(x).lt.0.2d0) then    !  Use series expansion.
 x2=x**2
 s15aff=x*(1.d0-A1*x2*(1.d0-A2*x2*(1.d0-A3*x2)))
else                     !  Use sampling theorem representation.
 xx=dabs(x)
 n0=2*nint(0.5d0*xx/H)
 xp=xx-dble(n0)*H
 e1=dexp(2.d0*xp*H)
 e2=e1**2
 dd1=dble(n0+1)
 dd2=dd1-2.d0
 sum=0.d0
 do i=1,NMAX
  sum=sum+c(i)*(e1/dd1+1.d0/(dd2*e1))
  dd1=dd1+2.d0
  dd2=dd2-2.d0
  e1=e2*e1
 enddo 
 s15aff=pisqin*dsign(dexp(-xp**2),x)*sum 
endif

return
end function s15aff

!______________________________________________________________________________

!FUNCTION erf(x)
FUNCTION s15aef(x)

! USES gammp
!Returns the error function erf(x).

IMPLICIT NONE
REAL*8 ::  s15aef,x,gammp

if(x.lt.0.d0)then
 s15aef=-gammp(0.5d0,x**2)
else
 s15aef= gammp(0.5d0,x**2)
endif

return
end function s15aef

!______________________________________________________________________________


FUNCTION gammp(a,x)

! USES gcf,gser
!Returns the incomplete gamma function P(a; x).

IMPLICIT NONE
REAL*8  :: a,gammp,x,gammcf,gamser,gln

if(x.lt.0.d0 .or.a.le.0.d0) print *, 'bad arguments in gammp'

if(x.lt.a+1.d0) then ! Use the series representation.
 call gser(gamser,a,x,gln)
 gammp=gamser
else ! Use the continued fraction representation
 call gcf(gammcf,a,x,gln)
 gammp=1.d0-gammcf ! and take its complement.
endif

 return
end function gammp

!______________________________________________________________________________

SUBROUTINE gcf(gammcf,a,x,gln)
IMPLICIT NONE

! USES gammln
!Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation as gammcf. Also returns as gln.
!Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accuracy;
!FPMIN is a number near the smallest representable floating-point number.

INTEGER, PARAMETER ::  ITMAX=1000
REAL*8,  PARAMETER ::  EPS=3.d-7,FPMIN=1.d-30
REAL*8  :: a,gammcf,gln,x,an,b,c,d,del,h,gammln
INTEGER :: i

gln=gammln(a)
b=x+1.d0-a   ! Set up for evaluating continued fraction by modied
c=1.d0/FPMIN !Lentz's method (x5.2) with b0 = 0.
d=1.d0/b
h=d

do i=1,ITMAX  ! Iterate to convergence.
 an=-dble(i)*(dble(i)-a)
 b=b+2.d0
 d=an*d+b
 if(dabs(d).lt.FPMIN) d=FPMIN
 c=b+an/c
 if(dabs(c).lt.FPMIN) c=FPMIN
 d=1.d0/d
 del=d*c
 h=h*del
 if(dabs(del-1.d0).lt.EPS) goto 1
enddo
 
print *, 'a too large, ITMAX too small in gcf',dabs(del-1.d0),EPS,del
print *,'cd=',an
stop 43957
1 gammcf=exp(-x+a*dlog(x)-gln)*h  ! Put factors in front.

return
end subroutine gcf

!______________________________________________________________________________

FUNCTION gammln(xx)
IMPLICIT NONE
!Returns the value ln[gamma(xx)] for xx > 0.
!Internal arithmetic will be done in double precision, 
!a nicety that you can omit if five-figure accuracy is good enough.

REAL*8  :: gammln,xx,ser,stp,tmp,x,y,cof(6)
INTEGER j

SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
             24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
              -.5395239384953d-5,2.5066282746310005d0/

x=xx
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*dlog(tmp)-tmp
ser=1.000000000190015d0

do j=1,6
 y=y+1.d0
 ser=ser+cof(j)/y
enddo 

gammln=tmp+dlog(stp*ser/x)
 
return

end function gammln

!______________________________________________________________________________


SUBROUTINE gser(gamser,a,x,gln)

! USES gammln
!Returns the incomplete gamma function P(a; x) evaluated by its series 
! representation as gamser. Also returns  as gln.

IMPLICIT NONE

REAL*8,  PARAMETER ::  EPS=3.d-7
INTEGER, PARAMETER ::  ITMAX=100

REAL*8  :: a,gamser,gln,x,ap,del,sum,gammln
INTEGER n

gln=gammln(a)

if(x.le.0.d0) then
 if(x.lt.0.d0) print *,'x < 0 in gser'
 gamser=0.d0
 return
endif

ap=a
sum=1.d0/a
del=sum

do n=1,ITMAX
 ap=ap+1.d0
 del=del*x/ap
 sum=sum+del
 if(dabs(del).lt.dabs(sum)*EPS) goto 1
enddo
 
print *, 'a too large, ITMAX too small in gser'

1 gamser=sum*dexp(-x+a*dlog(x)-gln)

return
end subroutine gser

!______________________________________________________________________________

FUNCTION s15adf(x) !FUNCTION erfc(x) 

! USES gammp,gammq 
! Returns thecomplementaryerror functionerfc(x). 

IMPLICIT NONE
REAL*8  :: s15adf,x,gammp,gammq 

if(x.lt.0.d0) then 
 s15adf=1.d0+gammp(0.5d0,x**2) 
else 
 s15adf=     gammq(0.5d0,x**2) 
endif 

return 
end function s15adf

!_______________________________________________________________________________

FUNCTION gammq(a,x) 

IMPLICIT NONE
REAL*8  :: a,gammq,x,gammcf,gamser,gln 

if(x.lt.0.d0 .or. a.le.0.d0) print *, 'bad arguments in gammq' 

if(x.lt.a+1.d0) then 
 call gser(gamser,a,x,gln) 
 gammq=1.d0-gamser ! and take its complement. 
else ! Use the continued fraction representation. 
 call gcf(gammcf,a,x,gln) 
 gammq=gammcf 
endif 

return 
end function gammq

!_______________________________________________________________________________
#endif
