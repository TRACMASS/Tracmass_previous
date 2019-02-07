#ifdef timeanalyt 
  !
  ! The subroutines related to the time-analytical option of TRCAMSS are located in this file below 
  ! The schemes are by Vries and Doos (2001), Doos et al. (2013) and Doos et al. (2017)
  !
  !         fn0im----------------------------------------fn0i0
  !          |F^{n}_{i-1}                        F^{n}_{i}|
  !          |                                            |
  !          |                                            |
  !          |                                            |
  !          |                                            |
  !          |                                            |
  !          |                                            |
  !          |                                        -->(r1)
  !          |                                    *       |
  !          |                                *           |
  !          |                             *              |
  !          |                          *                 |
  !          |                        *                   |
  !          |                      *                     |
  !          |                     *                      |
  !          |                    *                       |
  !          |                   *                        |
  !          |                   *                        |
  !          |F^{n-1}_{i-1}      *              F^{n-1}_{i}|
  !         fnmim --------------(r0)----------------------fnmi0
  !         ssiim                                         ssii
  !
subroutine cross_time(ijk,ia,ja,ka,rr0,sp,sn)

  ! subroutine to compute time (sp,sn) when trajectory 
  ! crosses face of box (ia,ja,ka) 
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
  !  intrpr   : Time interpolation constant between 0 and 1 
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
  ! nfn0imerical computation of root:
  ! let f(xi) = 0 be the root to be computed giving the time ss
  ! given some value x.neq.xi, it follows that
  ! 0 = f(xi) = (Taylor expansion) f(x) + (xi-x)*f'(x)
  ! with f' being the derivative of f (---> velocity)
  ! this gives the iteration method: xi_(n+1) = xi_n - f(xi_n)/f'(xi_n)
  ! choose initial value xi00 and use direction rijk for boundary condition 
  ! iteration variable xi_n is optimized during calculations for situations
  ! that would degrade convergency
  ! first, the accuracy parameter xxlim is used for minimizing f(x); then,
  ! xxlim is used to minimize the nfn0imerical error in the time variable ss 
  !
  !

  !
USE mod_precdef
USE mod_param
USE mod_vel
USE mod_traj, only: ntrac
USE mod_grid, only: dxyz,dzt,imt,km, nsp, nsm, dxdy
USE mod_time, only: ts, tt, tseas, intrpr, intrpg
USE mod_loopvars!, only: dsmin
USE mod_timeanalyt

IMPLICIT NONE

INTEGER :: iil,ijk,ia,ja,ka
REAL (DP)  :: sp,sn,rr0
REAL (QP)  :: dzs,dzu1,dzu2,alfa

!_______________________________________________________________________________
r0=rr0
sp=UNDEF ; sn=UNDEF

#ifdef twodim  
if(ijk==3) return
#endif 

if(dxyz<EPS) stop 8701
if(dsmin==1.0) stop 8702
s0=tt/dxyz
ss0=dble(idint(ts))*tseas/dxyz
if(s0<ss0) s0=ss0
!ss0=dble(idint(ts))*tseas/((intrpg * dzt(ia,ja,ka,nsp) + intrpr * dzt(ia,ja,ka,nsm))*dxdy(ia,ja))
!ss0=dble(idint(ts))*tseas/(dzt(ia,ja,ka,nsm)*dxdy(ia,ja))
!print *,tt,dble(idint(ts))*tseas,tt-dble(idint(ts))*tseas,

looop=0 ; rijk=0. ; ss=UNDEF ; f0=0.0 ; f1=0.0

if(ijk==1) then
 ii=ia
 iim=ia-1
 iil=iim
 if(iim==0) iil = imt
 fn0i0=uflux(ii ,ja,ka,nsm)
 fn0im=uflux(iil,ja,ka,nsm)
 fnmi0=uflux(ii ,ja,ka,nsp)
 fnmim=uflux(iil,ja,ka,nsp)
!#ifdef turb   
! if(r0/=dble(ja)) then
!  fn0i0=fn0i0+upr(7,2)  
!  fnmi0=fnmi0+upr(1,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fn0i0=fn0i0+upr(7,1) 
!  fnmi0=fnmi0+upr(1,1)  
! endif
! if(r0=/dble(ja-1)) then
!  fn0im=fn0im+upr(8,2)
!  fnmim=fnmim+upr(2,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr(8,1)  
!  fnmim=fnmim+upr(2,1)  
! endif
!#endif
elseif(ijk==2) then
 ii=ja
 iim=ja-1
 iil=iim
 fn0i0=vflux(ia,ii ,ka,nsm)
 fn0im=vflux(ia,iil,ka,nsm)
 fnmi0=vflux(ia,ii ,ka,nsp)
 fnmim=vflux(ia,iil,ka,nsp)
!#ifdef turb   
! if(r0/=dble(ii)) then
!  fnmi0=fnmi0+upr(3,2)  
!  fn0i0=fn0i0+upr(9,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fnmi0=fnmi0+upr(3,1)  
!  fn0i0=fn0i0+upr(9,1) 
! endif
! if(r0/=dble(iim)) then
!  fn0im=fn0im+upr( 4,2)
!  fnmim=fnmim+upr(10,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr( 4,1)  
!  fnmim=fnmim+upr(10,1)  
! endif
!#endif
elseif(ijk==3) then
 ii=ka
 iim=ka-1
 iil=iim
 fn0i0=wflux(ii ,nsm)
 fn0im=wflux(iil,nsm)
 fnmi0=wflux(ii ,nsp)
 fnmim=wflux(iil,nsp)
!#ifdef turb   
! if(r0/=dble(ii)) then
!  fnmi0=fnmi0+upr( 5,2)  
!  fn0i0=fn0i0+upr(11,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fnmi0=fnmi0+upr( 5,1)  
!  fn0i0=fn0i0+upr(11,1) 
! endif
! if(r0/=dble(iim)) then
!  fn0im=fn0im+upr( 6,2)
!  fnmim=fnmim+upr(12,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr( 6,1)  
!  fnmim=fnmim+upr(12,1)  
! endif
!#endif

!#ifdef zgrid3D 

  dzu1=dzt(ia,ja,ka,nsm) ! layer thickness at time step n-1
  dzu2=dzt(ia,ja,ka,nsp) ! layer thickness at time step n
!  dzs= intrpg*dzu1 + intrpr*dzu2   ! layer thickness at time interpolated for "present" ! (wrong?? 
  dzs= intrpr*dzu1 + intrpg*dzu2   ! layer thickness at time interpolated for "present" 
!  if(abs(dzu1)<eps) stop 8705
  if(dzu1==0.) stop 8705
!  if(abs(dzu2)<eps) stop 8706
  if(dzu2==0.) stop 8706
  f0=dzs/dzu1
  f1=dzs/dzu2
! print *,tt,ts, intrpg, intrpr,dzu1,dzu2,dzs,f0,f1
  fn0i0=fn0i0*f0
  fn0im=fn0im*f0
  fnmi0=fnmi0*f1
  fnmim=fnmim*f1
  if(abs(f0)<=eps) stop 8707
  if(abs(f1)<=eps) stop 8708
  f0=1.0/f0
  f1=1.0/f1
!#elif defined zgrid3D && defined freesurface
! if (ka==km) then
!  dzs= dz(km)+intrpg*hs(ia,ja,nsm)+intrpr*hs(ia,ja,nsp)
!  dzu1=dz(km)+hs(ia,ja,nsm)
!  dzu2=dz(km)+hs(ia,ja,nsp)
!  if(abs(dzu1)<eps) stop 8705
!  if(abs(dzu2)<eps) stop 8706
!  f0=dzs/dzu1
!  f1=dzs/dzu2
!  fn0i0=fn0i0*f0
!  fn0im=fn0im*f0
!  fnmi0=fnmi0*f1
!  fnmim=fnmim*f1
!  if(abs(f0)<eps) stop 8707
!  if(abs(f1)<eps) stop 8708
!  f0=1.0/f0
!  f1=1.0/f1
! endif 
!#endif /*zgrid3D*/


endif


alfa = -(fnmi0-fnmim-fn0i0+fn0im)
!print *,'alfa=',alfa,fnmi0,fnmim,fn0i0,fn0im
! ----- a > 0
!     if (alfa>EPS) then
     if (alfa>0.) then
      call apos !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)
! ----- a < 0      
!     elseif (alfa<-EPS) then
     elseif (alfa<0.) then
      call amin !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)
! ----- a = 0 
      else
      if (alfa/=0.) print *,'anil eps=',alfa
      call anil !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)
     endif
! translate direction and positions for old to new tracmass
!print *,'cross_time ',ijk,'u+',fn0i0,fnmi0,' u-',fn0im,fnmim,'ss',ss,ss0,' alfa',alfa
      if (rijk==dble(ii)) then
       sp=ss-s0
       sn=UNDEF
      elseif(rijk==dble(iim)) then
       sn=ss-s0
       sp=UNDEF
      else
       sp=UNDEF
       sn=UNDEF
      endif
      if(sp==0.0) sp=UNDEF
      if(sn==0.0) sn=UNDEF
      
!    if(sp/=undef .or. sn/=undef) print *,'cross_time',ss,s0,sp,sn
!      if(ijk==3) stop 49678
!stop 4856
return
end subroutine cross_time
!_______________________________________________________________________

subroutine pos_time(ijk,ia,ja,ka,rr0,r1)
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

USE mod_precdef
USE mod_param
USE mod_vel
USE mod_traj, only: ntrac
USE mod_grid, only: dxyz,dzt,imt,km, nsp, nsm
!USE mod_grid
USE mod_time!, only: tt,ts,tseas,intrpr, intrpg
USE mod_loopvars!, only: ds, dsmin
USE mod_precdef
USE mod_timeanalyt
IMPLICIT NONE

REAL (DP),  PARAMETER ::  xilim=3.0!,xxlim=1.d-9
INTEGER :: ijk,ia,ja,ka,iil
REAL (DP)  :: rr0,r1
REAL (DP)  :: xi,xi0,const,ga,erf0,daw0
REAL (DP)  :: dawson1, dawson2,s15adf,s15aef,errfun
REAL (QP)  :: alfa,beta,dzs,dzu1,dzu2


#ifdef twodim  
if(ijk==3) then
r1=r0
return
endif
#endif 

r0=rr0
s0=tt/dxyz-ds
if(s0<ss0) s0=ss0
!ss0=dble(idint(ts))*tseas/dxyz
ss=ts*tseas/dxyz

!print *,s0-ss0,ss-ss0,ss-s0




f0=0. ; f1=0.

if(ijk==1) then
 ii=ia
 iim=ia-1
 iil=iim
 if(iim==0) iil = imt
 fn0i0=uflux(ii ,ja,ka,nsm)
 fn0im=uflux(iil,ja,ka,nsm)
 fnmi0=uflux(ii ,ja,ka,nsp)
 fnmim=uflux(iil,ja,ka,nsp)
!#ifdef turb   
! if(r0/=dble(ja)) then
!  fn0i0=fn0i0+upr(7,2)  
!  fnmi0=fnmi0+upr(1,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fn0i0=fn0i0+upr(7,1) 
!  fnmi0=fnmi0+upr(1,1)  
! endif
! if(r0/=dble(ja-1)) then
!  fn0im=fn0im+upr(8,2)
!  fnmim=fnmim+upr(2,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr(8,1)  
!  fnmim=fnmim+upr(2,1)  
! endif
!#endif
elseif(ijk==2) then
 ii=ja
 iim=ja-1
 iil=iim
 fn0i0=vflux(ia,ii ,ka,nsm)
 fn0im=vflux(ia,iil,ka,nsm)
 fnmi0=vflux(ia,ii ,ka,nsp)
 fnmim=vflux(ia,iil,ka,nsp)
!#ifdef turb   
! if(r0/=dble(ii)) then
!  fnmi0=fnmi0+upr(3,2)  
!  fn0i0=fn0i0+upr(9,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fnmi0=fnmi0+upr(3,1)  
!  fn0i0=fn0i0+upr(9,1) 
! endif
! if(r0/=dble(iim)) then
!  fn0im=fn0im+upr( 4,2)
!  fnmim=fnmim+upr(10,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr( 4,1)  
!  fnmim=fnmim+upr(10,1)  
! endif
!#endif
elseif(ijk==3) then
 ii=ka
 iim=ka-1
 iil=iim
 fn0i0=wflux(ii ,nsm)
 fn0im=wflux(iil,nsm)
 fnmi0=wflux(ii ,nsp)
 fnmim=wflux(iil,nsp)
!#ifdef turb   
! if(r0/=dble(ii)) then
!  fnmi0=fnmi0+upr( 5,2)  
!  fn0i0=fn0i0+upr(11,2)  
! else  ! add u' from previous iterative time step if on box wall
!  fnmi0=fnmi0+upr( 5,1)  
!  fn0i0=fn0i0+upr(11,1) 
! endif
! if(r0/=dble(iim)) then
!  fn0im=fn0im+upr( 6,2)
!  fnmim=fnmim+upr(12,2)
! else  ! add u' from previous iterative time step if on box wall
!  fn0im=fn0im+upr( 6,1)  
!  fnmim=fnmim+upr(12,1)  
! endif
!#endif

!#ifdef zgrid3D 

  dzu1=dzt(ia,ja,ka,nsm) ! layer thickness at time step n-1
  dzu2=dzt(ia,ja,ka,nsp) ! layer thickness at time step n
!  dzs= intrpg*dzu1 + intrpr*dzu2   ! layer thickness at time interpolated for "present" ! (wrong?? 
  dzs= intrpr*dzu1 + intrpg*dzu2   ! layer thickness at time interpolated for "present" 
  if(abs(dzu1)<eps) stop 4966
  if(abs(dzu2)<eps) stop 4967
  f0=dzs/dzu1
  f1=dzs/dzu2
  fn0i0=fn0i0*f0
  fn0im=fn0im*f0
  fnmi0=fnmi0*f1
  fnmim=fnmim*f1
  if(abs(f0)<=eps) stop 4968
  if(abs(f1)<=eps) stop 4969
  f0=1.0/f0
  f1=1.0/f1
!  dzs= dz(km)+intrpg*hs(ia,ja,nsm)+intrpr*hs(ia,ja,nsp)



endif


!print *,ijk,'u+',fn0i0,fnmi0,' u-',fn0im,fnmim
!if(fn0i0==fnmi0 .and. fn0im==fnmim) then
! print *,'pos'
! print *,ijk,ia,ja,ka,r0,r1,ds,intrpr
! print *,fn0i0,fnmi0,fn0im,fnmim
! call pos_orgn(ijk,ia,ja,ka,r0,r1,ds,intrpr) 
!
! return
!endif

if(abs(dsmin)<=eps) stop 4981
alfa = -(fnmi0-fnmim-fn0i0+fn0im)/dsmin
beta = fn0im-fn0i0

if(alfa/=0.) then
 if(f0==0.) then
  ga = -iim + (fnmim-fn0im)/(fnmi0-fnmim-fn0i0+fn0im)
  const = (fn0im*fnmi0-fn0i0*fnmim)/(fnmi0-fnmim-fn0i0+fn0im)
 else
  ga = -iim + (f1*fnmim-f0*fn0im)/(fnmi0-fn0i0-fnmim+fn0im)
  const = (f0*fn0im*(fnmi0-fnmim)+f1*fnmim*(fn0im-fn0i0))/(fnmi0-fn0i0-fnmim+fn0im)
 endif
endif


if (alfa>0.) then
 const=2.*const/sqrt(2.*alfa)
 xi0=(beta+alfa*(s0-ss0))/sqrt(2.*alfa)
 xi =(beta+alfa*(ss-ss0))/sqrt(2.*alfa)
 daw0 = dawson2(xi0) 
 r1 = dawson1 (const,daw0,r0+ga,xi0,xi) -ga
!print *,'alfa>0',f0,r1,const,daw0,r0+ga,xi0,xi,ga,dawson1 (const,daw0,r0+ga,xi0,xi)
elseif (alfa<0.0) then
 const=const*sqrt(pi/(-2._qp*alfa))
 xi0=(beta+alfa*(s0-ss0))/sqrt(-2.*alfa)
 xi =(beta+alfa*(ss-ss0))/sqrt(-2.*alfa)
 if (xi0>xilim) then ! complementary error function
  erf0 = s15adf(xi0) !       erf0 = s15adf(xi0, ifai l)
 elseif(xi0<-xilim) then
  erf0 = - s15adf(-xi0) !       erf0 = - s15adf(-xi0, ifail)
 else ! error function
  erf0 = - s15aef(xi0)     !       erf0 = - s15aef(xi0, ifail)
 endif
 r1 = errfun (const,erf0,r0+ga,xi0,xi) -ga 
elseif (alfa==0.0) then
 if (f0==0.) then
  f0=1.
  f1=1.
 endif 
 alfa = fn0im-fn0i0
 beta =(f0*fn0im-f1*fnmim)/dsmin
 xi = ss-s0
 if (alfa==0.0) then
! if (abs(alfa)<EPS) then
  r1 = r0 - xi*(-f0*fn0im+0.5*beta*(s0-ss0+ss-ss0))
!  print *,'alfa,r1',alfa,r1
 else
  ga = (r0-iim) + (beta*(s0-ss0-1.0/alfa)-f0*fn0im)/alfa
  r1 = r0 + ga*(exp(-alfa*xi) - 1.0) -beta*xi/alfa
 endif 
endif

if (r1-dble(ii) >0.0 .and. r1-dble(ii) <EPS) r1=dble(ii )-EPS
if (dble(iim)-r1>0.0 .and. dble(iim)-r1<EPS) r1=dble(iim)+EPS

if(abs(r0-r1)>1.0) then
! print *,'warning time analytical solution outside the box',ntrac,ijk,r0,r1,r0-r1,ii,iim, alfa, beta
! print *,'ga',ga, -dble(iim) + (f1*fnmim-f0*fn0im)/(fnmi0-fn0i0-fnmim+fn0im), -iim + (f1*fnmim-f0*fn0im)/(fnmi0-fn0i0-fnmim+fn0im)
! print *,'xi,xi0',xi,xi0
! print *,'s0-ss0,ss-ss0',s0-ss0,ss-ss0,s0,ss,ss0
! print *,'ds',ds,dsmin
! print *,'fnmi0-fnmim-fn0i0+fn0im',fnmi0-fnmim-fn0i0+fn0im
! print *,'fnmi0,fnmim,fn0i0,fn0im',fnmi0,fnmim,fn0i0,fn0im
! print *,'dxyz',dxyz
! print *,'dzu1, dzs, dzu2', dzu1, dzs, dzu2
! print *,'fn0im------fn0i0', fn0im, fn0i0
! print *,'fnmim------fnmi0', fnmim, fnmi0
! stop 48956
 r1=r0
endif




! This needs to be included only if airborneError is commented out in loop
!if(ijk==3) then
! if(r1>dble(KM)) then
!  r1=dble(KM)-0.5
! endif
!endif

!print *,'time',ijk,r1,r0

return
end subroutine pos_time
!_______________________________________________________________________

subroutine apos !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)
! computation of time (ss) at crossing for a>0 

USE mod_precdef
USE mod_param
USE mod_vel
USE mod_traj, only: ntrac
USE mod_precdef
USE mod_loopvars , only: dsmin,ss0
USE mod_timeanalyt
IMPLICIT NONE

!REAL (DP),  PARAMETER ::  xxlim=1.d-7
INTEGER :: iconfig,i,loop
REAL (QP)  :: ssii,ssiim
REAL (QP)  :: alfa,beta
REAL (DP)  :: xi0,xin,xerr,xf1,xi,xi00,xia,xib,xibf
REAL (DP)  :: ga,const,daw0,xf,xf2,xiaf,dawson1,dawson2

!print *,'ii',ii

alfa = -(fnmi0-fnmim-fn0i0+fn0im)/dsmin                 ! same as alpha in paper
!if(alfa<EPS) stop 8771
if(alfa==0.) stop 8771
beta = fn0im-fn0i0                                ! not same as beta in paper
xi0=(beta+alfa*(s0-ss0))/sqrt(2.*alfa)  ! same as xi0 in paper where s=s0-ss0

! 1) ii or iim ---> land point 
if (fn0i0==0.0 .and. fnmi0==0.0) then
 xin =-log(dble(ii)-r0) 
 if (xi0<0.0 .and. xi0*xi0>=xin) then
  rijk=iim
  xin=-sqrt(xi0*xi0-xin)
  ss= 2.*(xin-xi0)/sqrt(2.*alfa) + s0
  return
 else ! no solution
  return 
 endif
elseif (fn0im==0.0 .and. fnmim==0.0 ) then
! xin =-log(r0-iim)
 xin =-log(r0-dble(iim))
 if (xi0<0.0 .and. xi0*xi0>=xin) then
  rijk=ii
  xin=-sqrt(xi0*xi0-xin)
  ss= 2.*(xin-xi0)/sqrt(2.*alfa) + s0
  return
 else ! no solution
  return
 endif
endif

if (f0==0.0) then
 ga = -dble(iim) + (fnmim-fn0im)/(fnmi0-fn0i0-fnmim+fn0im)   ! gamma/alfa
 const = (fn0im*fnmi0-fn0i0*fnmim)/(fnmi0-fn0i0-fnmim+fn0im)       ! (beta*gamma-alpha*delta)/alpha
 ssii=fn0i0-fnmi0
 if(ssii/=0.0) ssii = fn0i0/ssii
 ssiim=fn0im-fnmim
 if(ssiim/=0.0) ssiim = fn0im/ssiim
else
 ga = -dble(iim) + (f1*fnmim-f0*fn0im)/(fnmi0-fn0i0-fnmim+fn0im)
 const = (f0*fn0im*(fnmi0-fnmim)+f1*fnmim*(fn0im-fn0i0))/(fnmi0-fn0i0-fnmim+fn0im)
 ssii = fn0i0-fnmi0+fn0im*(f0-1.0)+fnmim*(1.0-f1)
 if(abs(ssii)>EPS) then
  ssii = (fn0i0+fn0im*(f0-1.0))/ssii
 else
  stop 8777
 endif
! ssiim= f0*fn0im/(f0*fn0im-f1*fnmim)
 ssiim= (f0*fn0im-f1*fnmim)
 if(abs(ssiim)>EPS) then
  ssiim= f0*fn0im/ssiim
 else
  ssiim= EPS ! This to correct when by accident f0*fn0im = f1*fnmim
  print *,'f0*fn0im = f1*fnmim',ssiim,f0,fn0im,f1,fnmim
  ssiim= f0*fn0im/ssiim
!  stop 8779
 endif
endif


const = 2.*const/sqrt(2.*alfa)
daw0 = dawson2(xi0)    
! 'velocity' at xi0
xf1 = const -(xi0+xi0)*(r0+ga)  

! 2) ii direction (four velocity configurations at edge)
if (ssii<=0.0 .or. ssii>=1.0) then
 if (fn0i0>0.0 .or. fnmi0>0.0) then
! 2a) configuration: +++
  iconfig = 1
  xi00=(fnmim-fnmi0)/sqrt(2.*alfa)
 else
! 2b) configuration: --- no solution, try rijk=iim
 ! print *,'2b configuration: --- no solution, try rijk=iim=',iim,ssii,fn0i0,fnmi0
  rijk=iim
  goto 3000
 endif 
else
 if (fnmi0>0.0) then
! 2c) configuration: -- 0 ++
  iconfig = 3   
  xi00=(fnmim-fnmi0)/sqrt(2.0*alfa)
 else
! 2d) configuration: ++ 0 --
  if (s0-ss0>=ssii*dsmin) goto 3000
  iconfig = 4
  xi00=(beta+alfa*ssii*dsmin)/sqrt(2.0*alfa) 
 endif 
endif

xf=dawson1(const,daw0,r0+ga,xi0,xi00)-ga-dble(ii)

if (xf<0.0) goto 3000
rijk=ii
xib =xi00
xibf=xf
 if (xf1>0.0) then
  if ( xf1*(xi00-xi0) > -(r0-rijk) ) then
   xi00=xi0
   xf=r0-rijk
  endif
 else
  if (ssiim>0.0 .and. ssiim<1.0 .and.fnmim>0.0) then
!  check crossing at 3c) configuration: -- 0 ++
   xf2=(beta+alfa*ssiim*dsmin)/sqrt(2.0*alfa)
   xf2=dawson1(const,daw0,r0+ga,xi0,xf2)-ga-dble(iim)
   if (xf2<0.0) goto 3000
  endif                
  if (iconfig==4) then
   xi00=(beta+alfa*ssiim*dsmin)/sqrt(2.0*alfa)
   xf=dawson1(const,daw0,r0+ga,xi0,xi00)-ga-dble(ii)
  endif
 endif     
 goto 1000

3000  continue


! 3) iim direction (four velocity configurations at edge)
     if (ssiim<=0.0 .or. ssiim>=1.0) then
      if (fn0im>0.0 .or. fnmim>0.0) then
! 3a) configuration: +++ no solution
       rijk=iim
       return
      else
! 3b) configuration: ---
       iconfig = 6
       xi00=(fnmim-fnmi0)/sqrt(2.0*alfa)
      endif 
     else
      if (fnmim>0.0) then
! 3c) configuration: -- 0 ++
       if (s0-ss0>=ssiim*dsmin) return
       iconfig = 7   
       xi00=(beta+alfa*ssiim*dsmin)/sqrt(2.0*alfa) 
      else
! 3d) configuration: ++ 0 --
       iconfig = 8
       xi00=(fnmim-fnmi0)/sqrt(2.0*alfa) 
      endif 
     endif
     xf=dawson1(const,daw0,r0+ga,xi0,xi00)-ga-dble(iim)
     if (xf>0.0) return
     rijk=iim
     xib =xi00
     xibf=xf
     if (xf1<0.0) then
      if (xf1*(xi00-xi0)<-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     elseif (iconfig==7) then
      xi00=(beta+alfa*ssii*dsmin)/sqrt(2.0*alfa)
      xf=dawson1(const,daw0,r0+ga,xi0,xi00)-ga-dble(iim)
     endif      

1000  continue
! calculation of root
     if (looop==100) then    ! why 100 ?
      xi = -0.05 + xi0
      do i=1,100              ! why 100 ?
       xi=xi+0.05
       xf2= dawson1 (const,daw0,r0+ga,xi0,xi)
       xf = xf2-ga-rijk
       xf1= const-(xi+xi)*xf2 
!       print *, ' looop', xi, xf, xf1
      enddo
      looop = 1000
     else

     xi = xi00
     xf1= const-(xi+xi)*(xf+rijk+ga)
     xia = xi0
     xiaf= r0-rijk 
     xerr = xxlim
     ssii = xxlim
100   continue
     looop=looop+1
     if(looop==90) then
      print *, ' looop stop1', xi0, xi00, xi, xf, ntrac, looop
      if (abs(xf)>100.*xxlim) looop=1000
!      print *, 'looopnr', abs(xf),xxlim,100*xxlim, looop
!       if (xerr==UNDEF) looop=1000
      goto 200
     endif
     if(abs(xf1)<=eps) stop 4992
     xin = xi - xf/xf1
     if ( abs(xf)>100.0 .or.(xf*xibf>0.0 .and.(xf1*xibf<0.0 .or.xin<xia)) ) then
      xib = xi
      xibf= xf
      xi = 0.5*(xi+xia)
     elseif (xf*xibf<0.0 .and. (xf1*xibf<0.0 .or.xin>xib) ) then
      xia = xi
      xiaf= xf
      xi = 0.5*(xi+xib)
     else
 !     ic=3   What is this????????? Should it be iconfig=3 ?????????????????????
       iconfig=3 ! ?????????????
! print *,'What is this?????',xi,xin
      xi = xin 
     endif

     xf2= dawson1 (const,daw0,r0+ga,xi0,xi)
     xf = xf2-ga-rijk
     xf1= const-(xi+xi)*xf2 
!     if(looop==1000) print *, ' looop', xi, xf, xf1
     if (abs(xf)>=xerr) then
      goto 100
     elseif ( abs(xf)>xxlim*abs(xf1*(xi-xi0)) .and. abs(xf)<ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
      xerr = UNDEF
      ssii = abs(xf)
      ssiim= xi
      goto 100
     endif
     if (xerr==UNDEF) xi = ssiim

     endif
200   continue
     ss= 2.0*(xi-xi0)/sqrt(2.0*alfa) + s0
!     nstat(iconfig,1,1)=nstat(iconfig,1,1)+1       ! removed statistics
!     nstat(iconfig,2,1)=nstat(iconfig,2,1)+looop
     if (xi/=xi0) then
      xin=abs( xf/((xi-xi0)*xf1) )
!      accu(1)=accu(1)+1.0
!      accu(2)=accu(2)+xin
!      if (xin>xxlim.and.ntrac<=100) print *,ntrac,xin,xf1,ss-s0
     endif
     if (looop==1000.or.ss-s0<0.0) then
      if (ss-s0<0.0) print *,' ss-s0 is negative '
      print *, '+ time cross =', ss-s0, ss, s0,looop
      print *, 'alfa',alfa,fn0i0,fn0im,fnmi0,fnmim
      print *, 'iconfig',iconfig,r0,rijk
       print *, 'xi',xi-xi0,xi,xi0,2.*(xi-xi0)/sqrt(2.*alfa),s0
      print *,'ii',ii,iim,abs(xf),xxlim,abs(xf1*(xi-xi0)),ssii,xerr
      STOP 7235
     endif

return
end subroutine apos
!_______________________________________________________________________

subroutine amin !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)

! computation of time (ss) at crossing for a<0 

USE mod_precdef
USE mod_param
USE mod_vel
USE mod_loopvars, only: dsmin,ss0
USE mod_timeanalyt
IMPLICIT NONE

REAL (DP),  PARAMETER ::  xilim=3.0!,xxlim=1.d-8
INTEGER    :: iconfig,i
REAL (QP)  :: ssii,ssiim
REAL (DP)  :: xi0,xin,xerr,xf1,xi,xi00,xia,xib,xibf
REAL (DP)  :: xiaf,xf,xf2
REAL (DP)  :: alfa,beta,ga,s15aef,s15adf,erf0,const,errfun

if(abs(dsmin)<eps) stop 8801
alfa = -(fnmi0-fnmim-fn0i0+fn0im)/dsmin
beta = fn0im-fn0i0
xi0=(beta+alfa*(s0-ss0))/sqrt(-2.0*alfa)
if(alfa>=EPS) then
! print *,alfa,fnmi0,fn0i0,fnmim,fn0im,dsmin,xi0
 stop 8802
endif
!print *,'alfabera=',alfa,beta,fnmi0,fn0i0,fnmim,fn0im,dsmin,xi0,f0


! 1) ii or iim ---> land point 
     if (fn0i0==0.0 .and. fnmi0==0.0) then
      rijk=iim
      xin =-log(dble(ii)-r0) 
      xin =-sqrt(xi0*xi0+xin)
      ss= 2.0*(xi0-xin)/sqrt(-2.0*alfa) + s0
      return
     elseif (fn0im==0.0 .and. fnmim==0.0 ) then
      rijk=ii
      xin =-log(r0-dble(iim))
      xin =-sqrt(xi0*xi0+xin)
      ss= 2.0*(xi0-xin)/sqrt(-2.0*alfa) + s0
      return
     endif 

     if (f0==0.0) then
      if(abs(fnmi0-fn0i0-fnmim+fn0im)<EPS) stop 8803
      ga = -dble(iim) + (fnmim-fn0im)/(fnmi0-fn0i0-fnmim+fn0im)
!      print *,'ga=',ga
      const = (fn0im*fnmi0-fn0i0*fnmim)/(fnmi0-fn0i0-fnmim+fn0im)
      ssii=fn0i0-fnmi0
      if(ssii/=0.0) ssii = fn0i0/ssii
      ssiim=fn0im-fnmim
      if(ssiim/=0.0) ssiim = fn0im/ssiim
     else
      if(fnmi0-fn0i0-fnmim+fn0im==0.0) stop 8806
      ga = -dble(iim) + (f1*fnmim-f0*fn0im)/(fnmi0-fn0i0-fnmim+fn0im)
      const = (f0*fn0im*(fnmi0-fnmim)+f1*fnmim*(fn0im-fn0i0))/(fnmi0-fn0i0-fnmim+fn0im)
      if(fn0i0-fnmi0+fn0im*(f0-1.0)+fnmim*(1.0-f1)==0.0) then
       fn0i0=fn0i0*1.000001 ! this might be a dirty fix in order tro avoid division by zero
     !  stop 8808
      endif
      ssii = (fn0i0+fn0im*(f0-1.0))/(fn0i0-fnmi0+fn0im*(f0-1.0)+fnmim*(1.0-f1))
      if(f0*fn0im-f1*fnmim==0.d0) then
       print *,'8810',f0*fn0im,f1*fnmim
       ssiim= f0*fn0im/eps
!      stop 8810
      else
       ssiim= f0*fn0im/(f0*fn0im-f1*fnmim)
      endif
     endif
     const = const*sqrt(pi/(-2.0*alfa))
     if (xi0>xilim) then
! complementary error function
      erf0 = s15adf(xi0)
     elseif(xi0<-xilim) then
      erf0 = -s15adf(-xi0)
     else
! error function
      erf0 = -s15aef(xi0)
     endif
! minus-velocity at xi0
     xf1 = xi0*(r0+ga)-(const/sqrt(pi)) 

! 2) ii direction (four velocity configurations at edge)
     if (ssii<=0.0 .or. ssii>=1.0) then
      if (fn0i0>0.0 .or. fnmi0>0.0) then
! 2a) configuration: +++
       iconfig = 1
       xi00=(fnmim-fnmi0)/sqrt(-2.0*alfa)
      else
! 2b) configuration: --- no solution, try rijk=iim
       rijk=iim
       goto 4000
      endif 
     else
      if (fnmi0>0.0) then
! 2c) configuration: -- 0 ++
       iconfig = 3   
       xi00=(fnmim-fnmi0)/sqrt(-2.0*alfa)
      else
! 2d) configuration: ++ 0 --
       if (s0-ss0>=ssii*dsmin) goto 4000
       iconfig = 4
       xi00=(beta+alfa*ssii*dsmin)/sqrt(-2.0*alfa) 
      endif 
     endif
     xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(ii)
     if (xf<0.0) goto 4000
     rijk=ii
     xib =xi00
     xibf=xf
     if (xf1<0.0) then
      if (xf1*(xi00-xi0)>-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     else
      if (ssiim>0.0 .and. ssiim<1.0 .and. fnmim>0.0) then
! check crossing at 3c) configuration: -- 0 ++
       xf2=(beta+alfa*ssiim*dsmin)/sqrt(-2.0*alfa)
       xf2=errfun(const,erf0,r0+ga,xi0,xf2)-ga-dble(iim)
       if (xf2<0.0) goto 4000
      endif
      if (iconfig==4) then
       xi00=(beta+alfa*ssiim*dsmin)/sqrt(-2.0*alfa)
       xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(ii)
      endif
     endif 
     goto 2000

4000  continue
!print *,'ssssss',ss, ssiim,fn0im,fnmim
! 3) iim direction (four velocity configurations at edge)
     if (ssiim<=0.0 .or. ssiim>=1.0) then
      if (fn0im>0.0 .or. fnmim>0.0) then
! 3a) configuration: +++ no solution
!     print *,'3a no solution'
       return
      else
! 3b) configuration: ---
       iconfig = 6
       xi00=(fnmim-fnmi0)/sqrt(-2.0*alfa)
!       print *,'xi00',xi00
      endif 
     else
      if (fnmim>0.0) then
! 3c) configuration: -- 0 ++
       if (s0-ss0>=ssiim*dsmin) return
       iconfig = 7   
       xi00=(beta+alfa*ssiim*dsmin)/sqrt(-2.0*alfa) 
      else
! 3d) configuration: ++ 0 --
       iconfig = 8
       xi00=(fnmim-fnmi0)/sqrt(-2.0*alfa) 
      endif 
     endif
     xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(iim)
     if (xf>0.0) return
     rijk=iim
     xib =xi00
     xibf=xf
     if (xf1>0.0) then
      if (xf1*(xi00-xi0)<-(r0-rijk)) then
       xi00=xi0
       xf=r0-rijk
      endif
     elseif (iconfig==7) then   
      xi00=(beta+alfa*ssii*dsmin)/sqrt(-2.0*alfa)
      xf=errfun(const,erf0,r0+ga,xi0,xi00)-ga-dble(iim)
     endif

2000  continue
! calculation of root
     if (looop==100) then
      xi = 0.005 + xi0
      do i=1,200
       xi=xi-0.005
       xf2= errfun (const,erf0,r0+ga,xi0,xi)
       xf = xf2-ga-rijk
       xf1= xi*xf2-(const/sqrt(pi)) 
       print *, ' looop', xi, xf, xf1
      enddo
      looop = 1000
     else

     xi = xi00
     xf1= xi*(xf+rijk+ga)-(const/sqrt(pi))
     xia = xi0
     xiaf= r0-rijk
     xerr = xxlim
     ssii = xxlim      
200   continue
     looop=looop+1
     if(looop==90) then
      print *, ' looop stop2', xi0, xi00, xi, xf
      if (abs(xf)>100.*xxlim) looop=1000
!       if (xerr==UNDEF) looop=1000
      goto 300
     endif
     xf1 = xf1+xf1
!     xin = xi - xf/xf1
     if(xf1==1.0) then
      xin=xi 
     else
      xin=xi-xf/xf1
     endif
     if ( abs(xf)>100.0 .or.(xf*xibf>0.0 .and. (xf1*xibf>0.0 .or.xin>xia)) )then
      xib = xi
      xibf= xf
      xi = 0.5*(xi+xia)
     elseif (xf*xibf<0.0 .and.  (xf1*xibf>0.0 .or.xin<xib) ) then
      xia = xi
      xiaf= xf
      xi = 0.5*(xi+xib)
     else
      xi = xin      
     endif
     xf2= errfun (const,erf0,r0+ga,xi0,xi)
     xf1= xi*xf2-(const/sqrt(pi))
     xf = xf2-ga-rijk
!      print *, ' looop', xi, xf, xf1
     if (abs(xf)>=xerr) then
      goto 200
     elseif ( abs(xf)>xxlim*abs(xf1*(xi-xi0)).and. abs(xf)<ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
      xerr = UNDEF
      ssii = abs(xf)
      ssiim= xi
      goto 200
     endif
     if (xerr==UNDEF) xi = ssiim

     endif
300   continue
     ss= 2.0*(xi0-xi)/sqrt(-2.0*alfa) + s0 
     if (xi/=xi0) xin=abs( xf/((xi-xi0)*xf1) )
     if (looop==1000.or.ss-s0<0.0) then
      if (ss-s0<0.0) then
       print *, ' ss-s0 is negative '
       stop 48967
      endif
     endif
     
     

return
end subroutine amin
!_______________________________________________________________________

subroutine anil !(r0,rijk,s0,ss,ss0,fn0i0,fn0im,fnmi0,fnmim,f0,f1)

! computation of time (ss) at crossing for alfa=0 

USE mod_precdef
USE mod_param
USE mod_vel
USE mod_loopvars, only: dsmin,ss0
USE mod_timeanalyt
IMPLICIT NONE

!REAL (DP),  PARAMETER :: xxlim=1.d-7
REAL (QP)  :: ssii,ssiim
REAL (DP)  :: alfa,beta,ga
REAL (DP)  :: xerr,xf,xf1,xf2,xi,xi0,xi00,xia,xiaf,xib,xibf,xin
INTEGER :: i,iconfig

if (f0==0.0) then
 f0=1.0
 f1=1.0
endif

alfa = fn0im-fn0i0
beta =(f0*fn0im-f1*fnmim)/dsmin 
!if(abs(beta)<EPS) then
! print *,'beta=',beta,f0*fn0im,f1*fnmim,dsmin
!endif

if (alfa==0.0) then
 if (beta==0.0) then
  ga = f0*fn0im
  if(ga>0.0) then
   rijk=ii
   ss=s0+(rijk-r0)/ga
  elseif (ga<0.0) then
   rijk=iim
   ss=s0+(rijk-r0)/ga
  endif
  return
 else
  if(abs(beta)<=eps) then
   print *,'beta=',beta,eps
   stop 4537
  endif
  ga = -(f0*fn0im+alfa*dble(iim))/beta -ss0
  if(s0+ga>=0.0) then
   if (beta>0.0) then
    rijk=iim
   else
    rijk=ii
   endif
   ss=-ga+sqrt( (s0+ga)**2-2.*(rijk-r0)/beta )
  else
   if (beta>0.0) then
    rijk=ii
   else
    rijk=iim
   endif         
   if ( (s0+ga)**2 >= 2.*(rijk-r0)/beta ) then
    ss=-ga-sqrt( (s0+ga)**2-2.*(rijk-r0)/beta )
   else
    if (beta>0.0) then
     rijk=iim
    else
     rijk=ii
    endif
    ss=-ga+sqrt( (s0+ga)**2-2.*(rijk-r0)/beta )
   endif
  endif         
 endif
 return
else
 if (beta==0.0) then  ! b=0: stationary velocity fields !!!
  if (fn0i0>0.0) then
   if(abs((fn0i0-fn0im*(f0-1.0)))<EPS) stop 8842
   xi=1.0+(r0-dble(ii))*(fn0i0-fn0im)/(fn0i0-fn0im*(f0-1.0))        
   if(xi>0.0) then
    if(xi<EPS) stop 8843
    rijk=ii
    if(fn0im/=fn0i0) then
     if(abs(fn0i0-fn0im)<eps) stop 8843
     ss= s0 - 1.0/(fn0i0-fn0im)*log(xi)
    else
     if(abs(fn0i0-fn0im*(f0-1.0))<eps) stop 8844
     ss= s0 + (dble(ii)-r0)/(fn0i0-fn0im*(f0-1.0)) 
    endif
    if(ss<eps) stop 8845
    if (ss<=0.0) ss=UNDEF
    return
   endif
  endif
  if (fn0im<0.0) then
   !if(f0*fn0im<eps) stop 8847
   xi=1.0+(r0-dble(iim))*(fn0i0-fn0im)/(f0*fn0im)        
   if(xi>0.0) then
    rijk=iim
    if (fn0im/=fn0i0) then
     if(abs(fn0i0-fn0im)<eps) stop 8849
     ss= s0 - 1.0/(fn0i0-fn0im)*log(xi)
    else
     if(abs(f0*fn0im)<eps) stop 8851
     ss= s0 + (dble(iim)-r0)/(f0*fn0im)
    endif
      if(ss<eps) stop 8853
   if (ss<=0.0) ss=UNDEF
   endif
   endif
!   ssii = dsmin+dsmin
!   ssiim= dsmin+dsmin
        return
       else
        ssii = (f0*fn0im-alfa)/beta
        ssiim= f0*fn0im/beta
       endif
       if(abs(alfa)<eps) stop 8855
       ga = (r0-dble(iim)) + (beta*(s0-ss0-1.0/alfa)-f0*fn0im )/alfa
       xi0 = 0.0
! ''velocity'' at xi0
       xf1 = f0*fn0im-beta*(s0-ss0)-alfa*(r0-dble(iim)) 

! 2) ii direction (four velocity configurations at edge)
       if (ssii<=0.0 .or. ssii>=dsmin) then
        if (fn0i0>0.0 .or. fnmi0>0.0) then
! 2a) configuration: +++
         iconfig = 1
         xi00=ss0+dsmin-s0
        else
! 2b) configuration: --- no solution, try rijk=iim
         rijk=iim
         goto 3000
        endif 
       else
        if (fnmi0>0.0) then
! 2c) configuration: -- 0 ++
         iconfig = 3   
         xi00=ss0+dsmin-s0
        else
! 2d) configuration: ++ 0 --
         if (s0-ss0>=ssii) goto 3000
         iconfig = 4
         xi00=ss0+ssii-s0 
        endif 
       endif
       xf=(r0-dble(ii))+ga*(exp(-alfa*xi00)-1.0)-beta*xi00/alfa
       if (xf<0.0) goto 3000
       rijk=ii
       xib =xi00
       xibf=xf
       if (xf1>0.0) then
        if (xf1*(xi00-xi0)>-(r0-rijk)) then
         xi00=xi0
         xf=r0-rijk
 !         if (beta/=0.) goto 1000
        endif
       else
        if (ssiim>0.0 .or. ssiim<dsmin .and.fnmim>0.0) then
! check crossing at 3c) configuration: -- 0 ++
         xf2=ss0+ssiim-s0
         xf2=(r0-dble(iim))+ga*(exp(-alfa*xf2)-1.0)-beta*xf2/alfa
         if (xf2<0.0) goto 3000
        endif                
       endif   
       goto 1000

3000    continue
! 3) iim direction (four velocity configurations at edge)
       if (ssiim<=0.0 .or. ssiim>=dsmin) then
        if (fn0im>0.0 .or. fnmim>0.0) then
! 3a) configuration: +++ no solution
         return
        else
! 3b) configuration: ---
         iconfig = 6
         xi00=ss0+dsmin-s0
        endif 
       else
        if (fnmim>0.0) then
! 3c) configuration: -- 0 ++
         if (s0-ss0>=ssiim) return
         iconfig = 7   
         xi00=ss0+ssiim-s0 
        else
! 3d) configuration: ++ 0 --
         iconfig = 8
         xi00=ss0+dsmin-s0 
        endif 
       endif
       xf=(r0-dble(iim))+ga*(exp(-alfa*xi00)-1.0)-beta*xi00/alfa
       if (xf>0.0) then
!        rijk=iim  ! nyinlagt
!        print *,'nytt prov'
!        stop 4968
        return
       endif
       rijk=iim
       xib =xi00
       xibf=xf
       if (xf1<0.0) then
        if (xf1*(xi00-xi0)<-(r0-rijk)) then
         xi00=xi0
         xf=r0-rijk
        endif
       endif

! determining root and crossing time ss        
1000    continue
       xi00 = abs(alfa)*xi00
       xib  = abs(alfa)*xib
       if (looop==100) then
        xi=xi0-0.05
        do i=1,100
         xi = xi+0.05
         xf = (r0-rijk)+ga*(exp(-sign(1.0_dp,alfa)*xi)-1.0_dp)-beta*xi/(alfa*abs(alfa)) 
         xf1= (-beta/alfa-alfa*ga*exp(-sign(1.0_dp,alfa)*xi))/abs(alfa)         
!         print *, ' looop', xi, xf, xf1
        enddo
        looop=1000
       else

       xi = xi00 
       xf1= (-beta/alfa-alfa*ga*exp(-sign(1.0_dp,alfa)*xi))/abs(alfa) 
       xia = xi0
       xiaf= r0-rijk
       xerr = xxlim
       ssii = xxlim
100     continue
!        print *, ' looop', xi, xf, xf1
       looop=looop+1
       if(looop==90) then
 !       print *, ' looop stop3', xi0, xi00, xi, xf
        if (abs(xf)>100.*xxlim) looop=1000
!         if (xerr==UNDEF) looop=1000
        goto 400
       endif
       
       if(abs(xf1)<eps) stop 8859
       xin = xi - xf/xf1
     if ( abs(xf)>100.0 .or.(xf*xibf>0.0 .and. (xf1*xibf<0.0 .or.xin<xia)) ) then
        xib = xi
        xibf= xf
        xi = 0.5*(xi+xia)
       elseif (xf*xibf<0.0 .and. (xf1*xibf<0.0 .or.xin>xib) ) then
        xia = xi
        xiaf= xf
        xi = 0.5*(xi+xib)
       else
        xi = xin      
       endif
       xf = (r0-rijk)+ga*(exp(-sign(1.0_dp,alfa)*xi)-1.0_dp)-beta*xi/(alfa*abs(alfa)) 
       xf1= (-beta/alfa-alfa*ga*exp(-sign(1.0_dp,alfa)*xi) )/abs(alfa) 
       if (abs(xf)>=xerr) then
        goto 100
       elseif ( abs(xf)>xxlim*abs(xf1*(xi-xi0)).and. abs(xf)<ssii ) then
! as long as the error can be reduced and the accuracy in the computed
! time (ss) is not the same as the one for xi: continue the iterations
        xerr = UNDEF
        ssii = abs(xf)
        ssiim= xi
        goto 100
       endif
       if (xerr==UNDEF) xi = ssiim

       endif
400     continue
       ss= s0 + xi/abs(alfa)
       if (xi/=0.0) then
        xin=abs( xf/((xi-xi0)*xf1) )
!        accu(1)=accu(1)+1.0
!        accu(2)=accu(2)+xin
!        if (xin>xxlim.and.ntrac<=10) print *,ntrac,xin,xf1,ss-s0
       endif
       if (looop==1000.or.ss-s0<0.0) then
        if (ss-s0<0.0) print *,' ss-s0 is negative '
        print *, '0 time cross =', ss, s0,looop
         print *, alfa,fn0i0,fn0im,fnmi0,fnmim
         print *, iconfig,r0,rijk,xi0,xi00,xib,xf1,xibf,xi,xf
         ss=s0
 !       STOP
       endif

      endif 

return
end subroutine anil
!_______________________________________________________________________

function dawson1 (const,daw0,r0,xi0,xi)
USE mod_precdef

 IMPLICIT NONE
 REAL (DP)  :: dawson1,const,daw0,r0,xi0,xi,daw,dawson2

 daw=dawson2(xi)

 dawson1=r0*exp(xi0**2-xi**2)+const*(daw-exp(xi0**2-xi**2)*daw0)
 
 return

end function dawson1
!_______________________________________________________________________

function errfun (const,erf0,r0,xi0,xi)

USE mod_param
USE mod_precdef
 IMPLICIT NONE
 REAL (dp),   PARAMETER :: xilim=3.0_dp
 REAL (dp)   :: errfun,const,erf0,r0,xi0,xi,erf,s15aef,s15adf,hh0,hh
 INTEGER :: i

if (xi0>xilim) then ! complementary error function
 erf  = s15adf(xi)  
elseif(xi0<-xilim) then
 erf  = -s15adf(-xi)
else
! if(xi/=0.0_dp) stop 3956
 erf  = -s15aef(xi) ! error function
endif
if (abs(xi)>25.) then 
 hh0=1.0_dp
 hh=1.0_dp
 do i=1,10
  hh0=-0.5_dp*hh0*dble(i+i-1)/(xi0*xi0)
  hh=hh+hh0
 enddo
 erf=exp(xi**2-xi0**2)*hh/xi0
 hh0=1.0_dp
 hh=1.0_dp
 do i=1,10
  hh0=-0.5_dp*hh0*dble(i+i-1)/(xi*xi)
  hh=hh+hh0
 enddo       
 erf = 1.0_dp/sqrt(pi)*((hh/xi)-erf) 
else
 erf = exp(xi**2)*(erf-erf0)
endif 
errfun = r0*exp(xi**2-xi0**2) + const*erf 

!print *,'errfun', errfun, const,erf0,r0,xi0,xi
return

end function errfun


!______________________________________________________________________________
!function dawson2(x)
!
!! Returns Dawson's integral for any real x. 
!
!USE mod_precdef
!USE mod_param
!IMPLICIT NONE
!INTEGER, PARAMETER ::  NMAX=6  ! Denna ska kollas och testas med högre värden
!REAL (DP),  PARAMETER ::  H=0.4,A1=2./3.,A2=0.4,A3=2./7.
!
!REAL (DP)        :: dawson2,x,dd1,dd2,e1,e2,sfn0im,x2,xp,xx,pisqin
!REAL (DP), SAVE  :: c(NMAX)
!
!INTEGER*4       :: i,n0
!INTEGER*4, SAVE :: init=0
!
!!SAVE init,c
!!DATA init/0/ !  Flag is 0 if we need to initialize, else 1.
!
!init=0
!pisqin=1./sqrt(pi)
!
!if(init==0) then
! init=1
! do i=1,NMAX
!  c(i)=exp(-((2.0_dp*dble(i)-1.)*H)**2)
! enddo 
!endif
!
!if(abs(x)<0.2_dp) then    !  Use series expansion.
! x2=x**2
! dawson2=x*(1.0_dp-A1*x2*(1.-A2*x2*(1.-A3*x2)))
!else                     !  Use sampling theorem representation.
! xx=abs(x)
! n0=2*nint(0.5*xx/H)
! xp=xx-dble(n0)*H
! e1=exp(2.0_dp*xp*H)
! e2=e1**2
! dd1=dble(n0+1)
! dd2=dd1-2.
! sfn0im=0.
! if(abs(dd1)<eps .or. abs(dd2*e1)<eps) stop 3434
! do i=1,NMAX
!  sfn0im=sfn0im+c(i)*(e1/dd1+1.0_dp/(dd2*e1))
!  dd1=dd1+2.0_dp
!  dd2=dd2-2.0_dp
!  e1=e2*e1
! enddo 
! dawson2=pisqin*dsign(exp(-xp**2),x)*sfn0im 
!endif
!
!return
!end function dawson2


FUNCTION dawson2(XX) RESULT(fn_val)
!FUNCTION dawson2(XX) RESULT(fn_val)
USE mod_precdef
!USE mod_param
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This function program evaluates Dawson's integral,

!                       2  / x   2
!                     -x   |    t
!             F(x) = e     |   e    dt
!                          |
!                          / 0

!   for a real argfn0iment x.

!   The calling sequence for this function is

!                   Y=DAW(X)

!   The main computation uses rational Chebyshev approximations published
!   in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek and Thacher.
!   This transmortable program is patterned after the machine-dependent
!   FUNPACK program DDAW(X), but cannot match that version for efficiency or
!   accuracy.  This version uses rational approximations that are
!   theoretically accurate to about 19 significant decimal digits.
!   The accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-dependent
!   constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XINF   = largest positive machine nfn0imber
!   XMIN   = the smallest positive machine nfn0imber.
!   EPS    = smallest positive nfn0imber such that 1+eps > 1.
!            Approximately  beta**(-p), where beta is the machine radix
!            and p is the nfn0imber of significant base-beta digits in a
!            floating-point nfn0imber.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XMAX   = absolute argfn0iment beyond which DAW(X) underflows.
!            XMAX = min(0.5/xmin, xinf).
!   XSMALL = absolute argfn0iment below DAW(X)  may be represented
!            by X.  We recommend XSMALL = sqrt(eps).
!   XLARGE = argfn0iment beyond which DAW(X) may be represented by
!            1/(2x).  We recommend XLARGE = 1/sqrt(eps).

!     Approximate values for some important machines are

!                        beta  p     eps     xmin       xinf

!  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
!  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307

!                         XSMALL     XLARGE     XMAX

!  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
!  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
!  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
!  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
!                (D.P.)  3.73E-09   2.68E+08  1.70E+38
!   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307

!*******************************************************************

! Error Returns

!  The program returns 0.0 for |X| > XMAX.

! Intrinsic functions required are:

!     ABS


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------

IMPLICIT NONE
!INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
!INTEGER, PARAMETER  :: dp = 16

!REAL (8), INTENT(IN)  :: xx
REAL (dp), INTENT(IN)  :: xx
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: i
REAL (dp)  :: frac, sfn0imp, sfn0imq, w2, x, y
!----------------------------------------------------------------------
!  Mathematical constants.
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         six25 = 6.25_dp, one225 = 12.25_dp, two5 = 25.0_dp
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
!REAL (dp), PARAMETER  :: EPS=1.d-10
REAL (dp), PARAMETER  :: XSMALL = 1.05D-08, XLARGE = 9.49D+07,   &
                         XMAX = 2.24D+307
!REAL (dp), PARAMETER  :: XSMALL = 3.73D-09, XLARGE = 2.68E+08,   &
!                         XMAX = 1.70E+38
!REAL (dp), PARAMETER  :: XSMALL = sqrt(eps), XLARGE = 1./sqrt(eps),   &
!                         XMAX = min(0.5/xmin, xinf)
!                         XMAX = 1.70E+38 !min(0.5/xmin, xinf)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P1(10) = (/  &
        -2.69020398788704782410D-12, 4.18572065374337710778D-10,  &
        -1.34848304455939419963D-08, 9.28264872583444852976D-07,  &
        -1.23877783329049120592D-05, 4.07205792429155826266D-04,  &
        -2.84388121441008500446D-03, 4.70139022887204722217D-02,  &
        -1.38868086253931995101D-01, 1.00000000000000000004D+00 /)
REAL (dp), PARAMETER  :: Q1(10) = (/  &
         1.71257170854690554214D-10, 1.19266846372297253797D-08,  &
         4.32287827678631772231D-07, 1.03867633767414421898D-05,  &
         1.78910965284246249340D-04, 2.26061077235076703171D-03,  &
         2.07422774641447644725D-02, 1.32212955897210128811D-01,  &
         5.27798580412734677256D-01, 1.00000000000000000000D+00 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [2.5, 3.5)
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P2(10) = (/  &
        -1.70953804700855494930D+00,-3.79258977271042880786D+01,  &
         2.61935631268825992835D+01, 1.25808703738951251885D+01,  &
        -2.27571829525075891337D+01, 4.56604250725163310122D+00,  &
        -7.33080089896402870750D+00, 4.65842087940015295573D+01,  &
        -1.73717177843672791149D+01, 5.00260183622027967838D-01 /)
REAL (dp), PARAMETER  :: Q2(9) = (/  &
         1.82180093313514478378D+00, 1.10067081034515532891D+03,  &
        -7.08465686676573000364D+00, 4.53642111102577727153D+02,  &
         4.06209742218935689922D+01, 3.02890110610122663923D+02,  &
         1.70641269745236227356D+02, 9.51190923960381458747D+02,  &
         2.06522691539642105009D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [3.5, 5.0]
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P3(10) = (/  &
        -4.55169503255094815112D+00,-1.86647123338493852582D+01,  &
        -7.36315669126830526754D+00,-6.68407240337696756838D+01,  &
         4.84507265081491452130D+01, 2.69790586735467649969D+01,  &
        -3.35044149820592449072D+01, 7.50964459838919612289D+00,  &
        -1.48432341823343965307D+00, 4.99999810924858824981D-01 /)
REAL (dp), PARAMETER  :: Q3(9) = (/  &
         4.47820908025971749852D+01, 9.98607198039452081913D+01,  &
         1.40238373126149385228D+01, 3.48817758822286353588D+03,  &
        -9.18871385293215873406D+00, 1.24018500009917163023D+03,  &
        -6.88024952504512254535D+01,-2.31251575385145143070D+00,  &
         2.50041492369922381761D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  |x| > 5.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P4(10) = (/  &
        -8.11753647558432685797D+00,-3.84043882477454453430D+01,  &
        -2.23787669028751886675D+01,-2.88301992467056105854D+01,  &
        -5.99085540418222002197D+00,-1.13867365736066102577D+01,  &
        -6.52828727526980741590D+00,-4.50002293000355585708D+00,  &
        -2.50000000088955834952D+00, 5.00000000000000488400D-01 /)
REAL (dp), PARAMETER  :: Q4(9) = (/  &
         2.69382300417238816428D+02, 5.04198958742465752861D+01,  &
         6.11539671480115846173D+01, 2.08210246935564547889D+02,  &
         1.97325365692316183531D+01,-1.22097010558934838708D+01,  &
        -6.99732735041547247161D+00,-2.49999970104184464568D+00,  &
         7.49999999999027092188D-01 /)
         
!----------------------------------------------------------------------
x = xx
IF (ABS(x) > xlarge) THEN
  IF (ABS(x) <= xmax) THEN
    fn_val = half / x
  ELSE
    fn_val = zero
  END IF
ELSE IF (ABS(x) < xsmall) THEN
  fn_val = x
ELSE
  y = x * x
  IF (y < six25) THEN
!----------------------------------------------------------------------
!  ABS(X) .LT. 2.5
!----------------------------------------------------------------------
    sfn0imp = p1(1)
    sfn0imq = q1(1)
    DO  i = 2, 10
      sfn0imp = sfn0imp * y + p1(i)
      sfn0imq = sfn0imq * y + q1(i)
    END DO
    fn_val = x * sfn0imp / sfn0imq
  ELSE IF (y < one225) THEN
!----------------------------------------------------------------------
!  2.5 .LE. ABS(X) .LT. 3.5
!----------------------------------------------------------------------
    frac = zero
    DO  i = 1, 9
      frac = q2(i) / (p2(i)+y+frac)
    END DO
    fn_val = (p2(10)+frac) / x
  ELSE IF (y < two5) THEN
!----------------------------------------------------------------------
!  3.5 .LE. ABS(X) .LT. 5.0
!---------------------------------------------------------------------
    frac = zero
    DO  i = 1, 9
      frac = q3(i) / (p3(i)+y+frac)
    END DO
    fn_val = (p3(10)+frac) / x
  ELSE
!----------------------------------------------------------------------
!  5.0 .LE. ABS(X) .LE. XLARGE
!------------------------------------------------------------------
    w2 = one / x / x
    frac = zero
    DO  i = 1, 9
      frac = q4(i) / (p4(i)+y+frac)
    END DO
    frac = p4(10) + frac
    fn_val = (half + half*w2*frac) / x
  END IF
END IF
RETURN
!---------- Last line of DAW ----------
END FUNCTION dawson2
!END FUNCTION dawson3





!______________________________________________________________________________

!FUNCTION erf(x)
FUNCTION s15aef(x)
USE mod_precdef

! USES gammp
!Returns the error function erf(x).

IMPLICIT NONE
REAL (DP) ::  s15aef,x,gammp
!REAL (dp)  :: s15aef,x,gammp

if(x<0.0_dp)then
 s15aef=-gammp(0.5_dp,x**2)
else
 s15aef= gammp(0.5_dp,x**2)
endif

return
end function s15aef

!______________________________________________________________________________


FUNCTION gammp(a,x)
USE mod_precdef

! USES gcf,gser
! Returns the incomplete gamma function P(a,x).

IMPLICIT NONE
REAL (DP)  :: a,gammp,x,gammcf,gamser,gln

if(x<0.0_dp .or. a<=0.0_dp) print *, 'bad argfn0iments in gammp'

if(x<a+1.0_dp) then ! Use the series representation.
 call gser(gamser,a,x,gln)
 gammp=gamser
else ! Use the continued fraction representation
 call gcf(gammcf,a,x,gln)
 gammp=1.0_dp-gammcf ! and take its complement.
endif

return
end function gammp

!______________________________________________________________________________

SUBROUTINE gcf(gammcf,a,x,gln)
USE mod_precdef
USE mod_param
IMPLICIT NONE

! USES gammln
!Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation as gammcf. Also returns as gln.
!Parameters: ITMAX is the maximfn0im allowed nfn0imber of iterations; EPS is the relative accuracy;
!FPMIN is a nfn0imber near the smallest representable floating-point nfn0imber.

INTEGER, PARAMETER ::  ITMAX=1000
REAL (DP),  PARAMETER ::  FPMIN=1.d-30
REAL (DP)  :: a,gammcf,gln,x,an,b,c,d,del,h,gammln
INTEGER :: i

gln=gammln(a)
b=x+1.0_dp-a   ! Set up for evaluating continued fraction by modied
c=1.0_dp/FPMIN !Lentz's method (x5.2) with b0 = 0.
d=1.0_dp/b
h=d

do i=1,ITMAX  ! Iterate to convergence.
 an=-dble(i)*(dble(i)-a)
 b=b+2.0_dp
 d=an*d+b
 if(abs(d)<FPMIN) d=FPMIN
 c=b+an/c
 if(abs(c)<FPMIN) c=FPMIN
 d=1.0_dp/d
 del=d*c
 h=h*del
 if(abs(del-1.0_dp)<EPS) goto 1
enddo
 
!print *, 'a too large, ITMAX too small in gcf',abs(del-1.0_dp),EPS,del
!print *,'cd=',an
!stop 8780

1 gammcf=exp(-x+a*log(x)-gln)*h  ! Put factors in front.

return
end subroutine gcf

!______________________________________________________________________________

FUNCTION gammln(xx)
USE mod_precdef
IMPLICIT NONE
!Returns the value ln[gamma(xx)] for xx > 0.
!Internal arithmetic will be done in double precision, 
!a nicety that you can omit if five-figure accuracy is good enough.

REAL (DP)  :: gammln,xx,ser,stp,tmp,x,y,cof(6)
INTEGER j

SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
             24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
              -.5395239384953d-5,2.5066282746310005d0/

x=xx
y=x
tmp=x+5.5_dp
tmp=(x+0.5_dp)*log(tmp)-tmp
ser=1.000000000190015

do j=1,6
 y=y+1.0_dp
 ser=ser+cof(j)/y
enddo 

if(x<1.d-10) stop 8890
gammln=tmp+log(stp*ser/x)
 
return

end function gammln

!______________________________________________________________________________


SUBROUTINE gser(gamser,a,x,gln)

! USES gammln
!Returns the incomplete gamma function P(a; x) evaluated by its series 
! representation as gamser. Also returns  as gln.

USE mod_precdef
USE mod_param
IMPLICIT NONE

INTEGER, PARAMETER ::  ITMAX=1000

REAL (DP)  :: a,gamser,gln,x,ap,del,sfn0im,gammln
INTEGER n

gln=gammln(a)

if(x<=0.0_dp) then
 if(x<0.0_dp) print *,'x < 0 in gser'
 gamser=0.0_dp
 return
endif

ap=a
sfn0im=1.0_dp/a
del=sfn0im

do n=1,ITMAX
 ap=ap+1.0_dp
 del=del*x/ap
 sfn0im=sfn0im+del
 if(abs(del)<abs(sfn0im)*EPS) goto 1
enddo
 
print *, 'a too large, ITMAX too small in gser'
stop 8860

1 gamser=sfn0im*exp(-x+a*log(x)-gln)

return
end subroutine gser

!______________________________________________________________________________

FUNCTION s15adf(x) !FUNCTION erfc(x) 

! USES gammp,gammq 
! Returns thecomplementaryerror functionerfc(x). 

USE mod_precdef
IMPLICIT NONE
REAL (DP)  :: s15adf,x,gammp,gammq

if(x<0.0_dp) then 
 s15adf=1.0_dp+gammp(0.5_dp,x**2) 
else 
 s15adf=       gammq(0.5_dp,x**2) 
endif 

return 
end function s15adf

!_______________________________________________________________________________

FUNCTION gammq(a,x) 

USE mod_precdef
IMPLICIT NONE
REAL (DP)  :: a,gammq,x,gammcf,gamser,gln 

if(x<0.0_dp .or. a<=0.0_dp) print *, 'bad argfn0iments in gammq' 

if(x<a+1.0_dp) then 
 call gser(gamser,a,x,gln) 
 gammq=1.0_dp-gamser ! and take its complement. 
else ! Use the continued fraction representation. 
 call gcf(gammcf,a,x,gln) 
 gammq=gammcf 
endif 

return 
end function gammq

!_______________________________________________________________________________
#endif
