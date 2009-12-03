#ifdef diffusion

!============================================================================
! Add a small displacement to a particle s.t. it is still in the model area.
! 
! Arguments
! x1, y1, z1 : Current position of the particle. Are updated by the subroutine.
! ib, jb, kb : Current box for the particle. Are updated by the subroutine.
! dt : The model time step
!
! FELKODER SÄTTS VID BEHOV, FÖRKLARA
! bla
!============================================================================
SUBROUTINE diffuse(x1, y1, z1, ib, jb, kb, dt)
	USE mod_coord
	USE mod_grid
	USE mod_param
	USE mod_precdef
	
	implicit none
 	INTEGER						:: ib,jb,kb		! Box indices
 	INTEGER						:: itno			! Number of iterations
	REAL						:: xd, yd, zd	! Displacement
	REAL						:: tmpX, tmpY, tmpZ		! Temporal position
	INTEGER						:: tmpi, tmpj, tmpk		! Temporal box indices
	REAL (KIND=DP), INTENT(OUT)	:: x1, y1, z1			! Final position
	REAL (KIND=DP), INTENT(IN)	:: dt			! Time step
	LOGICAL						:: tryAgain	! Tells whether to continue displace

	tryAgain = .FALSE.
	
	! Is particle within model area?
	if(ib>=1 .AND. ib<=IMT .AND. jb>=1 .AND. jb<=JMT .AND. KM+1-kmt(ib,jb)<=kb .AND. kb>=1 ) then
		tryAgain = .TRUE.
    else
     stop 86567
	end if
!	 print *,'ib,jb,kb',ib,jb,kb,kmt(ib,jb),x1,y1,z1
	if(.NOT. tryAgain) then
	 print *,'ib,jb,kb',kmt(ib,jb),ib,jb,kb,x1,y1,z1
		write(*,*)"========"
			write(*,*)"Particle outside model area. No diffusion added."
			write(*,*)"========"
			stop 3957
	end if
	
	itno=0
	do while(tryAgain)
	    itno=itno+1
		! find random displacement 
		CALL displacement(xd, yd, zd, ib, jb, kb, dt)
		! Convert displacement from meters to model coordinates
#ifdef rco
		xd = xd/(dx*cst(jb)*deg)  
		yd = yd/(dy*deg)          
		zd = zd/dz(kb)            !should be replaced for bottom box and include ssh
#elif orc
		xd = xd/dxv(ib,jb)  
		yd = yd/dyu(ib,jb)          
		zd = zd/dz(kb)            !should be replaced for bottom box and include ssh
#endif
		! Update position temporarily
		tmpX = x1 + xd
		tmpY = y1 + yd
		tmpZ = z1 + zd
		! Update box number temporarily
		tmpi = int(tmpX) + 1
		tmpj = int(tmpY) + 1
		tmpk = int(tmpZ) + 1

! Check that at the particle level, there is at least one adjacent open ocean velocity point.
#ifdef orc
! cyclic conditions for global Earth CCMs like ORCA or IFS
		if(tmpi>IMT) then 
		 tmpi=tmpi-IMT
		elseif(tmpi<1) then
		 tmpi=tmpi+IMT
		endif
        if(tmpX>dble(IMT)) then 
		 tmpX=tmpX-dble(IMT)
		elseif(tmpX<0.d0) then
		 tmpX=tmpX+dble(IMT)
		endif
		if(tmpj>JMT) tmpj=JMT  ! stop 34956 ! north fold for orca grids
		if(tmpY>dble(JMT)) tmpY=dble(JMT)-0.1d0  ! stop 34956 ! north fold for orca grids
#elif rco
		! Check if particle is on an open boundary
		if(tmpi==1 .AND. tmpj>=1 .AND. tmpj<=JMT .AND. KM+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
			tryAgain = .FALSE.
		end if
#endif
		! check that column is deep enough 	
		if( 1<=tmpi .AND. tmpi<=IMT .AND. 1<=tmpj .AND. tmpj<=JMT .AND. KM+1-kmt(tmpi,tmpj)<=tmpk .AND. tmpk>=1 ) then
            tryAgain = .FALSE. ! if false then a new position for the particle has been found and we exit the loop
!        print *,'hittat'
		end if 
		
		! If tryAgain is still true, the particle is outside model area. The
		! displacement is not saved, but we make a new try to displace.
		
		! "Infinite loop?"
		if(itno>=100000 .AND. tryAgain) then
			tryAgain = .FALSE.
			write(*,*)"Particle stuck in infinite diffusion loop. No diffusion added.",ib,jb,kb
			tmpX=x1 ; tmpY=y1 ; tmpZ=z1
			tmpi=ib ; tmpj=jb ; tmpk=kb 
		end if
		
	enddo
	
	! Update return position
	x1 = tmpX
	y1 = tmpY
	ib = tmpi
	jb = tmpj
#ifndef twodim   
	z1 = tmpZ
	kb = tmpk
#endif	

!print *,'slut',itno,kmt(ib,jb),ib,jb,kb,x1,y1,z1
	
END SUBROUTINE diffuse

!===============================================================================
! Calculate a random displacement
! (sqrt(-4Ah*dt*log(1-q1))*cos(2PIq2),
!  sqrt(-4Ah*dt*log(1-q1))*sin(2PIq2),
!  sqrt(-4Av*dt*log(1-q3))*cos(2PIq4))
! where Av and Ah are set in run.in, dt is the model time step and q1,q2,q3,q4 are random
! numbers between 0 and 1.
!
! Arguments :
! xd, yd, zd : Variables in which the displacement will be stored
! dt: Model time step
!===============================================================================
SUBROUTINE displacement(xd, yd, zd, ib, jb, kb, dt) 
	USE mod_precdef
	USE mod_diffusion
	USE mod_param
#ifdef anisodiffusion
	USE mod_grid
#endif
	IMPLICIT NONE
	
	REAL						:: q1, q2, q3, q4, R
	REAL, INTENT(OUT)			:: xd, yd, zd
	REAL (KIND=DP), INTENT(IN)	:: dt
!	REAL, PARAMETER				:: PI = 3.14159265358979323846
 	INTEGER						:: ib,jb,kb		! Box indices
#ifdef anisodiffusion 	
	REAL*8						:: Rx, Ry, grdx, grdy, grad, theta, xx, yy
 	INTEGER						:: ip,im,jp,jm
#endif
		
! random generated numbers between 0 and 1
	q1 = rand()	
	q2 = rand()
	q3 = rand()
	q4 = rand()

! Horizontal displacements in meters
	R = sqrt(-4*Ah*dt*log(1-q1))
	xd = R * cos(2*PI*q2)
	yd = R * sin(2*PI*q2)
		
	! Vertical displacement in meters
#ifndef twodim   
	R = sqrt(-4*Av*dt*log(1-q3))
	zd = R*cos(2*PI*q4)
#else
    zd = 0.
#endif

#ifdef anisodiffusion
!_______________________________________________________________________________
! The diffusion is here set on an ellipse instead of a circle
! so that the diffusion is higher along the isobaths and weaker
! in the perpendicular direction (up/downhill)

ip=ib+1
if(ip.eq.IMT+1) ip=1
im=ib-1
if(im.eq.0) im=IMT
jp=jb+1
if(jp.eq.JMT+1) jp=JMT
jm=jb-1
if(jm.eq.0) jm=1

grdx=float(kmt(ip,jb)-kmt(im,jb)) ! zonal      depth gradient (zw should be used)
grdy=float(kmt(ib,jp)-kmt(ib,jm)) ! meridional depth gradient
grad=dsqrt(grdx**2+grdy**2)       ! total depth gradient

! angle between the eastward direction and 
! the constant depth direction (=0 if only meridional slope)
if(grad .eq.0.) then
 theta=0
else
! theta=dasin(grdx/grad)
 theta=dasin(grdx/grad)
endif

! elliptic horizontal distribution of the diffusion
grad=dabs(grad)+1.d0
!grad=amin1(10.d0,grad) ! gives slightly higher relative dispersion
grad=amin1(5.d0,grad)
xx=xd*grad
yy=yd/grad

! coordinate transformation to put the diffusion on an 
! ellipse with the maxium diffusion along the isobaths
xd= xx*dcos(theta)-yy*dsin(theta)
yd=-xx*dsin(theta)+yy*dcos(theta)

!print *,ib,jb,grdx,grdy,grad,theta*180./pi,xd,xx,yd,yy


!if(jb.gt.400) stop 3096
!_______________________________________________________________________________
#endif

END SUBROUTINE displacement

#endif

