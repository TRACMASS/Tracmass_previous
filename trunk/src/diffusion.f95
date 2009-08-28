#ifdef diffusion

!===============================================================================
! Add a small displacement to a particle s.t. it is still in the model area.
! 
! Arguments
! x1, y1, z1 : Current position of the particle. Are updated by the subroutine.
! ib, jb, kb : Current box for the particle. Are updated by the subroutine.
! dt : The model time step
!
! FELKODER SÄTTS VID BEHOV, FÖRKLARA
!
!===============================================================================
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
		CALL displacement(xd, yd, zd, dt)
		! Convert to model coordinates
#ifdef rco
		xd = xd/(dx*cst(jb)*deg)  
		yd = yd/(dy*deg)          
		zd = zd/dz(kb)            !should be replaced for bottom box and include ssh
#elif orc
!if(xd.lt.0. .or. yd.lt.0.) print *,'------xdyd=',ib,jb,kb,xd,yd
!print *,'------xdyd=',xd,yd
		xd = xd/dxv(ib,jb)  
		yd = yd/dyu(ib,jb)          
!print *,'ibjbkb=',xd,yd
		zd = zd/dz(kb)            !should be replaced for bottom box and include ssh
#endif
		! Update position temporarily
		tmpX = x1 + xd
		tmpY = y1 + yd
!print *,'tmpX=',x1,xd,tmpX,y1,yd,tmpY
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
SUBROUTINE displacement(xd, yd, zd, dt) 
	USE mod_precdef
	USE mod_diffusion
	
!	REAL						:: Ah, Av, R
	REAL						:: q1, q2, q3, q4
	REAL, INTENT(OUT)			:: xd, yd, zd
	REAL (KIND=DP), INTENT(IN)	:: dt
	REAL, PARAMETER				:: PI = 3.14159265358979323846
		
! this is set in run.in now
!	Av = 0.0d0
!	Ah = 1.0d0  ! dgl0
!	Ah = 1.0d-2  ! dgl1
!	Ah = 1.0d2  ! dgl2

	q1 = rand()	
	q2 = rand()
	q3 = rand()
	q4 = rand()

	! Horizontal displacements
	R = sqrt(-4*Ah*dt*log(1-q1))
	xd = R * cos(2*PI*q2)
	yd = R * sin(2*PI*q2)
	
	! Vertical displacement
#ifndef twodim   
	R = sqrt(-4*Av*dt*log(1-q3))
	zd = R*cos(2*PI*q4)
#else
    zd = 0.
#endif
	
END SUBROUTINE displacement

#endif

