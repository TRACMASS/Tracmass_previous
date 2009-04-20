
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
	if(ib>1 .AND. ib<=IMT .AND. jb>1 .AND. jb<=JMT .AND. kb<=KM) then
		tryAgain = .TRUE.
	end if
	if(.NOT. tryAgain) then
		write(*,*)"========"
			write(*,*)"Particle outside model area. No diffusion added."
			write(*,*)"========"
	end if
	
	do while(tryAgain)
		! displace particle
		CALL displacement(xd, yd, zd, dt)
		! Convert to Earth coordinates(?)
		xd = xd/(dx*cst(jb)*deg)
		yd = yd/(dy*deg)
		zd = zd/dz(kb)
		! Update position temporarily
		tmpX = x1 + xd
		tmpY = y1 + yd
		tmpZ = z1 + zd
		! Update box number temporarily
		tmpi = int(tmpX) + 1
		tmpj = int(tmpY) + 1
		tmpk = KM - int(tmpZ)

		! Check that at the particle level, there is at least one adjacent
		! open ocean velocity point.
		if(tmpi>1 .AND. tmpi<=IMT .AND. tmpj>1 .AND. tmpj<=JMT .AND. tmpk>0) then
		! check that column is deep enough 	
			if(tmpk<=KM) then		! Borde vi ha med kmu också? För att kunna använda med andra projekt än tes???
				tryAgain = .FALSE.
			end if
		end if
		
		! Check if particle is on an open boundary
		if(tmpi==1 .AND. tmpj>=1 .AND. tmpj<=JMT .AND. tmpk>=1 .AND. &
		tmpk<=kmt(tmpi, tmpj)) then
			tryAgain = .FALSE.
		end if
		
		! Slå ihop de två ovanstående ifsatserna? kmu med i första men inte andra...
		
		! If tryAgain is still true, the particle is outside model area. The
		! displacement is not saved, but we make a new try to displace.
		
		! "Infinite loop?"
		if(itno>=10000 .AND. tryAgain) then
			tryAgain = .FALSE.
			write(*,*)"========"
			write(*,*)"Particle stuck in infinite diffusion loop. No diffusion added."
			write(*,*)"========"
		end if
		
	enddo
	
	! Update return position
	x1 = tmpX
	y1 = tmpY
	z1 = tmpZ
	ib = tmpi
	jb = tmpj
	kb = int(tmpZ) + 1
	
	! Check the vertical velocity if particle on boundary
	if(int(z1) == z1) then
		! Om partikeln är på väg nedåt (mot en låda med lägre index), ändra kb
		! till det lägre lådindexet.
	endif
	
	
END SUBROUTINE diffuse

!===============================================================================
! Calculate a random displacement
! (sqrt(-4Ah*dt*log(1-q1))*cos(2PIq2),
!  sqrt(-4Ah*dt*log(1-q1))*sin(2PIq2),
!  sqrt(-4Av*dt*log(1-q3))*cos(2PIq4))
! where
! Av=0.0001, Ah=200, dt is the model time step and q1,q2,q3,q4 are random
! numbers between 0 and 1.
!
! Arguments :
! xd, yd, zd : Variables in which the displacement will be stored
! dt: Model time step
!===============================================================================
SUBROUTINE displacement(xd, yd, zd, dt) 
	USE mod_precdef
	
	REAL						:: Ah, Av, R
	REAL						:: q1, q2, q3, q4
	REAL, INTENT(OUT)			:: xd, yd, zd
	REAL (KIND=DP), INTENT(IN)	:: dt
	REAL, PARAMETER				:: PI = 3.14159265358979323846
		
	Av = 1.0d-4
	Ah = 2.0d2
	Av = 1.0d-5
	Ah = 2.0d1
	Av = 1.0d-6
	Ah = 2.0d0
	Av = 1.0d-7
	Ah = 2.0d-1

	q1 = rand()		! Så?
	q2 = rand()
	q3 = rand()
	q4 = rand()

	! Horizontal displacements
	R = sqrt(-4*Ah*dt*log(1-q1))
	xd = R * cos(2*PI*q2)
	yd = R * sin(2*PI*q2)
	
	! Vertical displacement
	R = sqrt(-4*Av*dt*log(1-q3))
	zd = R*cos(2*PI*q4)
	
END SUBROUTINE displacement
