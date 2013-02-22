MODULE mod_seed
!!------------------------------------------------------------------------------!!
!!
!!       MODULE: mod_seed
!!
!!          Defines variables and matrices necessary for initializing
!!          particles.
!!
!!          Populates matrices trj and nrj, containing information of
!!          all particles. See init_par for allocation of them. 
!!             trj        - Dimension NTRACMAX x NTRJ
!!             nrj        - Dimension NTRACMAX x NNRJ
!!
!!
!!------------------------------------------------------------------------------
   USE mod_param
   USE mod_time
   USE mod_grid
   USE mod_buoyancy
   USE mod_vel
   USE mod_traj
   
   IMPLICIT NONE
  
   INTEGER                                    ::    nff,  isec,  idir,         &
                                               &   nqua,   num,   nsd, nsdMax, &
                                               & nsdTim,                       &
                                               & seedPos, seedTime, seedType,  &
                                               & seedAll, varSeedFile,         &
                                               & ist1, ist2, jst1, jst2,       &
                                               & kst1, kst2, tst1, tst2,       &
                                               & iist, ijst, ikst, jsd, jst,   &
                                               &  ijt,  ikt,  jjt, jkt, ntractot
   INTEGER*8                                  :: itim
   INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: seed_ijk, seed_set
   INTEGER, ALLOCATABLE, DIMENSION(:)         :: seed_tim
   REAL*8 , ALLOCATABLE, DIMENSION(:,:)       :: seed_xyz
   CHARACTER(LEN=200)                         :: seedDir, seedFile, timeFile
!!
!!------------------------------------------------------------------------------
  
CONTAINS


   SUBROUTINE seed (tt,ts)

!!------------------------------------------------------------------------------
      
      INTEGER                                  :: errCode,         &
                                                & ib, jb, kb, ibm, &
                                                &  i,  j,  k,   l, &
                                                &  m, ntrac
      REAL*8                                   :: temp,salt,dens
      REAL*8                                   :: tt, ts,     &
                                                & x1, y1, z1, &
                                                & vol, subvol

!!------------------------------------------------------------------------------

      ! --------------------------------------------
      ! --- Check if ntime is in vector seed_tim ---
      ! --------------------------------------------
      IF (seedTime == 2) THEN
         findTime: DO jsd=1,nsdTim
            IF (seed_tim(jsd) == ntime) THEN
               IF (seedAll == 1) THEN
                  itim = seed_tim(jsd)
               END IF
               EXIT findTime
            ELSE IF (seed_tim(jsd) /= ntime .AND. jsd == nsdTim) THEN
               RETURN
            END IF
         END DO findTime
      END IF
      
      ! ---------------------------------------
      ! --- Loop over the seed size, nsdMax ---
      !----------------------------------------
      startLoop: DO jsd=1,nsdMax
         iist  = seed_ijk (jsd,1)
         ijst  = seed_ijk (jsd,2)
         ikst  = seed_ijk (jsd,3)
         isec  = seed_set (jsd,1)
         idir  = seed_set (jsd,2)
         if (iist <= 1)   cycle startLoop
         if (ijst <= 1)   cycle startLoop
         if (iist >= imt) cycle startLoop
         if (ijst >= jmt) cycle startLoop
         if (ikst == km+1)  ikst = km
         if (ikst > km)   cycle startLoop
         if (kmt(iist,ijst) == 0) cycle startLoop



         IF (seedTime == 2 .AND. seedAll == 2) THEN
            itim  = seed_tim (jsd)
         END IF
         
#if defined baltix || defined rco
         ! -------------------------------------------------
         ! --- Test if it is time to launch the particle ---
         ! -------------------------------------------------
         IF ( (seedTime == 1 .AND. (ntime < tst1 .OR. ntime > tst2)) .OR. &
         &    (seedTime == 2 .AND. ntime /= itim) ) THEN
            CYCLE startLoop
         ELSE IF (seedTime /= 1 .AND. seedTime /= 2) THEN
            PRINT*,'timeStart =',seedTime,' is not a valid configuration!'
            STOP
         END IF
#endif         
         vol = 0    
         ib  = iist
         ibm = ib-1
         IF (ibm == 0) THEN
            ibm = IMT
         END IF
         jb  = ijst
         kb  = ikst
         
         ! -----------------------------------------------------------
         ! --- Determine the volume/mass flux through the grid box ---
         ! -----------------------------------------------------------
         SELECT CASE (isec)
         
            CASE (1)  ! Through eastern meridional-vertical surface
               vol = uflux (iist,ijst,ikst,1)
         
            CASE (2)  ! Through northern zonal-vertical surface
               vol = vflux (iist,ijst,ikst,1)
         
            CASE (3)  ! Through upper zonal-meridional surface
               CALL vertvel (1.d0,ib,ibm,jb,kb)
#if defined full_wflux || defined explicit_w
               vol=wflux(ib,jb,kb,1)
#elif twodim
               vol=1.
#else 
               vol=wflux(kb,1)
#endif
         
            CASE (4 ,5)   ! Total volume/mass transport through grid box
               IF (KM+1-kmt(iist,ijst) > kb) THEN
                  CYCLE startLoop
               ELSE
                  vol = uflux (ib, jb, kb, 1) + uflux (ibm, jb  , kb, 1) + & 
                  &     vflux (ib, jb, kb, 1) + vflux (ib , jb-1, kb, 1)
               ENDIF
               IF (vol == 0.d0) cycle startLoop
         
         END SELECT
      
         ! If the particle is forced to move in positive/negative direction
         IF ( (idir*ff*vol <= 0.d0 .AND. idir /= 0) .OR. (vol == 0.) ) THEN
            CYCLE startLoop
         ENDIF
      
         ! Volume/mass transport needs to be positive   
         vol = ABS (vol)
      
        
         ! Calculate transport of each individual trajectory
         IF (nqua == 3 .OR. isec > 4) THEN
#ifdef zgrid3Dt
            vol = dzt(ib,jb,kb,1)
#elif  zgrid3D
            vol = dzt(ib,jb,kb)
#elif  zgrid1D
            vol = dz(kb)
#endif /*zgrid*/
#ifdef varbottombox
            IF (kb == KM+1-kmt(ib,jb)) THEN
               vol = dztb (ib,jb,1)
            END IF
#endif /*varbottombox*/
#ifdef freesurface
            IF (kb == KM) THEN
               vol = vol+hs(ib,jb,1)
            END IF
            vol = vol*dxdy(ib,jb)
#endif /*freesurface*/
         END IF
          
         ! Number of trajectories for box (iist,ijst,ikst)
         SELECT CASE (nqua)
            CASE (1)
               num = partQuant
            CASE (2)
               num = vol/partQuant
            CASE (3)
               num = vol/partQuant
            CASE (4)
               num = partQuant
         END SELECT
         
         IF (num == 0 .AND. nqua /= 4) THEN
            num=1
         END IF
     
         ijt    = NINT (SQRT (FLOAT(num)) ) 
         ikt    = NINT (FLOAT (num) / FLOAT (ijt))
         subvol = vol / DBLE (ijt*ikt)
     
         IF (subvol == 0.d0) THEN
            PRINT*,' Transport of particle is zero!!!'
            PRINT*,' vol =',vol
            PRINT*,' subvol =',subvol
            STOP
         ENDIF
         
         ! --------------------------------------------------
         ! --- Determine start position for each particle ---
         ! --------------------------------------------------
         ijjLoop: DO jjt=1,ijt
            kkkLoop: DO jkt=1,ikt          
            
!            IF ( ib /= iist ) PRINT*,iist,ib
!            IF ( jb /= ijst ) PRINT*,ijst,jb
!            IF ( kb /= ikst ) PRINT*,ikst,kb
            ib = iist
            jb = ijst
            kb = ikst

            SELECT CASE (isec)
               CASE (1)   ! Meridional-vertical section
                  
                  x1 = DBLE (ib) 
                  y1 = DBLE (jb-1) + (DBLE (jjt) - 0.5d0) / DBLE (ijt) 
                  z1 = DBLE (kb-1) + (DBLE (jkt) - 0.5d0) / DBLE (ikt)
                  
                  IF (idir == 1) THEN
                     ib = iist+1
                  ELSE IF (idir == -1) THEN
                     ib=iist 
                  END IF
                  
               CASE (2)   ! Zonal-vertical section
                  
                  x1 = DBLE (ibm)  + (DBLE (jjt) - 0.5d0) / DBLE (ijt)
                  y1 = DBLE (jb)
                  z1 = DBLE (kb-1) + (DBLE (jkt) - 0.5d0) / DBLE (ikt) 
                  
                  IF (idir == 1) THEN
                     jb = ijst+1
                  ELSE IF (idir == -1) THEN
                     jb = ijst
                  END IF 
              
               CASE (3)   ! Horizontal section                  
                  x1 = DBLE (ibm)  + (DBLE (jjt) - 0.5d0) / DBLE (ijt)
                  y1 = DBLE (jb-1) + (DBLE (jkt) - 0.5d0) / DBLE (ikt) 
                  z1 = DBLE (kb)
                  
                  IF (idir == 1) THEN
                     kb = ikst+1
                  ELSE IF (idir == -1) THEN
                     kb = ikst
                  END IF
                  
               CASE (4)   ! Spread even inside box                  
                  x1 = DBLE (ibm)  + 0.25d0 * (DBLE(jjt) - 0.5d0) / DBLE(ijt)
                  y1 = DBLE (jb-1) + 0.25d0 * (DBLE(jkt) - 0.5d0) / DBLE(ikt)
                  z1 = DBLE (kb-1) + 0.5d0
                                 
               CASE (5)                  
                  x1 = seed_xyz (jsd,1)
                  y1 = seed_xyz (jsd,2) 
                  z1 = seed_xyz (jsd,3)
               
               END SELECT
           
           ! ------------------------------------------------------
           ! --- Check properties of water mass at initial time ---
           ! ------------------------------------------------------ 

#ifdef tempsalt 
               CALL interp (ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
               IF (temp < tmin0 .OR. temp > tmax0 .OR. &
               &   salt < smin0 .OR. salt > smax0 .OR. &
               &   dens < rmin0 .OR. dens > rmax0      ) THEN
                  CYCLE kkkLoop 
               END IF
#endif /*tempsalt*/

               ! Update trajectory numbers
               ntractot = ntractot+1
               ntrac = ntractot
           
#ifdef select
               ! Selects only one single trajectory
               if(ntrac.ne.57562) then 
                  nrj(ntrac,6)=1
                  cycle kkkLoop
               endif
#endif /*select*/
           
               ! ts - time, fractions of ints
               ! tt - time [s] rel to start
!               ts = ff * DBLE (ints-intstep) / tstep
               ts = DBLE (ints-1) 
               tt = ts * tseas
               
               ! ------------------------------------------------------------
               ! --- Put the new trajectory into the matrices trj and nrj ---
               ! ------------------------------------------------------------
               trj(ntrac,1:7) = [ x1, y1, z1, tt,    subvol, 0.d0, tt ]
               nrj(ntrac,1:5) = [ ib, jb, kb,  0, IDINT(ts)]
               nrj(ntrac,7)=1
           
            END DO kkkLoop
         END DO ijjLoop   
            
      END DO startLoop

   END SUBROUTINE seed

END MODULE mod_seed

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

