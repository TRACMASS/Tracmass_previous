!!----------------------------------------------------------------------
!!       Subroutine: INIT_SEED 
!!
!!          Populates the matrices seed_ijk, seed_xyz, seed_tim, seed_set
!!          which are then forwarded to Subroutine seed.
!!          
!!          seed_ijk   -  Indices of grid points to seed from.
!!          seed_xyz   -  Exact positions to seed from [optional]
!!          seed_set   -  isec (how to seed) and idir (direction to seed in)
!!          seed_tim   -  Time steps to seed at.
!!
!!----------------------------------------------------------------------

SUBROUTINE init_seed()
  
   USE mod_seed
   USE mod_time
   USE mod_param
   USE mod_write, only: rankstamp

   IMPLICIT NONE

   INTEGER                                    :: ji, jj, jk, filestat
   INTEGER                                    :: numsd, landsd=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)     :: seedMask
   CHARACTER(LEN=50)                          :: fileStamp
   CHARACTER(LEN=200)                         :: fullSeedFile
   CHARACTER(LEN=*), PARAMETER                :: xyzform = "(3f10.2,2i6,i12)"
   CHARACTER(LEN=*), PARAMETER                :: timform = "(42x      ,i12)"
#if ! defined baltix && ! defined rco
   CHARACTER(LEN=*), PARAMETER                :: ijkform = "(6i6)"
#else
   CHARACTER(LEN=*), PARAMETER                :: ijkform = "(3i10,  2i6,i12)"
#endif
   LOGICAL                                    :: fileexists

!-------------------------------------------------------------------------------
   print *, ' '
   SELECT CASE (seedType)
     
   CASE (1)      ! Seed particles within a given interval
                 ! defined by ist1, ist2, jst1, jst2, kst1, kst2 
      nsdMax = (kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1)
      ALLOCATE ( seed_ijk(nsdMax, 3), seed_set(nsdMax, 2) )
      seed_ijk = 0
      seed_set = 0
      numsd = 0
      landsd = 0
      IF (maxval(kmt) == 0) kmt = 1 
      DO ji=ist1,ist2
         DO jj=jst1,jst2
            IF (mask(ji,jj) .ne. 0) THEN 
               DO jk=kst1,kst2
                  numsd = numsd+1
                  seed_ijk (numsd,1:3) = [ ji, jj, jk ]
                  seed_set (numsd,1:2) = [ isec, idir ]
               END DO
            ELSE
               landsd = landsd+1
            END IF

         END DO
      END DO
      nsdMax = numsd
      print *,'Particles are seeded in a box defined as:'
      print '(A,I7,A,I7)', '        ist1 : ', ist1, '   ist2 : ', ist2
      print '(A,I7,A,I7)', '        jst1 : ', jst1, '   jst2 : ', jst2
      print '(A,I7,A,I7)', '        kst1 : ', kst1, '   kst2 : ', kst2

   
   CASE (2)      ! Seed particles according to indices given in a list
   
      IF (varSeedFile == 1) THEN
         fileStamp='/seed00000000.asc'
         WRITE (fileStamp(6:13),'(i8.8)') intstart/6
         fullSeedFile=trim(seedDir) // trim(fileStamp)
         PRINT *,'Particles are seeded from a dynamic listfile '
      ELSE
         fullSeedFile=trim(seedDir) // trim(seedFile)        
         PRINT *,'Particles are seeded from a given listfile  '
      END IF
      
      ! Test if file exists, and read it if it does
      INQUIRE (FILE = fullSeedFile, exist=fileexists)
      chFile: IF (fileexists) THEN
         
         nsdMax=0
         OPEN(unit=34,file=fullSeedFile, ACCESS = 'SEQUENTIAL', &
              FORM = 'FORMATTED', ACTION = 'READ')
         
         findRecl: DO
            READ (UNIT=34, fmt=timform,iostat=filestat)
            IF (filestat < 0) THEN
               EXIT findRecl
            END IF
            nsdMax = nsdMax+1
         END DO findRecl
         
         SELECT CASE (seedPos)
            CASE(1)
               ALLOCATE (seed_ijk(nsdMax,3))
            CASE(2)
               ALLOCATE (seed_ijk(nsdMax,3))
               ALLOCATE (seed_xyz(nsdMax,3))
         END SELECT
         
         nsdTim = nsdMax
         ALLOCATE (seed_set(nsdMax,2))
         
         REWIND (34)
         
         SELECT CASE (seedPos)
            CASE(1)
               DO jsd = 1,nsdMax
                  READ (unit=34, fmt=ijkform) seed_ijk(jsd,1), &
                                              seed_ijk(jsd,2), &
                                              seed_ijk(jsd,3), &
                                              seed_set(jsd,1), & 
                                              seed_set(jsd,2)
               END DO
            CASE(2)
               DO jsd = 1,nsdMax
                  READ (unit=34, fmt=xyzform) seed_xyz(jsd,1), &
                                              seed_xyz(jsd,2), &
                                              seed_xyz(jsd,3), &
                                              seed_set(jsd,1), &
                                              seed_set(jsd,2)
               END DO
               seed_ijk(:,:) = INT (seed_xyz(:,:)) + 1
         END SELECT
         
         CLOSE (34)
         
         PRINT *,'   File name   : '//trim(fullSeedFile)
      ELSE
         PRINT *,'-----------------------------------------------------'
         PRINT *,'*** ERROR!                                        ***'
         PRINT *,'*** Seed files does not exist                     ***' 
         PRINT *,'File name    : '//trim(fullSeedFile)
         PRINT *,'*** Run terminated.                               ***'
         STOP
      END IF chFile


   CASE (3)
   
      ALLOCATE ( seedMask(IMT, JMT, 1) )
      IF (varSeedFile == 1) then
         fullSeedFile=trim(seedFile) !// trim(inparg1)
         print *,'Variable file name not working!'
         stop
      ELSE
         fullSeedFile=trim(seedFile)
      END IF
      
      INQUIRE (FILE = fullSeedFile, exist=fileexists)
      
      chFile2d: IF (fileexists) THEN

         nsdMax=0
         OPEN (unit=34,file=fullSeedFile,form='unformatted', ACTION = 'READ')
            READ(unit=34) seedMask
         CLOSE (34)        
         nsdMax=0
         
         DO ji=1,imt
            DO jj=1,jmt
               DO jk=kst1,kst2
                  IF (seedMask(ji,jj,1) > 0) THEN
                     nsdMax = nsdMax + 1
                  END IF
               END DO
            END DO
         END DO
        
         ALLOCATE (seed_ijk (nsdMax,3), seed_set (nsdMax,2) )  

         jsd=0
         DO ji=1,imt
            DO jj=1,jmt
               IF (seedMask(ji,jj,1) > 0) then
                  DO jk=kst1,kst2
                     jsd = jsd + 1 
                     seed_ijk (jsd,:) = [ ji, jj, jk ]
                     seed_set (jsd,:) = [ isec, idir ]
                  END DO
               END IF
            END DO
         END DO
         print *,'Particles are seeded according to the listfile'
         print *,'   '//trim(fullSeedFile)
      ELSE
       stop 3857
      END IF chFile2d
       
   END SELECT
   print '(A,I9)','        Total number of cells : ', nsdMax + landsd
   print '(A,I9)','         Cells masked as land : ', landsd

   if (seedparts > 0) then
      if (seedpart_id > seedparts) then
         print *,"ERROR! seedpart_ID higher than number of groups"
         stop
      end if
      jj = nint(float(nsdMax)/seedparts) * (seedpart_id-1) + 1
      ji = min(jj+nsdMax/seedparts, nsdMax)
      nsdMax = ji - jj   
      seed_ijk(1:nsDmax,:) = seed_ijk(jj:ji,:)
      seed_set(1:nsDmax,:) = seed_set(jj:ji,:)
      print '(A,I2,A,I2,A)', '   seedpart is active with ', &
                            seedparts, ' groups, number ',      & 
                            seedpart_id, ' is seeded.'
      if (seedparts>0) write (rankstamp, '(A,i2.2)') '_r', seedpart_id  
   end if

   if (loneparticle>0) then
      print '(A,I7)','WARNING! Loneparticle is set to : ', loneparticle
   else
      print '(A,I9)','                   Cells used : ', nsdMax
   end if
   SELECT CASE (seedTime)
  
   CASE (1)
     
      PRINT *,'------------------------------------------------------'
      PRINT*,' Particles seeded in a time interval tst1 - tst2    '
      PRINT*,'    tst         : ',tst1,tst2
     
   CASE (2)
      
      fullSeedFile=trim(timeFile) 
      ! Test if file exists, and read it if it does
      INQUIRE (FILE = fullSeedFile, exist=fileexists)
      
      IF (fileexists) THEN
         OPEN(UNIT=34,FILE=fullSeedFile, ACCESS = 'SEQUENTIAL', &
         &    FORM = 'FORMATTED', ACTION = 'READ')
        
            ! Find number of time steps to seed    
            nsdTim = 0
            findNumberOfTimeSteps: DO
               READ (UNIT=34, fmt=timform,iostat=filestat)
               IF (filestat < 0) THEN
                  EXIT findNumberOfTimeSteps
               END IF
               nsdTim = nsdTim+1            
            END DO findNumberOfTimeSteps
          
            ALLOCATE( seed_tim(nsdTim) )
      
            REWIND (34)
         
            DO jsd=1,nsdTim
               READ (UNIT=34,FMT=timform,iostat=filestat) itim
               IF (filestat < 0) THEN
                  PRINT*,'Error in reading seed time steps'
                  PRINT*,'Seed time file: '//trim(timeFile)
                  STOP
               END IF
               seed_tim(jsd) = itim
            END DO
         
         CLOSE(34)
      
      END IF
   
      PRINT*,'------------------------------------------------------'
      PRINT*,' Particles are seeded at time steps given in list  '
      PRINT*,'File name: '//trim(timeFile)
      PRINT*,'Seed size: ',nsdTim
   
   END SELECT
   
   SELECT CASE (seedAll)
   
   CASE (1)
   
      PRINT*,'------------------------------------------------------'
      PRINT*,' All particles are seeded at each seeding time step   '
     
   CASE(2)
   
      PRINT*,'------------------------------------------------------'
      PRINT*,' Each seed position is used once                   '
   
   END SELECT
   
END SUBROUTINE init_seed
