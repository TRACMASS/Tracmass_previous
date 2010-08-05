
subroutine init_seed()
  
  use mod_seed
  use mod_time
  use mod_param

  IMPLICIT NONE

  CHARACTER(LEN=50)                          :: fileStamp
  CHARACTER(LEN=200)                         :: fullSeedFile
  CHARACTER(len=*), parameter                :: ijkform = "(6I6)"
  LOGICAL                                    :: fileexists
  INTEGER                                    :: i,j,k
  INTEGER                                    :: filestat
  INTEGER*4, ALLOCATABLE, DIMENSION(:,:)     :: seedMask
  select case (seedType)
     
  case (1) ! === ist, jst, kst method ===
     allocate (ijkst((kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1),6))
     ijk=1
     if (maxval(kmt) == 0) kmt = 1 
     do i=ist1,ist2
        do j=jst1,jst2
           if (kmt(i,j) .ne. 0) then 
              do k=kst1,kst2
                 ijkst(ijk,1)=i
                 ijkst(ijk,2)=j
                 ijkst(ijk,3)=k
                 ijkst(ijk,4)=idir
                 ijkst(ijk,5)=isec
                 ijkst(ijk,6)=-1
                 ijk=ijk+1
              end do
           end if
        end do
     end do
     ijkMax = ijk
     print *,'------------------------------------------------------'
     print *,'=== Particles are seeded using ist, jst, kst          '
     print *,'ist         : ', ist1,ist2
     print *,'jst         : ', jst1,jst2
     print *,'kst         : ', kst1,kst2
     print *,'Seed size   : ', ijkMax
  
  case (2) ! === seedlist method ===   
     print *,'------------------------------------------------------'
     if (varSeedFile == 1) then
        fileStamp='/seed00000000.asc'
        write (fileStamp(6:13),'(i8.8)') intstart/6
        fullSeedFile=trim(seedDir) // trim(fileStamp)
        print *,'===  Particles are seeded from a dynamic listfile  ==='
     else
        fullSeedFile=trim(seedFile)        
        print *,'===  Particles are seeded from a given listfile    ==='
     end if
     inquire(FILE = fullSeedFile, exist=fileexists)
     chFile: if (fileexists) then
        ijkMax=0
        open(unit=34,file=fullSeedFile, ACCESS = 'SEQUENTIAL', &
             FORM = 'FORMATTED', ACTION = 'READ')
        findRecl: do
           read (unit=34, fmt=ijkform,iostat=filestat)
           if (filestat < 0) exit findRecl
           ijkMax=ijkMax+1
        end do findRecl
        allocate (ijkst(ijkMax,6))
           
        rewind(34)
        read_ijkst: do ijk=1,ijkMax
           read (unit=34, fmt=ijkform) ijkst(ijk,1), ijkst(ijk,2), &
                ijkst(ijk,3),  ijkst(ijk,4),  ijkst(ijk,5), ijkst(ijk,6)


           ijkst(ijk,4)=idir
           ijkst(ijk,5)=isec
           ijkst(ijk,6)=-1
           
        end do read_ijkst
        
        print *,'File name    : '//trim(fullSeedFile)
        print *,'Seed size    : ', ijkMax
     else
        print *,'======================================================'
        print *,'*** ERROR!                                         ***'
        print *,'*** Seed files does not exisit                     ***' 
        print *,'File name    : '//trim(fullSeedFile)
        print *,'*** Run terminated.                                ***'
        stop
     end if chFile
     



  case (3) ! === 2-D Matrix method ===
     allocate ( seedMask(imt ,jmt) )
     if (varSeedFile == 1) then
        fullSeedFile=trim(seedFile) !// trim(inparg1)
        print *,'Variable file name not working!'
        stop
     else
        fullSeedFile=trim(seedFile)
     end if
     inquire(FILE = fullSeedFile, exist=fileexists)
     chFile2d: if (fileexists) then

        ijkMax=0
        open(unit=34,file=fullSeedFile,form='unformatted', ACTION = 'READ')
        read(unit=34) seedMask
        close(34)
        ijkMax=0

        do i=1,imt
         do j=1,jmt
          do k=kst1,kst2
!           if (seedMask(i,j) /= 0) ijkMax=ijkMax+1
           if (seedMask(i,j) /= 0) ijkMax=ijkMax+seedMask(i,j)
          end do
         end do
        end do
        
        allocate (ijkst(ijkMax,6))  

        ijk=0
        do i=1,imt
           do j=1,jmt
              if (seedMask(i,j) /= 0) then
                 do k=kst1,kst2
                    ijk=ijk+1
                    ijkst(ijk,1)=i
                    ijkst(ijk,2)=j
                    ijkst(ijk,3)=k
                    ijkst(ijk,4)=idir
                    ijkst(ijk,5)=isec
                    ijkst(ijk,6)=seedMask(i,j)
                 end do
              end if
           end do
        end do

        
        print *,'------------------------------------------------------'
        print *,'=== Particles are seeded according to a listfile   ==='
        print *,'File name    : '//trim(fullSeedFile)
        print *,'Seed size    : ', ijkMax
     end if chFile2d

#if defined rco
     do i=1,IMT
        do j=1,JMT
           if(seedMask(i,j).ne.0 .and. seedMask(i,j).le.4 .and. j.lt.215) seedMask(i,j)=-1  ! entire shallow Baltic south of 61N
        enddo
     enddo
#endif



!!$! === From loop.f95 ===
!!$#if defined for || sim 
!!$           if(mask(ist,jst).le.0) cycle jstLoop
!!$#ifdef twodim
!!$                 if(KM.ne.kst) cycle kstLoop
!!$#else
!!$                 if(KM+1-kmt(ist,jst).ne.kst) cycle kstLoop
!!$#endif
!!$                 luy=luy+1
!!$#endif
!!$#ifdef rco
!!$                 if(mask(ist,jst).ne.-1)  cycle jstLoop
!!$#endif
!!$#ifdef orc
!!$                 if(ist.le.0 .or. ist.gt.IMT .or. jst.le.0 .or. &
!!$                      jst.gt.JMT) print *,'hoooopsnasa',ist,jst
!!$                 if(mask(ist,jst).ne.-1)  cycle jstLoop
!!$#endif







  end select

end subroutine init_seed
