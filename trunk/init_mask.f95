
subroutine init_mask()
  
  use mod_seed
  CHARACTER(LEN=200)                         :: fullSeedFile
  CHARACTER(len=*), parameter                :: ijkform = "(6I6)"
  LOGICAL                                    :: fileexists
  INTEGER                                    :: filestat
  select case (seedType)
     
  case (1) ! === ist, jst, kst method ===
     allocate (ijkst((kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1),6))
     ijk=1
     do i=ist1,ist2
        do j=jst1,jst2
           do k=kst1,kst2
              ijkst(ijk,1)=i
              ijkst(ijk,2)=j
              ijkst(ijk,3)=k
              ijkst(ijk,4)=idir
              ijkst(ijk,5)=isec
              ijkst(ijk,6)=-1
              ijk=ijk+1
           end do
        end do
     end do
     ijkMax=(kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1)
  
  case (2) ! === seedlist method ===
     
     if (varSeedFileName == 1) then
        fullSeedFile=trim(seedFile) !// trim(inparg1)
        print *,'Variable file name not working!'
        stop
     else
        fullSeedFile=trim(seedFile)
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
        end do read_ijkst
        
        print *,'------------------------------------------------------'
        print *,'=== Particles are seeded according to a listfile   ==='
        print *,'File name : '//trim(fullSeedFile)
        print *,'Seed size : ', ijkMax
        stop
     end if chFile
        
        
  case (3) ! === Matrix method ===
     mask=-1.  ! define start section with ist1,ist2,jst1,jst2
#ifndef atm || ifs
     !open(21,file=directory//'topo/kmt',form='unformatted')
     !open(21,file=directory//'topo/maskust',form='unformatted')
     !read(21)mask
     !close(21)
#endif
     
#if defined rco
     do i=1,IMT
        do j=1,JMT
           if(mask(i,j).ne.0 .and. mask(i,j).le.4 .and. j.lt.215) mask(i,j)=-1  ! entire shallow Baltic south of 61N
        enddo
     enddo
#endif

  end select

end subroutine init_mask
