SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                 :: dates(62)*17
  INTEGER, SAVE                             :: nread,ndates
  CHARACTER(LEN=10), SAVE                   :: fstamp,rfilv,rfilh,rfilr
  logical around
  
  ! = Loop variables
  INTEGER                                   :: t ,i ,j ,k
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                       :: gridfile
  INTEGER                                   :: start1d, count1d
  INTEGER, DIMENSION(2)                     :: start2d, count2d
  INTEGER, DIMENSION(3)                     :: start3d, count3d
  INTEGER                                   :: ierr

  ! = ECCO Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)     :: gridDRC ,gridDRF
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDXC ,gridDXG
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDYC ,gridDYG
  REAL, DIMENSION(2)                        :: ttest1, ttest2
  
  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:,:) :: fieldx ,fieldy ,fieldw
 
  ! ===   ===   ===
  
  alloCondGrid: if(.not. allocated (gridDRC)) then
     allocate ( gridDRC(km)      ,gridDRF(km)      )
     allocate ( gridDXC(imt,jmt) ,gridDXG(imt,jmt) )
     allocate ( gridDYC(imt,jmt) ,gridDYG(imt,jmt) )
  end if alloCondGrid
  
  alloCondUVW: if(.not. allocated (fieldx)) then
     allocate ( fieldx(imt,jmt,km) ,fieldy(imt,jmt,km) )
     allocate ( fieldw(imt,jmt,km) )
  end if alloCondUVW
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  call datasetswap !Copy field(t+1) to field(t).
  
  ! === update the time counting ===
  iday=iday+5
  if(iday.gt.idmax(imon,iyear)) then
     iday=iday-idmax(imon,iyear)
     imon=imon+1
     if(imon.eq.13) then
        imon=1
        iyear=iyear+1
        ! if kan skrivas här om man vill börja om från iyear0
     endif
  endif
  ntime=10000*iyear+100*imon+iday

  ! === initialise ===
  !print *,'ints=',ints,intstart
  initCond: if(ints.eq.intstart) then
     ! call coordinat
     hs=0.
     u=0.
     v=0.
#ifdef tempsalt
     tem=0.
     sal=0.
     rho=0.
#endif
     ndates=0

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     start1d  =  1
     count1d  = 42
     gridfile = trim(directory) // '/GRID/' // 'DRC.data'
     gridDRC  = get1dfield()
     gridfile = trim(directory) // '/GRID/' // 'DRF.data'
     gridDRF  = get1dfield()
     
     start2d  = [   1,  1]
     count2d  = [2160,320]
     gridfile = trim(directory) // '/GRID/' // 'DXC.data'
     gridDXC  = get2dfield()
     gridfile = trim(directory) // '/GRID/' // 'DXG.data'
     gridDXG  = get2dfield()
     gridfile = trim(directory) // '/GRID/' // 'DYC.data'
     gridDYC  = get2dfield()
     gridfile = trim(directory) // '/GRID/' // 'DXG.data'
     gridDYG  = get2dfield() 
     
     dxdy=gridDXC*gridDYC
   
     !     dzt(:,:,:)=
     !     dzu(:,:,:)=fieldx(:,:,:,1)
     !     dzv(:,:,:)=fieldx(:,:,:,1)
     
  endif initCond   ! === End init section ===
  
  
  !  if(mod(ints,18).eq.1) then

  fstamp='0000000960'

  start3d  = [   1,  1, 1]
  count3d  = [2160,320,42]
  gridfile = trim(directory) // '/UVW/' // 'UVEL_5dy.' // fstamp // '.data'
  fieldx   = get3dfield()
  gridfile = trim(directory) // '/UVW/' // 'VVEL_5dy.' // fstamp // '.data'
  fieldy   = get3dfield()
  gridfile = trim(directory) // '/UVW/' // 'WVEL_5dy.' // fstamp // '.data'
  fieldw   = get3dfield()

  !Density not included
  !SSH not included
  
  iloop: do i=1,imt
     jloop: do j=1,jmt
        kloop: do k=1,km
           u(i,j,k,2)=fieldx(i,j,k)*gridDYC(i,j)*gridDRC(k)
           v(i,j,k,2)=fieldy(i,j,k)*gridDXC(i,j)*gridDRC(k)
           !w(i,j,k,2)=fieldz(i,j,k)*dxdy(i,j)
           rho(i,j,k,2)=0
        end do kloop
     enddo jloop
  enddo iloop
  
  return











  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  











contains
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  function get1dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfield
    INTEGER                                 :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count1d+start1d-1
    allocate ( get1dfield(d) )
    rl=d*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get1dfield
    close (3001)
  end function get1dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  function get2dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2dfield
    INTEGER,             DIMENSION(2)       :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count2d+start2d-1
    allocate ( get2dfield(d(1),d(2)) )
    rl=product((d))*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get2dfield
    close (3001)
  end function get2dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  function get3dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:,:)   :: get3dfield
    INTEGER,             DIMENSION(3)       :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count3d+start3d-1
    allocate ( get3dfield(d(1),d(2),d(3)) )
    rl=product((d))*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get3dfield
    close (3001)
  end function get3dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  !function get4dfield ()
  !  REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
  !  INTEGER,             DIMENSION(4)       :: d
  !end function get4dfield
  

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  subroutine datasetswap
    ! === swap between datasets ===
    hs(:,:,1)=hs(:,:,2)
    u(:,:,:,1)=u(:,:,:,2)
    v(:,:,:,1)=v(:,:,:,2)
#ifdef explicit_w
    w(:,:,:,1)=w(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)=tem(:,:,:,2)
    sal(:,:,:,1)=sal(:,:,:,2)
    rho(:,:,:,1)=rho(:,:,:,2)
#endif
  end subroutine datasetswap
  
end subroutine readfields
