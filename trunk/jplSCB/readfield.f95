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
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1
  INTEGER                                    :: yr2 ,mn2 ,dy2
  

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridfile ,getfile
  INTEGER, DIMENSION(1)                      :: start1d  ,count1d
  INTEGER, DIMENSION(4)                      :: start2d  ,count2d
  INTEGER, DIMENSION(4)                      :: start3d  ,count3d
  INTEGER, DIMENSION(4)                      :: start4d  ,count4d
  INTEGER                                    :: ierr
  CHARACTER (len=50)                         :: ncvar
  
  ! = ECCO Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)      :: valsz
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v ,e1t ,e2u ,e2t
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzu ,dzv ,dzt
  REAL, DIMENSION(2)                         :: ttest1, ttest2
  
  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh
  !REAL,       ALLOCATABLE, DIMENSION(:,:,:) :: uvel ,vvel 
  REAL,       ALLOCATABLE, DIMENSION(:,:,:)  :: fieldr
  ! ===   ===   ===
  
  alloCondGrid: if(.not. allocated (e1v)) then
     allocate ( valsz(km) )
     allocate ( e1v(imt+2,jmt)   ,e1t(imt+2,jmt) )
     allocate ( e2u(imt+2,jmt)   ,e2t(imt+2,jmt) )
     allocate ( dzu(imt+2,jmt,km) )
     allocate ( dzv(imt+2,jmt,km) )
     allocate ( dzt(imt+2,jmt,km) )
  end if alloCondGrid
  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
     !allocate ( uvel(imt+2,jmt,km) ,vvel(imt+2,jmt,km) )
     allocate ( fieldr(imt+2,jmt,km) )
  end if alloCondUVW
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  
  call datasetswap !Copy field(t+1) to field(t).

  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  ndates      = ints-intpart1
  dstamp      = 'scb_fcst_0000000015.nc'

  call  gdate (baseJD+intpart2 ,yr1 ,mn1 ,dy1)

  write (dstamp(10:17),'(i4i2.2i2.2)') yr1,mn1,dy1
  dataprefix  = trim(inDataDir) // '/ROMS/' // dstamp
  tpos        = intpart1+1
  
  ! === initialise ===
  initCond: if(ints.eq.intstart) then
     ! call coordinat
     hs    = 0.
     uflux = 0.
     vflux = 0.
#ifdef tempsalt
     tem   = 0.
     sal   = 0.
     rho   = 0.
#endif
     ndates=0
     
     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     getfile    = trim(dataprefix)
     start1d    = [ 1]
     count1d    = [km]
     ncvar      = 'depth'
     zw(0:km-1) = get1dfield()
     zw(km)     = 2500
     dz         = zw(km:1:-1)-zw(km-1:0:-1)
     
     ! ### FUSK! ###
     e2t=1153.4
     e2u=e2t
     e1t=2226.9
     e1v=e1t
     
     dxdy=e1t*e2t
     
     do k=1,km
        dzt(:,:,k)=dz(k)
     end do
     dzu=dzt
     dzv=dzt

  endif initCond   ! === End init section ===

  start2d   = [    1,   1, tpos,    1]
  count2d   = [imt  , jmt,    1,    1] 
  start3d   = [    1,   1,    1, tpos]
  count3d   = [imt  , jmt,   km,    1]
  getfile   = trim(dataprefix)
  ncvar     = 'u'
  uvel      = get3dfield()
  ncvar     = 'v'
  vvel      = get3dfield()
  ncvar    = 'zeta'
  ssh       = get2dfield()
  hs(:,:,2) = ssh
  
  where (uvel .eq. -9999) uvel=0
  where (vvel .eq. -9999) vvel=0
  where (ssh  .eq. -9999)  ssh=0

  dzu(1:imt-1,:,1) = dzu(1:imt-1,:,1)+0.5*(hs(1:imt-1,:,2)+hs(2:imt,:,2))
  dzv(:,1:jmt-1,1) = dzv(:,1:jmt-1,1)+0.5*(hs(:,1:jmt-1,2)+hs(:,2:jmt,2))
  dzu(imt,:,1)     = dzu(imt,:,1)    +hs(imt,:,2)
  dzv(:,jmt,1)     = dzv(:,jmt,1)    +hs(:,jmt,2)
  
  do k=1,km-1
     kk=km+1-k
     uflux(:,:,k,2)   = uvel(:,:,kk) * e2u * dzu(:,:,k)
     vflux(:,:,k,2)   = vvel(:,:,kk) * e1v * dzv(:,:,k)
     rho(:,:,k,2)     = fieldr(:,:,kk)
  end do
  
  if(ints.eq.intstart) then
     kmt=0
     do j=1,jmt
        do i=1,IMT
           do k=1,km
              kk=km+1-k
              if(uvel(i,j,k).ne.-9999) kmt(i,j)=k
              !kmt(i,j)=k
              !if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
           enddo
           if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dzt(i,j,kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        enddo
     enddo
  end if
   
 return










  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  










contains
  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  function get1dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfield
    INTEGER,             DIMENSION(1)       :: d    
    INTEGER                                 :: varid ,ncid
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count1d(1)+start1d(1)-1
    allocate ( get1dfield(d(1)) )
    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,getfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print *,'Error when trying to read the field   ',ncvar
       print *,'Error code: ' , ierr
       stop 3001
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get1dfield ,start1d ,count1d)
    fieldError: if(ierr.ne.0) then 
       print *,'Error when trying to read the field   ',ncvar
       print *, 'start1d =  ' ,start1d
       print *, 'count1d =  ' ,count1d
       print *,'Error code: ' ,ierr
       stop
    end if fieldError
    ierr=NF90_CLOSE(ncid)
  end function get1dfield

  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  function get2dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2dfield
    INTEGER,             DIMENSION(4)       :: d ,dimids ,r
    INTEGER                                 :: varid ,ncid ,i
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count2d+start2d-1
    allocate ( get2dfield(d(1),d(2)) )
    
    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print * ,'Error when trying to open the file'
       print * ,'       ' ,getfile
       print * ,'Error: ' ,NF90_STRERROR(ierr)
       stop
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print * ,'Error when trying to find the field   ',ncvar
       print * ,'Error: ' ,NF90_STRERROR(ierr)
       stop
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get2dfield ,start2d ,count2d)
    fieldError: if(ierr.ne.0) then 
       r=NF90_inquire_variable(ncid, varid, dimids = dimids)
       do i=1,4
          r=NF90_inquire_dimension(ncid, dimids(i), len=d(i))
       end do
       print * ,'Error when trying to read the field   ',ncvar
       print * ,'start2d =  ' ,start2d
       print * ,'count2d =  ' ,count2d
       print * ,'Dimensions: ' ,d
       print * ,'Error:      ' ,NF90_STRERROR(ierr)
       stop
    end if fieldError
    ierr=NF90_CLOSE(ncid)
  end function get2dfield


  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  function get3dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:,:)   :: get3dfield
    INTEGER,             DIMENSION(4)       :: d,dimids,r
    INTEGER                                 :: varid ,ncid ,i
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count3d+start3d-1
    allocate ( get3dfield(d(1),d(2),d(3)) )
    
    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print * ,'Error when trying to open the file'
       print * ,'   ' ,trim(getfile)
       print * ,'Error:      ' ,NF90_STRERROR(ierr)
       stop
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print * ,'Error when trying to find the field   ',ncvar
       print * ,'Error:      ' ,NF90_STRERROR(ierr)
       stop
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get3dfield ,start3d ,count3d)
    fieldError: if(ierr.ne.0) then 
       r=NF90_inquire_variable(ncid, varid, dimids = dimids)
       do i=1,4
          r=NF90_inquire_dimension(ncid, dimids(i), len=d(i))
       end do
       print * ,'Error when trying to read the field   ',ncvar
       print * ,'start3d =   ' ,start3d
       print * ,'count3d =   ' ,count3d
       print * ,'Dimensions: ' ,d
       print * ,'Error:      ' ,NF90_STRERROR(ierr)
       stop
    end if fieldError
    ierr=NF90_CLOSE(ncid)
  end function get3dfield

  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  function get4dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
    INTEGER,             DIMENSION(4)       :: d
    INTEGER                                 :: varid ,ncid
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count4d+start4d-1
    allocate ( get4dfield(d(1),d(2),d(3),d(4)) )

    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,getfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print *,'Error when trying to find the field   ',ncvar
       print *,'Error code: ' , ierr
       stop 3001
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get4dfield ,start1d ,count1d)
    if(ierr.ne.0) stop 3100
    ierr=NF90_CLOSE(ncid)
  end function get4dfield




  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !function get4dfield ()
  !  REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
  !  INTEGER,             DIMENSION(4)       :: d
  !end function get4dfield
  

  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  subroutine datasetswap
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    ! === swap between datasets ===
    hs(:,:,1)=hs(:,:,2)
    uflux(:,:,:,1)=uflux(:,:,:,2)
    vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1)=wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)=tem(:,:,:,2)
    sal(:,:,:,1)=sal(:,:,:,2)
    rho(:,:,:,1)=rho(:,:,:,2)
#endif
  end subroutine datasetswap


  !###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  subroutine printdiagnostics
    INTEGER                                 :: im
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    do j=2,jmt
       do i=1,IMT
          im=i-1
          if(im.eq.0) im=IMT
          do k=1,km
             if(uflux(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) then
                print *,'u',i,j,k,kmt(i,j),uflux(i,j,k,2)
             end if
             if(vflux(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) then
                print *,'v',i,j,k,kmt(i,j),vflux(i,j,k,2)
             end if
             if(uflux(im,j  ,k,2).ne.0..and. kmt(i,j).eq.0) then
                print *,'u',im,j,k,kmt(i,j),uflux(im,j,k,2)
             end if
             if(vflux(i ,j-1,k,2).ne.0..and. kmt(i,j).eq.0) then
                print *,'v',i,j-1,k,kmt(i,j),vflux(i,j-1,k,2)
             end if
          enddo
       enddo
    enddo
  end subroutine printdiagnostics

end subroutine readfields
