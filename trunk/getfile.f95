
MODULE mod_getfile
  INTEGER, DIMENSION(1)                      :: start1D  ,count1D
  INTEGER, DIMENSION(4)                      :: start2D  ,count2D ,map2D
  INTEGER, DIMENSION(4)                      :: start3D  ,count3D ,map3D
  INTEGER, DIMENSION(4)                      :: start4D  ,count4D ,map4D

  INTEGER                                 :: ierr  
  
CONTAINS

  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get1DfieldNC (fieldFile ,varName)
    USE netcdf
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfieldNC
    INTEGER,             DIMENSION(1)       :: d    
    INTEGER                                 :: varid ,ncid
  
    d=count1d(1)+start1d(1)-1
    allocate ( get1DfieldNC(d(1)) )

    ierr=NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printREadError(1)    
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printREadError(2)
    ierr=NF90_GET_VAR(ncid ,varid ,get1DfieldNC ,start1d ,count1d)
    if(ierr.ne.0) call printREadError(3)
    ierr=NF90_CLOSE(ncid)
  end function get1dfieldNC

  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get2DfieldNC (fieldFile ,varName)
    USE netcdf
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: tempField ,get2DfieldNC
    INTEGER,             DIMENSION(4)       :: d ,dimids ,r
    INTEGER                                 :: varid ,ncid ,i ,j
  
    d=count2d+start2d-1
    allocate ( tempField(d(map2D(1)),d(map2D(2))),get2DfieldNC(d(1),d(2)) )

    ierr=NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printReadError(1)
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printReadError(2)
    ierr=NF90_GET_VAR(ncid ,varid ,get2DfieldNC ,start2d ,count2d)
    print *,ierr 
    !if(ierr.ne.0) call printREadError(3)
    print *,IERR
    r=NF90_inquire_variable(ncid, varid, dimids = dimids)   
    do i=1,4
       r=NF90_inquire_dimension(ncid, dimids(i), len=d(i)) 
    end do
    print * ,'Error when trying to read the field   ',varName
    print * ,'start2d =  ' ,start2d 
    print * ,'count2d =  ' ,count2d 
    print * ,'Dimensions: ' ,d 
    print * ,'Error:      ' ,NF90_STRERROR(ierr)  
    ierr=NF90_CLOSE(ncid)

    print *,'kuk', count2D+start2D-1
    print *,shape(tempFile)
stop
    forall (i=1:imt,j=1:jmt)
       get2DfieldNC(i,j)=tempfield(map2D(1),map2D(2))
    end forall


  end function get2DfieldNC

  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get3DfieldNC (fieldFile ,varName)
    USE netcdf
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE,   DIMENSION(:,:,:)   :: get3DfieldNC
    INTEGER,             DIMENSION(4)       :: d,dimids,r
    INTEGER                                 :: varid ,ncid ,i

    d=count3d+start3d-1
    allocate ( get3DfieldNC(d(1),d(2),d(3)) )

    ierr=NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printREadError(1)
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printREadError(1)
    ierr=NF90_GET_VAR(ncid ,varid ,get3DfieldNC ,start3d ,count3d)
    if(ierr.ne.0) call printREadError(1)
    ierr=NF90_CLOSE(ncid)
  end function get3DfieldNC





  subroutine printReadError(sel)
    USE netcdf
    INTEGER                                 :: sel
    
    select case (sel)
    case (1)
       print *,'Error when trying to open the file'
       print *,'   ' ,getfile
       print *,'    Error: ' , NF90_STRERROR(ierr)
       stop 
    case (2)
       print *,'Error when trying to read the field   ',ncvar
       print *,'    Error: ' , NF90_STRERROR(ierr)
       stop 
    case (3)
       print *,'Error when trying to read the field   ',ncvar
       print *, 'start1d =  ' ,start1d
       print *, 'count1d =  ' ,count1d
       print *,'Error:      ' ,NF90_STRERROR(ierr)
       stop
       
       !r=NF90_inquire_variable(ncid, varid, dimids = dimids)
       !do i=1,4
       !   r=NF90_inquire_dimension(ncid, dimids(i), len=d(i))
       !end do
       !print * ,'Error when trying to read the field   ',varName
       !print * ,'start2d =  ' ,start2d
       !print * ,'count2d =  ' ,count2d
       !print * ,'Dimensions: ' ,d
       !print * ,'Error:      ' ,NF90_STRERROR(ierr)
       !stop

    end select
  end subroutine printReadError
end MODULE mod_getfile





!###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###     
!  function get4dfield ()
!    REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
!    INTEGER,             DIMENSION(4)       :: d
!    INTEGER                                 :: varid ,ncid
!  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   === 
!    d=count4d+start4d-1
!    allocate ( get4dfield(d(1),d(2),d(3),d(4)) )
!
!    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
!    fileError: if(ierr.ne.0) then
!       print *,'Error when trying to open the file'
!       print *,'   ' ,getfile
!       print *,'    Error code: ' , ierr
!       stop 3001
!    end if fileError
!    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
!    varError: if(ierr.ne.0) then
!       print *,'Error when trying to find the field   ',varName
!       print *,'Error code: ' , ierr
!       stop 3001
!    end if varError
!    ierr=NF90_GET_VAR(ncid ,varid ,get4dfield ,start1d ,count1d)
!    if(ierr.ne.0) stop 3100
!    ierr=NF90_CLOSE(ncid)
!  end function get4dfield
