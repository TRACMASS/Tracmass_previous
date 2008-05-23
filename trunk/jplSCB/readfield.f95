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
  
  ! === Variables for filename generation ===
  CHARACTER               :: dates(62)*17
  INTEGER, SAVE           :: nread,ndates
  CHARACTER(LEN=65), SAVE :: rfilu,rfilv,rfilh,rfilr

  ! === Loop variables ===
  integer i,ip,j,jp,k,kk,ints2,im
  
  ! === Variables used for netcdf procedures
  CHARACTER (len=200)               :: ncfile
  CHARACTER (len=30)                :: ncvar 
  REAL*4, DIMENSION(IMT+2,JMT,KM,1) :: ncfield
  INTEGER, DIMENSION(1)             :: start1d, count1d
  INTEGER, DIMENSION(2)             :: start2d, count2d
  INTEGER, DIMENSION(4)             :: start4d, count4d

  integer ncid !output ID index of netCDF file
  integer ierr !error 
  integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
  integer varid,varidx,varidy,varidz,varidt !output ID index of variable
  integer, dimension(4) :: start, count
  integer startA(1),startB(4),startC(2) !input index vector of position to start reading
  integer countA(1),countB(4),countC(2) !input lengths of 'volume' to be retrieved
  integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
  integer p, x1, y1, z1, t1 !?
  
  REAL*4, DIMENSION(KM) :: valsz
  !REAL*4, ALLOCATABLE, DIMENSION(:) :: ssh
  REAL*4, DIMENSION(IMT+2,JMT) :: ssh
  REAL*4, DIMENSION(IMT+2,JMT,KM,1) :: fieldx,fieldy,fieldr
  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: e1v,e1t,e2u,e2t
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dzu,dzv,dzt
  
  logical around
  
  if ( .not. allocated (e1v) ) then
     allocate ( e1v(IMT+2,JMT),e1t(IMT+2,JMT),e2u(IMT+2,JMT),e2t(IMT+2,JMT) )
     allocate ( dzu(IMT+2,JMT,KM),dzv(IMT+2,JMT,KM),dzt(IMT+2,JMT,KM) )
  end if
  
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

  !____________________________ initialise ___________________________
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
     ! ===  Set up the grid ===
     ! ======================================================
     ncfile  = trim(directory) // 'test.nc'
     ncvar   = 'lat'
     start1d = [1]
     count1d = [24]
     valsz   = get1dfield()
     
     ncvar   = 'lat'
     start2d = [1    ,1  ]
     count2d = [IMT+2,jmt]
     e1t     = get2dfield()
     
     ncvar   = 'lat'
     e1t     = get2dfield()
     print *, e1t

     ierr=NF90_CLOSE(ncid)
     if(ierr.ne.0) stop 3040
     
     do i=1,IMT
        do j=1,jmt
           dxdy(i,j)=e1t(i,j)*e2t(i,j)
        enddo
     enddo
     
     ! === Read ORCA grid horizontal ===

     ierr=NF90_OPEN(directory//'topo/mesh_zgr.nc',NF90_NOWRITE,ncid)
     if(ierr.ne.0) stop 4001
     ! dzt
     ierr=NF90_INQ_VARID(ncid,'e3t_ps',varid)
     if(ierr.ne.0) stop 4004
    
     startB=[1,1,1,nread]
     countB=[IMT+2,jmt,km,1]
     ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
     if(ierr.ne.0) stop 4002
     dzt(:,:,:)=fieldx(:,:,:,1)
     !print *,'dzt',(dzt(200,200,k),k=1,km)
     
     ! === dzu ===
     ierr=NF90_INQ_VARID(ncid,'e3u_ps',varid)
     if(ierr.ne.0) stop 4014
     startB=[1,1,1,nread]
     countB=[IMT+2,jmt,km,1]
     ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
     if(ierr.ne.0) stop 4012
     dzu(:,:,:)=fieldx(:,:,:,1)
     !print *,'dzu',(dzu(200,200,k),k=1,km)
     
     ! dzv
     ierr=NF90_INQ_VARID(ncid,'e3v_ps',varid)
     if(ierr.ne.0) stop 4024
 
     startB=[1,1,1,nread]
     countB=[IMT+2,jmt,km,1]
     ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
     if(ierr.ne.0) stop 4032
     dzv(:,:,:)=fieldx(:,:,:,1)
     !print *,'dzv',(dzv(200,200,k),k=1,km)
     
     ierr=NF90_CLOSE(ncid)
     if(ierr.ne.0) stop 4031
     
     !stop 4596
     
  endif initCond
  ! === End init ection ===
  
  if(mod(ints,18).eq.1) then
     
     ints2=ints
666  continue
     if(ints2.gt.intmax .or. ints2.lt.intmin) then
        ints2=ints2-intend+intstart-intstep
        goto 666
     endif
     
     if(iday0.eq.2-10 .and. imon0.eq.7 .and. iyear0.eq.1990) then
        ndates=ints2/18+1
     elseif(iday0.eq.1-10 .and. imon0.eq.10 .and. iyear0.eq.1992) then
        ndates=ints2/18+10
     else
        print *,iyear0,imon0,iday0,ndates
        print *,iyear,imon,iday,ndates
        stop 9567
     endif
     
     rfilu=directory//'gcm/KAB042j_5d_'//dates(ndates)//'_grid_U.nc'
     rfilv=directory//'gcm/KAB042j_5d_'//dates(ndates)//'_grid_V.nc'
     rfilr=directory//'gcm/KAB042j_5d_'//dates(ndates)//'_sigma.nc'
     rfilh=directory//'gcm/KAB042j_5d_'//dates(ndates)//'_SSH.nc'
     print *,'rfilu=',rfilu,ints,ints2,ndates,dates(ndates)
     inquire(file=rfilu,exist=around)
     if(.not.around) stop 4556
     inquire(file=rfilv,exist=around)
     if(.not.around) stop 4557
     inquire(file=rfilr,exist=around)
     if(.not.around) stop 4558
     inquire(file=rfilh,exist=around)
     if(.not.around) stop 4559
  endif
  !nread=mod(ints/5,18)+1
  nread=mod(ints,18)+1
  
  ! === zonal velocity ===
  ierr=NF90_OPEN(rfilu,NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3750
  
  
  ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) ! the main data fields
  if(ierr.ne.0) stop 3762
  
  startB=[1     ,1   ,1  ,nread]
  countB=[IMT+2 ,jmt ,km ,1    ]
  ierr=NF90_GET_VAR(ncid,varid,fieldx,startB,countB)
  if(ierr.ne.0) stop 3798
  
  ierr=NF90_CLOSE(ncid)
   
  ! === meridional velocity ===  
  ierr=NF90_OPEN(rfilv,NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3751
  
  ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! the main data fields
  if(ierr.ne.0) stop 3763
  
  startB=[1,1,1,nread]
  countB=[IMT+2,jmt,km,1]
  ierr=NF90_GET_VAR(ncid,varid,fieldy,startB,countB)
  if(ierr.ne.0) stop 3799
  ierr=NF90_CLOSE(ncid)
  
  ! === density ===
  ierr=NF90_OPEN(rfilr,NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3758
  
  ierr=NF90_INQ_VARID(ncid,'sigma',varid) ! the main data fields
  if(ierr.ne.0) stop 3767
  
  startB=[1,1,1,nread]
  countB=[IMT+2,jmt,km,1]
  ierr=NF90_GET_VAR(ncid,varid,fieldr,startB,countB)
  if(ierr.ne.0) stop 3799
  ierr=NF90_CLOSE(ncid)
  
  ! === sea surface height ===
  ierr=NF90_OPEN(rfilh,NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3759
  
  ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
  if(ierr.ne.0) stop 3763
  
  startB=[1,1,1,nread]
  countB=[IMT+2,jmt,1,1]
  ierr=NF90_GET_VAR(ncid,varid,ssh,startB,countB)
  if(ierr.ne.0) stop 3799
  ierr=NF90_CLOSE(ncid)
  
  do j=1,jmt
     do i=1,IMT
        hs(i,j,2)=0.01*ssh(i,j)
     enddo
  enddo
  
  do j=1,jmt
     jp=j+1
     if(jp.gt.jmt) jp=jmt
     do i=1,IMT
        ip=i+1
        if(ip.eq.IMT+1) ip=1
        do k=1,km-1
           kk=km+1-k
           u(i,j,k,2)=fieldx(i,j,kk,1)*e2u(i,j)*dzu(i,j,kk)
           v(i,j,k,2)=fieldy(i,j,kk,1)*e1v(i,j)*dzv(i,j,kk)
           rho(i,j,k,2)=fieldr(i,j,kk,1)
        enddo
        u(i,j,km,2)=fieldx(i,j,1,1)*e2u(i,j)*(dzu(i,j,1)+0.5*(hs(i,j,2)+hs(ip,j,2)))
        v(i,j,km,2)=fieldy(i,j,1,1)*e1v(i,j)*(dzv(i,j,1)+0.5*(hs(i,j,2)+hs(i,jp,2)))
        rho(i,j,km,2)=fieldr(i,j,1,1)
        !   if(i.eq.210 .and. j.eq.252) print *,'density=',rho(i,j,km,2)
     enddo
  enddo
  
  if(ints.eq.intstart) then
     do j=1,jmt
        do i=1,IMT
           do k=1,km
              kk=km+1-k
              if(fieldr(i,j,k,1).ne.0.) kmt(i,j)=k
              if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
           enddo
           if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dzt(i,j,kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        enddo
     enddo
     
     
679  format(i3,1x,400i1)
     
  endif
  
  do j=2,jmt
     do i=1,IMT
        im=i-1
        if(im.eq.0) im=IMT
        do k=1,km
           if(u(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'u',i,j,k,kmt(i,j),u(i,j,k,2)
           if(v(i ,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'v',i,j,k,kmt(i,j),v(i,j,k,2)
           if(u(im,j  ,k,2).ne.0..and. kmt(i,j).eq.0) print *,'u',im,j,k,kmt(i,j),u(im,j,k,2)
           if(v(i ,j-1,k,2).ne.0..and. kmt(i,j).eq.0) print *,'v',i,j-1,k,kmt(i,j),v(i,j-1,k,2)
        enddo
     enddo
  enddo
  return
  
!  #########         #########         #########         #########
!  #########         #########         #########         #########
!  #########         #########         #########         #########
!  #########         #########         #########         #########
!  #########         #########         #########         #########
  
contains
  function get1dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfield
    INTEGER,             DIMENSION(1)       :: d    
    d=count1d+start1d-1
    allocate ( get1dfield(d(1)) )
    
    ierr=NF90_OPEN(trim(ncfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,ncfile
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
    if(ierr.ne.0) stop 3100
    ierr=NF90_CLOSE(ncid)
  end function get1dfield

  function get2dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2dfield
    INTEGER,             DIMENSION(2)       :: d
    d=count2d+start2d-1
    allocate ( get2dfield(d(1),d(2)) )
    
    ierr=NF90_OPEN(trim(ncfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,ncfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print *,'Error when trying to read the field   ',ncvar
       print *,'Error code: ' , ierr
       stop 3001
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get2dfield ,start1d ,count1d)
    if(ierr.ne.0) stop 3100
    ierr=NF90_CLOSE(ncid)
  end function get2dfield
  
  function get4dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
    INTEGER,             DIMENSION(4)       :: d
    
    d=count4d+start4d-1
    allocate ( get4dfield(d(1),d(2),d(3),d(4)) )

    ierr=NF90_OPEN(trim(ncfile) ,NF90_NOWRITE ,ncid)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,ncfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    ierr=NF90_INQ_VARID(ncid ,ncvar ,varid)
    varError: if(ierr.ne.0) then
       print *,'Error when trying to read the field   ',ncvar
       print *,'Error code: ' , ierr
       stop 3001
    end if varError
    ierr=NF90_GET_VAR(ncid ,varid ,get4dfield ,start1d ,count1d)
    if(ierr.ne.0) stop 3100
    ierr=NF90_CLOSE(ncid)
  end function get4dfield
  





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
