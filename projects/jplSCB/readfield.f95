SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1,hr
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName
  
  ! = ECCO Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)      :: valsz
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v ,e1t ,e2u ,e2t
  REAL, DIMENSION(2)                         :: ttest1, ttest2
  
  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh
  ! ===   ===   ===
  
  alloCondGrid: if(.not. allocated (e1v)) then
     allocate ( valsz(km) )
     allocate ( e1v(imt+2,jmt)   ,e1t(imt+2,jmt) )
     allocate ( e2u(imt+2,jmt)   ,e2t(imt+2,jmt) )
  end if alloCondGrid
  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  
  call datasetswap !Copy field(t+1) to field(t).
  call updateClock

  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = 'scb_das_0000000000.nc'

  write (dstamp(9:18),'(i4i2.2i2.2i2.2)') & 
       currYear,currMon,currDay,currHour+3
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  print *,dataprefix
  ! === initialise ===
  initCond: if(ints.eq.intstart) then
     ! call coordinat
     ndates = 0
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
     
     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     zw(0:km-1) = get1DfieldNC(trim(dataprefix) ,'depth')
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

     ncTpos=1
     uvel = get3DfieldNC(trim(dataprefix) ,   'temp')

     kmt  = 0 
     do j=1,jmt
        do i=1,IMT
           do k=1,km
              kk=km+1-k
              if(uvel(i,j,k).ne.-9999) kmt(i,j)=k
           end do
           if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dzt(i,j,kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        end do
     end do


  endif initCond   ! === End init section ===

  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh       = get2dfieldNC(trim(dataprefix) ,'zeta')
  hs(:,:,2) = ssh

  where (uvel .eq. -9999) uvel=0
  where (vvel .eq. -9999) vvel=0
  where (ssh  .eq. -9999)  ssh=0

  dzu(1:imt-1,:,1) = dzu(1:imt-1,:,1)+0.5*(hs(1:imt-1,:,2)+hs(2:imt,:,2))
  dzv(:,1:jmt-1,1) = dzv(:,1:jmt-1,1)+0.5*(hs(:,1:jmt-1,2)+hs(:,2:jmt,2))
  dzu(imt,:,1)     = dzu(imt,:,1)    +hs(imt,:,2)
  dzv(:,jmt,1)     = dzv(:,jmt,1)    +hs(:,jmt,2)
  
  do k=1,km
     kk=km+1-k
     uflux(:,:,k,2)   = uvel(:,:,kk) * e2u * dzu(:,:,k)
     vflux(:,:,k,2)   = vvel(:,:,kk) * e1v * dzv(:,:,k)
  end do

 return

 !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

 
contains

  
  subroutine datasetswap
    hs(:,:,1)      = hs(:,:,2)
    uflux(:,:,:,1) = uflux(:,:,:,2)
    vflux(:,:,:,1) = vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1) = wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)   = tem(:,:,:,2)
    sal(:,:,:,1)   = sal(:,:,:,2)
    rho(:,:,:,1)   = rho(:,:,:,2)
#endif
  end subroutine datasetswap


end subroutine readfields
