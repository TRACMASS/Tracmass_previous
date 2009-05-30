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

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos

  integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
  integer varid,varidx,varidy,varidz,varidt !output ID index of variable
  integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
  integer p, x1, y1, z1, t1 !?

  ! = Variables for filename generation
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2 ,subYr
  INTEGER                                    :: filePos ,fileJD ,subYrD
  INTEGER                                    :: yr1 ,mn1 ,dy1
  INTEGER                                    :: yr2 ,mn2 ,dy2
  

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: ssh
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof
  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,e1t,e2u,e2t
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzu,dzv,dzt
  
  logical around
  
  alloCondGrid: if ( .not. allocated (e1v) ) then
     allocate ( e1v(IMT+2,JMT)    ,e1t(IMT+2,JMT) )
     allocate ( e2u(IMT+2,JMT)    ,e2t(IMT+2,JMT) )
     allocate ( dzu(IMT+2,JMT,KM) ,dzv(IMT+2,JMT,KM),dzt(IMT+2,JMT,KM) )
  end if alloCondGrid
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
     allocate ( rhof(IMT+2,JMT,KM) )
  end if alloCondUVW
  

  ! === swap between datasets ===
  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif
print *,' '
  do ints=1,200
     call updateClock
     ! === update the time counting ===
     filePos = mod((ints-1),fieldsPerFile)+1
     fileJD  = mod(floor(real((ints-1))/real(fieldsPerFile))+2,4)+1
     subYr   = mod(floor(real(ints-1)/real(fieldsPerFile)),4)+1
     dstamp  = 'KAB042j_5d_00000000_00000000'
     
     call  gdate (baseJD+fileJD ,yr1 ,mn1 ,dy1)
     call  gdate (baseJD+intpart2*90+intpart2+92 ,yr2 ,mn2 ,dy2)
  write (dstamp(12:19),'(i4i2.2i2.2)') currYear,mn1,dy1
  write (dstamp(21:28),'(i4i2.2i2.2)') currYear,mn2,dy2
  dataprefix  = trim(inDataDir) // '/' // dstamp
  tpos        = intpart1
 
  print *,ints, filePos ,fileJD ,currYear ,currMon ,currDay ,trim(dstamp)
enddo
 stop
 
!  iday=iday+5
!  if(iday.gt.idmax(imon,iyear)) then
!     iday=iday-idmax(imon,iyear)
!     imon=imon+1
!     if(imon.eq.13) then
!        imon=1
!        iyear=iyear+1
!        ! if kan skrivas här om man vill börja om från iyear0
!     endif
!  endif
!  ntime=10000*iyear+100*imon+iday

  start1d  = [ 1]
  count1d  = [km]
  start2d  = [  1   ,   1, tpos,    1]
  count2d  = [imt+2 , jmt,    1,    1]
  start3d  = [  1   ,   1,    1, tpos]
  count3d  = [imt+2 , jmt,   km,    1]
  

  ! === initialise ===
  initFieldcond: if(ints.eq.intstart) then
     ! call coordinat
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

     gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'
     zw(0)    = 0.d0
     zw(1:km) = get1DfieldNC(trim(gridFile) ,'nav_lev')

     do k=1,km
        kk=km+1-k
        dz(kk)=zw(k)-zw(k-1) 
        ! print *,k,zw(k),kk,dz(kk)
     end do

     e1t  = get2DfieldNC(trim(gridFile) ,'e1t')
     e2u  = get2DfieldNC(trim(gridFile) ,'e2u')
     e2t  = get2DfieldNC(trim(gridFile) ,'e2t')
     e1v  = get2DfieldNC(trim(gridFile) ,'e1v')
     dxdy = e1t * e2t
   
     gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
     dzt  = get3DfieldNC(trim(gridFile) ,'e3t_ps')
     dzu  = get3DfieldNC(trim(gridFile) ,'e3u_ps')
     dzv  = get3DfieldNC(trim(gridFile) ,'e3v_ps')

     fieldFile = trim(inDataDir)//dataprefix
     rhof = get3DfieldNC(trim(fieldFile)//'_sigma.nc'  ,'sigma')
     do j=1,jmt
        do i=1,IMT
           do k=1,km
              kk=km+1-k
              if(rhof(i,j,k).ne.0.) kmt(i,j)=k
              if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
           enddo
           if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dzt(i,j,kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        enddo
     enddo
     
  endif initFieldcond
  
  fieldFile = trim(inDataDir)//dataprefix
  
  uvel = get3DfieldNC(trim(fieldFile)//'_grid_U.nc' ,'vozocrtx')
  uvel = get3DfieldNC(trim(fieldFile)//'_grid_V.nc' ,'vomecrty')
  rhof = get3DfieldNC(trim(fieldFile)//'_sigma.nc'  ,'sigma')
  ssh  = get2DfieldNC(trim(fieldFile)//'_SSH.nc'    ,'sossheig')

  hs(:,:,2) = 0.01*ssh
  hs(imt+1,:,2) =hs(1,:,2)
  hs(:,jmt+1,2) =hs(:,1,2)

  do k=1,km-1
     kk=km+1-k
     uflux(:,:,k,2) = uvel(:,:,kk) * e2u(:,:) * dzu(:,:,kk)
     vflux(:,:,k,2) = vvel(:,:,kk) * e1v(:,:) * dzv(:,:,kk)
     rho  (:,:,k,2) = rhof(:,:,kk)
  end do

  uflux(:,:,km,2) = uvel(:,:,1) * e2u * (dzu(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(2:imt+1,:,2)))
  vflux(:,:,km,2) = vvel(:,:,1) * e1v * (dzv(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(:,2:jmt+1,2)))
  rho(:,:,km,2)   = rhof(:,:,1)
  
  return
end subroutine readfields
