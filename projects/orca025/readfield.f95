SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_traj
  USE mod_getfile
  use mod_seed
  USE mod_deformation
  USE mod_laplacian

#ifdef tempsalt
!  USE mod_dens
  USE mod_stat
#endif
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates

  ! = Loop variables
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: nread

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile, zfile

#ifdef orca1
  INTEGER, PARAMETER :: IMTG=???,JMTG=???
#elif orca025l75h6
  INTEGER, PARAMETER :: IMTG=1440,JMTG=1021,KMM=75
#endif

  REAL*4, ALLOCATABLE, DIMENSION(:,:)         	:: temp2d_simp
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:),SAVE    :: u_m, v_m
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)       	:: temp3d_simp
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
  REAL*4 										:: dd,dmult,uint,vint,zint
  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob, depthb,latb
#endif

 LOGICAL around

!---------------------------------------------------------------

#ifdef initxyt
  alloCondGrid: if ( .not. allocated (ntimask) ) then
     allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
  end if alloCondGrid
#endif


 alloCondUVW: if(.not. allocated (temp2d_simp)) then
   allocate ( temp3d_simp(IMT,JMT,KM), temp2d_simp(IMT,JMT)  )
#ifdef tempsalt
   allocate ( tempb(KM), saltb(KM), rhob(KM), depthb(KM), latb(KM))
#endif
 end if alloCondUVW

if(.not. allocated (zstot)) then
   allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
end if
 
  call datasetswap

! === Initialising fields ===
  initFieldcond: if(ints.eq.intstart) then
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif

 !currHour=startHour
 !currDay=startDay
 !currMon=startMon
 !currYear=startYear
 !if(ngcm.le.24) then
 ! currHour=24-ngcm
 ! currDay=startDay-5
 ! currMon=startMon
 !endif
 !    
 ! else

! ----------------------------------------------------------------

   ! === Update clockworks ===
 !     currDay=currDay+nff*ngcm/24
 !     
 ! if(currDay > idmax(currMon, currYear)) then ! why 1999 and not currYear?????
 !   currDay=currDay-idmax(currMon, currYear)
 !   currMon=currMon+1
 !   if(currMon == 13) then
 !      currMon=1
 !      currYear=currYear+1
 !    if(currYear.eq.yearmax+1) currYear=yearmin
 !   endif
 ! elseif(currDay <=0) then
 !   currMon=currMon-1
 !   if(currMon == 0) then
 !      currMon=12
 !      currYear=currYear-1
 !    if(currYear.eq.yearmin-1) currYear=yearmax
 !   endif
 !   currDay=currDay+idmax(currMon, currYear)
 !     endif
 !  
   endif initFieldcond

   call updateClock

   ! === Time number ===
   ntime=10000*currYear+100*currMon+currDay
   
   ! ------------------------------------------------------------
   if (ints == intstart) then
   !! Read mean
   if(.not. allocated (u_m)) then
   allocate ( u_m(imt,jmt,km), v_m(imt,jmt,km) )
   u_m(1:imt,1:jmt,1:km) = get3DfieldNC('/group_workspaces/jasmin2/aopp/joakim/ORCA025-N401/mean/ORCA025-N401_1978-2010m00U.nc',&
                            & 'vozocrtx')
   v_m(1:imt,1:jmt,1:km) = get3DfieldNC('/group_workspaces/jasmin2/aopp/joakim/ORCA025-N401/mean/ORCA025-N401_1978-2010m00V.nc',&
                            & 'vomecrty')   
   end if
   end if
   ! === Find the file for this timestep ===
   
   !dataprefix='xxxx/ORCA025-N112_xxxxxxxx'
   dataprefix='xxxx/ORCA025-N401_xxxxxxxx'
   write(dataprefix(1:4),'(i4)') currYear
   write(dataprefix(19:26),'(i4i2.2i2.2)') currYear,currMon,currDay
   !fieldFile = trim(inDataDir)//trim(dataprefix)//'d05'
   fieldFile = trim(inDataDir)//'/means/'//trim(dataprefix)//'d05'
 
#ifdef timestat
   fieldFile = trim(inDataDir)//'/fields/ORCA025-N112_1958to2001y01'
#endif 
   
   ! Read SSH                                                                           
   hs(:,     :, nsp) = get2DfieldNC(trim(fieldFile)//'T.nc', 'sossheig')
   hs(imt+1, :, nsp) = hs(1,:,nsp)
   
   do i=4,IMT !! is this north fold needed?
      ii=IMT+4-i
      hs(i,JMT+1,2)=hs(ii,JMT-3,2)  !  north fold 
   enddo

   ! Depth at U, V, T points as 2D arrays                                                                      
   allocate ( abyst(imt, jmt) , abysu(imt, jmt) , abysv(imt, jmt) )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   abysu = sum(dzu(:,:,:,1), dim=3)
   abysv = sum(dzv(:,:,:,1), dim=3)
   
   ! Calculate SSH/depth                                                                       
   where (abyst /= 0)
      zstot = hs(:imt,:jmt,nsp)/abyst + 1
   elsewhere
      zstot = 0.d0
   end where

   where (abysu /= 0)
      zstou = 0.5*(hs(:imt,:jmt,nsp)+hs(2:imt+1,:jmt,nsp))/abysu + 1
   elsewhere
      zstou = 0.d0
   end where

   where (abysv /= 0)
      zstov = 0.5*(hs(:imt,:jmt,nsp)+hs(:imt,2:jmt+1,nsp))/abysv + 1
   elsewhere
      zstov = 0.d0
   end where

   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
   end do
   end do
   end do


#ifdef tempsalt 
   ! Read temperature  
   temp3d_simp(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'votemper')
   tem(:,:,1:km,nsp) = temp3d_simp(:,:,km:1:-1)
   ! Read salinity                                                                                                                         
   temp3d_simp(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'vosaline')
   sal(:,:,1:km,nsp) = temp3d_simp(:,:,km:1:-1)
   ! Calculate potential density                                                                       
   do j=1,jmt
      latb=-80+0.25*float(j+subGridJmin-1)
      do i=1,IMT
         tempb(:) = tem(i,j,:,nsp)
         saltb(:) = sal(i,j,:,nsp)
         call statvd(tempb, saltb, rhob ,km ,depthb ,latb)
         rho(i,j,:,nsp)=rhob(:) - 1000.
      end do
   end do
#endif     
   
   dmult=1.  ! amplification of the velocity amplitude by simple multiplication
   uvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(trim(fieldFile)//'U.nc', 'vozocrtx')
   !print*,'uvel 1:5',uvel(867,664,1:5)
   !print*,'uvel 70:75',uvel(867,664,70:75)
   !print*,'u_m 1:5',u_m(867,664,1:5)
   !print*,'u_m 70:75',u_m(867,664,70:75)
   !uvel(1:imt,1:jmt,1:km) = uvel(1:imt,1:jmt,1:km) - u_m(1:imt,1:jmt,km:1:-1)
   !print*,'uvel 1:5',uvel(867,664,1:5)
   !print*,'uvel 70:75',uvel(867,664,70:75)
   !stop
   !! calculate volume fluxes                                                                                                 
             
   uflux(:,:,:,nsp) = 0.
   print*,u_m(648,118,3)
   do k = 1, km      
   do j = 1, jmt
   do i = 1, imt
      !if (i==648 .and. j==118) print*,'uvel',uvel(i,j,km+1-k),u_m(i,j,k)
      !print*,'uvel',uvel(i,j,km+1-k),u_m(i,j,k),i,j,k,size(u_m,1)
      uvel(i,j,km+1-k) = uvel(i,j,km+1-k) - u_m(i,j,k)                                                        
      !if (i==648 .and. j==118) print*,'uvel',uvel(i,j,km+1-k) 
      uflux(i,j,km+1-k,nsp) = uvel(i,j,km+1-k) * dyu(i,j) * dzu(i,j,km+1-k,1) * zstou(i,j)
   enddo
   enddo
   enddo

   !  north fold 
   do i=4,IMT
      ii=IMT+4-i
      ! vflux(i,JMT,:,2)=-vflux(ii,JMT-3,:,2)
   enddo
   
   vvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(trim(fieldFile)//'V.nc', 'vomecrty')
   
   !! calculate volume fluxes                                                                                                 
   vflux(:,:,:,nsp) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      vvel(i,j,km+1-k) = vvel(i,j,km+1-k) - v_m(i,j,k)
      vflux(i,j,km+1-k,nsp) = vvel(i,j,km+1-k) * dxv(i,j) * dzv(i,j,km+1-k,1) * zstov(i,j)
   enddo
   enddo
   enddo
  
   !! calculate laplacian of u,v                                   
   call laplacian 
   
   
#ifdef drifter
! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
kbot=59 ; ktop=60 ! number of surface layers to integrat over
do i=1,imt
 do j=1,jmt
 
  uint=0. ; vint=0. ; zint=0.
  if(ktop.eq.KM) zint=hs(i,j,2)
  do k=kbot,ktop
   uint=uint+uflux(i,j,k,2) ! integrated transport
   vint=vint+vflux(i,j,k,2)
   zint=zint+dz(k)          ! total depth of drougued drifter
  enddo
  ! weighted transport for each layer
  do k=kbot,KM
   if(k.ne.KM) then
    uflux(i,j,k,2)=uint*dz(k)/zint 
    vflux(i,j,k,2)=vint*dz(k)/zint
   else
    uflux(i,j,k,2)=uint*(hs(i,j,2)+dz(k))/zint
    vflux(i,j,k,2)=vint*(hs(i,j,2)+dz(k))/zint
   endif
  enddo

 enddo
enddo

#endif


 deallocate ( temp3d_simp, temp2d_simp )
     
#ifdef tempsalt
 deallocate ( tempb, saltb, rhob, depthb, latb )
#endif

  return
  
end subroutine readfields



