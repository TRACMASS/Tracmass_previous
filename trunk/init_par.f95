subroutine init_params
  USE mod_param
  USE mod_coord
  USE mod_grid
  USE mod_name
  USE mod_time 
  USE mod_domain
  USE mod_vel
  USE mod_dens
  USE mod_buoyancy
  USE mod_streamxy
  USE mod_streamv
  USE mod_streamr
  USE mod_tracer
#ifdef sediment
  USE mod_orbital
#endif
!  implicit none
  
  CHARACTER(LEN=23) ::  Project, Case
  namelist /INITGRIDGRID/ IMT, JMT, KM, JMAX, LBT,NEND
  namelist /INITGRIDNTRAC/ NTRACMAX
  namelist /INITGRIDDATE/ yearmin, yearmax
  namelist /INITGRIDTIME/ ngcm, iter, intmax

  namelist /INITRUNTIME/ intmin, intspin, intrun, intstep 
  namelist /INITRUNDATE/  ihour, iday, imon, iyear
  namelist /INITRUNWRITE/ ncoor, kriva, directory, name
  namelist /INITRUNSEED/ nff, isec, idir, nqua, partQuant, kst1, kst2, ist1, ist2, jst1, jst2, kst1, kst2

  Project = PROJECT_NAME
  Case    = CASE_NAME

  print *,trim(Project)//'/'//trim(Case)//'_run.in'
  stop
  open(8,file=trim(Project)//'/'//trim(Project)//'_grid.in',  &
       status='OLD', delim='APOSTROPHE')
  read(8,nml=INITGRIDGRID)
  read(8,nml=INITGRIDNTRAC)
  read(8,nml=INITGRIDTIME)
  read(8,nml=INITGRIDDATE)

  open(8,file=trim(Project)//'/'//trim(Case)//'_run.in',  &
       status='OLD', delim='APOSTROPHE')
  read(8,nml=INITRUNTIME)
  read(8,nml=INITRUNDATE)
  read(8,nml=INITRUNWRITE)  
  read(8,nml=INITRUNSEED)


  dstep=1.d0/dble(iter)
  dtmin=dtstep*tseas

  ! --ist -1 to imt
  if ( ist1 == -1) then 
     ist1=imt
  end if
  if ( ist2 == -1) then 
     ist2=imt
  end if
  ! --jst -1 to jmt
  if( jst1 == -1) then 
     jst1=jmt
  end if
  if ( jst2 == -1) then 
     jst2=jmt
  end if
  ! --ist -1 to imt  
  if ( kst1 == -1) then 
     kst1=km
  end if
  if ( kst2 == -1) then 
     kst2=km
  end if

  ! mod_coord
  allocate  ( csu (jmt), cst(jmt), zw(0:km) )       
  ! mod_grid
  allocate ( dxdy(imt,jmt), dztb(imt,jmt,kd) )   
  allocate (kmt(imt,jmt), dz(km) )
  ! mod_domain
  allocate ( ienw (nend),iene (nend) )
  allocate ( jens (nend),jenn (nend) )
  allocate ( mask (imt,jmt) )
  ! mod_vel
  allocate ( u(imt,0:jmax,km,nst), v(imt,0:jmax,km,nst) )
  allocate ( hs(imt,jmax,nst) )
  allocate ( w(0:km) )

  ! mod_dens
#ifdef tempsalt
  allocate ( tem(imt,jmax,km,nst) ,sal(imt,jmax,km,nst), rho(imt,jmax,km,nst) )
#endif
  ! mod_streamxy
#ifdef streamxy
  allocate ( stxyy(imt,jmt,lbt), stxyx(imt,jmt,lbt) )
#endif
  ! mod_streamv
#ifdef streamv
  allocate ( stxz(imt,km,lbt), styz(imt,km,lbt) )
#endif
  ! mod_streamr
#ifdef streamr
  allocate ( stxr(imt,mr,lbt,lov), styr(jmt,mr,lbt,lov) )
#endif
! mod_tracer
#ifdef tracer
  allocate ( tra(imt,jmt,km) )
#endif
  ! mod_sed
#ifdef sediment
  allocate (orb(km) )
#endif

end subroutine init_params

