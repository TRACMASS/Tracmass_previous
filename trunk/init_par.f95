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
  implicit none

  INTEGER :: argint1, argint2, dummy, factor, i,dtstep
  CHARACTER (LEN=30) :: inparg, argname
  CHARACTER (LEN=23) ::  Project, Case

  namelist /INITGRIDGRID/ IMT, JMT, KM, JMAX, LBT,NEND
  namelist /INITGRIDNTRAC/ NTRACMAX
  namelist /INITGRIDDATE/ yearmin, yearmax
  namelist /INITGRIDTIME/ ngcm, iter, intmax
  
  namelist /INITRUNTIME/ intmin, intspin, intrun, intstep 
  namelist /INITRUNDATE/  ihour, iday, imon, iyear
  namelist /INITRUNWRITE/ ncoor, kriva, directory, name
  namelist /INITRUNSEED/ nff, isec, idir, nqua, partQuant, kst1, kst2, ist1, ist2, jst1, jst2, kst1, kst2
#ifdef tempsalt
  namelist /INITRUNTEMPSALT/ tmin0, tmax0, smin0, smax0, rmin0, rmax0, &
                             tmine, tmaxe, smine, smaxe, rmine, rmaxe
#endif


!  allocate ( ienw (LBT),iene (LBT) )
!  allocate ( jens (LBT),jenn (LBT) )

  namelist /INITRUNEND/ ienw, iene, jens, jenn, timax


  Project = PROJECT_NAME
  Case    = CASE_NAME
  
  if ( (IARGC() .eq. 1 ) .or. (IARGC() .eq. 4 ) )  then
     call getarg(IARGC(),inparg)
     Case=inparg
  end if
  
 ! -- Check if there is a time argument and if so, use it.
  
  print *,trim(Project)//'/'//trim(Case)//'_run.in'
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
#ifdef tempsalt
  read(8,nml=INITRUNTEMPSALT)
#endif
  read(8,nml=INITRUNEND)
print *,'ienw',ienw
print *,'iene',iene
print *,'jens',jens
print *,'jenn',jenn
print *,'timax',timax
  timax=24.*3600.*timax ! convert time lengths from days to seconds

  dstep=1.d0/dble(iter)
  dtmin=dtstep*tseas

  if ((IARGC() > 1) .and. (IARGC() < 5) )  then
     call getarg(1,inparg)
     factor=1
     argint1=0
     do i=29,1,-1
        if (ichar(inparg(i:i)) .ne. 32) then
           argint1=argint1+(ichar(inparg(i:i))-48)*factor
           factor=factor*10
        end if
     end do
     ARG_INT1=argint1

     call getarg(2,inparg)
     factor=1
     argint2=0
     do i=29,1,-1
        if (ichar(inparg(i:i)) .ne. 32) then
           argint2=argint2+(ichar(inparg(i:i))-48)*factor
           factor=factor*10
        end if
     end do
     ARG_INT2=argint2

     call getarg(3,inparg)
     name=inparg
  end if
  
  
  print *,intmin,intrun,name

!stop

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
!  allocate ( ienw (LBT),iene (LBT) )
!  allocate ( jens (LBT),jenn (LBT) )
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
  allocate ( stxz(imt,km,lbt), styz(jmt,km,lbt) )
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

