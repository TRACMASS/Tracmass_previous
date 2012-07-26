!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
subroutine readfields

  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_stat
  
#ifdef tempsalt
  USE mod_dens
#endif

  IMPLICIT none
  
#ifdef tempsalt
  REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob
  integer kmm,got_tem,got_sal
#endif
  
!  integer NLEN,NSNAPS,NSNAPS_2D,NSNAPS_3D,nd,l_tke,nlength
  
!  parameter (NLEN=29557,NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,NSNAPS=1,nlength=2379)
!  parameter (l_tke=0)

  INTEGER, PARAMETER :: NLEN=29557,NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,NSNAPS=165
  INTEGER, PARAMETER :: nlength=2379,l_tke=0
  
  INTEGER :: itt0,year,month,day,hour,minute,second,ittstart,itt,dtts
  INTEGER :: imt0,jmt0,km0,nt0,NLEN0,NSNAPS0,nvars,got_hs,got_u,got_v
  INTEGER ::  i,j,k,m,kz,ii,ints2,kk,i0,ktop,kbot
!  INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:,:) :: kmu

  REAL*8 :: ird0,ird20,ird30,ird40,stlon,stlat
  REAL*4 :: snapd,totsec,ird,ird2,ird3,ird4
  REAL*4, ALLOCATABLE, DIMENSION(:)   :: rd1d_a, rd1d_b,zdzz,dzw,dxt,ispvar,isplev
  REAL*4, ALLOCATABLE, DIMENSION(:)   :: phit,yu,snap1d
  REAL*4, ALLOCATABLE, DIMENSION(:,:) :: rd2d
  REAL*4 :: uint,vint,zint

  REAL :: snap2d(imt,jmt) ! ??????????????????

  CHARACTER ofile*100,infile*100,zfile*193,rfile*99,a_exp1*3,a_exp2*2

  LOGICAL around
  
  REAL*8, SAVE :: dxa,dya

  alloCondGrid: if(.not. allocated (snap1d)) then
     allocate ( snap1d(NLEN),rd2d(IMT,JMT) )
     !allocate ( rd1d_a(NSNAPS),rd1d_b(NSNAPS) ,ispvar(NSNAPS), isplev(NSNAPS))
     allocate ( zdzz(KM),dzw(0:km),dxt(imt) ) 
     allocate ( phit(jmt),yu(jmt) )
#ifdef tempsalt
     allocate ( tempb(KM), saltb(KM), rhob(KM) )
#endif
  end if alloCondGrid

 if ( .not. allocated (kmu) ) then
     allocate ( kmu(imt,jmt) )
  end if

 
!_______________________ update the time counting ________________________________________
ihour=ihour+6
if(ihour.eq.24) then
 ihour=0
 iday=iday+1
 if(iday.gt.idmax(imon,iyear)) then
  iday=1
  imon=imon+1
!  if(imon.eq.9) imon=7
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
   if(iyear.gt.yearmax) iyear=yearmin ! recycle over gcm outputdata
  endif
 endif
endif
ntime=1000000*iyear+10000*imon+100*iday+ihour

!print *,'tiden=',ntime,iyear,imon,iday,ihour,iyear0
!____________________________ initialise ___________________________
if(ints.eq.intstart) then
hs=0.
uflux=0.
vflux=0.
#ifdef tempsalt
tem=0.
sal=0.
rho=0.
#endif
iyear=startyear
imon=startMon
iday=startDay
ihour=0

endif
ntime=1000000*iyear+10000*imon+100*iday+ihour

!_______________________________________________________________________

! swap between datasets

    ! === swap between datasets ===
    hs(:,:,1)=hs(:,:,2)
    uflux(:,:,:,1)=uflux(:,:,:,2)
    vflux(:,:,:,1)=vflux(:,:,:,2)

        !Save the velocities to file, this is another place to save vel. data, see alo in loop.f95 file.
        !open(unit=1122,file='velmatrix.dat')
        !do k=1,100
        !	write(1122,*) uflux(:,k,1,1),vflux(:,k,1,1)
        !enddo
        !close(1122)

#ifdef explicit_w
    wflux(:,:,:,1)=wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)=tem(:,:,:,2)
    sal(:,:,:,1)=sal(:,:,:,2)
    rho(:,:,:,1)=rho(:,:,:,2)
#endif


  !print *,'hhhhhhhhhhh',ints
                     
!     === Create filenames for the snap-files to be used ===
ofile='0000/d0000000000.snap1'
write(ofile(1:4), '(i4) ') iyear
write(ofile(7:16),'(i10)') ntime

infile=trim(inDataDir)//ofile

inquire(file=trim(infile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',infile,ntime
stop 4555
endif
zfile='gzip -c -d '//trim(infile)//'.gz > '//trim(inDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(inDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556



!     === open snap file ===
!open(unit=30,file=trim(rfile),status='old',form='unformatted',err=4000,convert='little_endian')
open(unit=30,file=trim(rfile),status='old',form='unformatted',err=4000,convert='big_endian')

!print *,'We start reading snap-file...'
! start read header:
read(30) ird
!Timestep counter
itt0 = ird
read(30) ird
!Number of depth levels
km0 = ird
!print *,'km0',ird
if(km0.ne.km) stop 2844
read(30) ird
!print *,'ird',ird
!Number of tracers,ignore 
nt0 = ird
read(30) ird
!This is y size(362)
imt0 = ird
!print *,'imt0',ird
if(imt0.ne.imt) stop 2845
read(30) ird
!print *,'ird',ird
!This is x size (320)
jmt0 = ird
if(jmt0.ne.jmt) stop 2846
read(30) ird
!Maximum number of active points on a single level.
!Can be used to allocate storage in this routine for snap1d
NLEN0 = ird
read(30) ird
!Number of fields in this file (counting each depth level as a field)
NSNAPS0 = ird
read(30) ird 
year = ird
read(30) ird
month = ird
read(30) ird
day = ird
read(30) ird
hour = ird
read(30) ird
minute = ird
read(30) ird
second = ird

!print *, 'Header is:itt0,km0,nt0,imt0,jmt0,NLEN0,NSNAPS0,year,month,day,hour,minute,second'
!print *, itt0,km0,nt0,imt0,jmt0,NLEN0,NSNAPS0,year,month,day,hour,minute,second

 if ( .not. allocated (rd1d_a) ) then
     allocate ( rd1d_a(NSNAPS0), rd1d_b(NSNAPS0) ,ispvar(NSNAPS0), isplev(NSNAPS0) )
  end if

read(30) ird0,ird20,ird30
!Timestep in second
dtts = ird0
!Number of seconds since referenc time (chosen in RCO at compile time)
totsec = ird20
!Time in days between snap-files
snapd = ird30
!print *,'dtts',ird0,ird20,ird30
read(30) ird0,ird20,ird30,ird40
!dx = ird0 !model grid in cm
!dy = ird20
dxdeg = ird30
dydeg = ird40
!print *,'dxdeg',ird0,ird20,ird30,ird40 
read(30,err=2000) ird0,ird20
! ird0  = stlon
! ird20 = stlat
! Lat and lon for the reference point of the grid.  This is half a grid
! square to the southwest of the first point in the t-grid or a full
! gridsquare to the southwest of the u-grid

read(30,err=2000) rd1d_a,rd1d_b

do i=1,NSNAPS0
        !print *,'i, rd1d_a(i),rd1d_b(i)=',i,rd1d_a(i),rd1d_b(i)
        !Parameter number of the field in this file
        ispvar(i) = rd1d_a(i)
        !Level of each field in this file
        !2d-fields have level = 1
        isplev(i) = rd1d_b(i)
enddo

read(30,err=2000) rd2d

if(ints.eq.intstart) then

stlon1 = ird0
stlat1 = ird20
!They are longitude and lattitude values
!dya=0.005*dydeg
!dxa=0.005*dxdeg

!kmt = number of depth levels in each grid cell on the t-grid.
do i=1,imt
 do j=1,jmt
  kmt(i,j) = rd2d(i,j)
 enddo
enddo
!print*,kmt


!open(unit=747,file='/Users/doos/data/rco/topo/kmt-new',form='unformatted')
!write(unit=747) kmt
!close(747)

!kmu = number of depth levels in each grid cell on the u-grid.
 do i=1,imt-1
  do j=1,jmt-1
   kmu(i,j)=min(kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
  enddo
 enddo

call coordinat

dya=0.005*dy*deg  
dxa=0.005*dx*deg

endif

!print *,'dx=',dxa,dya,dx,dy,deg

nvars = 0
got_hs = 0
got_u = 0
got_v = 0
#ifdef tempsalt
got_tem = 0
got_sal = 0
#endif

!Read the field until we found them all.
!For each field the number of active (wet) points is stored first.
!If this number is nonzero, then the active points is stored as a 1d-array.

do while (got_hs < 1 .or. got_u < km .or. got_v < km &
#ifdef tempsalt
        .or. got_tem < km .or. got_sal < km &
#endif
        )
        nvars = nvars + 1
        !print *, 'ispvar(nvars), # of variables',ispvar(nvars),nvars
        if (ispvar(nvars) == 1) then !Free surface height
        ! ssh - sea surface height
        !print *,'ssh'
        got_hs = 1 

        read(30,err=2000) ird

        i0=nint(ird)
        !print *,'i0',ird,i0
        !print *,'snap1d',snap1d
        if (i0 == 0) stop 'Empty ssh field'

        read(30,err=2000) snap1d(1:i0)


        ii=0
                do j=1,jmt
                        do i=1,imt
                          hs(i,j,2)=0. !Missing value (dry point)
                          if(kmt(i,j).ge.1) then
                           ii=ii+1
                           hs(i,j,2)=0.01*snap1d(ii)
                          !print *,'ssh,ii=', 0.01*snap1d(ii),ii
      !   if(kmt(i,j).ne.41) print *,ii,i,j,kmt(i,j),hs(1,j,2)
                        endif
                        enddo
                enddo

!hs=0.

#ifdef tempsalt
   else if (ispvar(nvars) == 4) then  ! Potential temperature
      got_tem = got_tem + 1

      !print *,'Temperature'
      kk=km+1-isplev(nvars)
      read(30,err=2001) ird 
      i0=nint(ird)
      if(i0.gt.0) then
         read(30,err=2001) snap1d(1:i0)


         ii=0
         do j=1,jmt
            do i=1,imt
               tem(i,j,kk,2)=0.  ! Missing value (dry point)
               if(kmt(i,j).ge.isplev(nvars)) then
                  ii=ii+1
                  tem(i,j,kk,2)=snap1d(ii)
                  !print *,'tem,ii=',snap1d(ii),ii
               endif
            enddo
         enddo
      else
         tem(:,:,kk,2) = 0.0  ! Missing value (dry point/level)
      endif

   else if (ispvar(nvars) == 5) then  ! Salinity
      got_sal = got_sal + 1

      !print *,'Salinity'
      kk=km+1-isplev(nvars)
      read(30,err=2002) ird
      i0=nint(ird)
      if(i0.gt.0) then
         read(30,err=2002) snap1d(1:i0)

         ii=0
         do j=1,jmt
            do i=1,imt
               sal(i,j,kk,2)=0.   ! Missing value (dry point)
               if(kmt(i,j).ge.isplev(nvars)) then
                  ii=ii+1
                  sal(i,j,kk,2)=snap1d(ii)
                  !print *,'sal,ii=', snap1d(ii),ii
               endif
            enddo
         enddo
      else
         sal(:,:,kk,2) = 0.0  ! Missing value (dry point/level)
      endif
#endif

 else if (ispvar(nvars) == 8) then ! u-velocity
      !Reads zonal velocity values for each depth level
      ! zonal velocity

      got_u = got_u + 1
        !print *,' U-velocity'
      kk=km+1-isplev(nvars)
      read(30,err=2003) ird
      i0=nint(ird)
      if(i0.gt.0) then
         read(30,err=2003) snap1d(1:i0)
         ii=0
         do j=1,jmt
            do i=1,imt
               snap2d(i,j)=0.     ! Missing value (dry point)
               if(kmu(i,j).ge.isplev(nvars)) then
                  ii=ii+1
                  snap2d(i,j)=snap1d(ii)
                  !print *,'u-vel,ii=', snap1d(ii),ii
               endif
            enddo
         enddo

         ! u -> transports
         do j=2,jmt
            do i=1,imt-1
               if(kk.ne.km) then
                  uflux(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))*dz(kk)
               else
                  uflux(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))* &
                       ( dz(kk)+0.5*(hs(i,j,2)+hs(i+1,j,2)) )
               endif
            enddo
         enddo
      else
         uflux(:,:,kk,2) = 0.0    ! Missing value (dry point/level)
      endif

   else if (ispvar(nvars) == 9) then ! v-velocity
      !Reads meridional velocity values for each depth level
      ! meridional velocity

      got_v = got_v + 1
      !print *,'V-velocity'
      kk=km+1-isplev(nvars) !kk 41...1, isplev(nvars)=1..41
      read(30,err=2004) ird
      i0=nint(ird)
      if(i0.gt.0) then
         read(30,err=2004) snap1d(1:i0)
         ii=0
         do j=1,jmt
            do i=1,imt
               snap2d(i,j)=0.    ! Missing value (dry point)
               if(kmu(i,j).ge.isplev(nvars)) then
                  ii=ii+1
                  snap2d(i,j)=dxa*csu(j)*snap1d(ii)
                  !print *,'v-vel,ii,km,kk,i0',dxa*csu(j)*snap1d(ii),ii,km,kk,i0
               endif
            enddo
         enddo
         ! v -> transport
         do j=1,jmt-1
            do i=2,imt
               if(kk.ne.km) then
                  vflux(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*dz(kk)
               else
                  vflux(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*&
                       ( dz(kk)+0.5*(hs(i,j,2)+hs(i,j+1,2)) )
!                       if(dz(kk)-0.5*(hs(i,j,2)+hs(i,j+1,2)).lt.0.) print *,dz(kk),0.5*(hs(i,j,2)+hs(i+1,j,2))
               endif
            enddo
         enddo
      else
         vflux(:,:,kk,2) = 0.0    ! Missing value (dry point/level)
      endif

   else
      ! Skip field

      !print *,'! End of file here means some fields are missing in this file.'
      read(30,err=2000) ird
      i0=nint(ird)
      if (i0 /= 0) then
         read(30,err=2000) snap1d(1:i0)
      endif
   endif

enddo


goto 9000

2001  continue
2002  continue
2003  continue
2004  continue
    
4000  continue
2000  print*,'error 2000'
stop 2000
9000 continue


close(30)  


#ifdef tempsalt
! the density
do i=1,imt
 do j=1,jmt
  if(kmt(i,j).ne.0) then
   kmm=kmt(i,j)
   do k=1,kmm
    kk=km+1-k
    tempb(k)=tem(i,j,kk,2)
!    saltb(k)=(sal(i,j,kk,2)-35.)/1000.
    saltb(k)=sal(i,j,kk,2)
    if(saltb(k).lt.0.) saltb(k)=0.
    enddo
    call statv(tempb,saltb,rhob,kmm)
    do k=1,kmm
     kk=km+1-k
     rho(i,j,kk,2)=rhob(k)
    enddo
   endif
  enddo
 enddo
#endif

#ifdef drifter
! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
kbot=36 ; ktop=37 ! number of surface layers to integrat over

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
  do k=KM-kk,KM
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

uflux(:,:,:,2) = uflux(:,:,:,2)
vflux(:,:,:,2) = vflux(:,:,:,2)

!deallocate ( snap1d, rd2d )
!deallocate ( rd1d_a, rd1d_b )
!deallocate ( zdzz,dzw,dxt )
!deallocate ( phit, yu )
!deallocate ( tempb, saltb, rhob )
!deallocate ( kmu )


!print *,'readfield slut',ints
!if(ints.gt.1) stop 49678

!print *,'...end of reading snap-file.'
return
end subroutine readfields

!________________________________________________________________________________________
      
