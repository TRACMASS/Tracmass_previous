!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine readfields
USE HDF5
USE OC5
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

INTEGER ints2,i,ip,im,j,jp,k,jm,kmu(IMT,JMT),kk,kmm
CHARACTER day*4

#ifdef fiveday
CHARACTER day2*4,dar*10
INTEGER iday2
#endif

REAL fd1(IMT,JMT),saltb(km),tempb(km),rhob(km)
REAL tp(KM,IMT),sp(KM,IMT),vp(KM,IMT),up(KM,IMT),ssh(IMT,JMT),um(KM,IMT),dus,dun,dvw,dve

CHARACTER(LEN=80)  :: filename1
INTEGER (HID_T)    :: file_in(2), gr_in(2), data_in(2)
INTEGER (HID_T)    :: sdsidt_1, sdsids_1, sdsidu_1, sdsidv_1
INTEGER   :: ivalues(8), ifmode, iftype, idomain
REAL    ::  vmask(4)
LOGICAL:: around
save kmu,fd1

! ta bort detta kanske och ersätt med längre ner
hs(:,:,1)=hs(:,:,2)
u(:,:,:,1)=u(:,:,:,2)
v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
tem(:,:,:,1)=tem(:,:,:,2)
sal(:,:,:,1)=sal(:,:,:,2)
rho(:,:,:,1)=rho(:,:,:,2)
#endif

!____________________________ initialise ___________________________
if(ints.eq.intstart) then
 hs=0.
 u=0.
 v=0.
#ifdef tempsalt
 tem=0.
 sal=0.
 rho=0.
#endif
open(45,file='../../tmp/'//name)
 98 format(i4)
 99 format(4a)

 open(86,file='/Volumes/sjo2/data/occam/p250/topo/kmt1',form='unformatted')
 read(86) kmt
 close(86)
 open(86,file='/Volumes/sjo2/data/occam/p250/topo/kmu1',form='unformatted')
 read(86) kmu
 close(86)
 open(86,file='/Volumes/sjo2/data/occam/p250/topo/fd1',form='unformatted')
 read(86) fd1
 close(86)

 do j=1,JMT
  do i=1,IMT
   dztb(i,j,1)=0.d0
   if(fd1(i,j).eq.100000.) fd1(i,j)=0.
  enddo
 enddo

!print *,zw

! bottom grid box thickness on T-points
 do j=1,JMT
  jm = j-1
  if(jm.eq.0) jm = 1   
  do i=1,IMT
   im = i-1
   if(im.eq.0) im = IMT
   if(kmt(i,j).gt.1) dztb(i,j,1)=0.25*( &
                            fd1(i,j )-zw(kmu(i,j )-1) + fd1(im,j )-zw(kmu(im,j )-1) &
                           +fd1(i,jm)-zw(kmu(i,jm)-1) + fd1(im,jm)-zw(kmu(im,jm)-1) )
   if(dztb(i,j,1).lt.0.) then
    print *,i,im,j,jm,dztb(i,j,1),fd1(i,j ),zw(kmu(i,j )-1)
    stop 2875
   endif
  enddo
 enddo

endif

!______________________  Time or stat   _________________________
!#if defined stat     
!open(87,file='/Volumes/sjo2/data/occam/huvry/huvr.9102',form='unformatted')
!read(87) huvr
!close(87)
!#endif     

#ifdef fiveday     
      ints2=ints
  666  continue
      if(ints2.gt.intmax .or. ints2.lt.intmin) then
       ints2=ints2-intend+intstart-intstep
       goto 666
      endif

rewind(45)
write(45,98) ints2
rewind(45)
read (45,99)  day
iday2=ints2+5
rewind(45)
write(45,98) iday2
rewind(45)
read (45,99)  day2      
dar=day//'to'//day2

call oc5init
filename1='/Volumes/sjo2/data/occam/p25k66r202/5day/d'//dar//'.h5m1'
!write(6,*) filename1
inquire(file = filename1, exist = around )
if(.not.around) then
 write(6,*) filename1, ': file does not exist'
 stop 4965
endif

call oc5open_rd(filename1,file_in(1),gr_in(1), ifmode, iftype, idomain, vmask )
call oc5select('SEA SURFACE HEIGHT (MEAN)',gr_in(1),data_in(1))
call oc5getksb(data_in(1),0,1,IMT,1,JMT,ssh)
call oc5endacc(data_in(1))
do j=1,JMAX
 do i=1,IMT
  hs(i,j,2)=0.01*ssh(i,j) ! sea surface from cm to m
 enddo
enddo
!c  read temp, salt and vel
call oc5select('POTENTIAL TEMPERATURE (MEAN)',gr_in(1),sdsidt_1)
call oc5select('SALINITY (MEAN)',             gr_in(1),sdsids_1)
call oc5select('U VELOCITY (MEAN)',           gr_in(1),sdsidu_1)
call oc5select('V VELOCITY (MEAN)',           gr_in(1),sdsidv_1)

!      if(iday.eq.idaystart) then
!c depth
!call oc5getax(file_in(1), sdsidt_1, 1,1,KM,zw)
!print *,zw

do 1000 j=1,JMAX
um=up ! swap from north to south row for u
jp = j+1
if(jp.gt.JMAX) jp = JMAX
jm = j-1
if(jm.eq.0) jm = 1   ! borde inte detta vara ngt annat?
call oc5getjsb(sdsidt_1,j,1,KM,1,IMT,tp)
call oc5getjsb(sdsids_1,j,1,KM,1,IMT,sp)
call oc5getjsb(sdsidu_1,j,1,KM,1,IMT,up) ! in cm/s
call oc5getjsb(sdsidv_1,j,1,KM,1,IMT,vp) ! in cm/s
do 900 i=1,IMT
 im = i-1
 if(im.eq.0) im = IMT
 ip = i+1
 if(ip.eq.IMT+1) ip = 1
 kmm=kmt(i,j)
 do 800 k=1,kmm
  kk=km+1-k
  if(k.eq.1) then
   dun=dz(km)+0.5*(hs(i,j,2)+hs(ip,j ,2))
   dve=dz(km)+0.5*(hs(i,j,2)+hs(i ,jp,2))
   dus=dun
   dvw=dve
  else
   if (k.eq.kmu(i,j)) then    
    dun = fd1(i,j) - zw(k-1)
   else
    dun = dz(kk)
   endif
   dve = dun
   if (k.eq.kmu(im,j)) then    
    dvw = fd1(im,j) - zw(k-1)
   else
    dvw = dz(kk)
   endif
   if (k.eq.kmu(i,jm)) then    
    dus = fd1(i,jm) - zw(k-1)
   else
    dus = dz(kk)
   endif
  endif
if(dve.le.0.) then
print *,dve,i,j,jp,k,kk
print *,dz(km),hs(i,j,2),hs(i ,jp,2)
stop 3865
endif

  u(i,j,kk,NST)=0.005*dy*deg*       (dun*up(k,i)+dus*um(k,i ))  ! u transport in m3/s
  v(i,j,kk,NST)=0.005*dx*deg*csu(j)*(dve*vp(k,i)+dvw*vp(k,im))  ! v transport in m3/s
#ifdef tempsalt 
  tem(i,j,kk,2)=tp(k,i)
  sal(i,j,kk,2)=35.+1000.*sp(k,i)
  tempb(k)=tp(k,i)
  saltb(k)=(sp(k,i)-35.)/1000.
  saltb(k)=sp(k,i)
#endif

800 continue

!stop 9457
! temperature, salinity and density
! kmm=kmt(i,j)
! if(kmm.ne.0) then
!  do k=1,kmm
!   kk=km+1-k
!   tem(i,j,kk,1)=tem(i,j,kk,2)
!   sal(i,j,kk,1)=sal(i,j,kk,2)
!   rho(i,j,kk,1)=rho(i,j,kk,2)

!   tem(i,j,kk,2)=tp(k,i)
!   sal(i,j,kk,2)=35.+1000.*sp(k,i)
!   tempb(k)=tp(k,i)
!   saltb(k)=(sp(k,i)-35.)/1000.
!   saltb(k)=sp(k,i)

!   if(tempb(k).lt.-10. .or. tempb(k).gt.35.) then
!    print *,'hoppsan'
!    print *,kmm,i,j,k,tp(k,i),sp(k,i),up(k,i),vp(k,im)
!   endif
!  enddo
 if(kmm.ne.0) then
  call statv(tempb,saltb,rhob,kmm)
  do k=1,kmm
   kk=km+1-k
   rho(i,j,kk,2)=rhob(k)
  enddo

!  print *,kmm,i,j
!  print *,(tempb(k),k=1,kmm)
!  print *,(saltb(k),k=1,kmm)
!  print *,(rhob(k),k=1,kmm)
!         stop 2356
endif




  900 continue
 1000 continue



! End access to the file
call oc5endacc(sdsidt_1)
call oc5endacc(sdsids_1)
call oc5endacc(sdsidu_1)
call oc5endacc(sdsidv_1)
call oc5close(gr_in(1), file_in(1))



if(mod(ints,365).eq.0) then
 print *,'psi written for ints=',ints,ints2
 call writepsi
 endif
#endif     



return
end subroutine readfields      
