

program main

  integer                                  :: ierr ,recPosRun
  integer                                  :: ntrac ,ints
  

  open(unit=70 ,file=trim(fullWritePref)//'_run.bin' &  ! Trajectory path  
       ,access='direct' ,form='unformatted' ,recl=20 ,status='replace') 

  open(unit=71 ,file=trim(fullWritePref)//'_run_//.bin' &
       ,access='direct' ,form='unformatted' ,recl=20 ,status='replace') 


  do while (ierr .ne. 0)
     recPosRun = recPosRun+1

     if (recPosRun/1e+6 == floor(recPosRun/1e+6)) then
        close (unit=71)



     end if
     
     read  (unit=70 ,rec=recPosRun) ntrac,ints,x14,y14,z14
     write (unit=71 ,rec=recPosRun) ntrac,ints,x14,y14,z14

  end do
