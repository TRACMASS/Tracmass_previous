
module mod_write

  USE mod_time, only: intstart,ints

  IMPLICIT NONE
  INTEGER                                    :: intminInOutFile
  CHARACTER(LEN=200)                         :: outDataDir, outDataFile
  INTEGER                                    :: twritetype = 0
  INTEGER                                    :: fileseq = 0

CONTAINS


  subroutine open_outfiles

    IMPLICIT NONE
    CHARACTER(LEN=200)                         :: fullWritePref
    CHARACTER(LEN=20)                          :: WriteStamp

    if (intminInOutFile.eq.2) then
       writeStamp='00000000_000000'
       write (writeStamp(1:8),'(i8.8)') intstart !,ints-intstart
       write (writeStamp(10:15),'(i6.6)') max(ints-intstart,0)+1
       fullWritePref =  trim(outDataDir)//trim(outDataFile)//trim(writeStamp)
    elseif (intminInOutFile.eq.1) then
       writeStamp='00000000'
       write (writeStamp,'(i8.8)') intstart
       fullWritePref =  trim(outDataDir)//trim(outDataFile)//trim(writeStamp)
    else
       fullWritePref =  trim(outDataDir)//trim(outDataFile)
    end if

#if defined textwrite
    open(56,file=trim(fullWritePref)//'_run.asc')    
    open(57,file=trim(fullWritePref)//'_out.asc')  
    open(58,file=trim(fullWritePref)//'_ini.asc')   
    open(59,file=trim(fullWritePref)//'_err.asc')
#endif
#if defined binwrite
    open(unit=75 ,file=trim(fullWritePref)//'_out.bin' &  
         ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
    open(unit=76 ,file=trim(fullWritePref)//'_run.bin' &  
         ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
    open(unit=77 ,file=trim(fullWritePref)//'_kll.bin' &
         ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
    open(unit=78 ,file=trim(fullWritePref)//'_ini.bin' &  
         ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
    open(unit=79 ,file=trim(fullWritePref)//'_err.bin' &  
         ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
#endif
#if defined csvwrite
    open(unit=85 ,file=trim(fullWritePref)//'_out.csv', status='replace')
    open(unit=86 ,file=trim(fullWritePref)//'_run.csv', status='replace')
    open(unit=87 ,file=trim(fullWritePref)//'_kll.csv', status='replace')
    open(unit=88 ,file=trim(fullWritePref)//'_ini.csv',  status='replace')
    open(unit=89 ,file=trim(fullWritePref)//'_err.csv', status='replace')
#endif
  end subroutine open_outfiles

  subroutine close_outfiles
#if defined textwrite
    close(56)
    close(57)
    close(58)
    close(59)
#endif
#if defined binwrite
    close(75)
    close(76)
    close(77)
    close(78)
    close(79)
#endif
#if defined csvwrite
    close(75)
    close(76)
    close(77)
    close(78)
    close(79)
#endif
  end subroutine close_outfiles

  subroutine writedata(sel)
    
    USE mod_time
    USE mod_grid
    USE mod_pos
    USE mod_traj
    USE mod_loopvars
    
    IMPLICIT NONE

    REAL                                 :: vort
    INTEGER                              :: sel ,xf ,yf ,zf ,n
    INTEGER*8, SAVE                      :: recPosIn=0  ,recPosOut=0
    INTEGER*8, SAVE                      :: recPosRun=0 ,recPosErr=0
    INTEGER*8, SAVE                      :: recPosKll=0
    REAL                                 :: x14 ,y14 ,z14
    REAL*8                               :: twrite

#if defined for || sim 
566 format(i8,i7,f7.2,f7.2,f7.1,f10.2,f10.2 &
         ,f10.1,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined rco || baltix 
566 format(i8,i7,f7.2,f7.2,f7.1,2f12.4 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined tes 
566 format(i8,i7,f8.3,f8.3,f7.3,2f10.2 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined ifs 
566 format(i8,i7,f7.2,f7.2,f7.2,f10.2,f10.2 &
         ,f15.0,f8.2,f8.2,f8.2,f6.0,8e8.1 )
#elif defined orc
    !566 format(i8,i7,2f8.2,f6.2,2f10.2 &
    !         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#else
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
    !566 format(i7,i7,f7.2,f7.2,f7.1,f10.4,f10.4 &
    !         ,f13.4,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#endif
    
    xf   = floor(x1)
    yf   = floor(y1)
    zf   = floor(z1)
    
    !if ((sel .ne. 19) .and. (sel.ne.40)) then
       ! this requires too much memory
       !       vort = (vvel(xf+1,yf,zf)-vvel(xf-1,yf,zf))/4000 - &
       !            (uvel(xf,yf+1,zf)-uvel(xf,yf-1,zf))/4000   
    !end if

#if defined textwrite 
    select case (sel)
    case (10)
       write(58,566) ntrac,niter,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
    case (11)
       if(  (kriva == 1 .AND. nrj(ntrac,4) == niter-1      ) .or. &
            (kriva == 2 .AND. scrivi                       ) .or. &
            (kriva == 3                                    ) .or. &
            (kriva == 4 .AND. niter == 1                   ) .or. &
            (kriva == 5 .AND.                                  &
          &  MOD((REAL(tt)-REAL(t0))*REAL(NGCM)/REAL(ITER), 3600.) == 0.d0 ) .or. &
            (kriva == 6 .AND. .not.scrivi                  ) ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
#if defined biol
          write(56,566) ntrac,ints,x1,y1,z1,tt/3600.,t0/3600.
#else
          write(56,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
#endif        
       endif
    case (13)
       ! === write sed pos ===
       write(57,566) ntrac,niter,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens 
    case (14)
       write(56,566) ntrac,ints,x1,y1,z1, &
            tt/60.,t0/3600.,subvol,temp,salt,dens
    case (15)
       write(57,566) ntrac,ints,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens
    case (16)
       if(kriva.ne.0 ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
          write(56,566) ntrac,ints,x1,y1,z1, &
               tt/tday,t0/tday,subvol,temp,salt,dens
       end if
    case (17)
       write(57,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (19)
       ! === write last sedimentation positions ===
       open(34,file=trim(outDataDir)//trim(outDataFile)//'_sed.asc') 
       do n=1,ntracmax
        if(nrj(n,1).ne.0) then
         write(34,566) n,nrj(n,4),trj(n,1),trj(n,2),trj(n,3),trj(n,4)/tday,trj(n,7)/tday
      endif
       enddo
       close(34)
    case (40)
       write(59,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (99) !switch
       
    end select
#endif 
   
#if defined binwrite 
    x14=real(x1,kind=4)
    y14=real(y1,kind=4)
    z14=real(z1,kind=4)
    if (twritetype==1) then
       twrite = tt
    else if (twritetype==2) then
       call updateclock
       twrite = currJDtot
    else
       twrite = real(ints,kind=8)
    end if
    select case (sel)       
    case (10) !in
       recPosIn = recPosIn + 1
       write(unit=78 ,rec=recPosIn) ntrac,twrite,x14,y14,z14
       return
    case (11)
       if(  (kriva == 1 .and. nrj(ntrac,4)  ==  niter-1    ) .or. &
            (kriva == 2 .and. scrivi                       ) .or. &
            (kriva == 3                                    ) .or. &
            (kriva == 4 .and. niter == 1                   ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi                  ) ) then
#if defined tempsalt
          call interp(ib,jb,kb,x1,y1,z1,temp, salt,  dens,1)
          call interp(ib,jb,kb,x1,y1,z1,temp2,salt2, dens2,2)
          !z14=real(salt*rb+salt2*(1-rb),kind=4)
#endif
          recPosRun = recPosRun+1
          write(unit=76 ,rec=recPosRun) ntrac,twrite,x14,y14,z14
       end if
    case (13)
       recPosKll = recPosKll + 1
       write(unit=77 ,rec=recPosKll) ntrac,twrite,x14,y14,z14   
    case (15)
       recPosRun = recPosRun + 1
       write(unit=76 ,rec=recPosRun) ntrac,twrite,x14,y14,z14   
    case (17) !out
       recPosOut = recPosOut + 1
       write(unit=77 ,rec=recPosOut) ntrac,twrite,x14,y14,z14   
    case (19) !end
       recPosOut = recPosOut + 1
       write(unit=75 ,rec=recPosOut) ntrac,twrite,x14,y14,z14
    case (40) !error
       recPosErr=recPosErr + 1    
       write(unit=79 ,rec=recPosErr) ntrac,twrite,x14,y14,z14   
    case (99) !switch
       if ((recPosRun > 50000000).and.(intminInOutFile.eq.2)) then
          call close_outfiles
          call open_outfiles
          recPosRun = 0
          recPosIn  = 0
          recPosOut = 0
          recPosErr = 0
          print *, "Switched run file" 
       end if
    end select
#endif    

#if defined csvwrite 
    x14=real(x1,kind=4)
    y14=real(y1,kind=4)
    z14=real(z1,kind=4)
    if (twritetype==1) then
       twrite = tt
    else if (twritetype==2) then
       call updateclock
       twrite = currJDtot
    else
       twrite = real(ints,kind=8)
    end if
    select case (sel)       
    case (10)
       write(88,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
       return
    case (11)
       if(  (kriva == 1 .and. nrj(ntrac,4)  ==  niter-1    ) .or. &
            (kriva == 2 .and. scrivi                       ) .or. &
            (kriva == 3                                    ) .or. &
            (kriva == 4 .and. niter == 1                   ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi                  )  ) then
          !!!! CALL FIELD-INTERP !!!!
          write(86,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
       end if
    case (13)
       write(87,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (15)
       write(86,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (17)
       write(87,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (19)
       write(85,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (40)
       write(89,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (99) !switch
       
    end select
#endif   
  end subroutine writedata

end module mod_write
