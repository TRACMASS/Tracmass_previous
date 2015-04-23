PROGRAM TRACMASS
 
  USE mod_seed, only:  nqua, nff, num
  USE mod_grid, only:  dyu, dxv
  USE mod_time, only:  intmin, intstart, intend, intspin, intrun, intmax, &
       partQuant, minvelints, voltr
  USE mod_write, only: open_outfiles, close_outfiles
  USE mod_print, only: print_header_main, writesetup_main
#ifdef diffusion
  USE mod_diffusion
#endif
  
  IMPLICIT none
  INTEGER                                    :: i,j,n

  call print_header_main
  call init_params
  call coordinat
  call writesetup_main
  
  modrundirCond: if(nff == 1) then ! forward 
     intstart =  intmin          
     intend   =  intmax
  elseif(nff == 2) then ! backward
     intstart =  intmin+intrun
     minvelints = minvelints + intrun
     intend   =  intmin
     intspin  = -intspin
     intrun   = -intrun
     nff      =  -1    
  end if modrundirCond

  call setupgrid

  if (minval(dxv) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dxv contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if
  if (minval(dyu) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dyu contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if

  call init_seed

  if(nqua.eq.1) then ! number of trajectories (per time resolution)
     ! num=NTRACMAX
     num = partQuant
  elseif(nqua.eq.2) then 
     voltr = partQuant 
  elseif(nqua.eq.3) then 
     voltr = partQuant
  end if
  
  call open_outfiles
  call loop
  call close_outfiles

end PROGRAM TRACMASS


