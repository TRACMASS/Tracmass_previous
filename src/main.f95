PROGRAM main
 
  USE mod_param
  USE mod_seed
  USE mod_name
  USE mod_time 
  USE mod_write
  USE mod_domain
  USE mod_buoyancy
  USE mod_print
#ifdef diffusion
  USE mod_diffusion
#endif
  
  IMPLICIT none
  INTEGER                                    :: i,j,n

  call print_header_main
  call init_params
  call coordinat
  call writesetup_main

  modrundirCond: if(intstep.gt.0) then ! forward 
!  modrundirCond: if(nff == 1) then ! forward 
     intstart =  intmin          
     intend   =  intmax
!  elseif(nff == -1) then ! backward
  elseif(intstep.lt.0) then ! backward
     intstart =  intmin+intrun
     minvelints = minvelints + intrun
     intend   =  intmin
     intspin  = -intspin
     intrun   = -intrun    
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
     num=partQuant
  elseif(nqua.eq.2) then 
     voltr=partQuant 
  elseif(nqua.eq.3) then 
     voltr=partQuant
  endif
  
  call open_outfiles
  call loop
  call close_outfiles

  return


end PROGRAM main


