  subroutine writetracer
  
  ! === writes the simulated tracer ===


#ifdef tracer
    USE mod_param
    USE mod_name
    USE mod_time
    USE mod_tracer
    IMPLICIT none
    
    INTEGER nwx
    
    nwx=IMT*JMT*KM*4
    open(53,file=trim(outDataDir)//trim(outDataFile)//'tracer.bin',form='unformatted', &
         access='DIRECT',recl=nwx)
    write(53,rec=ints+1) tra
close(53)
return
#endif /*tracer*/
end subroutine writetracer


