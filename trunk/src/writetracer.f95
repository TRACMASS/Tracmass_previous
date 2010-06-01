  ! write the simulated tracer
  subroutine writetracer
#ifdef tracer
    USE mod_param
    USE mod_name
    USE mod_time
    USE mod_tracer
    IMPLICIT none
    
    INTEGER nwx
    
    nwx=IMT*JMT*KM*4
    open(53,file=directory//'orm/tracer.'//name,form='unformatted', &
         access='DIRECT',recl=nwx)
    write(53,rec=ints+1) tra
close(53)
return
#endif /*tracer*/
end subroutine writetracer


