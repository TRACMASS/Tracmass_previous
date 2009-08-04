!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
! write the simulated tracer
#ifdef tracer
subroutine writetracer
USE mod_param
USE mod_name
USE mod_time
USE mod_tracer
IMPLICIT none

INTEGER nwx

nwx=IMT*JMT*KM*4
open(53,file=directory//'orm/tracer.'//name,form='unformatted',access='DIRECT',recl=nwx)
write(53,rec=ints+1) tra
close(53)

return
end subroutine writetracer

#endif
