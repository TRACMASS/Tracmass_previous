!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
! write the Lagrangian stream functions & a simulated tracer
subroutine writepsi
USE mod_param
USE mod_name
USE mod_time
#ifdef streamxy
USE mod_streamxy
#endif
#ifdef streamv
USE mod_streamv
#endif
#ifdef streamr
USE mod_streamr
#endif
#ifdef stream_thermohaline
  USE mod_stream_thermohaline
#endif
IMPLICIT none


print *,trim(outDataDir)//'stxy.'//trim(outDataFile)

#ifdef streamxy
open(51,file=trim(outDataDir)//'stxy.'//trim(outDataFile),form='unformatted')
write(51) stxyy
write(51) stxyx
close(51)
#endif
#if defined streamv
open(52,file=trim(outDataDir)//'stv.'//trim(outDataFile),form='unformatted')
write(52) styz
write(52) stxz
close(52)
#endif
#if defined streamr 
open(53,file=trim(outDataDir)//'str.'//trim(outDataFile),form='unformatted')
write(53) stxr
write(53) styr
close(53)
#endif
#ifdef stream_thermohaline
open(54,file=trim(outDataDir)//'str_thermohaline.'//trim(outDataFile),form='unformatted')
write(54) psi_ts
close(54)
#endif


return
end subroutine writepsi

