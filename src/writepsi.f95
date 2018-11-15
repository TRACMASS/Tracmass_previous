subroutine write_streamfunctions

  ! === writes the Lagrangian stream functions ===

  USE mod_param
  USE mod_name
  USE mod_time
  USE mod_streamfunctions

IMPLICIT none

!call open_outfiles

#ifdef streamxy
rewind(51)
write(51) stxyy
write(51) stxyx
#endif

#if defined streamv
rewind(52)
write(52) styz
write(52) stxz
#endif

#if defined streamr 
rewind(53)
write(53) stxr
write(53) styr
write(53) stzr
#endif

#ifdef stream_thermohaline
rewind(54)
write(54) psi_ts
#endif

#if defined streamxy || streamv || streamr || stream_thermohaline
print *,'stream function written'
#endif

#ifdef tracer_convergence
rewind(55)
write(55) uct
write(55) vct
write(55) wct
#endif

#ifdef tracer_convergence
rewind(50)
write(50) ucs
write(50) vcs
write(50) wcs
#endif

#if defined tracer_convergence
PRINT *, 'Tracer convergence written'
#endif


return
end subroutine write_streamfunctions

