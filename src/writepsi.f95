subroutine writepsi

  ! === writes the Lagrangian stream functions ===

  USE mod_param
  USE mod_name
  USE mod_time
  USE mod_streamfunctions

IMPLICIT none

#ifdef streamxy
open(51,file=trim(outDataDir)//trim(outDataFile)//'_psi_xy_yx.bin',form='unformatted')
write(51) stxyy
write(51) stxyx
close(51)
#endif
#if defined streamv
open(52,file=trim(outDataDir)//trim(outDataFile)//'_psi_yz_xz.bin',form='unformatted')
write(52) styz
write(52) stxz
close(52)
#endif
#if defined streamr 
open(53,file=trim(outDataDir)//trim(outDataFile)//'_psi_xr_yr.bin',form='unformatted')
write(53) stxr
write(53) styr
close(53)
#endif
#ifdef stream_thermohaline
open(54,file=trim(outDataDir)//trim(outDataFile)//'_psi_ts.bin',form='unformatted')
write(54) psi_ts
close(54)
#endif

print *,'stream functions written'

return
end subroutine writepsi

