module mod_psi

  USE mod_loopvars
  USE mod_grid
#ifdef streamxy
  USE mod_streamxy
#endif
#ifdef streamr
  USE mod_streamr
#endif
#ifdef streamv
  USE mod_streamv
#endif
#ifdef stream_thermohaline
  USE mod_stream_thermohaline
#endif

  CONTAINS

subroutine savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,xy,dir,flux)

IMPLICIT NONE

! Subroutine that writes the stream functions
! call save_psi(i,j,k,lbas,xy,dir,flux)
! i,j,k - coordinate where the psi is to be written
! xy - zonal/meridional, dir - north or south
! flux = subvol*ff

INTEGER             :: ia,ja,ka   !where to write
REAL*8              :: x1,y1,z1
INTEGER             :: xy, dir !1 - zonal, 2 - meridional, 3 - vertical
REAL                :: flux
REAL                :: temp,salt,dens
INTEGER             :: mrb,mtb,msb
INTEGER             :: mta,msa,m


#ifdef stream_thermohaline
		  do m=mta,mtb-1
           psi_ts(m,msb,1,lbas) = psi_ts(m,msb,1,lbas) + flux
          enddo
		  do m=mtb,mta-1
           psi_ts(m,msb,1,lbas) = psi_ts(m,msb,1,lbas) - flux
          enddo
		  do m=msa,msb-1
           psi_ts(mtb,m,2,lbas) = psi_ts(mtb,m,2,lbas) + flux
          enddo
		  do m=msb,msa-1
           psi_ts(mtb,m,2,lbas) = psi_ts(mtb,m,2,lbas) - flux
          enddo
#endif

flux = flux*real(dir)
select case(xy)
     ! === Zonal component ===
     case(1)
          ! === Barotropic stream function ===
#ifdef streamxy
          stxyx(ia,ja,lbas) = stxyx(ia,ja,lbas) + flux
#endif
          ! === Overturning stream function ===
#ifdef streamv
          stxz(ia,ka,lbas) = stxz(ia,ka,lbas) + flux
#endif
          ! === Overturning dens/temp/salt stream function ===
#ifdef streamr
          stxr(ia,mrb,lbas,1) = stxr(ia,mrb,lbas,1) + flux
#endif
#ifdef streamts
          stxr(ia,mtb,lbas,2) = stxr(ia,mtb,lbas,2) + flux
          stxr(ia,msb,lbas,3) = stxr(ia,msb,lbas,3) + flux
#endif

    
     ! === Meridional component ===
     case(2)
          ! === Barotropic stream function ===
#ifdef streamxy
          stxyy(ia,ja,lbas) = stxyy(ia,ja,lbas) + flux
#endif
          ! === Overturning stream function ===
#ifdef streamv
          styz(ja,ka,lbas) = styz(ja,ka,lbas) + flux
#endif
          ! === Overturning dens/temp/salt stream function ===
#ifdef streamr
          styr(ja,mrb,lbas,1) = styr(ja,mrb,lbas,1) + flux
#endif
#ifdef streamts
          styr(ja,mtb,lbas,2) = styr(ja,mtb,lbas,2) + flux
          styr(ja,msb,lbas,3) = styr(ja,msb,lbas,3) + flux 
#endif

!     ! === Vertical component ===
!     case(3)


end select

end subroutine

end module mod_psi