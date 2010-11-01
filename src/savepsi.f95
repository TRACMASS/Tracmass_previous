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

  CONTAINS

subroutine savepsi(ia,ja,ka,mra,mta,msa,xy,dir,flux)

IMPLICIT NONE

! Subroutine that writes the stream functions
! call save_psi(i,j,k,lbas,xy,dir,flux)
! i,j,k - coordinate where the psi is to be written
! xy - zonal/meridional, dir - north or south
! flux = subvol*ff

INTEGER             :: ia,ja,ka   !where to write
REAL*8              :: x1,y1,z1
INTEGER             :: xy, dir !1 - zonal, 2 - meridional
REAL                :: flux
REAL                :: temp,salt,dens
INTEGER             :: mra,mta,msa


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
          stxr(ia,mra,lbas,1) = stxr(ia,mra,lbas,1) + flux
#endif
#ifdef streamts
          stxr(ia,mta,lbas,2) = stxr(ia,mta,lbas,2) + flux
          stxr(ia,msa,lbas,3) = stxr(ia,msa,lbas,3) + flux
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
          styr(ja,mra,lbas,1) = styr(ja,mra,lbas,1) + flux
#endif
#ifdef streamts
          styr(ja,mta,lbas,2) = styr(ja,mta,lbas,2) + flux
          styr(ja,msa,lbas,3) = styr(ja,msa,lbas,3) + flux 
#endif

end select

end subroutine

end module mod_psi