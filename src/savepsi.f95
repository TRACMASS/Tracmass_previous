module mod_psi

  USE mod_loopvars
  USE mod_grid
  USE mod_streamfunctions
  USE mod_tempsalt
  USE mod_traj

  CONTAINS

subroutine savepsi(ia,ja,ka,mrb,mta,mtb,msa,msb,xy,dir,flux)

IMPLICIT NONE

! Lagt till mod_tempsalt (för att få temp och salt)
! och mod_traj för ib,jb,kb för att använda interp2.

! Subroutine that writes the stream functions
! call save_psi(i,j,k,lbas,xy,dir,flux)
! i,j,k - coordinate where the psi is to be written
! xy - zonal/meridional, dir - north or south
! flux = subvol*ff

INTEGER             :: ia,ja,ka   !where to write
!INTEGER             :: ib,jb,kb
REAL*8              :: x1,y1,z1
INTEGER             :: xy, dir !1 - zonal, 2 - meridional, 3 - vertical
REAL                :: flux
REAL                :: tempa,salta,densa,tempb,saltb,densb
INTEGER             :: mrb,mtb,msb
INTEGER             :: mta,msa,m


#ifdef stream_thermohaline
!		  do m=mta,mtb-1
!           psi_ts(m,msb,1,lbas) = psi_ts(m,msb,1,lbas) + flux
!          enddo
!		  do m=mtb,mta-1
!           psi_ts(m,msb,1,lbas) = psi_ts(m,msb,1,lbas) - flux
!          enddo
!		  do m=msa,msb-1
!           psi_ts(mtb,m,2,lbas) = psi_ts(mtb,m,2,lbas) + flux
!          enddo
!		  do m=msb,msa-1
!           psi_ts(mtb,m,2,lbas) = psi_ts(mtb,m,2,lbas) - flux
!          enddo
          
          
DO m=mta,mtb-1
   psi_ts(m,(msa+msb)/2.,1,lbas) = psi_ts(m,(msa+msb)/2.,1,lbas) + flux
!PRINT *, '1', ia, ja, ka, mta, mtb, msa, msb
ENDDO
DO m=mtb,mta-1
   psi_ts(m,(msa+msb)/2.,1,lbas) = psi_ts(m,(msa+msb)/2.,1,lbas) - flux
!PRINT *, '2', ia, ja, ka, mta, mtb, msa, msb
ENDDO
DO m=msa,msb-1
   psi_ts((mta+mtb)/2.,m,2,lbas) = psi_ts((mta+mtb)/2.,m,2,lbas) + flux
!PRINT *, '3', ia, ja, ka, mta, mtb, msa, msb
ENDDO
DO m=msb,msa-1
   psi_ts((mta+mtb)/2.,m,2,lbas) = psi_ts((mta+mtb)/2.,m,2,lbas) - flux
!PRINT *, '4', ia, ja, ka, mta, mtb, msa, msb
ENDDO


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
! LAGT TILL DETTA !Sara
     case(3)

          ! === Overturning dens/temp/salt stream function with depth-coordinates ===
#ifdef streamr
          stzr(ka,mrb,lbas,1) = stzr(ka,mrb,lbas,1) + flux
#endif

#ifdef streamts
          stzr(ka,mtb,lbas,2) = stzr(ka,mtb,lbas,2) + flux
          stzr(ka,msb,lbas,3) = stzr(ka,msb,lbas,3) + flux 
#endif
! ===========

end select

!==== trajecory convergence of heat and salinity, can be used for other tracers in future (?)
#ifdef tracer_convergence

select case(xy)

! === Zonal Component ===
     case(1)
#ifdef tempsalt
        CALL interp2(ia+1,ja,ka,tempb,saltb,densb)
#endif
        CALL interp2(ia,ja,ka,tempa,salta,densa)
        !IF (mlh(ia,ja,nsp) < dep(42+1-ka)) THEN
           !PRINT *, ka, 42+1-ka, dep(42+1-ka), mlh(ia,ja,nsp)
           uct(ia,ja,ka,lbas) = uct(ia,ja,ka,lbas) + flux*0.5*(tempb+tempa)
           ucs(ia,ja,ka,lbas) = ucs(ia,ja,ka,lbas) + flux*0.5*(saltb+salta)

           IF (tempb == 0) THEN
              PRINT *, '1', ia, ia+1, ja, ka
           ENDIF
        !ENDIF

! === Meridional Component ===
     case(2)
#ifdef tempsalt
        CALL interp2(ia,ja+1,ka,tempb,saltb,densb)
#endif
        CALL interp2(ia,ja,ka,tempa,salta,densa)
       ! IF (mlh(ia,ja,nsp) < dep(42+1-ka)) THEN
           vct(ia,ja,ka,lbas) = vct(ia,ja,ka,lbas) + flux*0.5*(tempb+tempa)
           vcs(ia,ja,ka,lbas) = vcs(ia,ja,ka,lbas) + flux*0.5*(saltb+salta)

           IF (tempb == 0) THEN
              PRINT *, '2', ia, ja, ja+1, ka
           ENDIF
       ! ENDIF

! === Vertical Component ===
     case(3)
#ifdef tempsalt
        CALL interp2(ia,ja,ka+1,tempb,saltb,densb)
#endif
        CALL interp2(ia,ja,ka,tempa,salta,densa)
        !IF ( mlh(ia,ja,nsp) < dep(42+1-ka)) THEN
           IF (ka /= KM) THEN
              wct(ia,ja,ka,lbas) = wct(ia,ja,ka,lbas) + flux*0.5*(tempb+tempa)
              wcs(ia,ja,ka,lbas) = wcs(ia,ja,ka,lbas) + flux*0.5*(saltb+salta)
           ENDIF

           IF (tempb == 0) THEN
              PRINT *, '3', ia, ja, ka, ka+1
           ENDIF
        !ENDIF

end select
#endif



end subroutine

end module mod_psi
