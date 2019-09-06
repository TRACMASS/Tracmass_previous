module mod_psi_new

  USE mod_precdef 
  USE mod_loopvars
  USE mod_grid
  USE mod_streamfunctions
  USE mod_tempsalt, only: n3Dtracers

  CONTAINS

  subroutine savepsi_new(ia,ja,ka,mtra,mtrb,xy,dir,flux)

    implicit none 
    
    ! Subroutine that writes the stream functions
    ! call save_psi(i,j,k,lbas,xy,dir,flux)
    ! i,j,k - coordinate where the psi is to be written
    ! xy - zonal/meridional, dir - north or south
    ! flux = subvol*ff
    
    integer, intent(in)             :: ia,ja,ka   !where to write
    integer                         :: jt, m
    real(DP)            :: x1,y1,z1
    integer, intent(in)             :: xy, dir !1 - zonal, 2 - meridional, 3 - vertical
    real, intent(in)                :: flux
    real(PP)                :: temp,salt,dens, flx
    integer, dimension(n3Dtracers), intent(in) :: mtra, mtrb


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
          
          
         do m=mtra(2),mtrb(2)-1
           psi_ts(m,mtra(1),1,lbas) = psi_ts(m,mtra(1),1,lbas) + flux
         enddo
         do m=mtrb(2)+1,mtra(2)
           psi_ts(m,mtra(1),1,lbas) = psi_ts(m,mtra(1),1,lbas) - flux
         enddo
         do m=mtra(1),mtrb(1)-1
           psi_ts(mtra(2),m,2,lbas) = psi_ts(mtra(2),m,2,lbas) + flux
         enddo
         do m=mtrb(1)+1,mtra(1)
           psi_ts(mtra(2),m,2,lbas) = psi_ts(mtra(2),m,2,lbas) - flux
         enddo
#endif

flx = flux*real(dir)

select case(xy)
     ! === Zonal component ===
     case(1)
          ! === Barotropic stream function ===
#ifdef streamxy
          stxyx(ia,ja,lbas) = stxyx(ia,ja,lbas) + flx
#endif
          ! === Overturning stream function ===
#ifdef streamv
          stxz(ia,ka,lbas) = stxz(ia,ka,lbas) + flx
#endif
          ! === Overturning dens/temp/salt stream function ===
          do jt = 1, n3Dtracers
             stxr(ia,mtrb(jt),lbas,jt) = stxr(ia,mtrb(jt),lbas,jt) + flx
          end do 

    
     ! === Meridional component ===
     case(2)
          ! === Barotropic stream function ===
#ifdef streamxy
          stxyy(ia,ja,lbas) = stxyy(ia,ja,lbas) + flx
#endif
          ! === Overturning stream function ===
#ifdef streamv
          styz(ja,ka,lbas) = styz(ja,ka,lbas) + flx
#endif
          ! === Overturning dens/temp/salt stream function ===
          do jt = 1, n3Dtracers
             styr(ja,mtrb(jt),lbas,jt) = styr(ja,mtrb(jt),lbas,jt) + flx
          end do 

     ! === Vertical component ===
       case(3)
          ! === Depth dens/temp/salt stream function ===          
          do jt = 1, n3Dtracers
             stzr(ka,mtrb(jt),lbas,jt) = stzr(ka,mtrb(jt),lbas,jt) + flx
          end do 


end select

!==== trajecory convergence of heat and maybe other tracers
#ifdef tracer_convergence
select case(xy)
     case(1)
   	  converg(ia,ja,ka,lbas,1) = converg(ia,ja,ka,lbas,1) + flx
     case(2)
   	  converg(ia,ja,ka,lbas,2) = converg(ia,ja,ka,lbas,2) + flx
     case(3)
   	  converg(ia,ja,ka,lbas,3) = converg(ia,ja,ka,lbas,3) + flx
end select
#endif

end subroutine

end module mod_psi_new
