
module mod_laplacian
  USE mod_param
  USE mod_grid
  USE mod_vel
  USE mod_loopvars
  USE mod_time
  USE mod_streamfunctions
  USE mod_deformation  
  
  IMPLICIT none
  
contains

subroutine laplacian
   !!
   !! Calculate Laplacian of velocity field 
   !! Note that 
   !! nabla^2 u = d/dx div - d/dy vort
   !! nabla^2 v = d/dy div + d/dx vort
   !! so we calculate divergence and vorticity as well
   !! in case anyone needs it in the future
   !!
   
   integer :: ji,jj,jk,jn
   
   do jn = 1,2
      do jk = 1,km
         do jj = 2,jmt-1
            do ji = 2,imt-1
               !vort(ji,jj,jk,2) = (vvel(ji+1,jj  ,jk) - vvel(ji  ,jj  ,jk)) / dxv(ji,jj) - &
               !                 & (uvel(ji  ,jj+1,jk) - uvel(ji  ,jj  ,jk)) / dyu(ji,jj)
               !hdiv(ji,jj,jk,2) = (uvel(ji  ,jj  ,jk) - uvel(ji-1,jj  ,jk)) / dxv(ji,jj) + &
               !                 & (vvel(ji  ,jj  ,jk) - vvel(ji  ,jj-1,jk)) / dyu(ji,jj)
               vort(ji,jj,jk,jn) = (vflux(ji+1,jj  ,jk, jn) - vflux(ji  ,jj  ,jk, jn)) / dxv(ji,jj) - &       
                                 & (uflux(ji  ,jj+1,jk, jn) - uflux(ji  ,jj  ,jk, jn)) / dyu(ji,jj)         
               hdiv(ji,jj,jk,jn) = (uflux(ji  ,jj  ,jk, jn) - uflux(ji-1,jj  ,jk, jn)) / dxv(ji,jj) + &     
                                 & (vflux(ji  ,jj  ,jk, jn) - vflux(ji  ,jj-1,jk, jn)) / dyu(ji,jj)
            end do
         end do
         
         do jj = 2,jmt-1
            do ji = 2,imt-1
               !lapu(ji,jj,jk,2) = (hdiv(ji+1,jj,jk,2) - hdiv(ji,jj,jk,2)) / dxv(ji,jj) - &
               !                 & (vort(ji,jj,jk,2) - vort(ji,jj-1,jk,2)) / dyu(ji,jj)
               !lapv(ji,jj,jk,2) = (hdiv(ji,jj+1,jk,2) - hdiv(ji,jj,jk,2)) / dyu(ji,jj) - &
               !                 & (vort(ji,jj,jk,2) - vort(ji-1,jj,jk,2)) / dxv(ji,jj)
               lapu(ji,jj,jk,jn) = (hdiv(ji+1,jj  ,jk,jn) - hdiv(ji  ,jj  ,jk,jn)) / dxv(ji,jj) - &
                                 & (vort(ji  ,jj  ,jk,jn) - vort(ji  ,jj-1,jk,jn)) / dyu(ji,jj)
               lapv(ji,jj,jk,jn) = (hdiv(ji  ,jj+1,jk,jn) - hdiv(ji  ,jj  ,jk,jn)) / dyu(ji,jj) - &
                                 & (vort(ji  ,jj  ,jk,jn) - vort(ji-1,jj  ,jk,jn)) / dxv(ji,jj)
            end do
         end do
      end do
   end do
end subroutine

end module mod_laplacian 
