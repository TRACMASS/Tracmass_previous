SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------

  ! === Init local variables for the subroutine ===
  INTEGER                                     :: i ,j ,k, ip, jp, kk, kku, kkv, umin, vmin
  REAL*4, DIMENSION(2)                        :: dz_uvec, dz_vvec
  REAL*4, DIMENSION(3)                        :: km_uvec, km_vvec  
  REAL*4, ALLOCATABLE, DIMENSION(:,:)         :: e1t,e2t,e3t_ps
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)       :: dzt0
  CHARACTER (len=200)                         :: gridFile

  ALLOCATE ( e1t(IMT,JMT) , e2t(IMT,JMT) , e3t_ps(IMT,JMT) )
  ALLOCATE ( dzt0(IMT,JMT,KM), dzu(imt,jmt,km), dzv(imt,jmt,km)) 
  ALLOCATE ( kmu(IMT,JMT) , kmv(IMT,JMT) )

  ! === Variable mapping order for reading in netCDF files ===
  ! Order is  t  k  i  j 
  map2d    = [3, 4, 1, 2]
  map3d    = [2, 3, 4, 1]

  ncTpos = 1
  ! === Open mesh file ===
  gridFile = trim(inDataDir)//'GL2V1_mesh_mask_new.nc'

  ! === Read dx and dy for T points ===
  e1t  = get2DfieldNC(gridFile, 'e1t')
  e2t  = get2DfieldNC(gridFile, 'e2t')
  ! === Compute area of grid box    ===
  dxdy = e1t * e2t  

  ! === Read dy for U points ===
  dyu  = get2DfieldNC(gridFile, 'e2u')
  ! === Read dx for V points ===
  dxv  = get2DfieldNC(gridFile, 'e1v')
 
  ! === Read dz for T points ===
  dz   = get1DfieldNC(gridFile, 'e3t_0')
  ! mbathy indicates which is the bottom level at a given (i,j) t-point 
  kmt  = get2DfieldNC(gridFile, 'mbathy') + 1 !! SWITCH FROM 0 (GLORYS) to 1 (FORTRAN) base indexing - NOT IN ORC
  ! === e3t_ps indicates dz for the bottom level ===
  ! === - equivalent to 'botbox' variable in previous NEMO methods ===
  e3t_ps = get2DfieldNC(gridFile,'e3t_ps')

  ! === Assign place holder for grid dzt pre-ssh ===
  dzt0 = 0  
  ! === Propagate grid dzt array ===
  do j=1,JMT
     do  i=1,IMT
         dzt0(i,j,:)  = dz
         kk           = kmt(i,j) !! kk = index for partial cell 
         dzt0(i,j,kk) = e3t_ps(i,j)
         if (kk.lt.KM) then
            do k = kk+1, KM
               dzt0(i,j,k) = 0
            enddo
         endif
     enddo
  enddo

  kmu = 0 
  kmv = 0
  dzu = 0
  dzv = 0

  ! Identify land mask
  mask = 1.0
  where(dzt0(:,:,1).eq.0)
       mask = 0
  end where


  do j=1,jmt
     jp=j+1
     if(jp.eq.jmt+1) jp=jmt
     do i=1,imt
        ip=i+1
        !if(ip.eq.IMT+1) ip=1  !! <-FOR GLOBAL GRID DOMAIN THAT CONNECTS LON=0 AND LON=360
        if(ip.eq.IMT+1) ip=imt
        km_uvec = (/kmt(i,j), kmt(ip,j), KM/) ! KM justification: vertical subdomain may be < kmt 
        km_vvec = (/kmt(i,j), kmt(i,jp), KM/) 

        dz_uvec = [e3t_ps(i,j),e3t_ps(ip,j)]
        dz_vvec = [e3t_ps(i,j),e3t_ps(i,jp)]
        
        umin = minloc(km_uvec,1)
        vmin = minloc(km_vvec,1)
        kmu(i,j) = km_uvec(umin) !! SAME AS: kmu(i,j)=min(kmt(i,j), kmt(ip,j))
        kmv(i,j) = km_vvec(vmin) !! SAME AS: kmv(i,j)=min(kmt(i,j), kmt(i,jp))
        kku = kmu(i,j)

        dzu(i,j,:)  = dz

        if (km_uvec(1).eq.km_uvec(2)) then
           dzu(i,j,kku) = min(dz_uvec(1),dz_uvec(2))
        else
           dzu(i,j,kku) = dz_uvec(umin) 
        endif

        if (kku.lt.KM) then
            do k = kku+1, KM
               dzu(i,j,k) = 0
            enddo
        endif

        kkv = kmv(i,j)
        dzv(i,j,:)  = dz
        if (km_vvec(1).eq.km_vvec(2)) then
           dzv(i,j,kkv) = min(dz_vvec(1),dz_vvec(2))
        else
           dzv(i,j,kkv) = dz_vvec(vmin)
        endif

        if (kkv.lt.KM) then
            do k = kkv+1, KM
               dzv(i,j,k) = 0
            enddo
        endif
     enddo
  enddo

 ! Flipping dzt in the vertical direction to match tracmass convention
 do k=1,km
    dzt(:,:,k,1) = dzt0(:,:,km+1-k)
    dzt(:,:,k,2) = dzt0(:,:,km+1-k)
 end do
 
end SUBROUTINE setupgrid
