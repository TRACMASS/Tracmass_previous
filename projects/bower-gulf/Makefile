# Makfile for orca1
#========================================================================
  INPUT_INT1 = intmin       # Variables set from command line
  INPUT_INT2 = intrun       # Use 'dummy' if not used.   
#------------------------------------------------------------------------
#   fl01   = -Dstationary    # Stationary velocity fields with analytical solution
   fl01   = -Dtimestep      # Time steps with analytical stationary scheme differential Eqs.
#  fl01   = -Dtimeanalyt    #  Analytical time scheme used to solve the
#  fl02   = -Dregulardt     # Regular time steps to be used with -Dtimestep
#  fl05   = -Dmonthly       # Monthly GCM data
#  fl06   = -Dseasonal      # More than two data sets in the central memory 
#------------------------------------------------------------------------
#   fl07   = -Dtempsalt      # Include temp and salt
   fl08   = -Dtextwrite     # Write results to textfile
#  fl09   = -Dbinwrite      # Write results to binaryfile
#  fl10   = -Dmysqlwrite    # Write results to mysql
#------------------------------------------------------------------------
#   fl11  = -Dstreamxy       # Calculates the barotropic stream function.
#   fl12  = -Dstreamr        #    vertical stream function, z=density
#   fl13  = -Dstreamts       #    vertical stream function, z=temp&salt
#   fl14  = -Dstreamv        #    vertical stream function, z=depth
#   fl15  = -Dstream_thermohaline        #    psi(T,S)
#   fl16  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#  fl17  = -Drst_rest       # Restictions on what water masses are alowed
#------------------------------------------------------------------------
  fl21  = -Dtwodim         # Turn off vertical velocities.
#  fl22  = -Dfull_wflux     # Use a full 3D wflux field.
#  fl23  = -Dexplicit_w     # Use a given vertical velocity.
#------------------------------------------------------------------------
#   fl25  = -Dvarbottombox   # Variable bottom box to match actual depth
#   fl26  = -Dfreesurface    # Variable bottom box to match actual depth
   fl27  = -Dzgrid1D        # Cell depths defined as vector (for z-coord?)
#  fl28  = -Dzgrid3D        # Cell depths defined as 3D grid (for sigma)
#  fl29  = -Dzgrid3Dt       # Cell depths 3D and time interp. (for atm)
#------------------------------------------------------------------------
#  fl31  = -Dselect         # Select only one trajectory (for debugging)
#  fl32  = -Dtracer         # Stores a simulated tracer
#  fl33  = -Dsediment       # Sediment code developed for RCO
#------------------------------------------------------------------------
#  fl35  = -Dturb           # Adds subgrid turbulent velocities 
#  fl36  = -Ddiffusion      # Adds a diffusion on trajectory
#  fl37  = -Danisodiffusion # Adds an anisotropic diffusion on trajectory
#------------------------------------------------------------------------
#   fl41   = -Ddrifter       # surface drifter depth average of uflux/vflux
#------------------------------------------------------------------------
#   fl42   = -Dzipped        # zipped files to be unzipped
#========================================================================

ORM_FLAGS= -D$(PROJECT) \
$(fl01)$(fl02)$(fl03)$(fl04)$(fl05)$(fl06)$(fl07)$(fl08)$(fl09)$(fl10)\
$(fl11)$(fl12)$(fl13)$(fl14)$(fl15)$(fl16)$(fl17)$(fl18)$(fl19)$(fl20)\
$(fl21)$(fl22)$(fl23)$(fl24)$(fl25)$(fl26)$(fl27)$(fl28)$(fl29)$(fl30)\
$(fl31)$(fl32)$(fl33)$(fl34)$(fl35)$(fl36)$(fl37)$(fl38)$(fl39)$(fl40)\
$(fl41)$(fl42)$(fl43)$(fl44)$(fl45)$(fl46)$(fl47)$(fl48)$(fl49)$(fl50)\

