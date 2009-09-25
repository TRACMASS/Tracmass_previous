
#========================================================================
  INPUT_INT1        = intmin            # Variables set from command line
  INPUT_INT2        = intrun            # Use 'dummy' if not used.   
#------------------------------------------------------------------------
   fl01   = -Dtimestep       # Analytical time scheme used to solve the differential Eqs.
#  fl01   = -Dtimeanalyt    # Time steps with analytical stationary scheme 
   fl02   = -Dregulardt      # Regular time steps to be used with -Dtimestep
#------------------------------------------------------------------------
#  fl03   = -Dtempsalt       # Include temp and salt
   fl04   = -Dtextwrite      # Write results to textfile
#  fl05   = -Dbinwrite       # Write results to binaryfile
#  fl06   = -Dmysqlwrite     # Write results to mysql
#------------------------------------------------------------------------
#  fl07  = -Dstreamxy       # Calculates the barotropic stream function.
#  fl08  = -Dstreamr        #    vertical stream function, z=density
#  fl09  = -Dstreamts       #    vertical stream function, z=???
#  fl10  = -Dstreamv        #    vertical stream function, z=depth
#  fl11  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#------------------------------------------------------------------------
   fl12  = -Dtwodim         # Turn off vertical velocities.
#  fl13  = -Dfull_wflux     # Use a full 3D wflux field.
#  fl14  = -Dexplicit_w     # Use a given vertical velocity.
#------------------------------------------------------------------------
#  fl15  = -Dselect         # Select only one trajectory (for debugging)
#  fl16  = -Dtracer         # Stores a simulated tracer
#  fl17  = -Dsediment       # Sediment code developed for RCO
#   fl18  = -Dturb           # Adds subgrid turbulent velocities 
  fl18  = -Ddiffusion      # Adds a diffusion on trajectory
  fl19  = -Dellipticdiffusion      # Adds a diffusion on trajectory
#========================================================================

ORM_FLAGS= -D$(PROJECT) $(fl01)$(fl02)$(fl03)$(fl04)$(fl05)$(fl06)\
			$(fl07)$(fl08)$(fl09)$(fl10)$(fl11)$(fl12)\
			$(fl13)$(fl14)$(fl15)$(fl16)$(fl17)$(fl18)\
			$(fl19)$(fl20)$(fl21)$(fl22)$(fl23)$(fl24)\
