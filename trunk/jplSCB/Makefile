
#========================================================================
  INPUT_INT1        = intmin            # Variables set from command line
  INPUT_INT2        = intrun            # Use 'dummy' if not used.   
#------------------------------------------------------------------------
  fl01   = -Dmean           #
  fl02   = -Dtempsalt       # Include temp and salt
#  fl03   = -Dtextwrite      # Write results to textfile
  fl04   = -Dbinwrite       # Write results to binaryfile
#  fl05   = -Dmysqlwrite     # Write results to mysql
#------------------------------------------------------------------------
#  fl06  = -Dstreamxy       # Calculates the barotropic stream function.
#  fl07  = -Dstreamr        #    vertical stream function, z=density
#  fl08  = -Dstreamts       #    vertical stream function, z=???
#  fl09  = -Dstreamv        #    vertical stream function, z=depth
#  fl10  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#------------------------------------------------------------------------
#  fl11   = -Dselect         # Select only one trajectory (for debugging)
#  fl12  = -Dtracer         # Stores a simulated tracer
#  fl13  = -Dsediment       # Sediment code developed for RCO
#  fl14  = -Dtwodim         # Uses a two-dimensional velocity field. 
#========================================================================

ORM_FLAGS= -D$(PROJECT) $(fl01)$(fl02)$(fl03)$(fl04)$(fl05)$(fl06)\
			$(fl07)$(fl08)$(fl09)$(fl10)$(fl11)$(fl12)\
			$(fl13)$(fl14)$(fl15)$(fl16)$(fl17)$(fl18)
