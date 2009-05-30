
#========================================================================
  fl01   = -Dmean           #
  fl02   = -Dtempsalt       # Include temp and salt
  fl03   = -Dmysqlwrite     # Write results to mysql
#------------------------------------------------------------------------
#  fl04  = -Dstreamxy       # Calculates the barotropic stream function.
#  fl05  = -Dstreamr        #    vertical stream function, z=density
#  fl06  = -Dstreamts       #    vertical stream function, z=???
#  fl07  = -Dstreamv        #    vertical stream function, z=depth
#  fl08  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#------------------------------------------------------------------------
#  fl09   = -Dselect         # Select only one trajectory (for debugging)
#  fl10  = -Dtracer         # Stores a simulated tracer
#  fl11  = -Dsediment       # Sediment code developed for RCO
#  fl12  = -Dtwodim         # Uses a two-dimensional velocity field. 
#========================================================================

ORM_FLAGS= -D$(PROJECT) $(fl01)$(fl02)$(fl03)$(fl04)$(fl05)$(fl06)\
			$(fl07)$(fl08)$(fl09)$(fl10)$(fl11)$(fl02)
