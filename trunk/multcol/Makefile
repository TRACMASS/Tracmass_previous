
#========================================================================
  fl1   = -Dmean           #
  fl2   = -Dtempsalt       # Include temp and salt
  fl3   = -Dmysqlwrite     # Write results to mysql
  fl4   = -Dtime
#------------------------------------------------------------------------
#  fl4   = -Dstreamxy       # Calculates the barotropic stream function.
#  fl5   = -Dstreamr        #    vertical stream function, z=density
#  fl6   = -Dstreamts       #    vertical stream function, z=???
#  fl7   = -Dstreamv        #    vertical stream function, z=depth
#  fl11  = -Drerun          # Stores the Lagrangian stream functions as a
                            # function of the end positions that has been
                            # calculated in an identical previous run.    
#------------------------------------------------------------------------
#  fl8   = -Dselect         # Select only one trajectory (for debugging)
#  fl9   = -Dtracer         # Stores a simulated tracer
#  fl10  = -Dsediment       # Sediment code developed for RCO
#  fl12  = -Dtwodim         # Uses a two-dimensional velocity field. 
#========================================================================

ORM_FLAGS= -D$(PROJECT) $(fl1) $(fl2) $(fl3) $(fl4) $(fl5) $(fl6)\
			$(fl7) $(fl8) $(fl9) $(fl11) $(fl12)
