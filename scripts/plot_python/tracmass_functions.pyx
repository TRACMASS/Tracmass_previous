import numpy as np
cimport numpy as np

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPE32 = np.int32
ctypedef np.int32_t DTYPE32_t
FTYPE = np.float
ctypedef np.float_t FTYPE_t
FTYPE32 = np.float32
ctypedef np.float32_t FTYPE32_t

cpdef find_numsteps(np.ndarray[np.float32_t, ndim=2] rawdata, int ntracmax=100000, float fill_value=-999.):
   
   cdef int npoints = rawdata.shape[0]
   cdef int ntrac, jj
   cdef int nsteps
   cdef np.ndarray[np.int32_t, ndim=1] steps = np.zeros([ntracmax], dtype=np.int32)
   nsteps = 0
   #print(nsteps)
   for jj in xrange(0,npoints):
      ntrac = <int>rawdata[jj,0]
      steps[ntrac] += 1
      #print(ntrac,steps[ntrac])
      if np.mod(jj,10000)==0: 
         print(jj,npoints)
   
   
   return steps
   

cpdef sort_particles(np.ndarray[np.float32_t, ndim=2] rawdata, int ntracmax=100000, int nsteps=100, float fill_value=-999.):
   """
   Takes data 
   rawdata(ntrac,data)
   where 
   rawdata(ntrac,data) is data position for trac at time t
   
   Returns
   data(trac,time,data)
   """
   
   cdef int ii, jj
   cdef int ntrac
   cdef int npoints = rawdata.shape[0]
   cdef int ndata   = rawdata.shape[1]
   cdef np.ndarray[np.float32_t, ndim=3] data = np.ones([ntracmax,nsteps,ndata], dtype=np.float32) 
   cdef np.ndarray[np.int32_t, ndim=1] steps = np.zeros([ntracmax], dtype=np.int32)
   
   data = data * fill_value
   
   for jj in xrange(0,npoints):
      ntrac = <int>rawdata[jj,0]
      step = steps[ntrac]
      for ii in range(0,ndata):
         data[ntrac,step,ii] = rawdata[jj,ii]
      steps[ntrac] += 1
      
   return data
