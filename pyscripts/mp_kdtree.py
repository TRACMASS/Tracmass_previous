import os, sys, ctypes
import multiprocessing as mp

import numpy as np
import scipy
import scipy.spatial

class Scheduler(object):
    
    def __init__(self, ndata, nprocs, chunk=None,schedule='guided'):
        if not schedule in ['guided','dynamic','static']:
            raise ValueError, 'unknown scheduling strategy'
        self._ndata = mp.RawValue(ctypes.c_int,ndata)
        self._start = mp.RawValue(ctypes.c_int,0)
        self._lock = mp.Lock()
        self._schedule = schedule
        self._nprocs = nprocs
        if schedule == 'guided' or schedule == 'dynamic':
            min_chunk = ndata // (10*nprocs)
            if chunk:
                min_chunk = chunk
            min_chunk = 1 if min_chunk < 1 else min_chunk
            self._chunk = min_chunk
        elif schedule == 'static':
            min_chunk = ndata // nprocs
            if chunk:
                min_chunk = chunk if chunk > min_chunk else min_chunk
            min_chunk = 1 if min_chunk < 1 else min_chunk
            self._chunk = min_chunk

    def __iter__(self):
        return self

    def next(self):
        self._lock.acquire()
        ndata = self._ndata.value
        nprocs = self._nprocs
        start = self._start.value
        if self._schedule == 'guided':
            _chunk = ndata // nprocs
            chunk = max(self._chunk, _chunk)
        else:
            chunk = self._chunk
        if ndata:
            if chunk > ndata:
                s0 = start
                s1 = start + ndata
                self._ndata.value = 0
            else:
                s0 = start
                s1 = start + chunk
                self._ndata.value = ndata - chunk
                self._start.value = start + chunk
            self._lock.release()
            return slice(s0, s1)
        else:
            self._lock.release()
            raise StopIteration

def tonumpyarray(mp_arr, dtype=np.float64):
    return np.frombuffer(mp_arr.get_obj(), dtype)

                
class mpKDTree( scipy.spatial.cKDTree ):
    """ Multiprocessing cKDTree subclass, shared memory """

    def __init__(self, data, leafsize=10):
        """ Init cKDTree class with shared memory
        
        Same as cKDTree.__init__ except that an internal copy of
        data to shared memory is made.
        """
        n,m = data.shape
        # Allocate shared memory for data
        self.shmem_data = mp.Array(ctypes.c_double, n*m)
        # View shared memory as ndarray, and copy over the data.
        # The RawArray objects have information about
        #the dtype and # buffer size.
        _data = tonumpyarray(self.shmem_data,np.float64).reshape((n,m))
        _data[:,:] = data
        super(mpKDTree, self).__init__(_data, leafsize=leafsize)

    def parallel_query(self, x, k=1, eps=0, p=2,
                       distance_upper_bound=np.inf, chunk=None,
                       schedule='guided'):
        """Query kdtree with multiple processes and shared memory
        
        Same as cKDTree.query except parallelized with multiple processes
        and shared memory.

        Extra keyword arguments:
        chunk:    Minimum chunk size for the load balancer.
        schedule: Strategy for balancing work load
                  ('static', 'dynamic' or 'guided').

        """
        # allocate shared memory for x and result
        nx = x.shape[0]
        shmem_x = mp.Array(ctypes.c_double, nx*self.m)
        shmem_d = mp.Array(ctypes.c_double, nx*k)
        shmem_i = mp.Array(ctypes.c_int,    nx*k)
        # view shared memory as ndarrays
        _x = tonumpyarray(shmem_x, np.float64).reshape((nx, self.m))
        _d = tonumpyarray(shmem_d, np.float64).reshape((nx, k))
        _i = tonumpyarray(shmem_i, np.int32).reshape((nx, k))
        # copy x to shared memory
        _x[:] = x
        # set up a scheduler to load balance the query
        nprocs = mp.cpu_count()
        scheduler = Scheduler(nx, nprocs, chunk=chunk, schedule=schedule)
        # return status in shared memory
        # access to these values are serialized automatically
        ierr = mp.Value(ctypes.c_int, 0)
        err_msg = mp.Array(ctypes.c_char, 1024)
        # query with multiple processes
        query_args = (scheduler,
                      self.shmem_data, self.n, self.m, self.leafsize,
                      shmem_x, nx, shmem_d, shmem_i,
                      k, eps, p, distance_upper_bound,
                      ierr)
        query_fun = _parallel_query
        pool = [mp.Process(target=query_fun, args=query_args)
                for n in range(nprocs)]
        for p in pool: p.start()
        for p in pool: p.join()
        if ierr.value != 0:
            mess =  ('%d errors in worker processes. Last one reported:\n%s'
                     % (ierr.value, err_msg.value))
            raise RuntimeError, mess
        # return results (private memory)
        return np.squeeze(_d.copy()), np.squeeze(_i.copy())

# This is executed in an external process:
def _parallel_query(scheduler,
                    data, ndata, ndim, leafsize,
                    x, nx, d, i,
                    k, eps, p, dub,
                    ierr):
                    # scheduler for load balancing
                    # data needed to reconstruct the kd-tree
                    # query data and results
                    # auxillary query parameters
                    # return values (0 on success)

    try:
        # View shared memory as ndarrays.
        _data = tonumpyarray(data, np.float64).reshape((ndata,ndim))
        _x    = tonumpyarray(x, np.float64).reshape((nx,ndim))
        _d    = np.squeeze(tonumpyarray(d, np.float64).reshape((nx,k)))
        _i    = np.squeeze(tonumpyarray(i, np.int32).reshape((nx,k)))
        # Reconstruct the kd-tree from the data.
        # This is relatively inexpensive.
        kdtree = scipy.spatial.cKDTree(_data, leafsize=leafsize)
        # Query for nearest neighbours, using slice ranges,
        # from the load balancer.
        for s in scheduler:
            _d[s],_i[s] = kdtree.query(_x[s,:], k=k, eps=eps, p=p,
                                           distance_upper_bound=dub)
    except IOError:
        ierr.value += 1
