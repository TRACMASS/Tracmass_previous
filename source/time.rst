
Time Conventions
================

TRACMASS uses a number of different time definitions, which can easily be confusing. Thre are three main clocks that are more or less run in parallell: switching velocity fields (:f:var:`modtime/ints`), the interpolated advetion of particles (:f:var:`modtime/tt`), and an inverse time used tofind the solution associated with sub-time interpolation (:f:var:`modtime/dt`). There are also absolute "real time" variables available to simplify reading velocity fields and post processing results.

The discrete GCM timestep (ints)
--------------------------------
:f:var:`modtime/ints` is an integer starting at the beginning of the run and increased with one each time a new velocity field is read. It is possible to convert ints to time units using :f:var:`modtime/ngcm`. This timestep also controls seeding since particles are potentially seeded with each new velocity field. The extent of seeding is also set in this unit.

The internal interpolated timestep (tt)
---------------------------------------
TRACMASS instrnal clock is mainly used as the time reference for interpolated particle postions between ints steps. tt is defined as seconds 


Absolute time
-------------
Absolute time can be acces both as Julian Dates (JD) and as calendar dates (year, month, day, hour, seconds). The associated avriables are not updated automatically since that would require a substanitial CPU overhead. Call the function :f:subr:`mod_time/updateclock` to set the variables to current model time. There are two sets of absolute time variables:

 :f:var:`modtime/curryear`, :f:var:`modtime/currmon`,   
 :f:var:`modtime/currday`, :f:var:`modtime/currmin`, 
 :f:var:`modtime/currsec`, :f:var:`modtime/currtime`

These variables follow the simulated model time with the help of  :f:var:`modtime/tt`,  :f:var:`modtime/ints`, and  :f:var:`modtime/ngcm`. They continue to increase even if the velocity fields are looped.

 :f:var:`modtime/loopyear`, :f:var:`modtime/loopmon`,   
 :f:var:`modtime/loopday`, :f:var:`modtime/loopmin`, 
 :f:var:`modtime/loopsec`, :f:var:`modtime/looptime`

These variables are are reset each time intmax or maxvelJD are reached. Use them to create the correct filename in readfield.f95 if you want to have the vability to loop velocity fields. This is useful so that particles can be advect for longer times than available GCM runs. :f:var:`modtime/loopints` is provided as an alternative to :f:var:`modtime/ints` for old code. Use the absolute time variables when setting up new projects. 
