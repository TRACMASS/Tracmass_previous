
Environmental Variables
=======================

TRACMASS uses a couple of envronmental variables. This is paricularly useful when you run the same cases on different computers where velocity fields and output files are stored at different locations. This way, you can define paths locally without needing to change your code each time.


TRACMASS uses the following environmantal variables:


 .. envvar:: TRMPROJDIR

    Path to the project directory. Use this variable if you want to
    place your project outside the tracmass directory.

 .. envvar:: TRMDIR

    Path to the tracmass directory. Useful if you want to run
    tracmass using a script that is located outside the tracmass 
    directory, or if you want to move runtraj somewhere else.

 .. envvar:: TRMINDATADIR

    Path to where velocity fields are saved. the project name is
    automatiacally appended so that TRACMASS expects files to be
    stored at /data/gcm/ if TRMINDATADIR=/data and the project is
    gcm.

 .. envvar:: TRMOUTDATADIR

    Path to where output is saved. the project name is 
    automatiacally appended so that TRACMASS will save files at 
    /output/gcm/ if TRMOUTDATADIR=/output and the project is
    gcm.

These instructions are for :command:`bash`. You need to modify the setup if you use :command:`zsh` or, god have mercy on your soul, :command:`tcsh`. In Mac OS X and Linux, you can set the environment variables in one of the following files::

  ~/.bashrc
  ~/.bash_profile
  ~/.profile


By default, Mac OS X does not has any of these files so you need to create one of them manually. Homwebrew or macports might have created them though.

Just add the following lines::

  export TRMDIR=/path/to/tracmass
  export TRMINDATADIR=/path/to/velocity/fields
  export TRMOUTDATADIR=/path/to/output/dir
 
to any of those files. For example::

  cd ~
  touch .bash_profile
  nano .bash_profile

Copy and paste the files, change paths to your liking, and close the editor (ctrl-x y).
 
Restart the terminal (close the window and open a new terminal) and run :command:`env` to make sure that the variables are set.
