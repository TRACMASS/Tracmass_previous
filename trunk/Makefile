
PROJECT	          = rco
CASE              = $(PROJECT)
INPUT_INT1        = intmin		
INPUT_INT2        = intrun		#Use 'dummy' if not used.



PROJECT_FLAG      = -DPROJECT_NAME=\'$(PROJECT)\'
CASE_FLAG         = -DCASE_NAME=\'$(CASE)\'
ARG_FLAGS         = -DARG_INT1=$(INPUT_INT1) -DARG_INT2=$(INPUT_INT2)

#MYCFG            = /usr/local/mysql/bin/mysql_config
#MYI_FLAGS        = `$(MYCFG) --cflags` 
#MYL_FLAGS        = `$(MYCFG) --libs` 

LIB_DIR           = -L/sw/lib -L/sw/lib/netcdf-g95/lib
INC_DIR           = -I/sw/include -I/sw/lib/netcdf-g95/include \
                    -I/usr/local/mysql/include

ORM_FLAGS         = -D$(PROJECT) -Dmean -Dstreamxy  -Dstreamr -Dstreamv  -Dtracer -Dtime  -Dtempsalt -Dmysqlwrite

F90_FLAGS         =
LNK_FLAGS         = -lnetcdf -lSystemStubs

GFORTRAN_FLAGS    = -c -integrated-cpp -fno-underscoring  
G95_FLAGS         = -c -cpp -fendian=big -fno-underscoring 
FF_FLAGS          = $(GFORTRAN_FLAGS)

###FF             = g95 $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
FF                = gfortran $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
CC                = gcc -O  $(INC_DIR)

objects           = modules.o sw_stat.o loop.o vertvel.o coord.o \
                    cross.o init_par.o interp2.o pos.o arclength.o \
                    writepsi.o writetracer.o turb.o main.o
#jacket.o

runtraj : $(objects) readfield.o
	$(FF)  $(MYI_FLAGS) -o runtraj $(objects) readfield.o $(LNK_FLAGS) $(MYL_FLAGS)

%.o : %.f95
	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT_FLAG) $(CASE_FLAG) $(ARG_FLAGS)  $< -o $@

$(objects) : 

readfield.o:  $(PROJECT)/readfield.f95
	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT)/readfield.f95


#stat.o:  $(PROJECT)/stat.f95
#	$(FF) $(FF_FLAGS) $(ORM_FLAGS) $(PROJECT)/stat.f95

jacket.o : ../mysql/jacket.c
	$(CC)  -c ../mysql/jacket.c

#main.o : main.f95 
#	$(FF) $(FF_FLAGS) $(ORM_FLAGS) main.f95

.PHONY : clean
clean :
	-rm runtraj $(objects) *.mod

