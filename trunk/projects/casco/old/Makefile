

MYCFG=/usr/local/mysql/bin/mysql_config
MYI_FLAGS=`$(MYCFG) --cflags` 
MYL_FLAGS=`$(MYCFG) --libs` 
LIB_DIR=-L/usr/local/netcdf/lib
INC_DIR=-I/usr/local/netcdf/include -I/usr/local/mysql/include

ORM_FLAGS=-Dmaine -Dmean -Dstreamxy  -Dstreamr -Dstreamv -Dtracer \
          -Dtime  -Dtempsalt -ftrace=full -Dmysqlwrite

F90_FLAGS=-O3 -C  -g
LNK_FLAGS=-lnetcdf -lSystemStubs

FF = /sw/bin/g95 $(LIB_DIR) $(INC_DIR) $(F90_FLAGS) $(ORM_FLAGS)
CC = gcc -O  $(INC_DIR)

objects = main.o jacket.o

runtraj : $(objects)
	$(FF)  $(MYI_FLAGS) -o runtraj $(objects) $(LNK_FLAGS) $(MYL_FLAGS)

main.o : main.f95 
	$(FF)  -fno-underscoring -c -cpp $(ORM_FLAGS) main.f95

jacket.o : ../mysql/jacket.c
	$(CC) -c ../mysql/jacket.c

.PHONY : clean
clean :
	-rm runtraj $(objects)