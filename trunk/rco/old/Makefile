#name=../../rco/run/fi.00.a0
#name=hejmake
#	rm -f $(name)

INC_DIR=-I/Users/doos/data/rco12/tmp/

ORM_FLAGS=-Drco -Dmean -Dtempsalt -Dfrsurf -Dstreamxy -Dstreamr \
          -Dstreamv -Dstreamts -Dtime -Dturb 

F90_FLAGS=-O3 -C -g -ftrace=frame -fendian=big
LNK_FLAGS= -lSystemStubs

FF = /Applications/g95-install/bin/i386-apple-darwin8.8.1-g95 $(F90_FLAGS) $(ORM_FLAGS) $(INC_DIR)
CC = gcc -O $(INC_DIR)

objects = main.o 

runtraj : $(objects)
	$(FF) -o runtraj $(objects) $(LNK_FLAGS)

main.o : main.f95 
	$(FF) -fno-underscoring -c -cpp $(ORM_FLAGS) main.f95

.PHONY : clean
clean :
	-rm runtraj $(objects)

