IFLAGS=-I./ -I/d/libs_extra_ame/GSL/ -I/d/libs_extra_ame/qpOASES/lib/utils/include -I./include
LFLAGS=-L./ -L/d/libs_extra_ame/lib -L/d/libs_extra_ame/qpOASES/win32


SRC=src/helper.c src/model.c src/kalman.c src/mpccontrol.c 
INCL=./include
OBJ=$(SRC:%.c=%.o) #tells make we are converting all source files to object files
BIN=./bin
LIB=./lib/libkalman.a
output=mpcprog

$(BIN)/$(output): src/$(output).o $(LIB)
	gcc src/$(output).o $(LIB) $(LFLAGS) -o $@ -lm -lgsl -lgslcblas -lOASES_win32-gcc
	
	
$(LIB):$(OBJ)
	ar rcs $@ $(OBJ)ls
%.o: %.c
	gcc -c  $< -o $@ $(IFLAGS) 

clean:
	rm src/*.o
	

	
