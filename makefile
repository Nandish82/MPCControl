FINAL_DIR=/d/libs_extra_ame_delta

IFLAGS=-I./ -I/d/libs_extra_ame_delta/GSL -I/d/libs_extra_ame_delta/qpOASES/lib/utils/include -I./include
LFLAGS=-L./ -L/d/libs_extra_ame_delta/lib -L/d/libs_extra_ame_delta/qpOASES/win32


SRC=src/helper.c src/model.c src/kalman.c src/mpccontrol.c 
INCL=./include
OBJ=$(SRC:%.c=%.o) #tells make we are converting all source files to object files
BIN=./bin
LIB=./lib/libdelta.a
output=mpcprog



$(FINAL_DIR)/lib/libdelta.a:$(LIB)
	cp $(LIB) $(FINAL_DIR)/lib
	cp $(INCL)/*.h $(FINAL_DIR)/include


$(BIN)/$(output): src/$(output).o $(LIB)
	gcc src/$(output).o $(LIB) $(LFLAGS) -o $@ -lm -lgsl -lgslcblas -lOASES_win32-gcc
	
	
$(LIB):$(OBJ)
	ar rcs $@ $(OBJ)
%.o: %.c
	gcc -c  $< -o $@ $(IFLAGS) 

clean:
	rm src/*.o
	

	
